"""Build the SBGN-AF layout, in either ``plain`` or ``dot`` mode.

For each AF activity we build a fresh ``BiologicalActivityLayout`` /
``PhenotypeLayout`` with a typed unit-of-information sublayout and a label.
Influence arcs are drawn afresh between the activity positions (PD arcs ran
entity->process, so they cannot be copied). The ``layout_model_mapping`` is
built per the catalogue in ``momapy/sbgn/af/__init__.py`` (singleton keys for
activities and units of information; a ``frozenset`` key ``{arc, source,
target}`` anchored on the arc for each influence).

Two layout modes:

* ``plain`` -- reuse the curated SBGN-PD geometry: each activity is placed at
  its input element's layout position. An input element drawn with several
  layouts (a cloned entity pool) yields one activity layout per input layout,
  and influence arcs fan out across the source/target layouts per
  ``--influence-pairing`` (``cross`` -- one arc per pair -- or ``nearest`` -- a
  single arc between the closest pair). Activities/compartments whose input has
  no layout are skipped.
* ``dot`` -- build every element at a placeholder position with the default
  size, then hand the whole layout to ``pd2af.utils.make_auto_layout``
  (graphviz) in ``build.py`` for repositioning. Required by the merged
  ``normal`` / ``normal-no-complex`` modes, where a merged activity has no
  single input layout.

``overlay`` for SBGN-AF is unsupported (see ``pd2af.core``).
"""

import momapy.builder
import momapy.core.elements
import momapy.core.mapping
import momapy.geometry
import momapy.sbgn.af
import momapy.sbgn.io.sbgnml._reading_layout
import momapy.sbgn.layout

import pd2af.utils


_UNIT_OF_INFORMATION_CLASS_TO_LAYOUT_CLASS = {
    momapy.sbgn.af.MacromoleculeUnitOfInformation: momapy.sbgn.af.MacromoleculeUnitOfInformationLayout,
    momapy.sbgn.af.NucleicAcidFeatureUnitOfInformation: momapy.sbgn.af.NucleicAcidFeatureUnitOfInformationLayout,
    momapy.sbgn.af.SimpleChemicalUnitOfInformation: momapy.sbgn.af.SimpleChemicalUnitOfInformationLayout,
    momapy.sbgn.af.ComplexUnitOfInformation: momapy.sbgn.af.ComplexUnitOfInformationLayout,
    momapy.sbgn.af.UnspecifiedEntityUnitOfInformation: momapy.sbgn.af.UnspecifiedEntityUnitOfInformationLayout,
    momapy.sbgn.af.PerturbationUnitOfInformation: momapy.sbgn.af.PerturbationUnitOfInformationLayout,
}

_INFLUENCE_CLASS_TO_LAYOUT_CLASS = {
    momapy.sbgn.af.PositiveInfluence: momapy.sbgn.af.PositiveInfluenceLayout,
    momapy.sbgn.af.NegativeInfluence: momapy.sbgn.af.NegativeInfluenceLayout,
    momapy.sbgn.af.UnknownInfluence: momapy.sbgn.af.UnknownInfluenceLayout,
    momapy.sbgn.af.NecessaryStimulation: momapy.sbgn.af.NecessaryStimulationLayout,
}

_OPERATOR_CLASS_TO_LAYOUT_CLASS = {
    momapy.sbgn.af.AndOperator: momapy.sbgn.af.AndOperatorLayout,
    momapy.sbgn.af.OrOperator: momapy.sbgn.af.OrOperatorLayout,
    momapy.sbgn.af.NotOperator: momapy.sbgn.af.NotOperatorLayout,
    momapy.sbgn.af.DelayOperator: momapy.sbgn.af.DelayOperatorLayout,
}

# The operator glyph layout classes, for `isinstance_or_builder` checks (the
# layouts are builders during construction).
_OPERATOR_LAYOUT_CLASSES = tuple(_OPERATOR_CLASS_TO_LAYOUT_CLASS.values())

# Nodes (activities, phenotypes, units of information) use momapy's default
# layout-class sizes, so we never pass width/height when building them. The one
# exception is a compartment in `plain` mode, which reuses its input layout size
# so it still encloses its members at their curated positions.
#
# dot-mode placeholder position; graphviz repositions everything afterwards.
_PLACEHOLDER_POSITION = momapy.geometry.Point(0.0, 0.0)

# How far right of an activity's north-west corner the unit of information sits.
_UNIT_OF_INFORMATION_X_OFFSET = 15.0


def _builder(layout_class, **kwargs):
    return momapy.builder.get_or_make_builder_cls(layout_class)(**kwargs)


def make_and_add_layout(context):
    if context.layout_mode not in ("plain", "dot"):
        raise NotImplementedError(
            "SBGN-AF output supports the 'plain' and 'dot' layout modes "
            f"(got {context.layout_mode!r}); 'overlay' is unsupported."
        )
    context.layout = _builder(momapy.sbgn.af.SBGNAFLayout)
    context.layout_model_mapping = momapy.core.mapping.LayoutModelMappingBuilder()

    for compartment in context.model.compartments:
        _make_and_add_compartment_layout(context, compartment)
    for activity, input_element in context.activity_emissions:
        _make_and_add_activity_layout(context, activity, input_element)
    # Operators after activities (their logic arcs target activity layouts) and
    # before influences (an operator-sourced influence resolves its source
    # through the operator layout registered here).
    for operator, input_operator in context.operator_emissions:
        _make_and_add_operator_layout(context, operator, input_operator)
    for influence in context.model.influences:
        _make_and_add_influence_layout(context, influence)

    pd2af.utils.harmonize_root_layout(context.layout)


def _get_input_layouts(context, input_element):
    """Return the input map's layout elements for a model element, as a tuple
    (empty when the element has none). A cloned entity pool maps to several."""
    input_layouts = context.input_map.layout_model_mapping.get_mapping(
        input_element
    )
    return tuple(input_layouts) if input_layouts else ()


def _make_and_add_compartment_layout(context, compartment):
    if context.layout_mode == "dot":
        # Default size; graphviz fits the cluster around its members afterwards.
        compartment_layout = _builder(
            momapy.sbgn.af.CompartmentLayout, position=_PLACEHOLDER_POSITION
        )
    else:
        input_compartment = _input_compartment_for(context, compartment)
        if input_compartment is None:
            return
        input_layouts = _get_input_layouts(context, input_compartment)
        if not input_layouts:
            return
        input_layout = input_layouts[0]
        # Reuse the input layout's size so the compartment still encloses its
        # members at their curated positions.
        compartment_layout = _builder(
            momapy.sbgn.af.CompartmentLayout,
            position=input_layout.position,
            width=input_layout.width,
            height=input_layout.height,
        )
    compartment_layout.label = (
        momapy.sbgn.io.sbgnml._reading_layout.make_text_layout(
            compartment.label, compartment_layout.position
        )
    )
    context.layout.layout_elements.append(compartment_layout)
    context.layout_model_mapping.add_mapping(compartment_layout, compartment)
    context.model_element_to_layout_elements[id(compartment)] = (compartment_layout,)


def _input_compartment_for(context, af_compartment):
    return context.af_compartment_to_input_compartment.get(id(af_compartment))


def _make_and_add_activity_layout(context, activity, input_element):
    # The activity uses momapy's default size; only its position(s) differ by
    # mode. In dot a single placeholder graphviz repositions; in plain the
    # curated input layout positions -- one activity layout per input layout, so
    # a cloned entity pool keeps each of its placements.
    if context.layout_mode == "dot":
        positions = (_PLACEHOLDER_POSITION,)
    else:
        input_layouts = _get_input_layouts(context, input_element)
        if not input_layouts:
            return
        positions = tuple(input_layout.position for input_layout in input_layouts)
    activity_layout_class = (
        momapy.sbgn.af.PhenotypeLayout
        if isinstance(activity, momapy.sbgn.af.Phenotype)
        else momapy.sbgn.af.BiologicalActivityLayout
    )
    activity_layouts = []
    for position in positions:
        activity_layout = _builder(activity_layout_class, position=position)
        activity_layout.label = (
            momapy.sbgn.io.sbgnml._reading_layout.make_text_layout(
                activity.label, position
            )
        )
        context.layout.layout_elements.append(activity_layout)
        context.layout_model_mapping.add_mapping(activity_layout, activity)
        for unit_of_information in getattr(
            activity, "units_of_information", frozenset()
        ):
            _make_and_add_unit_of_information_layout(
                context, activity_layout, unit_of_information
            )
        activity_layouts.append(activity_layout)
    context.model_element_to_layout_elements[id(activity)] = tuple(activity_layouts)


def _make_and_add_unit_of_information_layout(
    context, activity_layout, unit_of_information
):
    layout_class = _UNIT_OF_INFORMATION_CLASS_TO_LAYOUT_CLASS.get(
        type(unit_of_information)
    )
    if layout_class is None:
        return
    # Default size; place the UoI on the activity's top border (the SBGN
    # convention), offset right of the north-west corner.
    north_west = activity_layout.north_west()
    unit_layout = _builder(
        layout_class,
        position=momapy.geometry.Point(north_west.x + _UNIT_OF_INFORMATION_X_OFFSET, north_west.y),
    )
    # In the merged modes the entity's unit-of-information block is carried on
    # the glyph (e.g. ``[ct:mRNA]``) rather than inlined in the activity label;
    # render it so it is written out and round-trips. Bare glyphs stay bare.
    if unit_of_information.label:
        unit_layout.label = (
            momapy.sbgn.io.sbgnml._reading_layout.make_text_layout(
                unit_of_information.label,
                unit_layout.position,
                # Match the smaller auxiliary-unit font the momapy sbgn reader
                # uses for units of information (not the default node font).
                font_size=momapy.sbgn.layout.DEFAULT_AUXILIARY_UNIT_FONT_SIZE,
            )
        )
    activity_layout.layout_elements.append(unit_layout)
    context.layout_model_mapping.add_mapping(unit_layout, unit_of_information)


def _make_and_add_operator_layout(context, operator, input_operator):
    """Build an operator glyph and its logic arcs.

    The glyph is built fresh (like every AF activity) at the input operator's
    curated position (``plain``) or a placeholder (``dot``). One
    ``LogicArcLayout`` runs from the operator to each input activity's layout
    (operator -> input, the SBGN convention). The operator maps to a frozenset
    of {glyph, logic arcs, input activity layouts} anchored on the glyph -- the
    catalogue the SBGN-AF writer expects -- and is registered in
    ``model_element_to_layout_elements`` so the influence pass can resolve an
    operator-sourced influence."""
    operator_layout_class = _OPERATOR_CLASS_TO_LAYOUT_CLASS.get(type(operator))
    if operator_layout_class is None:
        return
    if context.layout_mode == "dot":
        operator_layout = _builder(
            operator_layout_class, position=_PLACEHOLDER_POSITION
        )
        # make_auto_layout reverses the logic-arc dot edges, so graphviz ranks
        # the inputs above the operator and the target below it; a vertical,
        # left-to-right operator then points its input connector up (toward the
        # inputs) and its output connector down (toward the target).
        operator_layout.orientation = momapy.core.elements.Orientation.VERTICAL
        operator_layout.left_to_right = True
    else:
        input_glyph = _input_operator_glyph(context, input_operator)
        if input_glyph is None:
            return
        operator_layout = _builder(
            operator_layout_class, position=input_glyph.position
        )
        # plain: inherit the curated input operator's connector geometry, so the
        # arcs meet the same connectors (e.g. a vertical operator with ports
        # up/down) the input map drew.
        operator_layout.orientation = getattr(
            input_glyph, "orientation", operator_layout.orientation
        )
        operator_layout.left_to_right = getattr(
            input_glyph, "left_to_right", operator_layout.left_to_right
        )
        operator_layout.left_connector_length = getattr(
            input_glyph, "left_connector_length", operator_layout.left_connector_length
        )
        operator_layout.right_connector_length = getattr(
            input_glyph,
            "right_connector_length",
            operator_layout.right_connector_length,
        )
    context.layout.layout_elements.append(operator_layout)
    logic_arcs = []
    input_layouts = []
    for operator_input in operator.inputs:
        input_activity_layouts = context.model_element_to_layout_elements.get(
            id(operator_input.referred_element)
        )
        if not input_activity_layouts:
            continue
        input_layout = input_activity_layouts[0]
        arc = _make_logic_arc(operator_layout, input_layout)
        context.layout.layout_elements.append(arc)
        logic_arcs.append(arc)
        input_layouts.append(input_layout)
    frozenset_key = frozenset([operator_layout, *logic_arcs, *input_layouts])
    context.layout_model_mapping.add_mapping(
        frozenset_key, operator, anchor=operator_layout
    )
    context.model_element_to_layout_elements[id(operator)] = (operator_layout,)


def _input_operator_glyph(context, input_operator):
    """Return the input SBGN-PD operator's glyph layout (the node anchoring its
    frozenset mapping), or ``None`` when the input operator has no layout."""
    mapping = context.input_map.layout_model_mapping.get_mapping(input_operator)
    if not mapping:
        return None
    for layout_key in mapping:
        if isinstance(layout_key, frozenset):
            for element in layout_key:
                if (
                    context.input_map.layout_model_mapping.get_mapping(element)
                    is input_operator
                ):
                    return element
        else:
            return layout_key
    return None


def _operator_connector_segment(operator_layout, other_layout, is_logic_arc):
    """Segment between a logical operator and one of its arc endpoints, attached
    to the operator's *connector tip* rather than its circle border.

    A logic arc (operator -> input) meets the **input** connector; an influence
    arc (operator -> target) leaves the **output** connector -- the opposite
    side. Which physical side that is depends on the operator's ``left_to_right``
    (and ``direction``), mirroring momapy's own connector selection in
    ``momapy.sbgn.utils.set_arcs_to_borders``. The operator is the arc's source
    in both cases, so its connector tip is the segment's start point."""
    if is_logic_arc:
        tip = (
            operator_layout.left_connector_tip()
            if operator_layout.left_to_right
            else operator_layout.right_connector_tip()
        )
    else:
        tip = (
            operator_layout.right_connector_tip()
            if operator_layout.left_to_right
            else operator_layout.left_connector_tip()
        )
    other_point = other_layout.border(tip)
    if other_point is None:
        other_point = other_layout.center()
    return momapy.geometry.Segment(tip, other_point)


def resolve_operator_arc_segments(arc, source_builder, target_builder):
    """Per-arc hook for :func:`pd2af.utils.make_auto_layout`'s arc-geometry
    step: when ``arc`` is sourced by an operator glyph, return its
    connector-attached segments recomputed from the graphviz-repositioned
    geometry; otherwise return ``None`` so the caller keeps its normal border
    geometry. This reuses the same connector helper as plain-mode arc creation,
    so the auto snap is folded into the recompute that already iterates every
    arc."""
    if not momapy.builder.isinstance_or_builder(
        source_builder, _OPERATOR_LAYOUT_CLASSES
    ):
        return None
    is_logic_arc = momapy.builder.isinstance_or_builder(
        arc, momapy.sbgn.af.LogicArcLayout
    )
    return [
        _operator_connector_segment(source_builder, target_builder, is_logic_arc)
    ]


def _make_logic_arc(operator_layout, input_layout):
    segment = _operator_connector_segment(
        operator_layout, input_layout, is_logic_arc=True
    )
    return _builder(
        momapy.sbgn.af.LogicArcLayout,
        source=operator_layout,
        target=input_layout,
        segments=(segment,),
    )


def _make_and_add_influence_layout(context, influence):
    source_layouts = context.model_element_to_layout_elements.get(
        id(influence.source)
    )
    target_layouts = context.model_element_to_layout_elements.get(
        id(influence.target)
    )
    if not source_layouts or not target_layouts:
        return
    prefer_nearest = (
        context.influence_pairing == "nearest" and context.layout_mode == "plain"
    )
    for source_layout, target_layout in pd2af.utils.influence_layout_pairs(
        source_layouts, target_layouts, prefer_nearest
    ):
        arc = _make_influence_arc(influence, source_layout, target_layout)
        context.layout.layout_elements.append(arc)
        context.layout_model_mapping.add_mapping(
            frozenset([arc, source_layout, target_layout]), influence, anchor=arc
        )


def _make_influence_arc(influence, source_layout, target_layout):
    arc_class = _INFLUENCE_CLASS_TO_LAYOUT_CLASS[type(influence)]
    # An operator-sourced influence leaves the operator's output connector; every
    # other influence runs plain border-to-border (untouched).
    if momapy.builder.isinstance_or_builder(source_layout, _OPERATOR_LAYOUT_CLASSES):
        segments = (
            _operator_connector_segment(
                source_layout, target_layout, is_logic_arc=False
            ),
        )
    else:
        segments = pd2af.utils.make_arc_segments_from_source_and_target(
            source_layout, target_layout
        )
    return _builder(
        arc_class, source=source_layout, target=target_layout, segments=tuple(segments)
    )
