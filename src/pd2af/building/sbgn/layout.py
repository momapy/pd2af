"""Build the SBGN-AF layout, in either ``plain`` or ``dot`` mode.

For each AF activity we build a fresh ``BiologicalActivityLayout`` /
``PhenotypeLayout`` with a typed unit-of-information sublayout and a label.
Influence arcs are drawn afresh between the activity positions (PD arcs ran
entity->process, so they cannot be copied). Every glyph is paired with the
model element it draws, and ``momapy.core.make_map`` computes the keys per the
catalogue in ``momapy/sbgn/af/__init__.py``.

Two layout modes:

* ``plain`` -- reuse the curated SBGN-PD geometry: each activity is placed at
  its input element's layout position. An input element drawn with several
  layouts (a cloned entity pool) yields one activity layout per input layout,
  and influence arcs fan out across the source/target layouts per
  ``--influence-pairing`` (``cross`` -- one arc per pair -- or ``nearest`` -- a
  single arc between the closest pair). Activities/compartments whose input has
  no layout are skipped.
* ``dot`` -- build every element at its own placeholder position with the
  default size, then hand the whole layout to
  ``pd2af.building.layout.make_auto_layout`` (graphviz) for repositioning. The
  placeholder positions differ so that no two glyphs are equal: ``make_map``
  would merge them. Required by the merged
  reading (``keep_species`` off), where a merged activity has no
  single input layout.

``overlay`` for SBGN-AF is unsupported (see ``pd2af.core``).
"""

import dataclasses
import typing

import momapy.builder
import momapy.core.elements
import momapy.geometry
import momapy.sbgn.af
import momapy.sbgn.io.sbgnml._reading_layout
import momapy.sbgn.layout

import pd2af.building.context
import pd2af.building.layout

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
# layouts are builders during the dot repositioning).
_OPERATOR_LAYOUT_CLASSES = tuple(_OPERATOR_CLASS_TO_LAYOUT_CLASS.values())

# Nodes (activities, phenotypes, units of information) use momapy's default
# layout-class sizes, so we never pass width/height when building them. The one
# exception is a compartment in `plain` mode, which reuses its input layout size
# so it still encloses its members at their curated positions.

# How far right of an activity's north-west corner the unit of information sits.
_UNIT_OF_INFORMATION_X_OFFSET = 15.0


def _get_placeholder_position(
    context: pd2af.building.context.BuilderContext,
) -> momapy.geometry.Point:
    """The next dot-mode placeholder position; graphviz repositions it afterwards."""
    position = momapy.geometry.Point(float(context.synthetic_index), 0.0)
    context.synthetic_index += 1
    return position


def make_and_add_layout(context: pd2af.building.context.BuilderContext):
    """Fill ``context.layout_elements`` and ``context.layout_model_pairs`` (pass 2)."""
    if context.layout_mode not in ("plain", "dot"):
        raise NotImplementedError(
            "SBGN-AF output supports the 'plain' and 'dot' layout modes "
            f"(got {context.layout_mode!r}); 'overlay' is unsupported."
        )
    for compartment in _get_model_elements_of_class(
        context, momapy.sbgn.af.Compartment
    ):
        _make_and_add_compartment_layout(context, compartment)
    for activity, input_element in context.activity_emissions:
        _make_and_add_activity_layout(context, activity, input_element)
    # Operators after activities (their logic arcs target activity layouts) and
    # before influences (an operator-sourced influence resolves its source
    # through the operator layout registered here).
    for operator, input_operator in context.operator_emissions:
        _make_and_add_operator_layout(context, operator, input_operator)
    for influence in _get_model_elements_of_class(context, momapy.sbgn.af.Influence):
        _make_and_add_influence_layout(context, influence)


def _get_model_elements_of_class(
    context: pd2af.building.context.BuilderContext, model_element_class: type
) -> list:
    return [
        model_element
        for model_element in context.model_elements
        if isinstance(model_element, model_element_class)
    ]


def _get_input_layouts(
    context: pd2af.building.context.BuilderContext, input_element: typing.Any
) -> typing.Any:
    """Return the input map's layout elements for a model element, as a tuple.

    Empty when the element has none; a cloned entity pool maps to several.
    """
    return tuple(
        pd2af.building.layout.get_layout_elements_sorted_by_id(
            context.input_map.layout_model_mapping, input_element
        )
    )


def _make_and_add_compartment_layout(
    context: pd2af.building.context.BuilderContext, compartment: typing.Any
):
    if context.layout_mode == "dot":
        # Default size; graphviz fits the cluster around its members afterwards.
        position = _get_placeholder_position(context)
        size = {}
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
        position = input_layout.position
        size = {"width": input_layout.width, "height": input_layout.height}
    compartment_layout = momapy.sbgn.af.CompartmentLayout(
        position=position,
        label=momapy.sbgn.io.sbgnml._reading_layout.make_text_layout(
            compartment.label, position
        ),
        **size,
    )
    context.layout_elements.append(compartment_layout)
    context.layout_model_pairs.append((compartment_layout, compartment))
    context.model_element_to_layout_elements[id(compartment)] = (compartment_layout,)


def _input_compartment_for(
    context: pd2af.building.context.BuilderContext, af_compartment: typing.Any
) -> typing.Any:
    return context.af_compartment_to_input_compartment.get(id(af_compartment))


def _make_and_add_activity_layout(
    context: pd2af.building.context.BuilderContext,
    activity: typing.Any,
    input_element: typing.Any,
):
    # The activity uses momapy's default size; only its position(s) differ by
    # mode. In dot a single placeholder graphviz repositions; in plain the
    # curated input layout positions -- one activity layout per input layout, so
    # a cloned entity pool keeps each of its placements.
    if context.layout_mode == "dot":
        positions = (_get_placeholder_position(context),)
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
        activity_layout = activity_layout_class(
            position=position,
            label=momapy.sbgn.io.sbgnml._reading_layout.make_text_layout(
                activity.label, position
            ),
        )
        unit_layout_pairs = []
        for unit_of_information in getattr(
            activity, "units_of_information", frozenset()
        ):
            unit_layout = _make_unit_of_information_layout(
                activity_layout, unit_of_information
            )
            if unit_layout is not None:
                unit_layout_pairs.append((unit_layout, unit_of_information))
        activity_layout = dataclasses.replace(
            activity_layout,
            layout_elements=tuple(unit_layout for unit_layout, _ in unit_layout_pairs),
        )
        context.layout_elements.append(activity_layout)
        context.layout_model_pairs.append((activity_layout, activity))
        context.layout_model_pairs.extend(unit_layout_pairs)
        activity_layouts.append(activity_layout)
    context.model_element_to_layout_elements[id(activity)] = tuple(activity_layouts)


def _make_unit_of_information_layout(
    activity_layout: typing.Any, unit_of_information: typing.Any
) -> typing.Any:
    layout_class = _UNIT_OF_INFORMATION_CLASS_TO_LAYOUT_CLASS.get(
        type(unit_of_information)
    )
    if layout_class is None:
        return None
    # Default size; place the UoI on the activity's top border (the SBGN
    # convention), offset right of the north-west corner.
    north_west = activity_layout.north_west()
    position = momapy.geometry.Point(
        north_west.x + _UNIT_OF_INFORMATION_X_OFFSET, north_west.y
    )
    # In the merged modes the entity's unit-of-information block is carried on
    # the glyph (e.g. ``[ct:mRNA]``) rather than inlined in the activity label;
    # render it so it is written out and round-trips. Bare glyphs stay bare.
    label = None
    if unit_of_information.label:
        label = momapy.sbgn.io.sbgnml._reading_layout.make_text_layout(
            unit_of_information.label,
            position,
            # Match the smaller auxiliary-unit font the momapy sbgn reader
            # uses for units of information (not the default node font).
            font_size=momapy.sbgn.layout.DEFAULT_AUXILIARY_UNIT_FONT_SIZE,
        )
    return layout_class(position=position, label=label)


def _make_and_add_operator_layout(
    context: pd2af.building.context.BuilderContext,
    operator: typing.Any,
    input_operator: typing.Any,
):
    """Build an operator glyph and its logic arcs.

    The glyph is built fresh (like every AF activity) at the input operator's
    curated position (``plain``) or a placeholder (``dot``). One
    ``LogicArcLayout`` runs from the operator to each input activity's layout
    (operator -> input, the SBGN convention) and is paired with its operator
    input, so ``momapy.core.make_map`` puts the logic arcs and the input
    activity layouts in the operator's key. The glyph is registered in
    ``model_element_to_layout_elements`` so the influence pass can resolve an
    operator-sourced influence.
    """
    operator_layout_class = _OPERATOR_CLASS_TO_LAYOUT_CLASS.get(type(operator))
    if operator_layout_class is None:
        return
    if context.layout_mode == "dot":
        # make_auto_layout reverses the logic-arc dot edges, so graphviz ranks
        # the inputs above the operator and the target below it; a vertical,
        # left-to-right operator then points its input connector up (toward the
        # inputs) and its output connector down (toward the target).
        operator_layout = operator_layout_class(
            position=_get_placeholder_position(context),
            orientation=momapy.core.elements.Orientation.VERTICAL,
            left_to_right=True,
        )
    else:
        input_glyph = _input_operator_glyph(context, input_operator)
        if input_glyph is None:
            return
        operator_layout = operator_layout_class(position=input_glyph.position)
        # plain: inherit the curated input operator's connector geometry, so the
        # arcs meet the same connectors (e.g. a vertical operator with ports
        # up/down) the input map drew.
        operator_layout = dataclasses.replace(
            operator_layout,
            orientation=getattr(
                input_glyph, "orientation", operator_layout.orientation
            ),
            left_to_right=getattr(
                input_glyph, "left_to_right", operator_layout.left_to_right
            ),
            left_connector_length=getattr(
                input_glyph,
                "left_connector_length",
                operator_layout.left_connector_length,
            ),
            right_connector_length=getattr(
                input_glyph,
                "right_connector_length",
                operator_layout.right_connector_length,
            ),
        )
    context.layout_elements.append(operator_layout)
    context.layout_model_pairs.append((operator_layout, operator))
    for operator_input in operator.inputs:
        input_activity_layouts = context.model_element_to_layout_elements.get(
            id(operator_input.referred_element)
        )
        if not input_activity_layouts:
            continue
        arc = _make_logic_arc(context, operator_layout, input_activity_layouts[0])
        context.layout_elements.append(arc)
        context.layout_model_pairs.append((arc, operator_input))
    context.model_element_to_layout_elements[id(operator)] = (operator_layout,)


def _input_operator_glyph(
    context: pd2af.building.context.BuilderContext, input_operator: typing.Any
) -> typing.Any:
    """Return the input SBGN-PD operator's glyph layout, or ``None``.

    When the operator is drawn several times, the glyph with the smallest
    ``id_`` is returned.
    """
    input_operator_layouts = pd2af.building.layout.get_layout_elements_sorted_by_id(
        context.input_map.layout_model_mapping, input_operator
    )
    if not input_operator_layouts:
        return None
    return input_operator_layouts[0]


def _operator_connector_segment(
    operator_layout: typing.Any, other_layout: typing.Any, is_logic_arc: bool
) -> typing.Any:
    """Segment between a logical operator and one of its arc endpoints.

    It attaches to the operator's *connector tip* rather than its circle border.

    A logic arc (operator -> input) meets the **input** connector; an influence
    arc (operator -> target) leaves the **output** connector -- the opposite
    side. Which physical side that is depends on the operator's ``left_to_right``
    (and ``direction``), mirroring momapy's own connector selection in
    ``momapy.sbgn.utils.set_arcs_to_borders``. The operator is the arc's source
    in both cases, so its connector tip is the segment's start point.
    """
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


def resolve_operator_arc_segments(
    arc: typing.Any, source_builder: typing.Any, target_builder: typing.Any
) -> typing.Any:
    """Per-arc hook for :func:`pd2af.building.layout.make_auto_layout`'s arc-geometry step.

    When ``arc`` is sourced by an operator glyph, return its connector-attached
    segments recomputed from the graphviz-repositioned geometry; otherwise
    return ``None`` so the caller keeps its normal border geometry. This reuses
    the same connector helper as plain-mode arc creation, so the auto snap is
    folded into the recompute that already iterates every arc.
    """
    if not momapy.builder.isinstance_or_builder(
        source_builder, _OPERATOR_LAYOUT_CLASSES
    ):
        return None
    is_logic_arc = momapy.builder.isinstance_or_builder(
        arc, momapy.sbgn.af.LogicArcLayout
    )
    return [_operator_connector_segment(source_builder, target_builder, is_logic_arc)]


def _make_logic_arc(
    context: pd2af.building.context.BuilderContext,
    operator_layout: typing.Any,
    input_layout: typing.Any,
) -> typing.Any:
    if context.layout_mode == "dot":
        segments = pd2af.building.layout.PLACEHOLDER_ARC_SEGMENTS
    else:
        segments = (
            _operator_connector_segment(
                operator_layout, input_layout, is_logic_arc=True
            ),
        )
    return momapy.sbgn.af.LogicArcLayout(
        source=operator_layout,
        target=input_layout,
        segments=segments,
    )


def _make_and_add_influence_layout(
    context: pd2af.building.context.BuilderContext, influence: typing.Any
):
    source_layouts = context.model_element_to_layout_elements.get(id(influence.source))
    target_layouts = context.model_element_to_layout_elements.get(id(influence.target))
    if not source_layouts or not target_layouts:
        return
    prefer_nearest = (
        context.influence_pairing == "nearest" and context.layout_mode == "plain"
    )
    for source_layout, target_layout in pd2af.building.layout.influence_layout_pairs(
        source_layouts, target_layouts, prefer_nearest
    ):
        arc = _make_influence_arc(context, influence, source_layout, target_layout)
        context.layout_elements.append(arc)
        context.layout_model_pairs.append((arc, influence))


def _make_influence_arc(
    context: pd2af.building.context.BuilderContext,
    influence: typing.Any,
    source_layout: typing.Any,
    target_layout: typing.Any,
) -> typing.Any:
    arc_class = _INFLUENCE_CLASS_TO_LAYOUT_CLASS[type(influence)]
    # An operator-sourced influence leaves the operator's output connector; every
    # other influence runs plain border-to-border (untouched).
    if context.layout_mode == "dot":
        segments = pd2af.building.layout.PLACEHOLDER_ARC_SEGMENTS
    elif momapy.builder.isinstance_or_builder(source_layout, _OPERATOR_LAYOUT_CLASSES):
        segments = (
            _operator_connector_segment(
                source_layout, target_layout, is_logic_arc=False
            ),
        )
    else:
        segments = pd2af.building.layout.make_arc_segments_from_source_and_target(
            source_layout, target_layout
        )
    return arc_class(
        source=source_layout, target=target_layout, segments=tuple(segments)
    )
