"""Build the SBGN-AF layout, in either ``plain`` or ``auto`` mode.

For each AF activity we build a fresh ``BiologicalActivityLayout`` /
``PhenotypeLayout`` with a typed unit-of-information sublayout and a label.
Influence arcs are drawn afresh between the activity positions (PD arcs ran
entity->process, so they cannot be copied). The ``layout_model_mapping`` is
built per the catalogue in ``momapy/sbgn/af/__init__.py`` (singleton keys for
activities and units of information; a ``frozenset`` key ``{arc, source,
target}`` anchored on the arc for each influence).

Two layout modes:

* ``plain`` -- reuse the curated SBGN-PD geometry: each activity is placed at
  its input PD entity's glyph position and size. Activities/compartments whose
  input has no glyph are skipped. Used for the keep-species modes (every
  activity maps back to exactly one input glyph).
* ``auto`` -- build every element at a placeholder position with the default
  size, then hand the whole layout to ``pd2af.utils.auto_layout`` (graphviz) in
  ``build.py`` for repositioning. Required by the merged ``normal`` /
  ``no-complex`` modes, where a merged activity has no single input glyph.

``overlay`` for SBGN-AF is unsupported (see ``pd2af.core``).
"""

import momapy.builder
import momapy.core.mapping
import momapy.geometry
import momapy.sbgn.af
import momapy.sbgn.io.sbgnml._reading_layout

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

# Nodes (activities, phenotypes, units of information) use momapy's default
# layout-class sizes, so we never pass width/height when building them. The one
# exception is a compartment in `plain` mode, which reuses its input glyph size
# so it still encloses its members at their curated positions.
#
# auto-mode placeholder position; graphviz repositions everything afterwards.
_PLACEHOLDER_POSITION = momapy.geometry.Point(0.0, 0.0)

# How far right of an activity's north-west corner the unit of information sits.
_UNIT_OF_INFORMATION_X_OFFSET = 15.0


def _builder(layout_class, **kwargs):
    return momapy.builder.get_or_make_builder_cls(layout_class)(**kwargs)


def make_and_add_layout(context):
    if context.layout_mode not in ("plain", "auto"):
        raise NotImplementedError(
            "SBGN-AF output supports the 'plain' and 'auto' layout modes "
            f"(got {context.layout_mode!r}); 'overlay' is unsupported."
        )
    context.layout = _builder(momapy.sbgn.af.SBGNAFLayout)
    context.layout_model_mapping = momapy.core.mapping.LayoutModelMappingBuilder()

    for compartment in context.model.compartments:
        _make_and_add_compartment_layout(context, compartment)
    for activity, input_element in context.activity_emissions:
        _make_and_add_activity_layout(context, activity, input_element)
    for influence in context.model.influences:
        _make_and_add_influence_layout(context, influence)

    pd2af.utils.harmonize_root_layout(context.layout)


def _input_glyph(context, input_element):
    glyphs = context.input_map.layout_model_mapping.get_mapping(input_element)
    if not glyphs:
        return None
    return glyphs[0]


def _make_and_add_compartment_layout(context, compartment):
    if context.layout_mode == "auto":
        # Default size; graphviz fits the cluster around its members afterwards.
        compartment_layout = _builder(
            momapy.sbgn.af.CompartmentLayout, position=_PLACEHOLDER_POSITION
        )
    else:
        input_compartment = _input_compartment_for(context, compartment)
        if input_compartment is None:
            return
        input_glyph = _input_glyph(context, input_compartment)
        if input_glyph is None:
            return
        # Reuse the input glyph size so the compartment still encloses its
        # members at their curated positions.
        compartment_layout = _builder(
            momapy.sbgn.af.CompartmentLayout,
            position=input_glyph.position,
            width=input_glyph.width,
            height=input_glyph.height,
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
    # The activity uses momapy's default size; only its position differs by mode
    # (curated input glyph position in plain, placeholder in auto).
    if context.layout_mode == "auto":
        position = _PLACEHOLDER_POSITION
    else:
        input_glyph = _input_glyph(context, input_element)
        if input_glyph is None:
            return
        position = input_glyph.position
    activity_layout_class = (
        momapy.sbgn.af.PhenotypeLayout
        if isinstance(activity, momapy.sbgn.af.Phenotype)
        else momapy.sbgn.af.BiologicalActivityLayout
    )
    activity_layout = _builder(activity_layout_class, position=position)
    activity_layout.label = (
        momapy.sbgn.io.sbgnml._reading_layout.make_text_layout(
            activity.label, position
        )
    )
    context.layout.layout_elements.append(activity_layout)
    context.layout_model_mapping.add_mapping(activity_layout, activity)
    context.model_element_to_layout_elements[id(activity)] = (activity_layout,)

    for unit_of_information in getattr(activity, "units_of_information", frozenset()):
        _make_and_add_unit_of_information_layout(
            context, activity_layout, unit_of_information
        )


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
    activity_layout.layout_elements.append(unit_layout)
    context.layout_model_mapping.add_mapping(unit_layout, unit_of_information)


def _make_and_add_influence_layout(context, influence):
    source_layouts = context.model_element_to_layout_elements.get(
        id(influence.source)
    )
    target_layouts = context.model_element_to_layout_elements.get(
        id(influence.target)
    )
    if not source_layouts or not target_layouts:
        return
    source_layout = source_layouts[0]
    target_layout = target_layouts[0]
    arc = _make_influence_arc(influence, source_layout, target_layout)
    context.layout.layout_elements.append(arc)
    context.layout_model_mapping.add_mapping(
        frozenset([arc, source_layout, target_layout]), influence, anchor=arc
    )


def _make_influence_arc(influence, source_layout, target_layout):
    arc_class = _INFLUENCE_CLASS_TO_LAYOUT_CLASS[type(influence)]
    if source_layout is target_layout:
        start_point = source_layout.anchor_point("north_north_west")
        end_point = source_layout.anchor_point("north_north_east")
    else:
        start_point = source_layout.own_border(target_layout.center())
        end_point = target_layout.own_border(source_layout.center())
        if start_point is None:
            start_point = source_layout.north_west()
        if end_point is None:
            end_point = target_layout.north_east()
    segment = momapy.geometry.Segment(start_point, end_point)
    return _builder(
        arc_class, source=source_layout, target=target_layout, segments=(segment,)
    )
