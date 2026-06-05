"""Build the SBGN-AF layout by reusing the curated SBGN-PD geometry (plain).

For each AF activity we look up the input PD entity's glyph and build a fresh
``BiologicalActivityLayout`` / ``PhenotypeLayout`` at the *same* position and
size, with a typed unit-of-information sublayout and a label. Influence arcs are
drawn afresh between the reused activity positions (PD arcs ran entity->process,
so they cannot be copied). The ``layout_model_mapping`` is built per the
catalogue in ``momapy/sbgn/af/__init__.py`` (singleton keys for activities and
units of information; a ``frozenset`` key ``{arc, source, target}`` anchored on
the arc for each influence).

Only ``plain`` layout is implemented; ``auto``/``overlay`` for SBGN-AF are
future work (see the plan / ``pd2af.utils.auto_layout``).
"""

import momapy.builder
import momapy.core.layout
import momapy.core.mapping
import momapy.geometry
import momapy.sbgn.af

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

_UNIT_OF_INFORMATION_SIZE = 12.0


def _builder(layout_class, **kwargs):
    return momapy.builder.get_or_make_builder_cls(layout_class)(**kwargs)


def make_and_add_layout(context):
    if context.layout_mode != "plain":
        raise NotImplementedError(
            "SBGN-AF output currently supports only the 'plain' layout mode "
            f"(got {context.layout_mode!r}); 'auto'/'overlay' are future work."
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
    input_compartment = _input_compartment_for(context, compartment)
    if input_compartment is None:
        return
    input_glyph = _input_glyph(context, input_compartment)
    if input_glyph is None:
        return
    compartment_layout = _builder(
        momapy.sbgn.af.CompartmentLayout,
        position=input_glyph.position,
        width=input_glyph.width,
        height=input_glyph.height,
    )
    context.layout.layout_elements.append(compartment_layout)
    context.layout_model_mapping.add_mapping(compartment_layout, compartment)
    context.model_element_to_layout_elements[id(compartment)] = (compartment_layout,)


def _input_compartment_for(context, af_compartment):
    for input_id, canonical in context.input_compartment_to_af_compartment.items():
        if canonical is af_compartment:
            return context.clingo_id_to_model_element.get(input_id)
    return None


def _make_and_add_activity_layout(context, activity, input_element):
    input_glyph = _input_glyph(context, input_element)
    if input_glyph is None:
        return
    position = input_glyph.position
    width = input_glyph.width
    height = input_glyph.height
    activity_layout_class = (
        momapy.sbgn.af.PhenotypeLayout
        if isinstance(activity, momapy.sbgn.af.Phenotype)
        else momapy.sbgn.af.BiologicalActivityLayout
    )
    activity_layout = _builder(
        activity_layout_class, position=position, width=width, height=height
    )
    activity_layout.label = momapy.core.layout.TextLayout(
        text=activity.label or "", position=position
    )
    context.layout.layout_elements.append(activity_layout)
    context.layout_model_mapping.add_mapping(activity_layout, activity)
    context.model_element_to_layout_elements[id(activity)] = (activity_layout,)

    for unit_of_information in getattr(activity, "units_of_information", frozenset()):
        _make_and_add_unit_of_information_layout(
            context, activity_layout, unit_of_information, position, width, height
        )


def _make_and_add_unit_of_information_layout(
    context, activity_layout, unit_of_information, position, width, height
):
    layout_class = _UNIT_OF_INFORMATION_CLASS_TO_LAYOUT_CLASS.get(
        type(unit_of_information)
    )
    if layout_class is None:
        return
    unit_position = momapy.geometry.Point(
        position.x - width / 2 + _UNIT_OF_INFORMATION_SIZE,
        position.y - height / 2,
    )
    unit_layout = _builder(
        layout_class,
        position=unit_position,
        width=_UNIT_OF_INFORMATION_SIZE,
        height=_UNIT_OF_INFORMATION_SIZE,
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
