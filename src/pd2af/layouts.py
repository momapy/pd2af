"""Leaf primitives for layout construction.

The orchestration lives in :mod:`pd2af.build`; this module only exposes
the small pure-ish helpers it calls.
"""

import momapy.builder
import momapy.celldesigner
import momapy.core.layout
import momapy.core.mapping
import momapy.geometry

import pd2af.predicates


_SPECIES_CLASS_TO_LAYOUT_CLASS = {
    momapy.celldesigner.GenericProtein: momapy.celldesigner.GenericProteinLayout,
    momapy.celldesigner.TruncatedProtein: momapy.celldesigner.TruncatedProteinLayout,
    momapy.celldesigner.Receptor: momapy.celldesigner.ReceptorLayout,
    momapy.celldesigner.IonChannel: momapy.celldesigner.IonChannelLayout,
    momapy.celldesigner.Gene: momapy.celldesigner.GeneLayout,
    momapy.celldesigner.RNA: momapy.celldesigner.RNALayout,
    momapy.celldesigner.AntisenseRNA: momapy.celldesigner.AntisenseRNALayout,
    momapy.celldesigner.Phenotype: momapy.celldesigner.PhenotypeLayout,
    momapy.celldesigner.Ion: momapy.celldesigner.IonLayout,
    momapy.celldesigner.SimpleMolecule: momapy.celldesigner.SimpleMoleculeLayout,
    momapy.celldesigner.Drug: momapy.celldesigner.DrugLayout,
    momapy.celldesigner.Unknown: momapy.celldesigner.UnknownLayout,
    momapy.celldesigner.Complex: momapy.celldesigner.ComplexLayout,
}


def new_layout_and_mapping_builders():
    layout_builder_class = momapy.builder.get_or_make_builder_cls(
        momapy.core.layout.Layout
    )
    return (
        layout_builder_class(),
        momapy.core.mapping.LayoutModelMappingBuilder(),
    )


def clone_input_layout_and_mapping(
    map_,
    object_to_builder,
    input_model_element_to_canonical_model_element=None,
):
    # Round-tripping the mapping through `LayoutModelMappingBuilder.from_object`
    # would round-trip the model-side values too, producing fresh clones
    # with different `id()` that the writer's identity lookup can't find
    # (empty <listOfSpeciesAliases>). Build the mapping by hand: convert
    # layout-side keys to the cloned builders, keep model-side values as
    # the input frozen objects.
    layout_builder = momapy.builder.builder_from_object(
        map_.layout, object_to_builder=object_to_builder
    )
    mapping_builder = momapy.core.mapping.LayoutModelMappingBuilder()
    input_layout_model_mapping = map_.layout_model_mapping
    anchor_for_key_id = {}
    for input_anchor, input_key in input_layout_model_mapping._singleton_to_key.items():
        anchor_for_key_id[id(input_key)] = input_anchor
    for input_layout_key, input_model_value in input_layout_model_mapping.items():
        new_layout_key = _convert_layout_side(input_layout_key, object_to_builder)
        if input_model_element_to_canonical_model_element is not None:
            input_model_value = input_model_element_to_canonical_model_element.get(
                id(input_model_value), input_model_value
            )
        input_anchor = anchor_for_key_id.get(id(input_layout_key))
        if input_anchor is not None:
            new_anchor = object_to_builder.get(id(input_anchor), input_anchor)
            mapping_builder.add_mapping(
                new_layout_key, input_model_value, anchor=new_anchor
            )
        else:
            mapping_builder[new_layout_key] = input_model_value
    return layout_builder, mapping_builder


def _convert_layout_side(layout_key, object_to_builder):
    if isinstance(layout_key, frozenset):
        return frozenset(
            object_to_builder.get(id(element), element) for element in layout_key
        )
    return object_to_builder.get(id(layout_key), layout_key)


def make_synthetic_layout(species, index):
    """Build a placeholder Node for ``species``.

    ``index`` seeds the position so two synthetic layouts for
    content-equal species (same class, same name) are themselves
    content-distinct. Layout dataclasses have ``compare=False`` on
    ``id_``, so identical content collapses to one mapping key — and
    the layout-model mapping then drops every synthetic but the first.
    The position itself is throwaway: ``auto_layout`` repositions every
    node before render.
    """
    layout_class = _SPECIES_CLASS_TO_LAYOUT_CLASS.get(type(species))
    if layout_class is None:
        raise ValueError(
            f"no default layout class registered for species type "
            f"{type(species).__name__}"
        )
    position = momapy.geometry.Point(float(index), 0.0)
    label = momapy.core.layout.TextLayout(
        text=getattr(species, "name", "") or "",
        position=position,
    )
    return layout_class(position=position, label=label)


def make_modulation_arc(modulation, source_layout, target_layout):
    arc_class = pd2af.predicates.model_element_class_to_layout_element_class[
        type(modulation)
    ]
    if source_layout is target_layout:
        start_point = source_layout.anchor_point("north_north_west")
        end_point = source_layout.anchor_point("north_north_east")
    else:
        start_point = source_layout.border(target_layout.center())
        end_point = target_layout.border(source_layout.center())
        if start_point is None:
            start_point = source_layout.north_west()
        if end_point is None:
            end_point = target_layout.north_east()
    segment = momapy.geometry.Segment(start_point, end_point)
    return arc_class(
        source=source_layout,
        target=target_layout,
        segments=(segment,),
    )


def make_overlay_modulation_arc(
    context, modulation, source_layout, target_layout
):
    # Build the arc against cloned source/target builders. The arc's
    # internal Segment/Point entries get throwaway addresses; isolate
    # them in a per-arc cache so they don't poison object_to_builder.
    per_arc_cache = dict(context.object_to_builder)
    arc = make_modulation_arc(modulation, source_layout, target_layout)
    arc_builder = momapy.builder.builder_from_object(
        arc, object_to_builder=per_arc_cache
    )
    source_builder = momapy.builder.builder_from_object(
        source_layout, object_to_builder=context.object_to_builder
    )
    target_builder = momapy.builder.builder_from_object(
        target_layout, object_to_builder=context.object_to_builder
    )
    return arc_builder, source_builder, target_builder


def add_mappings_for_layout_and_descendants(
    input_layout_model_mapping,
    layout_element,
    mapping_builder,
    input_model_element_to_canonical_model_element=None,
):
    elements = [layout_element] + list(layout_element.descendants())
    for element in elements:
        if element in input_layout_model_mapping:
            model_value = input_layout_model_mapping[element]
            if input_model_element_to_canonical_model_element is not None:
                model_value = input_model_element_to_canonical_model_element.get(
                    id(model_value), model_value
                )
            mapping_builder.add_mapping(element, model_value)


def add_modulation_mapping(
    mapping_builder, arc, source_layout, target_layout, modulation
):
    source_key = mapping_builder._singleton_to_key.get(source_layout)
    source_cluster = (
        source_key if source_key is not None else frozenset([source_layout])
    )
    target_key = mapping_builder._singleton_to_key.get(target_layout)
    target_cluster = (
        target_key if target_key is not None else frozenset([target_layout])
    )
    mapping_builder.add_mapping(
        frozenset([arc]) | source_cluster | target_cluster,
        modulation,
        anchor=arc,
    )
