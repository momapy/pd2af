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


def clone_layout_pruning_foreground(
    input_layout_element, foreground_ids, object_to_builder
):
    """Clone an input layout element as a builder for use as dimmed
    background, omitting any descendant subtree whose input original is in
    ``foreground_ids`` (already drawn in the foreground). Returns ``None``
    when the element itself belongs to the foreground.

    Clones share ``object_to_builder`` so structure shared across the input
    layout stays shared once cloned, and they are builders -- not the input
    frozen objects -- so dimming the background never mutates the input
    map's layout.
    """
    if id(input_layout_element) in foreground_ids:
        return None
    clone = momapy.builder.builder_from_object(
        input_layout_element, object_to_builder=object_to_builder
    )
    _prune_foreground_from_clone(
        input_layout_element, clone, foreground_ids, object_to_builder
    )
    return clone


def _prune_foreground_from_clone(
    input_layout_element, clone, foreground_ids, object_to_builder
):
    # The only foreground glyphs that can sit *inside* a background clone are
    # subunits nested via `layout_elements` (e.g. promoted subunits of a
    # dissolved complex). `builder_from_object` already cloned the whole
    # subtree into `object_to_builder`; rebuild `layout_elements` keeping only
    # the clones whose input original is not in the foreground, recursing so
    # deeper nestings are pruned too.
    input_subunits = getattr(input_layout_element, "layout_elements", None)
    if not input_subunits:
        return
    surviving_clones = []
    for input_subunit in input_subunits:
        if id(input_subunit) in foreground_ids:
            continue
        subunit_clone = object_to_builder[id(input_subunit)]
        _prune_foreground_from_clone(
            input_subunit, subunit_clone, foreground_ids, object_to_builder
        )
        surviving_clones.append(subunit_clone)
    clone.layout_elements = surviving_clones


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
