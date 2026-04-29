import momapy.builder
import momapy.celldesigner
import momapy.coloring
import momapy.core.layout
import momapy.core.mapping
import momapy.geometry
import momapy.positioning

import pd2af.predicates
import pd2af.solver
import pd2af.utils


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
    momapy.celldesigner.Degraded: momapy.celldesigner.DegradedLayout,
}


# TODO: delete when issue is solved in momapy
def _get_layout_elements_for_model_element(map_, model_element):
    layout_elements = map_.layout_model_mapping.get_mapping(model_element)
    if layout_elements is None:
        for key in map_.layout_model_mapping.inverse.keys():
            if isinstance(key, tuple) and key[0] == model_element:
                layout_elements = map_.layout_model_mapping.get_mapping(key)
                break
    return layout_elements


def _make_synthetic_species_layout(species):
    layout_cls = _SPECIES_CLASS_TO_LAYOUT_CLASS.get(type(species))
    if layout_cls is None:
        raise ValueError(
            f"no default layout class registered for species type "
            f"{type(species).__name__}"
        )
    label = momapy.core.layout.TextLayout(
        text=species.name or "",
        position=momapy.geometry.Point(0.0, 0.0),
    )
    return layout_cls(
        position=momapy.geometry.Point(0.0, 0.0),
        label=label,
    )


def _resolve_species_layouts(cd_map, new_cd_model):
    """Resolve `(species, layout, synthetic)` for each species in the new model.

    `synthetic=True` means no original layout could be found (or it's a
    merged-proteoform synthesized species) and a fresh layout was created.
    """
    resolved = []
    for species in new_cd_model.species:
        if pd2af.solver.is_synthesized_species(species):
            layout_elements = None
        else:
            layout_elements = _get_layout_elements_for_model_element(
                cd_map, species
            )
        if layout_elements:
            resolved.append((species, layout_elements[0], False))
        else:
            resolved.append((species, _make_synthetic_species_layout(species), True))
    return resolved


def _partition_top_level_species_layouts(resolved):
    """Split resolved species layouts into top-level vs nested.

    A species layout is "nested" if it appears as a descendant of another
    kept species's layout — i.e. a subunit of a kept complex. Such a
    species must not be promoted to a top-level layout element: its
    position is governed by its parent complex layout. CellDesigner
    layouts form a tree, so this partition is well-defined.
    """
    descendant_ids = set()
    for _, layout, _ in resolved:
        for descendant in layout.descendants():
            descendant_ids.add(id(descendant))
    top_level = []
    nested = []
    for entry in resolved:
        _, layout, synthetic = entry
        if not synthetic and id(layout) in descendant_ids:
            nested.append(entry)
        else:
            top_level.append(entry)
    return top_level, nested


def _make_modulation_arc_tuple(modulation, source_layout, target_layout):
    cls = pd2af.predicates.model_element_class_to_layout_element_class[
        type(modulation)
    ]
    start_point = source_layout.border(target_layout.center())
    end_point = target_layout.border(source_layout.center())
    if start_point is None:
        start_point = source_layout.north_west()
    if end_point is None:
        end_point = target_layout.north_east()
    segment = momapy.geometry.Segment(start_point, end_point)
    arc = cls(
        source=source_layout,
        target=target_layout,
        segments=(segment,),
    )
    return arc, source_layout, target_layout


def _iter_modulation_arc_tuples(cd_map, modulation):
    sources = _get_layout_elements_for_model_element(cd_map, modulation.source)
    targets = _get_layout_elements_for_model_element(cd_map, modulation.target)
    for source in sources:
        for target in targets:
            yield _make_modulation_arc_tuple(modulation, source, target)


def _copy_subtree_singleton_mappings(source_mapping, layout_element, target_builder):
    elements = [layout_element] + list(layout_element.descendants())
    for element in elements:
        if element in source_mapping:
            target_builder.add_mapping(element, source_mapping[element])


def _add_modulation_mapping(
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


def _finalize(new_cd_model, layout_builder, mapping_builder):
    builder_to_object = {}
    layout = momapy.builder.object_from_builder(
        layout_builder, builder_to_object=builder_to_object
    )
    mapping = momapy.builder.object_from_builder(
        mapping_builder, builder_to_object=builder_to_object
    )
    return momapy.celldesigner.CellDesignerMap(
        model=new_cd_model,
        layout=layout,
        layout_model_mapping=mapping,
    )


def make_overlay(cd_map, new_cd_model):
    object_to_builder = {}
    layout_builder = momapy.builder.builder_from_object(
        cd_map.layout, object_to_builder=object_to_builder
    )
    mapping_builder = momapy.core.mapping.LayoutModelMappingBuilder.from_object(
        cd_map.layout_model_mapping,
        object_to_builder=object_to_builder,
    )
    to_keep = set()
    for compartment in new_cd_model.compartments:
        compartment_layouts = cd_map.layout_model_mapping.get_mapping(compartment)
        if compartment_layouts is None:
            continue
        for compartment_layout in compartment_layouts:
            builder = object_to_builder.get(id(compartment_layout))
            if builder is not None:
                to_keep.add(builder)
    for species in new_cd_model.species:
        species_layouts = _get_layout_elements_for_model_element(cd_map, species)
        for species_layout in species_layouts:
            builder = object_to_builder.get(id(species_layout))
            if builder is not None:
                to_keep.add(builder)
    for modulation in new_cd_model.modulations:
        for arc, source_layout, target_layout in _iter_modulation_arc_tuples(
            cd_map, modulation
        ):
            # The shared `object_to_builder` is keyed by `id()`. The freshly
            # built `Segment`/`Point` instances inside `arc` are short-lived
            # and would have their addresses recycled across iterations,
            # causing later `builder_from_object` calls to hit stale entries.
            # A per-arc cache snapshot keeps the original-layout entries
            # (so `arc.source`/`arc.target` reuse the cached builders) while
            # confining the throwaway segment/point entries to this iteration.
            per_arc_cache = dict(object_to_builder)
            arc_builder = momapy.builder.builder_from_object(
                arc, object_to_builder=per_arc_cache
            )
            source_builder = momapy.builder.builder_from_object(
                source_layout, object_to_builder=object_to_builder
            )
            target_builder = momapy.builder.builder_from_object(
                target_layout, object_to_builder=object_to_builder
            )
            layout_builder.layout_elements.append(arc_builder)
            _add_modulation_mapping(
                mapping_builder,
                arc_builder,
                source_builder,
                target_builder,
                modulation,
            )
            to_keep.update([arc_builder, source_builder, target_builder])
    layout_builder = pd2af.utils.highlight_layout_elements(to_keep, layout_builder)
    return _finalize(new_cd_model, layout_builder, mapping_builder)


def make_auto(cd_map, new_cd_model):
    layout_builder_cls = momapy.builder.get_or_make_builder_cls(
        momapy.core.layout.Layout
    )
    layout_builder = layout_builder_cls()
    mapping_builder = momapy.core.mapping.LayoutModelMappingBuilder()
    for compartment in new_cd_model.compartments:
        compartment_layouts = cd_map.layout_model_mapping.get_mapping(compartment)
        if compartment_layouts is None:
            continue
        for compartment_layout in compartment_layouts:
            layout_builder.layout_elements.append(compartment_layout)
            _copy_subtree_singleton_mappings(
                cd_map.layout_model_mapping, compartment_layout, mapping_builder
            )
            break
    resolved = _resolve_species_layouts(cd_map, new_cd_model)
    top_level, _ = _partition_top_level_species_layouts(resolved)
    species_to_layout_element = {
        species: layout for species, layout, _ in resolved
    }
    for species, species_layout, synthetic in top_level:
        layout_builder.layout_elements.append(species_layout)
        if synthetic:
            mapping_builder.add_mapping(species_layout, species)
        else:
            _copy_subtree_singleton_mappings(
                cd_map.layout_model_mapping, species_layout, mapping_builder
            )
    for modulation in new_cd_model.modulations:
        source_layout = species_to_layout_element[modulation.source]
        target_layout = species_to_layout_element[modulation.target]
        arc, source_layout, target_layout = _make_modulation_arc_tuple(
            modulation, source_layout, target_layout
        )
        layout_builder.layout_elements.append(arc)
        _add_modulation_mapping(
            mapping_builder, arc, source_layout, target_layout, modulation
        )
    new_map = _finalize(new_cd_model, layout_builder, mapping_builder)
    return pd2af.utils.auto_layout(new_map)


def make_plain(cd_map, new_cd_model):
    layout_builder_cls = momapy.builder.get_or_make_builder_cls(
        momapy.core.layout.Layout
    )
    layout_builder = layout_builder_cls()
    mapping_builder = momapy.core.mapping.LayoutModelMappingBuilder()
    for compartment in new_cd_model.compartments:
        compartment_layouts = cd_map.layout_model_mapping.get_mapping(compartment)
        if compartment_layouts is None:
            continue
        for compartment_layout in compartment_layouts:
            layout_builder.layout_elements.append(compartment_layout)
            _copy_subtree_singleton_mappings(
                cd_map.layout_model_mapping, compartment_layout, mapping_builder
            )
    resolved = _resolve_species_layouts(cd_map, new_cd_model)
    top_level, _ = _partition_top_level_species_layouts(resolved)
    for species, species_layout, synthetic in top_level:
        layout_builder.layout_elements.append(species_layout)
        if synthetic:
            mapping_builder.add_mapping(species_layout, species)
        else:
            _copy_subtree_singleton_mappings(
                cd_map.layout_model_mapping, species_layout, mapping_builder
            )
    for modulation in new_cd_model.modulations:
        for arc, source_layout, target_layout in _iter_modulation_arc_tuples(
            cd_map, modulation
        ):
            layout_builder.layout_elements.append(arc)
            _add_modulation_mapping(
                mapping_builder, arc, source_layout, target_layout, modulation
            )
    layout_builder.fill = momapy.coloring.white
    momapy.positioning.set_fit(
        layout_builder,
        layout_builder.layout_elements,
        xsep=10.0,
        ysep=10.0,
    )
    return _finalize(new_cd_model, layout_builder, mapping_builder)


def make_no_layout(new_cd_model):
    return momapy.celldesigner.CellDesignerMap(model=new_cd_model)


_BUILDERS = {
    "plain": make_plain,
    "overlay": make_overlay,
    "auto": make_auto,
}


def build_map(cd_map, new_cd_model, layout_mode):
    if layout_mode is None:
        return make_no_layout(new_cd_model)
    builder = _BUILDERS.get(layout_mode)
    if builder is None:
        raise ValueError(f"unsupported layout mode {layout_mode!r}")
    return builder(cd_map, new_cd_model)
