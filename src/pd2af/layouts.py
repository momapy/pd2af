import momapy.builder
import momapy.celldesigner
import momapy.coloring
import momapy.core.layout
import momapy.core.mapping
import momapy.geometry
import momapy.positioning

import pd2af.predicates
import pd2af.utils


# TODO: delete when issue is solved in momapy
def _get_layout_elements_for_model_element(map_, model_element):
    layout_elements = map_.layout_model_mapping.get_mapping(model_element)
    if layout_elements is None:
        for key in map_.layout_model_mapping.inverse.keys():
            if isinstance(key, tuple) and key[0] == model_element:
                layout_elements = map_.layout_model_mapping.get_mapping(key)
                break
    return layout_elements


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
            arc_builder = momapy.builder.builder_from_object(
                arc, object_to_builder=object_to_builder
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
    species_to_layout_element = {}
    for species in new_cd_model.species:
        species_layouts = _get_layout_elements_for_model_element(cd_map, species)
        for species_layout in species_layouts:
            layout_builder.layout_elements.append(species_layout)
            _copy_subtree_singleton_mappings(
                cd_map.layout_model_mapping, species_layout, mapping_builder
            )
            species_to_layout_element[species] = species_layout
            break
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
    for species in new_cd_model.species:
        species_layouts = _get_layout_elements_for_model_element(cd_map, species)
        for species_layout in species_layouts:
            layout_builder.layout_elements.append(species_layout)
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
