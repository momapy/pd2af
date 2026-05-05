"""Helpers shared by placers."""

import momapy.geometry

import pd2af.predicates


def get_input_layouts_for(cd_map, input_element):
    layouts = cd_map.layout_model_mapping.get_mapping(input_element)
    if layouts is not None:
        return layouts
    for key in cd_map.layout_model_mapping.inverse.keys():
        if isinstance(key, tuple) and key[0] == input_element:
            return cd_map.layout_model_mapping.get_mapping(key)
    return None


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


def copy_subtree_singleton_mappings(source_mapping, layout_element, target_builder):
    elements = [layout_element] + list(layout_element.descendants())
    for element in elements:
        if element in source_mapping:
            target_builder.add_mapping(element, source_mapping[element])
