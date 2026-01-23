import collections

import momapy.core
import momapy.styling
import momapy.coloring
import momapy.drawing
import momapy.builder
import momapy.geometry
import momapy.positioning
import momapy.celldesigner.core
import pydot

_POINTS_PER_INCH = 96


def highlight_layout_elements(layout_elements, layout):
    all_layout_elements = []
    for layout_element in layout_elements:
        all_layout_elements.append(layout_element)
        all_layout_elements += layout_element.descendants()
    id_selectors = [
        momapy.styling.IdSelector(layout_element.id_)
        for layout_element in all_layout_elements
    ]
    not_selector = momapy.styling.NotSelector(tuple(id_selectors))
    layout_element_selector = momapy.styling.CompoundSelector(
        tuple([momapy.styling.ClassSelector("LayoutElement"), not_selector])
    )
    production_layout_selector = momapy.styling.CompoundSelector(
        tuple([momapy.styling.TypeSelector("ProductionLayout"), not_selector])
    )
    compartment_layout_selector = momapy.styling.CompoundSelector(
        tuple([momapy.styling.TypeSelector("RectangleCompartmentLayout"), not_selector])
    )
    reaction_layout_selector = momapy.styling.CompoundSelector(
        tuple(
            [
                momapy.styling.OrSelector(
                    tuple(
                        [
                            momapy.styling.TypeSelector(class_name)
                            for class_name in [
                                "StateTransitionLayout",
                                "HeterodimerAssociationLayout",
                                "KnownTransitionOmittedLayout",
                                "UnknownTransitionLayout",
                                "TransportLayout",
                                "TranslationLayout",
                                "TranscriptionLayout",
                            ]
                        ]
                    )
                ),
                not_selector,
            ]
        )
    )
    style_sheet = momapy.styling.StyleSheet(
        {
            layout_element_selector: momapy.styling.StyleCollection(
                {
                    "stroke": None,
                    "fill": momapy.drawing.NoneValue,
                    "path_stroke": None,
                    "end_arrowhead_stroke": None,
                    "start_arrowhead_stroke": None,
                    "arrowhead_stroke": None,
                    "reaction_node_stroke": None,
                    "active_stroke": None,
                    "inner_stroke": None,
                    "group_stroke": momapy.coloring.lightgray,
                }
            ),
            production_layout_selector: momapy.styling.StyleCollection(
                {
                    "arrowhead_stroke": momapy.coloring.lightgray,
                    "arrowhead_fill": momapy.coloring.lightgray,
                }
            ),
            compartment_layout_selector: momapy.styling.StyleCollection(
                {
                    "fill": momapy.coloring.lightgray,
                }
            ),
            reaction_layout_selector: momapy.styling.StyleCollection(
                {
                    "end_arrowhead_stroke": momapy.coloring.lightgray,
                    "end_arrowhead_fill": momapy.coloring.lightgray,
                    "start_arrowhead_stroke": momapy.coloring.lightgray,
                    "start_arrowhead_fill": momapy.coloring.lightgray,
                    "_font_fill": momapy.coloring.lightgray,
                }
            ),
        }
    )
    layout = momapy.styling.apply_style_sheet(layout, style_sheet, strict=False)
    return layout


def _translate_layout_element(layout_element, tx, ty):
    layout_element.position = momapy.geometry.Point(
        layout_element.position.x + tx, layout_element.position.y + ty
    )
    for sub_layout_element in layout_element.children():
        _translate_layout_element(sub_layout_element, tx, ty)


def _get_coordinates_from_pydot_node(dot_node):
    x, y = [
        float(coordinate) for coordinate in dot_node.get("pos").strip('"').split(",")
    ]
    width = float(dot_node.get("width").strip('"'))
    height = float(dot_node.get("height").strip('"'))
    x = x + width / 2
    y = y + height / 2
    return x, y


def _get_flatten_dot_nodes(dot_graph):
    dot_nodes = dot_graph.get_nodes()
    for dot_subgraph in dot_graph.get_subgraphs():
        dot_nodes += _get_flatten_dot_nodes(dot_subgraph)
    return dot_nodes


def auto_layout(cd_map):
    new_map_builder = momapy.builder.builder_from_object(cd_map)
    new_layout_builder = new_map_builder.layout
    dot_graph = pydot.Dot(graph_type="digraph")
    compartment_to_dot_cluster = {}
    compartment_layout_element_to_included_layout_elements = collections.defaultdict(
        list
    )
    for compartment in new_map_builder.model.compartments:
        compartment_layout_elements = new_map_builder.layout_model_mapping.get_mapping(
            compartment
        )
        if compartment_layout_elements is not None:
            compartment_layout_element = compartment_layout_elements[0]
            dot_cluster = pydot.Cluster(compartment.id_)
            compartment_to_dot_cluster[compartment] = dot_cluster
            outside_compartment = compartment.outside
            outside_compartment_layout_elements = (
                new_map_builder.layout_model_mapping.get_mapping(outside_compartment)
            )
            if outside_compartment_layout_elements is not None:
                outside_compartment_layout_element = (
                    outside_compartment_layout_elements[0]
                )
                compartment_layout_element_to_included_layout_elements[
                    outside_compartment_layout_element
                ].append(compartment_layout_element)
    for compartment, compartment_dot_cluster in compartment_to_dot_cluster.items():
        outside_compartment = compartment.outside
        if outside_compartment is not None:
            outside_dot_cluster = compartment_to_dot_cluster.get(outside_compartment)
            if outside_dot_cluster is not None:
                outside_dot_cluster.add_subgraph(compartment_dot_cluster)
            else:
                dot_graph.add_subgraph(compartment_dot_cluster)
    id_to_layout_element = {}
    for layout_element_builder in new_layout_builder.layout_elements:
        if momapy.builder.isinstance_or_builder(
            layout_element_builder, momapy.core.Node
        ) and not momapy.builder.isinstance_or_builder(
            layout_element_builder,
            (
                momapy.celldesigner.core.RectangleCompartmentLayout,
                momapy.celldesigner.core.OvalCompartmentLayout,
            ),
        ):
            dot_node = pydot.Node(layout_element_builder.id_)
            dot_node.set("width", layout_element_builder.width / _POINTS_PER_INCH)
            dot_node.set("height", layout_element_builder.width / _POINTS_PER_INCH)
            model_element = new_map_builder.layout_model_mapping.get_mapping(
                layout_element_builder
            )
            compartment = model_element.compartment
            if compartment is not None:
                compartment_dot_cluster = compartment_to_dot_cluster.get(
                    model_element.compartment
                )
                if compartment_dot_cluster is not None:
                    compartment_dot_cluster.add_node(dot_node)
                    compartment_layout_element = (
                        new_map_builder.layout_model_mapping.get_mapping(compartment)[0]
                    )
                    compartment_layout_element_to_included_layout_elements[
                        compartment_layout_element
                    ].append(layout_element_builder)
                else:
                    dot_graph.add_node(dot_node)
            else:
                dot_graph.add_node(dot_node)
        elif momapy.builder.isinstance_or_builder(
            layout_element_builder, momapy.core.Arc
        ):
            dot_graph.add_edge(
                pydot.Edge(
                    layout_element_builder.source.id_, layout_element_builder.target.id_
                )
            )
        id_to_layout_element[layout_element_builder.id_] = layout_element_builder
    dot_graph.set("ranksep", 3.0)
    dot_graph.set("rankdir", "BT")
    dot = dot_graph.create_dot(prog="dot").decode("utf-8")
    dot_graph = pydot.graph_from_dot_data(dot)[0]
    id_to_new_layout_element_builder = {}
    dot_nodes = _get_flatten_dot_nodes(dot_graph)
    for dot_node in dot_nodes:
        dot_node_id = dot_node.get_name().strip('"')
        if dot_node_id not in ["graph", "node"]:
            x, y = _get_coordinates_from_pydot_node(dot_node)
            layout_element_builder = id_to_layout_element[dot_node_id]
            _translate_layout_element(
                layout_element_builder,
                x - layout_element_builder.x,
                y - layout_element_builder.y,
            )
            id_to_new_layout_element_builder[layout_element_builder.id_] = (
                layout_element_builder
            )
    for layout_element_builder in new_layout_builder.layout_elements:
        if momapy.builder.isinstance_or_builder(
            layout_element_builder, momapy.core.Arc
        ):
            source_layout_element_builder = id_to_new_layout_element_builder[
                layout_element_builder.source.id_
            ]
            target_layout_element_builder = id_to_new_layout_element_builder[
                layout_element_builder.target.id_
            ]
            start_point = source_layout_element_builder.border(
                target_layout_element_builder.center()
            )
            end_point = target_layout_element_builder.border(
                source_layout_element_builder.center()
            )
            layout_element_builder.segments = momapy.core.TupleBuilder(
                [momapy.geometry.Segment(start_point, end_point)]
            )
    for (
        compartment_layout_element,
        included_layout_elements,
    ) in compartment_layout_element_to_included_layout_elements.items():
        momapy.positioning.set_fit(
            compartment_layout_element, included_layout_elements, xsep=10.0, ysep=10.0
        )
        compartment_layout_element.label.position = compartment_layout_element.position
    momapy.positioning.set_fit(new_layout_builder, new_layout_builder.layout_elements)
    new_map = momapy.builder.object_from_builder(new_map_builder)
    return new_map
