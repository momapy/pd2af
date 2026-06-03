import collections
import dataclasses
import math

import momapy.core.layout
import momapy.styling
import momapy.coloring
import momapy.drawing
import momapy.builder
import momapy.geometry
import momapy.positioning
import momapy.celldesigner
import pydot


@dataclasses.dataclass(frozen=True)
class _NotInIdSetSelector(momapy.styling.Selector):
    """Selects elements whose `id_` is not in `keep_ids` (or which have none).

    Replaces a NotSelector wrapping one IdSelector per kept element: that
    construction made selection O(N) per visited element and quadratic over
    the whole layout. A frozenset lookup keeps it O(1)."""

    keep_ids: frozenset

    def select(self, obj, ancestors):
        obj_id = getattr(obj, "id_", None)
        return obj_id is None or obj_id not in self.keep_ids


@dataclasses.dataclass(frozen=True)
class _ClassNameSuffixSelector(momapy.styling.Selector):
    """Selects elements whose class name (ignoring a trailing 'Builder')
    ends with `suffix`.

    Used to target CellDesigner active-border layouts (`*ActiveLayout`). The
    'Builder' strip mirrors `TypeSelector`: `apply_style_sheet` visits the
    layout as builders (e.g. `GenericProteinActiveLayoutBuilder`)."""

    suffix: str

    def select(self, obj, ancestors):
        class_name = type(obj).__name__
        if class_name.endswith("Builder"):
            class_name = class_name[: -len("Builder")]
        return class_name.endswith(self.suffix)


_POINTS_PER_INCH = 96
_BEZIER_OFFSET = 30.0

_ROOT_LAYOUT_SEP = 15.0


def harmonize_root_layout(layout_builder):
    """Set the root layout's fill to white and fit it tightly around its
    children. Applied uniformly across all layout modes so the rendered
    canvas has consistent background and padding."""
    layout_builder.fill = momapy.coloring.white
    momapy.positioning.set_fit(
        layout_builder,
        layout_builder.layout_elements,
        xsep=_ROOT_LAYOUT_SEP,
        ysep=_ROOT_LAYOUT_SEP,
    )


def highlight_layout_elements(layout_elements, layout):
    keep_ids = set()
    for layout_element in layout_elements:
        keep_ids.add(layout_element.id_)
        for descendant in layout_element.descendants():
            keep_ids.add(descendant.id_)
    not_selector = _NotInIdSetSelector(frozenset(keep_ids))
    layout_element_selector = momapy.styling.CompoundSelector(
        tuple([momapy.styling.ClassSelector("LayoutElement"), not_selector])
    )
    active_border_selector = momapy.styling.CompoundSelector(
        tuple([_ClassNameSuffixSelector("ActiveLayout"), not_selector])
    )
    text_layout_selector = momapy.styling.CompoundSelector(
        tuple([momapy.styling.TypeSelector("TextLayout"), not_selector])
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
                    "fill": momapy.coloring.white,
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
            # Active-border sublayouts (`*ActiveLayout`) are transparent by
            # default and drawn larger, on top of the species body. The blanket
            # rule above sets `fill: white`, which makes the border opaque and
            # paints over the body; restore transparency so the dimmed body
            # shows through. Must follow `layout_element_selector`: matching
            # selectors apply in dict order with no specificity, last wins.
            active_border_selector: momapy.styling.StyleCollection(
                {
                    "fill": momapy.drawing.NoneValue,
                }
            ),
            text_layout_selector: momapy.styling.StyleCollection(
                {
                    "stroke": momapy.drawing.NoneValue,
                    "fill": momapy.coloring.lightgray,
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
    directed_pairs = set()
    # Map every descendant id to its top-level Node ancestor (the one that
    # gets added as a pydot node). Arc endpoints can reference descendants
    # (e.g. a subunit inside a complex); those must be redirected to the
    # top-level node so dot doesn't auto-create phantom nodes.
    descendant_id_to_top_level_id = {}
    for layout_element_builder in new_layout_builder.layout_elements:
        if momapy.builder.isinstance_or_builder(
            layout_element_builder, momapy.core.layout.Node
        ) and not momapy.builder.isinstance_or_builder(
            layout_element_builder,
            (
                momapy.celldesigner.RectangleCompartmentLayout,
                momapy.celldesigner.OvalCompartmentLayout,
            ),
        ):
            descendant_id_to_top_level_id[layout_element_builder.id_] = (
                layout_element_builder.id_
            )
            for descendant in layout_element_builder.descendants():
                descendant_id = getattr(descendant, "id_", None)
                if descendant_id is not None:
                    descendant_id_to_top_level_id.setdefault(
                        descendant_id, layout_element_builder.id_
                    )
    for layout_element_builder in new_layout_builder.layout_elements:
        if momapy.builder.isinstance_or_builder(
            layout_element_builder, momapy.core.layout.Node
        ) and not momapy.builder.isinstance_or_builder(
            layout_element_builder,
            (
                momapy.celldesigner.RectangleCompartmentLayout,
                momapy.celldesigner.OvalCompartmentLayout,
            ),
        ):
            dot_node = pydot.Node(layout_element_builder.id_)
            dot_node.set("width", layout_element_builder.width / _POINTS_PER_INCH)
            dot_node.set("height", layout_element_builder.height / _POINTS_PER_INCH)
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
            layout_element_builder, momapy.core.layout.Arc
        ):
            source_id = descendant_id_to_top_level_id.get(
                layout_element_builder.source.id_,
                layout_element_builder.source.id_,
            )
            target_id = descendant_id_to_top_level_id.get(
                layout_element_builder.target.id_,
                layout_element_builder.target.id_,
            )
            dot_graph.add_edge(pydot.Edge(source_id, target_id))
            directed_pairs.add((source_id, target_id))
        id_to_layout_element[layout_element_builder.id_] = layout_element_builder
    dot_graph.set("ranksep", 1.0)
    dot_graph.set("nodesep", 0.5)
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
            # Descendants (e.g. subunit layouts inside a complex layout) are
            # auto-translated by `_translate_layout_element` along with their
            # parent. Register them so arcs whose source/target is a nested
            # subunit can still resolve.
            for descendant in layout_element_builder.descendants():
                descendant_id = getattr(descendant, "id_", None)
                if descendant_id is not None:
                    id_to_new_layout_element_builder.setdefault(
                        descendant_id, descendant
                    )
    for layout_element_builder in new_layout_builder.layout_elements:
        if momapy.builder.isinstance_or_builder(
            layout_element_builder, momapy.core.layout.Arc
        ):
            source_layout_element_builder = id_to_new_layout_element_builder[
                layout_element_builder.source.id_
            ]
            target_layout_element_builder = id_to_new_layout_element_builder[
                layout_element_builder.target.id_
            ]
            source_top_level_id = descendant_id_to_top_level_id.get(
                layout_element_builder.source.id_,
                layout_element_builder.source.id_,
            )
            target_top_level_id = descendant_id_to_top_level_id.get(
                layout_element_builder.target.id_,
                layout_element_builder.target.id_,
            )
            is_self_loop = source_top_level_id == target_top_level_id
            is_bidirectional = not is_self_loop and (
                target_top_level_id,
                source_top_level_id,
            ) in directed_pairs
            if is_self_loop:
                start_point = source_layout_element_builder.own_angle(120)
                end_point = source_layout_element_builder.own_angle(60)
                if start_point is None:
                    start_point = source_layout_element_builder.north_west()
                if end_point is None:
                    end_point = source_layout_element_builder.north_east()
                center = source_layout_element_builder.center()
                start_delta_x = start_point.x - center.x
                start_delta_y = start_point.y - center.y
                end_delta_x = end_point.x - center.x
                end_delta_y = end_point.y - center.y
                start_length = math.hypot(start_delta_x, start_delta_y) or 1.0
                end_length = math.hypot(end_delta_x, end_delta_y) or 1.0
                start_control_point = momapy.geometry.Point(
                    start_point.x + start_delta_x / start_length * _BEZIER_OFFSET,
                    start_point.y + start_delta_y / start_length * _BEZIER_OFFSET,
                )
                end_control_point = momapy.geometry.Point(
                    end_point.x + end_delta_x / end_length * _BEZIER_OFFSET,
                    end_point.y + end_delta_y / end_length * _BEZIER_OFFSET,
                )
                # Express the loop as a polyline through the two control
                # points so the CellDesigner writer can recover them as edit
                # points (it only sees segment endpoints). Without this, a
                # single Bezier collapses to start/end on the same node and
                # the reader's modulation-geometry call hits a None border.
                layout_element_builder.segments = [
                    momapy.geometry.Segment(start_point, start_control_point),
                    momapy.geometry.Segment(start_control_point, end_control_point),
                    momapy.geometry.Segment(end_control_point, end_point),
                ]
            elif is_bidirectional:
                source_center = source_layout_element_builder.center()
                target_center = target_layout_element_builder.center()
                delta_x = target_center.x - source_center.x
                delta_y = target_center.y - source_center.y
                length = math.hypot(delta_x, delta_y)
                if length == 0:
                    start_point = source_layout_element_builder.own_border(
                        target_center
                    )
                    end_point = target_layout_element_builder.own_border(
                        source_center
                    )
                    if start_point is None:
                        start_point = source_layout_element_builder.north_west()
                    if end_point is None:
                        end_point = target_layout_element_builder.north_east()
                    layout_element_builder.segments = [
                        momapy.geometry.Segment(start_point, end_point)
                    ]
                else:
                    normal_x = -delta_y / length
                    normal_y = delta_x / length
                    # Deterministic side rule: A→B and B→A get opposite offsets,
                    # so the two curves bow away from each other.
                    if (source_top_level_id, target_top_level_id) > (
                        target_top_level_id,
                        source_top_level_id,
                    ):
                        normal_x = -normal_x
                        normal_y = -normal_y
                    middle_x = (source_center.x + target_center.x) / 2
                    middle_y = (source_center.y + target_center.y) / 2
                    control_point = momapy.geometry.Point(
                        middle_x + _BEZIER_OFFSET * normal_x,
                        middle_y + _BEZIER_OFFSET * normal_y,
                    )
                    start_point = source_layout_element_builder.own_border(
                        control_point
                    )
                    end_point = target_layout_element_builder.own_border(
                        control_point
                    )
                    if start_point is None:
                        start_point = source_layout_element_builder.north_west()
                    if end_point is None:
                        end_point = target_layout_element_builder.north_east()
                    # Polyline through the control point so the writer
                    # serializes it as an edit point (see self-loop note).
                    layout_element_builder.segments = [
                        momapy.geometry.Segment(start_point, control_point),
                        momapy.geometry.Segment(control_point, end_point),
                    ]
            else:
                start_point = source_layout_element_builder.border(
                    target_layout_element_builder.center()
                )
                end_point = target_layout_element_builder.border(
                    source_layout_element_builder.center()
                )
                if start_point is None:
                    start_point = source_layout_element_builder.north_west()
                if end_point is None:
                    end_point = target_layout_element_builder.north_east()
                layout_element_builder.segments = [
                    momapy.geometry.Segment(start_point, end_point)
                ]
    for (
        compartment_layout_element,
        included_layout_elements,
    ) in compartment_layout_element_to_included_layout_elements.items():
        momapy.positioning.set_fit(
            compartment_layout_element, included_layout_elements, xsep=10.0, ysep=10.0
        )
        compartment_layout_element.label.position = compartment_layout_element.position
    harmonize_root_layout(new_layout_builder)
    new_map = momapy.builder.object_from_builder(new_map_builder)
    return new_map
