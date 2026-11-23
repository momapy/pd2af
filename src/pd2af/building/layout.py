"""Layout helpers shared by both language builders: styling, arcs and dot."""

import dataclasses
import itertools
import math
import typing

import momapy.builder
import momapy.celldesigner
import momapy.coloring
import momapy.core.layout
import momapy.drawing
import momapy.geometry
import momapy.positioning
import momapy.styling
import pydot


@dataclasses.dataclass(frozen=True)
class _ClassNameSuffixSelector(momapy.styling.Selector):
    """Selects elements whose class name ends with `suffix`.

    A trailing 'Builder' is ignored. Used to target CellDesigner active-border layouts (`*ActiveLayout`). The
    'Builder' strip mirrors `TypeSelector`: `apply_style_sheet` visits the
    layout as builders (e.g. `GenericProteinActiveLayoutBuilder`).
    """

    suffix: str

    def select(self, obj: typing.Any, ancestors: typing.Any) -> bool:
        class_name = type(obj).__name__
        class_name = class_name.removesuffix("Builder")
        return class_name.endswith(self.suffix)


_POINTS_PER_INCH = 96
_BEZIER_OFFSET = 30.0

# Where the first self-loop of a node sits (degrees, 90 being straight up) and
# how much of the border one loop spans.
_SELF_LOOP_CENTER_ANGLE = 90.0
_SELF_LOOP_SPAN = 60.0

_ROOT_LAYOUT_SEP = 15.0


# Padding (in points) dot leaves between a compartment's contents and its
# cluster bounding box. graphviz defaults to 8, which is too tight; this is the
# auto-layout equivalent of the old set_fit xsep/ysep.
_DOT_CLUSTER_SEP = 40.0


def get_layout_elements_sorted_by_id(
    layout_model_mapping: typing.Any, model_element: typing.Any
) -> list:
    """Return the glyphs that draw ``model_element``, sorted by ``id_``."""
    return sorted(
        layout_model_mapping.get_representative_layout_elements(model_element),
        key=lambda layout_element: layout_element.id_,
    )


def get_nearest_layout_element(
    layout_elements: typing.Any, position: typing.Any
) -> typing.Any:
    """Return the layout element whose position is nearest ``position``."""
    return min(
        layout_elements,
        key=lambda layout_element: math.dist(
            layout_element.position.to_tuple(), position.to_tuple()
        ),
    )


def nearest_layout_pair(
    source_layouts: typing.Any, target_layouts: typing.Any
) -> tuple:
    """Return the (source_layout, target_layout) pair closest in Euclidean distance.

    ``min`` over ``itertools.product`` is deterministic: on a distance tie it
    keeps the first pair in product order (source-major), so the choice mirrors
    the order the cross product would have emitted.
    """
    return min(
        itertools.product(source_layouts, target_layouts),
        key=lambda pair: math.dist(
            pair[0].position.to_tuple(), pair[1].position.to_tuple()
        ),
    )


def influence_layout_pairs(
    source_layouts: typing.Any, target_layouts: typing.Any, prefer_nearest: bool
) -> list[tuple]:
    """Return the (source_layout, target_layout) pairs an arc is drawn between.

    By default the full ``source x target`` cross product. ``prefer_nearest``
    (``-p nearest``) collapses the fan-out to the single closest pair, but only
    where there *is* fan-out (more than one layout on a side) and every layout
    has a real ``position`` -- in ``dot`` mode positions are throwaway
    placeholders that graphviz overwrites, so the cross product is kept.
    """
    use_nearest = (
        prefer_nearest
        and (len(source_layouts) > 1 or len(target_layouts) > 1)
        and all(
            layout.position is not None
            for layout in itertools.chain(source_layouts, target_layouts)
        )
    )
    if use_nearest:
        return [nearest_layout_pair(source_layouts, target_layouts)]
    return list(itertools.product(source_layouts, target_layouts))


# In ``dot`` mode every node is built at a placeholder position and graphviz
# repositions it afterwards, so any arc geometry built before that is replaced
# by :func:`_rebuild_arc_geometry`. Arcs are built with this instead.
PLACEHOLDER_ARC_SEGMENTS = (
    momapy.geometry.Segment(
        momapy.geometry.Point(0.0, 0.0), momapy.geometry.Point(0.0, 0.0)
    ),
)


def make_arc_segments_from_source_and_target(
    source_layout: typing.Any, target_layout: typing.Any
) -> list:
    """Segments for a modulation / influence arc connecting two node layouts.

    A self-loop (source is target) bows out into a visible loop via
    :func:`_make_self_loop_segments`; otherwise it is drawn straight by
    :func:`_make_straight_segments`. Returns a list of segments.
    """
    if source_layout is target_layout:
        start_angle, end_angle = _make_self_loop_angles(0, 1)
        return _make_self_loop_segments(source_layout, start_angle, end_angle)
    return _make_straight_segments(source_layout, target_layout)


def harmonize_root_layout(layout: typing.Any) -> typing.Any:
    """Return ``layout`` with a white fill, fitted around its children.

    Applied uniformly across all layout modes so the rendered canvas has
    consistent background and padding.
    """
    if not layout.layout_elements:
        # An empty output map (e.g. every activity suppressed via
        # `--set-all-inactive`) has nothing to fit around; `fit` raises
        # on an empty element list. Give the root a minimal, valid canvas so
        # the map still has numeric width/height/position for the writer.
        side = 2 * _ROOT_LAYOUT_SEP
        return dataclasses.replace(
            layout,
            fill=momapy.coloring.white,
            width=side,
            height=side,
            position=momapy.geometry.Point(_ROOT_LAYOUT_SEP, _ROOT_LAYOUT_SEP),
        )
    bounding_box = momapy.positioning.fit(
        layout.layout_elements, xsep=_ROOT_LAYOUT_SEP, ysep=_ROOT_LAYOUT_SEP
    )
    return dataclasses.replace(
        layout,
        fill=momapy.coloring.white,
        width=bounding_box.width,
        height=bounding_box.height,
        position=bounding_box.position,
    )


def dim_layout_elements(layout_element_builders: typing.Any) -> None:
    """Wash out ``layout_element_builders`` in place: the overlay background."""
    layout_element_selector = momapy.styling.ClassSelector("LayoutElement")
    active_border_selector = _ClassNameSuffixSelector("ActiveLayout")
    text_layout_selector = momapy.styling.TypeSelector("TextLayout")
    production_layout_selector = momapy.styling.TypeSelector("ProductionLayout")
    compartment_layout_selector = momapy.styling.TypeSelector(
        "RectangleCompartmentLayout"
    )
    reaction_layout_selector = momapy.styling.OrSelector(
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
    for layout_element_builder in layout_element_builders:
        momapy.styling.apply_style_sheet(
            layout_element_builder, style_sheet, strict=False
        )


def _translate_layout_element(
    layout_element: typing.Any, translation_x: float, translation_y: float
):
    layout_element.position = momapy.geometry.Point(
        layout_element.position.x + translation_x,
        layout_element.position.y + translation_y,
    )
    for sub_layout_element in layout_element.children():
        _translate_layout_element(sub_layout_element, translation_x, translation_y)


def _get_coordinates_from_pydot_node(dot_node: typing.Any) -> tuple[float, float]:
    x, y = [
        float(coordinate) for coordinate in dot_node.get("pos").strip('"').split(",")
    ]
    width = float(dot_node.get("width").strip('"'))
    height = float(dot_node.get("height").strip('"'))
    x = x + width / 2
    y = y + height / 2
    return x, y


def _collect_dot_nodes_recursively(dot_graph: typing.Any) -> list:
    dot_nodes = dot_graph.get_nodes()
    for dot_subgraph in dot_graph.get_subgraphs():
        dot_nodes += _collect_dot_nodes_recursively(dot_subgraph)
    return dot_nodes


def _build_dot_graph(
    new_map_builder: typing.Any,
    new_layout_builder: typing.Any,
    compartment_layout_classes: typing.Any,
    reversed_arc_classes: typing.Any = (),
) -> typing.Any:
    """Build the pydot graph from the built layout.

    Compartments become dot clusters, Node layout elements become dot nodes
    (placed in their compartment's cluster when there is one), and Arc layout
    elements become dot edges. Returns the graph plus the bookkeeping the
    repositioning and arc-geometry phases need.

    ``reversed_arc_classes``: arc layout classes whose dot edge is added with
    source and target **swapped**, so graphviz ranks the arc's target upstream
    of its source. The layout arc object is unchanged; only the ranking flips.
    Used for SBGN logic arcs (stored operator -> input) so an operator's inputs
    rank above it and its output target below it.
    """
    dot_graph = pydot.Dot(graph_type="digraph")
    compartment_to_dot_cluster = {}
    # dot cluster name -> compartment layout element, so the repositioning phase
    # can copy each cluster's dot-computed bounding box onto its compartment.
    dot_cluster_name_to_compartment_layout_element = {}
    for compartment in new_map_builder.model.compartments:
        compartment_layout_elements = get_layout_elements_sorted_by_id(
            new_map_builder.layout_model_mapping, compartment
        )
        if compartment_layout_elements:
            compartment_layout_element = compartment_layout_elements[0]
            dot_cluster = pydot.Cluster(compartment.id_)
            dot_cluster.set("margin", _DOT_CLUSTER_SEP)
            compartment_to_dot_cluster[compartment] = dot_cluster
            # pydot quotes cluster names whose ids need it (uuid-like ids);
            # strip so the name matches the post-dot subgraph's stripped one.
            dot_cluster_name_to_compartment_layout_element[
                dot_cluster.get_name().strip('"')
            ] = compartment_layout_element
    for compartment, compartment_dot_cluster in compartment_to_dot_cluster.items():
        outside_compartment = getattr(compartment, "outside", None)
        outside_dot_cluster = (
            compartment_to_dot_cluster.get(outside_compartment)
            if outside_compartment is not None
            else None
        )
        if outside_dot_cluster is not None:
            outside_dot_cluster.add_subgraph(compartment_dot_cluster)
        else:
            dot_graph.add_subgraph(compartment_dot_cluster)

    def is_node_layout_element(layout_element_builder: typing.Any) -> bool:
        return momapy.builder.isinstance_or_builder(
            layout_element_builder, momapy.core.layout.Node
        ) and not momapy.builder.isinstance_or_builder(
            layout_element_builder, compartment_layout_classes
        )

    id_to_layout_element = {}
    # Map every descendant id to its top-level Node ancestor (the one that
    # gets added as a pydot node). Arc endpoints can reference descendants
    # (e.g. a subunit inside a complex); those must be redirected to the
    # top-level node so dot doesn't auto-create phantom nodes.
    descendant_id_to_top_level_id = {}
    for layout_element_builder in new_layout_builder.layout_elements:
        if is_node_layout_element(layout_element_builder):
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
        if is_node_layout_element(layout_element_builder):
            dot_node = pydot.Node(layout_element_builder.id_)
            dot_node.set("width", layout_element_builder.width / _POINTS_PER_INCH)
            dot_node.set("height", layout_element_builder.height / _POINTS_PER_INCH)
            model_element = new_map_builder.layout_model_mapping.get_mapping(
                layout_element_builder
            )
            # A logical operator (CellDesigner gate / SBGN-AF operator) has no
            # `compartment` field -- it is never placed inside a compartment
            # cluster -- so resolve defensively.
            compartment = getattr(model_element, "compartment", None)
            if compartment is not None:
                compartment_dot_cluster = compartment_to_dot_cluster.get(compartment)
                if compartment_dot_cluster is not None:
                    compartment_dot_cluster.add_node(dot_node)
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
            if momapy.builder.isinstance_or_builder(
                layout_element_builder, reversed_arc_classes
            ):
                source_id, target_id = target_id, source_id
            dot_graph.add_edge(pydot.Edge(source_id, target_id))
        id_to_layout_element[layout_element_builder.id_] = layout_element_builder
    dot_graph.set("ranksep", 1.0)
    dot_graph.set("nodesep", 0.5)
    dot_graph.set("rankdir", "BT")
    return (
        dot_graph,
        id_to_layout_element,
        descendant_id_to_top_level_id,
        dot_cluster_name_to_compartment_layout_element,
    )


def _apply_dot_cluster_bounding_boxes_to_compartments(
    dot_graph: typing.Any, dot_cluster_name_to_compartment_layout_element: typing.Any
):
    """Copy each dot cluster's bounding box onto its compartment layout element.

    dot lays clusters out already nested and non-overlapping; we use its `bb`
    directly rather than refitting compartments around their members.

    pydot exposes a cluster's `bb` only through a synthetic node named ``graph``
    inside the subgraph (``subgraph.get("bb")`` returns ``None``). `bb` is
    ``llx,lly,urx,ury`` in points -- the same coordinate space as the node
    positions consumed elsewhere. An empty cluster has no `bb`; its compartment
    keeps its built geometry.
    """
    for dot_subgraph in dot_graph.get_subgraphs():
        compartment_layout_element = dot_cluster_name_to_compartment_layout_element.get(
            dot_subgraph.get_name().strip('"')
        )
        if compartment_layout_element is not None:
            bounding_box = None
            for dot_node in dot_subgraph.get_nodes():
                if dot_node.get_name().strip('"') == "graph":
                    bounding_box = dot_node.get("bb")
            if bounding_box is not None:
                lower_left_x, lower_left_y, upper_right_x, upper_right_y = [
                    float(coordinate)
                    for coordinate in bounding_box.strip('"').split(",")
                ]
                compartment_layout_element.position = momapy.geometry.Point(
                    (lower_left_x + upper_right_x) / 2,
                    (lower_left_y + upper_right_y) / 2,
                )
                compartment_layout_element.width = upper_right_x - lower_left_x
                compartment_layout_element.height = upper_right_y - lower_left_y
                compartment_layout_element.label.position = (
                    compartment_layout_element.position
                )
        _apply_dot_cluster_bounding_boxes_to_compartments(
            dot_subgraph, dot_cluster_name_to_compartment_layout_element
        )


def _reposition_from_dot(
    dot_graph: typing.Any,
    id_to_layout_element: typing.Any,
    dot_cluster_name_to_compartment_layout_element: typing.Any,
) -> typing.Any:
    """Run graphviz `dot` and apply the geometry it computes to the layout.

    Every layout element is translated to the position dot computed for it, and
    each compartment sized to its dot cluster bounding box. Returns id -> layout
    element builder (including nested descendants) so the arc-geometry phase can
    resolve endpoints.
    """
    dot = dot_graph.create_dot(prog="dot").decode("utf-8")
    dot_graph = pydot.graph_from_dot_data(dot)[0]
    _apply_dot_cluster_bounding_boxes_to_compartments(
        dot_graph, dot_cluster_name_to_compartment_layout_element
    )
    id_to_new_layout_element_builder = {}
    dot_nodes = _collect_dot_nodes_recursively(dot_graph)
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
    return id_to_new_layout_element_builder


def _make_self_loop_angles(index: int, count: int) -> tuple[float, float]:
    """The (start, end) border angles of one of several parallel self-loops.

    The ``count`` loops on a node are spread evenly around it, each spanning
    :data:`_SELF_LOOP_SPAN` degrees, the first one on top; this returns the
    ``index``-th.
    """
    center_angle = (_SELF_LOOP_CENTER_ANGLE + index * 360 / count) % 360
    return (
        (center_angle + _SELF_LOOP_SPAN / 2) % 360,
        (center_angle - _SELF_LOOP_SPAN / 2) % 360,
    )


def _make_self_loop_segments(
    layout_element: typing.Any, start_angle: float, end_angle: float
) -> list:
    """Segments for an arc whose source and target resolve to the same node.

    The loop leaves the border at ``start_angle`` and returns at ``end_angle``,
    bowing out through two control points, expressed as a polyline so the
    CellDesigner writer can recover them as edit points (it only sees segment
    endpoints). Without this, a single Bezier collapses to start/end on the same
    node and the reader's modulation-geometry call hits a None border.
    """
    start_point = layout_element.own_angle(start_angle)
    end_point = layout_element.own_angle(end_angle)
    center = layout_element.center()
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
    return [
        momapy.geometry.Segment(start_point, start_control_point),
        momapy.geometry.Segment(start_control_point, end_control_point),
        momapy.geometry.Segment(end_control_point, end_point),
    ]


def _make_offset_ladder(count: int) -> list[float]:
    """The perpendicular offsets spreading ``count`` arcs between the same two nodes.

    They form a ladder symmetric about the straight line joining the nodes, in
    steps of :data:`_BEZIER_OFFSET`. An odd count puts one arc on the line
    (offset 0); an even count straddles it, so no arc is drawn straight.
    """
    if count % 2 == 1:
        return [(index - count // 2) * _BEZIER_OFFSET for index in range(count)]
    half_count = count // 2
    return [
        (index - half_count) * _BEZIER_OFFSET
        if index < half_count
        else (index - half_count + 1) * _BEZIER_OFFSET
        for index in range(count)
    ]


def _make_offset_segments(
    source_layout_element: typing.Any,
    target_layout_element: typing.Any,
    frame_start_center: typing.Any,
    frame_end_center: typing.Any,
    offset: float,
) -> typing.Any:
    """Segments for one arc of a group connecting the same two nodes.

    The curve bows ``offset`` points off the straight line through a single
    control point, serialized as a polyline so the writer keeps that control
    point (see the self-loop note).

    The normal is taken in the group's canonical frame -- from
    ``frame_start_center`` to ``frame_end_center``, the centers of the two nodes
    in sorted-id order -- and not from the arc's own direction, so the sign of
    ``offset`` names the same side of the line for every arc of the group
    whichever way it points. A zero offset is a straight segment.
    """
    delta_x = frame_end_center.x - frame_start_center.x
    delta_y = frame_end_center.y - frame_start_center.y
    length = math.hypot(delta_x, delta_y)
    if offset == 0 or length == 0:
        return _make_straight_segments(source_layout_element, target_layout_element)
    normal_x = -delta_y / length
    normal_y = delta_x / length
    middle_x = (frame_start_center.x + frame_end_center.x) / 2
    middle_y = (frame_start_center.y + frame_end_center.y) / 2
    control_point = momapy.geometry.Point(
        middle_x + offset * normal_x,
        middle_y + offset * normal_y,
    )
    start_point = source_layout_element.own_border(control_point)
    end_point = target_layout_element.own_border(control_point)
    return [
        momapy.geometry.Segment(start_point, control_point),
        momapy.geometry.Segment(control_point, end_point),
    ]


def _make_straight_segments(
    source_layout_element: typing.Any, target_layout_element: typing.Any
) -> list:
    """Segments for an arc drawn straight: one segment between the node borders.

    Used for the only arc between its two nodes, or the middle one of an odd
    group.
    """
    start_point = source_layout_element.own_border(target_layout_element.center())
    end_point = target_layout_element.own_border(source_layout_element.center())
    return [momapy.geometry.Segment(start_point, end_point)]


@dataclasses.dataclass
class _ArcPlacement:
    """One arc of a node-pair group, with everything the geometry needs.

    That is the arc itself, the layout elements its endpoints resolve to (a
    subunit when the arc attaches to one) and the top-level nodes those sit in.
    """

    arc: object
    source_layout_element: object
    target_layout_element: object
    source_top_level_id: str
    target_top_level_id: str


def _group_arcs_by_node_pair(
    new_layout_builder: typing.Any,
    id_to_new_layout_element_builder: typing.Any,
    descendant_id_to_top_level_id: typing.Any,
    operator_arc_resolver: typing.Any = None,
) -> typing.Any:
    """Group the layout's arcs by the *unordered* pair of nodes they connect.

    Both directions between two nodes thus land in one group. Returns
    ``{(smaller_id, greater_id): [placement, ...]}``, each group ordered
    deterministically. Arcs the ``operator_arc_resolver`` claims get their
    segments set here and are left out of the groups.
    """
    arcs_by_node_pair = {}
    for layout_element_builder in new_layout_builder.layout_elements:
        if not momapy.builder.isinstance_or_builder(
            layout_element_builder, momapy.core.layout.Arc
        ):
            continue
        source_layout_element_builder = id_to_new_layout_element_builder[
            layout_element_builder.source.id_
        ]
        target_layout_element_builder = id_to_new_layout_element_builder[
            layout_element_builder.target.id_
        ]
        if operator_arc_resolver is not None:
            resolved_segments = operator_arc_resolver(
                layout_element_builder,
                source_layout_element_builder,
                target_layout_element_builder,
            )
            if resolved_segments is not None:
                layout_element_builder.segments = resolved_segments
                continue
        source_top_level_id = descendant_id_to_top_level_id.get(
            layout_element_builder.source.id_,
            layout_element_builder.source.id_,
        )
        target_top_level_id = descendant_id_to_top_level_id.get(
            layout_element_builder.target.id_,
            layout_element_builder.target.id_,
        )
        node_pair = tuple(sorted((source_top_level_id, target_top_level_id)))
        arcs_by_node_pair.setdefault(node_pair, []).append(
            _ArcPlacement(
                arc=layout_element_builder,
                source_layout_element=source_layout_element_builder,
                target_layout_element=target_layout_element_builder,
                source_top_level_id=source_top_level_id,
                target_top_level_id=target_top_level_id,
            )
        )
    for placements in arcs_by_node_pair.values():
        placements.sort(
            key=lambda placement: (
                placement.source_top_level_id,
                placement.target_top_level_id,
                placement.arc.id_,
            )
        )
    return arcs_by_node_pair


def _rebuild_arc_geometry(
    new_layout_builder: typing.Any,
    id_to_new_layout_element_builder: typing.Any,
    descendant_id_to_top_level_id: typing.Any,
    operator_arc_resolver: typing.Any = None,
):
    """Rebuild each arc's segments from the repositioned node geometry.

    Arcs are handled a node pair at a time rather than one by one: every arc
    between the same two top-level nodes -- in either direction -- takes its own
    slot in the symmetric offset ladder of :func:`_make_offset_ladder`, measured
    in the pair's canonical frame, so no two arcs of a pair share a curve. A
    group whose two nodes are the same node is a bundle of self-loops instead,
    fanned around it.

    ``operator_arc_resolver``: an optional ``(arc, source_builder,
    target_builder) -> segments | None`` callback consulted first; when it
    returns segments they are used as-is (and the arc skips the default
    dispatch). It lets a caller attach an arc to special geometry -- e.g. an
    SBGN operator's connector tips -- without this generic routine knowing about
    those classes. ``None`` returned (or no callback) keeps the default.
    """
    arcs_by_node_pair = _group_arcs_by_node_pair(
        new_layout_builder,
        id_to_new_layout_element_builder,
        descendant_id_to_top_level_id,
        operator_arc_resolver=operator_arc_resolver,
    )
    for node_pair, placements in arcs_by_node_pair.items():
        first_node_id, second_node_id = node_pair
        if first_node_id == second_node_id:
            for index, placement in enumerate(placements):
                start_angle, end_angle = _make_self_loop_angles(index, len(placements))
                placement.arc.segments = _make_self_loop_segments(
                    placement.source_layout_element, start_angle, end_angle
                )
            continue
        frame_start_center = id_to_new_layout_element_builder[first_node_id].center()
        frame_end_center = id_to_new_layout_element_builder[second_node_id].center()
        offsets = _make_offset_ladder(len(placements))
        for offset, placement in zip(offsets, placements):
            placement.arc.segments = _make_offset_segments(
                placement.source_layout_element,
                placement.target_layout_element,
                frame_start_center,
                frame_end_center,
                offset,
            )


def make_auto_layout(
    cd_map: typing.Any,
    compartment_layout_classes: typing.Any = (
        momapy.celldesigner.RectangleCompartmentLayout,
        momapy.celldesigner.OvalCompartmentLayout,
    ),
    reversed_arc_classes: typing.Any = (),
    operator_arc_resolver: typing.Any = None,
) -> typing.Any:
    """Reposition an already-built layout with graphviz (dot).

    Generic over the map language except for the compartment-layout classes,
    which differ per language: pass ``compartment_layout_classes`` to identify
    the compartment containers (so they become dot clusters rather than nodes).
    Defaults to the CellDesigner compartment-layout classes.

    ``reversed_arc_classes`` flips the dot-edge direction of the given arc
    classes for ranking only (see :func:`_build_dot_graph`);
    ``operator_arc_resolver`` overrides the rebuilt segments of selected arcs
    (see :func:`_rebuild_arc_geometry`). Both default to inert, so callers
    that pass neither (e.g. CellDesigner) are unaffected.

    The root layout of the returned map is fitted around the repositioned
    elements (see :func:`harmonize_root_layout`).
    """
    new_map_builder = momapy.builder.builder_from_object(cd_map)
    new_layout_builder = new_map_builder.layout
    (
        dot_graph,
        id_to_layout_element,
        descendant_id_to_top_level_id,
        dot_cluster_name_to_compartment_layout_element,
    ) = _build_dot_graph(
        new_map_builder,
        new_layout_builder,
        compartment_layout_classes,
        reversed_arc_classes=reversed_arc_classes,
    )
    id_to_new_layout_element_builder = _reposition_from_dot(
        dot_graph,
        id_to_layout_element,
        dot_cluster_name_to_compartment_layout_element,
    )
    _rebuild_arc_geometry(
        new_layout_builder,
        id_to_new_layout_element_builder,
        descendant_id_to_top_level_id,
        operator_arc_resolver=operator_arc_resolver,
    )
    new_map = momapy.builder.object_from_builder(new_map_builder)
    return dataclasses.replace(new_map, layout=harmonize_root_layout(new_map.layout))
