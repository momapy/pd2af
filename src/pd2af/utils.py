import collections
import dataclasses
import itertools
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


def register_or_reuse(element, cache):
    """Intern ``element`` by content in ``cache``. First-registered wins:
    if a content-equal element is already cached, return it; otherwise
    record ``element`` as the canonical instance and return it.
    """
    existing = cache.get(element)
    if existing is not None:
        return existing
    cache[element] = element
    return element


def add_model_element_if_new(collection, model_element, seen_identities):
    """Append ``model_element`` to ``collection`` unless an object with the
    same identity was already appended (tracked in ``seen_identities``).

    ``model_element`` is assumed to already be the canonical instance (e.g. the
    result of :func:`register_or_reuse`); this only guards against adding the
    same identity twice. Returns ``True`` if it was added this call, ``False``
    if it was a duplicate.
    """
    if id(model_element) in seen_identities:
        return False
    seen_identities.add(id(model_element))
    collection.add(model_element)
    return True


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


def nearest_layout_pair(source_layouts, target_layouts):
    """Return the (source_layout, target_layout) pair whose ``.position``s are
    closest in Euclidean distance.

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


def influence_layout_pairs(source_layouts, target_layouts, prefer_nearest):
    """Return the (source_layout, target_layout) pairs to draw an influence /
    modulation arc between.

    By default the full ``source x target`` cross product. ``prefer_nearest``
    (``-p nearest``) collapses the fan-out to the single closest pair, but only
    where there *is* fan-out (more than one layout on a side) and every layout
    has a real ``position`` -- in ``auto`` mode positions are throwaway
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


def make_arc_segment_from_source_and_target(source_layout, target_layout):
    """Straight Segment connecting two node layouts for a modulation /
    influence arc. A self-loop (source is target) anchors to the node's
    north edge; otherwise the segment runs border-to-border, falling back
    to corner anchors when a border point is undefined."""
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
    return momapy.geometry.Segment(start_point, end_point)


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


def _build_dot_graph(
    new_map_builder,
    new_layout_builder,
    compartment_layout_classes,
    reversed_arc_classes=(),
):
    """Build the pydot graph from the built layout: compartments become dot
    clusters, Node layout elements become dot nodes (placed in their
    compartment's cluster when there is one), and Arc layout elements become
    dot edges. Returns the graph plus the bookkeeping the repositioning and
    arc-geometry phases need.

    ``reversed_arc_classes``: arc layout classes whose dot edge is added with
    source and target **swapped**, so graphviz ranks the arc's target upstream
    of its source. The layout arc object is unchanged; only the ranking flips.
    Used for SBGN logic arcs (stored operator -> input) so an operator's inputs
    rank above it and its output target below it."""
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
            # `outside` is a CellDesigner-only relation (compartment nesting);
            # SBGN compartments have no such concept, so don't assume it exists.
            outside_compartment = getattr(compartment, "outside", None)
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
        outside_compartment = getattr(compartment, "outside", None)
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
            compartment_layout_classes,
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
            compartment_layout_classes,
        ):
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
                compartment_dot_cluster = compartment_to_dot_cluster.get(
                    compartment
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
            if momapy.builder.isinstance_or_builder(
                layout_element_builder, reversed_arc_classes
            ):
                source_id, target_id = target_id, source_id
            dot_graph.add_edge(pydot.Edge(source_id, target_id))
            directed_pairs.add((source_id, target_id))
        id_to_layout_element[layout_element_builder.id_] = layout_element_builder
    dot_graph.set("ranksep", 1.0)
    dot_graph.set("nodesep", 0.5)
    dot_graph.set("rankdir", "BT")
    return (
        dot_graph,
        id_to_layout_element,
        descendant_id_to_top_level_id,
        directed_pairs,
        compartment_layout_element_to_included_layout_elements,
    )


def _reposition_from_dot(dot_graph, id_to_layout_element):
    """Run graphviz `dot`, then translate every layout element to the position
    dot computed for it. Returns id -> layout element builder (including nested
    descendants) so the arc-geometry phase can resolve endpoints."""
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
    return id_to_new_layout_element_builder


def _self_loop_segments(layout_element):
    """Segments for an arc whose source and target resolve to the same node: a
    loop bowing out through two control points, expressed as a polyline so the
    CellDesigner writer can recover them as edit points (it only sees segment
    endpoints). Without this, a single Bezier collapses to start/end on the same
    node and the reader's modulation-geometry call hits a None border."""
    start_point = layout_element.own_angle(120)
    end_point = layout_element.own_angle(60)
    if start_point is None:
        start_point = layout_element.north_west()
    if end_point is None:
        end_point = layout_element.north_east()
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


def _bidirectional_segments(
    source_layout_element,
    target_layout_element,
    source_top_level_id,
    target_top_level_id,
):
    """Segments for one arc of a bidirectional pair (A->B and B->A both exist):
    a curve bowing away from its reverse via a single offset control point,
    serialized as a polyline so the writer keeps the control point (see the
    self-loop note)."""
    source_center = source_layout_element.center()
    target_center = target_layout_element.center()
    delta_x = target_center.x - source_center.x
    delta_y = target_center.y - source_center.y
    length = math.hypot(delta_x, delta_y)
    if length == 0:
        start_point = source_layout_element.own_border(target_center)
        end_point = target_layout_element.own_border(source_center)
        if start_point is None:
            start_point = source_layout_element.north_west()
        if end_point is None:
            end_point = target_layout_element.north_east()
        return [momapy.geometry.Segment(start_point, end_point)]
    normal_x = -delta_y / length
    normal_y = delta_x / length
    # Deterministic side rule: A→B and B→A get opposite offsets, so the two
    # curves bow away from each other.
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
    start_point = source_layout_element.own_border(control_point)
    end_point = target_layout_element.own_border(control_point)
    if start_point is None:
        start_point = source_layout_element.north_west()
    if end_point is None:
        end_point = target_layout_element.north_east()
    return [
        momapy.geometry.Segment(start_point, control_point),
        momapy.geometry.Segment(control_point, end_point),
    ]


def _simple_segments(source_layout_element, target_layout_element):
    """Segments for a one-directional arc: a single straight segment between
    the two node borders."""
    start_point = source_layout_element.border(target_layout_element.center())
    end_point = target_layout_element.border(source_layout_element.center())
    if start_point is None:
        start_point = source_layout_element.north_west()
    if end_point is None:
        end_point = target_layout_element.north_east()
    return [momapy.geometry.Segment(start_point, end_point)]


def _arc_geometry(
    new_layout_builder,
    id_to_new_layout_element_builder,
    descendant_id_to_top_level_id,
    directed_pairs,
    operator_arc_resolver=None,
):
    """Rebuild each arc's segments from the repositioned node geometry,
    dispatching on whether the arc is a self-loop, one of a bidirectional pair,
    or a plain one-directional arc.

    ``operator_arc_resolver``: an optional ``(arc, source_builder,
    target_builder) -> segments | None`` callback consulted first; when it
    returns segments they are used as-is (and the arc skips the default
    dispatch). It lets a caller attach an arc to special geometry -- e.g. an
    SBGN operator's connector tips -- without this generic routine knowing about
    those classes. ``None`` returned (or no callback) keeps the default."""
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
        is_self_loop = source_top_level_id == target_top_level_id
        is_bidirectional = not is_self_loop and (
            target_top_level_id,
            source_top_level_id,
        ) in directed_pairs
        if is_self_loop:
            layout_element_builder.segments = _self_loop_segments(
                source_layout_element_builder
            )
        elif is_bidirectional:
            layout_element_builder.segments = _bidirectional_segments(
                source_layout_element_builder,
                target_layout_element_builder,
                source_top_level_id,
                target_top_level_id,
            )
        else:
            layout_element_builder.segments = _simple_segments(
                source_layout_element_builder,
                target_layout_element_builder,
            )


def auto_layout(
    cd_map,
    compartment_layout_classes=(
        momapy.celldesigner.RectangleCompartmentLayout,
        momapy.celldesigner.OvalCompartmentLayout,
    ),
    reversed_arc_classes=(),
    operator_arc_resolver=None,
):
    """Reposition an already-built layout with graphviz (dot).

    Generic over the map language except for the compartment-layout classes,
    which differ per language: pass ``compartment_layout_classes`` to identify
    the compartment containers (so they become dot clusters rather than nodes).
    Defaults to the CellDesigner compartment-layout classes.

    ``reversed_arc_classes`` flips the dot-edge direction of the given arc
    classes for ranking only (see :func:`_build_dot_graph`);
    ``operator_arc_resolver`` overrides the rebuilt segments of selected arcs
    (see :func:`_arc_geometry`). Both default to inert, so callers that pass
    neither (e.g. CellDesigner) are unaffected."""
    new_map_builder = momapy.builder.builder_from_object(cd_map)
    new_layout_builder = new_map_builder.layout
    (
        dot_graph,
        id_to_layout_element,
        descendant_id_to_top_level_id,
        directed_pairs,
        compartment_layout_element_to_included_layout_elements,
    ) = _build_dot_graph(
        new_map_builder,
        new_layout_builder,
        compartment_layout_classes,
        reversed_arc_classes=reversed_arc_classes,
    )
    id_to_new_layout_element_builder = _reposition_from_dot(
        dot_graph, id_to_layout_element
    )
    _arc_geometry(
        new_layout_builder,
        id_to_new_layout_element_builder,
        descendant_id_to_top_level_id,
        directed_pairs,
        operator_arc_resolver=operator_arc_resolver,
    )
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
