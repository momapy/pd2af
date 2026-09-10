import os
import tempfile
import types

import pytest

import momapy.celldesigner
import momapy.core.layout
import momapy.core.mapping
import momapy.geometry
import momapy.io.core
import momapy.sbgn.pd

import pd2af
import pd2af.asp.predicates
import pd2af.building.celldesigner.layout
import pd2af.building.layout

from tests._helpers import MAPS_DIR, has_dot_binary, read_cd_map


def _node(position):
    return momapy.celldesigner.GenericProteinLayout(position=position)


class TestModulationLayoutMap:
    def test_modulation_class_to_layout_class(self):
        mapping = pd2af.building.celldesigner.layout._MODULATION_CLASS_TO_LAYOUT_CLASS
        celldesigner = momapy.celldesigner
        # NegativeInfluence / UnknownNegativeInfluence have no own *Layout
        # class; they reuse the inhibition arc layouts (mirrors the reader).
        assert mapping == {
            celldesigner.PositiveInfluence: celldesigner.PositiveInfluenceLayout,
            celldesigner.NegativeInfluence: celldesigner.InhibitionLayout,
            celldesigner.Modulation: celldesigner.ModulationLayout,
            celldesigner.Triggering: celldesigner.TriggeringLayout,
            celldesigner.UnknownPositiveInfluence: celldesigner.UnknownPositiveInfluenceLayout,
            celldesigner.UnknownNegativeInfluence: celldesigner.UnknownInhibitionLayout,
            celldesigner.UnknownModulation: celldesigner.UnknownModulationLayout,
            celldesigner.UnknownTriggering: celldesigner.UnknownTriggeringLayout,
        }

    def test_every_output_influence_class_has_a_layout(self):
        # Cross-module invariant: every model class an output influence
        # predicate maps to must have an arc layout class to draw it.
        output_classes = pd2af.asp.predicates.predicate_to_model_element_class.values()
        layout_map = (
            pd2af.building.celldesigner.layout._MODULATION_CLASS_TO_LAYOUT_CLASS
        )
        for model_class in output_classes:
            assert model_class in layout_map, model_class


class TestNearestInfluencePairing:
    def test_picks_the_closest_source_target_pair(self):
        target = _node(momapy.geometry.Point(0.0, 0.0))
        near_source = _node(momapy.geometry.Point(1.0, 0.0))
        far_source = _node(momapy.geometry.Point(10.0, 0.0))
        source_layout, target_layout = pd2af.building.layout.nearest_layout_pair(
            (far_source, near_source), (target,)
        )
        assert source_layout is near_source
        assert target_layout is target

    def test_tie_break_keeps_first_pair_in_product_order(self):
        target = _node(momapy.geometry.Point(0.0, 0.0))
        first_source = _node(momapy.geometry.Point(0.0, 5.0))
        second_source = _node(momapy.geometry.Point(5.0, 0.0))
        source_layout, _ = pd2af.building.layout.nearest_layout_pair(
            (first_source, second_source), (target,)
        )
        assert source_layout is first_source


class TestArcOffsetLadder:
    def test_a_lone_arc_stays_on_the_straight_line(self):
        assert pd2af.building.layout._make_offset_ladder(1) == [0.0]

    def test_an_even_count_straddles_the_line(self):
        # No arc of an even group is straight: they pair off on either side.
        offset = pd2af.building.layout._BEZIER_OFFSET
        assert pd2af.building.layout._make_offset_ladder(2) == [-offset, offset]
        assert pd2af.building.layout._make_offset_ladder(4) == [
            -2 * offset,
            -offset,
            offset,
            2 * offset,
        ]

    def test_an_odd_count_keeps_a_middle_arc_straight(self):
        offset = pd2af.building.layout._BEZIER_OFFSET
        assert pd2af.building.layout._make_offset_ladder(3) == [-offset, 0.0, offset]

    def test_offsets_are_all_distinct(self):
        for count in range(1, 8):
            offsets = pd2af.building.layout._make_offset_ladder(count)
            assert len(offsets) == count
            assert len(set(offsets)) == count


class TestArcOffsetGeometry:
    def _control_point(self, source_layout, target_layout, frame, offset):
        frame_start_layout, frame_end_layout = frame
        segments = pd2af.building.layout._make_offset_segments(
            source_layout,
            target_layout,
            frame_start_layout.center(),
            frame_end_layout.center(),
            offset,
        )
        return segments[0].p2

    def test_opposite_directions_bow_to_opposite_sides(self):
        # The canonical frame is what makes this hold: computing the normal
        # from each arc's own direction would negate it a second time and put
        # both arcs on the same curve.
        source = _node(momapy.geometry.Point(0.0, 0.0))
        target = _node(momapy.geometry.Point(100.0, 0.0))
        frame = (source, target)
        forward_offset, backward_offset = pd2af.building.layout._make_offset_ladder(2)
        forward_control_point = self._control_point(
            source, target, frame, forward_offset
        )
        backward_control_point = self._control_point(
            target, source, frame, backward_offset
        )
        assert forward_control_point.x == backward_control_point.x
        assert forward_control_point.y == -backward_control_point.y
        assert forward_control_point.y != 0.0

    def test_same_direction_arcs_get_distinct_geometries(self):
        source = _node(momapy.geometry.Point(0.0, 0.0))
        target = _node(momapy.geometry.Point(100.0, 0.0))
        control_points = {
            self._control_point(source, target, (source, target), offset)
            for offset in pd2af.building.layout._make_offset_ladder(4)
        }
        assert len(control_points) == 4

    def test_a_zero_offset_is_a_straight_segment(self):
        source = _node(momapy.geometry.Point(0.0, 0.0))
        target = _node(momapy.geometry.Point(100.0, 0.0))
        segments = pd2af.building.layout._make_offset_segments(
            source, target, source.center(), target.center(), 0.0
        )
        assert len(segments) == 1


class TestSelfLoopFan:
    def test_a_lone_self_loop_sits_on_top_of_its_node(self):
        assert pd2af.building.layout._make_self_loop_angles(0, 1) == (120.0, 60.0)

    def test_parallel_self_loops_are_spread_around_the_node(self):
        angles = [
            pd2af.building.layout._make_self_loop_angles(index, 3) for index in range(3)
        ]
        assert len(set(angles)) == 3
        # Each loop keeps the same span; only its center moves.
        for start_angle, end_angle in angles:
            assert (
                start_angle - end_angle
            ) % 360 == pd2af.building.layout._SELF_LOOP_SPAN


def _compartment_layout_element(x_offset=0.0):
    return momapy.celldesigner.RectangleCompartmentLayout(
        position=momapy.geometry.Point(x_offset, 0.0), width=100.0, height=100.0
    )


class TestDotGraphCompartmentClusters:
    """Every compartment cluster must end up attached: to its outside
    compartment's cluster when it has one, to the root graph otherwise. An
    unattached cluster is invisible to graphviz, so its member nodes drop out
    of the graph and the compartment never gets computed bounds."""

    @staticmethod
    def _build_dot_graph(compartments, compartment_layouts, member_nodes=()):
        mapping_builder = momapy.core.mapping.LayoutModelMappingBuilder()
        for compartment, compartment_layout in zip(compartments, compartment_layouts):
            if compartment_layout is not None:
                mapping_builder.add_mapping(compartment_layout, compartment)
        for node, member_of in member_nodes:
            mapping_builder.add_mapping(
                node, types.SimpleNamespace(compartment=member_of)
            )
        map_builder = types.SimpleNamespace(
            model=types.SimpleNamespace(compartments=list(compartments)),
            layout_model_mapping=mapping_builder,
        )
        layout_builder = types.SimpleNamespace(
            layout_elements=[node for node, _ in member_nodes]
        )
        dot_graph, *_ = pd2af.building.layout._build_dot_graph(
            map_builder, layout_builder, (), ()
        )
        return dot_graph

    def test_root_celldesigner_compartment_cluster_is_attached(self):
        compartment = momapy.celldesigner.Compartment(id_="c")
        dot_graph = self._build_dot_graph(
            [compartment], [_compartment_layout_element()]
        )
        assert [cluster.get_name() for cluster in dot_graph.get_subgraphs()] == [
            "cluster_c"
        ]

    def test_sbgn_compartment_cluster_is_attached(self):
        # An SBGN compartment has no `outside` attribute at all.
        compartment = momapy.sbgn.pd.Compartment(id_="c", label=None)
        dot_graph = self._build_dot_graph(
            [compartment], [_compartment_layout_element()]
        )
        assert [cluster.get_name() for cluster in dot_graph.get_subgraphs()] == [
            "cluster_c"
        ]

    def test_nested_cluster_attaches_to_its_parent(self):
        parent = momapy.celldesigner.Compartment(id_="parent")
        child = momapy.celldesigner.Compartment(id_="child", outside=parent)
        # Child before parent: attachment must not depend on iteration order.
        dot_graph = self._build_dot_graph(
            [child, parent],
            [_compartment_layout_element(200.0), _compartment_layout_element()],
        )
        (root_cluster,) = dot_graph.get_subgraphs()
        assert root_cluster.get_name() == "cluster_parent"
        assert [cluster.get_name() for cluster in root_cluster.get_subgraphs()] == [
            "cluster_child"
        ]

    def test_nested_cluster_attaches_to_root_when_parent_has_no_layout(self):
        parent = momapy.celldesigner.Compartment(id_="parent")
        child = momapy.celldesigner.Compartment(id_="child", outside=parent)
        dot_graph = self._build_dot_graph(
            [child, parent], [_compartment_layout_element(200.0), None]
        )
        assert [cluster.get_name() for cluster in dot_graph.get_subgraphs()] == [
            "cluster_child"
        ]

    def test_member_node_lands_inside_its_compartment_cluster(self):
        # A member without any influence arc (no dot edges) must still be
        # handed to graphviz inside its compartment's cluster.
        compartment = momapy.celldesigner.Compartment(id_="c")
        node = _node(momapy.geometry.Point(10.0, 10.0))
        dot_graph = self._build_dot_graph(
            [compartment], [_compartment_layout_element()], [(node, compartment)]
        )
        (cluster,) = dot_graph.get_subgraphs()
        assert [node_layout.get_name() for node_layout in cluster.get_nodes()] == [
            node.id_
        ]
        assert dot_graph.get_nodes() == []


_CELLDESIGNER_COMPARTMENT_LAYOUT_CLASSES = tuple(
    layout_class
    for layout_class in vars(momapy.celldesigner).values()
    if isinstance(layout_class, type)
    and layout_class.__name__.endswith("CompartmentLayout")
)

_SNCA_EXPRESSION_MAP_PATH = os.path.join(MAPS_DIR, "SNCA_expression.xml")


def _is_inside(compartment_layout, node_layout):
    return (
        compartment_layout.position.x - compartment_layout.width / 2
        <= node_layout.position.x
        <= compartment_layout.position.x + compartment_layout.width / 2
        and compartment_layout.position.y - compartment_layout.height / 2
        <= node_layout.position.y
        <= compartment_layout.position.y + compartment_layout.height / 2
    )


def _celldesigner_compartment_layouts(layout):
    return [
        element
        for element in layout.layout_elements
        if isinstance(element, _CELLDESIGNER_COMPARTMENT_LAYOUT_CLASSES)
    ]


class TestCelldesignerAutoLayoutCompartments:
    """End to end, auto layout: each compartment's dot-computed bounds enclose
    its members. Before the fix, the clusters never reached graphviz, so the
    compartments kept their curated input geometry while members moved away."""

    @pytest.mark.parametrize("keep_species", (True, False))
    def test_compartments_enclose_their_members(self, keep_species):
        if not has_dot_binary():
            pytest.skip("graphviz `dot` binary not on PATH")
        out = pd2af.transform(
            read_cd_map(_SNCA_EXPRESSION_MAP_PATH),
            mode="normal",
            keep_species=keep_species,
            layout_mode="auto",
        ).obj
        mapping = out.layout_model_mapping
        compartment_layouts = _celldesigner_compartment_layouts(out.layout)
        assert compartment_layouts
        for compartment_layout in compartment_layouts:
            compartment = mapping.get_mapping(compartment_layout)
            members = [
                element
                for element in out.layout.layout_elements
                if isinstance(element, momapy.core.layout.Node)
                and getattr(mapping.get_mapping(element), "compartment", None)
                is compartment
            ]
            assert members
            for member in members:
                assert _is_inside(compartment_layout, member)

    def test_compartment_map_round_trips(self):
        if not has_dot_binary():
            pytest.skip("graphviz `dot` binary not on PATH")
        out = pd2af.transform(
            read_cd_map(_SNCA_EXPRESSION_MAP_PATH), mode="normal", layout_mode="auto"
        ).obj
        path = os.path.join(tempfile.gettempdir(), "pd2af_test_compartments.xml")
        momapy.io.core.write(out, path, writer="celldesigner")
        back = momapy.io.core.read(path).obj
        assert len(back.model.compartments) == len(out.model.compartments)
        assert sorted(species.name for species in back.model.species) == sorted(
            species.name for species in out.model.species
        )
        assert len(back.model.modulations) == len(out.model.modulations)
