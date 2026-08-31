import momapy.celldesigner
import momapy.geometry

import pd2af.predicates
import pd2af.celldesigner.building_layout
import pd2af.utils


def _node(position):
    return momapy.celldesigner.GenericProteinLayout(position=position)


class TestModulationLayoutMap:
    def test_modulation_class_to_layout_class(self):
        mapping = pd2af.celldesigner.building_layout._MODULATION_CLASS_TO_LAYOUT_CLASS
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
        output_classes = pd2af.predicates.predicate_to_model_element_class.values()
        layout_map = (
            pd2af.celldesigner.building_layout._MODULATION_CLASS_TO_LAYOUT_CLASS
        )
        for model_class in output_classes:
            assert model_class in layout_map, model_class


class TestNearestInfluencePairing:
    def test_picks_the_closest_source_target_pair(self):
        target = _node(momapy.geometry.Point(0.0, 0.0))
        near_source = _node(momapy.geometry.Point(1.0, 0.0))
        far_source = _node(momapy.geometry.Point(10.0, 0.0))
        source_layout, target_layout = pd2af.utils.nearest_layout_pair(
            (far_source, near_source), (target,)
        )
        assert source_layout is near_source
        assert target_layout is target

    def test_tie_break_keeps_first_pair_in_product_order(self):
        target = _node(momapy.geometry.Point(0.0, 0.0))
        first_source = _node(momapy.geometry.Point(0.0, 5.0))
        second_source = _node(momapy.geometry.Point(5.0, 0.0))
        source_layout, _ = pd2af.utils.nearest_layout_pair(
            (first_source, second_source), (target,)
        )
        assert source_layout is first_source


class TestArcOffsetLadder:
    def test_a_lone_arc_stays_on_the_straight_line(self):
        assert pd2af.utils._make_offset_ladder(1) == [0.0]

    def test_an_even_count_straddles_the_line(self):
        # No arc of an even group is straight: they pair off on either side.
        offset = pd2af.utils._BEZIER_OFFSET
        assert pd2af.utils._make_offset_ladder(2) == [-offset, offset]
        assert pd2af.utils._make_offset_ladder(4) == [
            -2 * offset,
            -offset,
            offset,
            2 * offset,
        ]

    def test_an_odd_count_keeps_a_middle_arc_straight(self):
        offset = pd2af.utils._BEZIER_OFFSET
        assert pd2af.utils._make_offset_ladder(3) == [-offset, 0.0, offset]

    def test_offsets_are_all_distinct(self):
        for count in range(1, 8):
            offsets = pd2af.utils._make_offset_ladder(count)
            assert len(offsets) == count
            assert len(set(offsets)) == count


class TestArcOffsetGeometry:
    def _control_point(self, source_layout, target_layout, frame, offset):
        frame_start_layout, frame_end_layout = frame
        segments = pd2af.utils._make_offset_segments(
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
        forward_offset, backward_offset = pd2af.utils._make_offset_ladder(2)
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
            for offset in pd2af.utils._make_offset_ladder(4)
        }
        assert len(control_points) == 4

    def test_a_zero_offset_is_a_straight_segment(self):
        source = _node(momapy.geometry.Point(0.0, 0.0))
        target = _node(momapy.geometry.Point(100.0, 0.0))
        segments = pd2af.utils._make_offset_segments(
            source, target, source.center(), target.center(), 0.0
        )
        assert len(segments) == 1


class TestSelfLoopFan:
    def test_a_lone_self_loop_sits_on_top_of_its_node(self):
        assert pd2af.utils._make_self_loop_angles(0, 1) == (120.0, 60.0)

    def test_parallel_self_loops_are_spread_around_the_node(self):
        angles = [pd2af.utils._make_self_loop_angles(index, 3) for index in range(3)]
        assert len(set(angles)) == 3
        # Each loop keeps the same span; only its center moves.
        for start_angle, end_angle in angles:
            assert (start_angle - end_angle) % 360 == pd2af.utils._SELF_LOOP_SPAN
