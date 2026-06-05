import momapy.celldesigner
import momapy.geometry

import pd2af.predicates
import pd2af.celldesigner.building_layout


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
        layout_map = pd2af.celldesigner.building_layout._MODULATION_CLASS_TO_LAYOUT_CLASS
        for model_class in output_classes:
            assert model_class in layout_map, model_class


class TestNearestInfluencePairing:
    def test_picks_the_closest_source_target_pair(self):
        target = _node(momapy.geometry.Point(0.0, 0.0))
        near_source = _node(momapy.geometry.Point(1.0, 0.0))
        far_source = _node(momapy.geometry.Point(10.0, 0.0))
        source_layout, target_layout = pd2af.celldesigner.building_layout._nearest_layout_pair(
            (far_source, near_source), (target,)
        )
        assert source_layout is near_source
        assert target_layout is target

    def test_tie_break_keeps_first_pair_in_product_order(self):
        target = _node(momapy.geometry.Point(0.0, 0.0))
        first_source = _node(momapy.geometry.Point(0.0, 5.0))
        second_source = _node(momapy.geometry.Point(5.0, 0.0))
        source_layout, _ = pd2af.celldesigner.building_layout._nearest_layout_pair(
            (first_source, second_source), (target,)
        )
        assert source_layout is first_source
