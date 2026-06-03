import momapy.celldesigner

import pd2af.predicates
import pd2af._building_layout


class TestModulationLayoutMap:
    def test_modulation_class_to_layout_class(self):
        mapping = pd2af._building_layout._MODULATION_CLASS_TO_LAYOUT_CLASS
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
        layout_map = pd2af._building_layout._MODULATION_CLASS_TO_LAYOUT_CLASS
        for model_class in output_classes:
            assert model_class in layout_map, model_class
