import clorm

import momapy.celldesigner

import pd2af.predicates


class TestPredicateClass:
    def test_activity_is_clorm_predicate(self):
        assert issubclass(pd2af.predicates.activity, clorm.Predicate)

    def test_positively_influences_is_clorm_predicate(self):
        assert issubclass(pd2af.predicates.positivelyInfluences, clorm.Predicate)

    def test_negatively_influences_is_clorm_predicate(self):
        assert issubclass(pd2af.predicates.negativelyInfluences, clorm.Predicate)

    def test_new_is_clorm_predicate(self):
        assert issubclass(pd2af.predicates.new, clorm.Predicate)

    def test_activity_can_be_constructed(self):
        atom = pd2af.predicates.activity(name="some_id")
        assert atom.name == "some_id"

    def test_positively_influences_can_be_constructed(self):
        atom = pd2af.predicates.positivelyInfluences(source="a", target="b")
        assert atom.source == "a"
        assert atom.target == "b"

    def test_negatively_influences_can_be_constructed(self):
        atom = pd2af.predicates.negativelyInfluences(source="a", target="b")
        assert atom.source == "a"
        assert atom.target == "b"


class TestMappingDicts:
    def test_predicate_to_model_element_class(self):
        mapping = pd2af.predicates.predicate_to_model_element_class
        assert (
            mapping[pd2af.predicates.positivelyInfluences]
            is momapy.celldesigner.PositiveInfluence
        )
        assert (
            mapping[pd2af.predicates.negativelyInfluences]
            is momapy.celldesigner.Inhibition
        )

    def test_model_element_class_to_layout_element_class(self):
        mapping = pd2af.predicates.model_element_class_to_layout_element_class
        assert (
            mapping[momapy.celldesigner.PositiveInfluence]
            is momapy.celldesigner.PositiveInfluenceLayout
        )
        assert (
            mapping[momapy.celldesigner.Inhibition]
            is momapy.celldesigner.InhibitionLayout
        )
