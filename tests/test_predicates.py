import clorm

import momapy.celldesigner

import pd2af.predicates


class TestPredicateClass:
    def test_kept_species_key_is_clorm_predicate(self):
        assert issubclass(pd2af.predicates.keptSpeciesKey, clorm.Predicate)

    def test_promoted_subunit_key_is_clorm_predicate(self):
        assert issubclass(pd2af.predicates.promotedSubunitKey, clorm.Predicate)

    def test_activity_is_clorm_predicate(self):
        assert issubclass(pd2af.predicates.activity, clorm.Predicate)

    def test_positively_influences_is_clorm_predicate(self):
        assert issubclass(pd2af.predicates.positivelyInfluences, clorm.Predicate)

    def test_negatively_influences_is_clorm_predicate(self):
        assert issubclass(pd2af.predicates.negativelyInfluences, clorm.Predicate)

    def test_new_is_clorm_predicate(self):
        assert issubclass(pd2af.predicates.new, clorm.Predicate)

    def test_kept_species_key_can_be_constructed(self):
        atom = pd2af.predicates.keptSpeciesKey(species="some_id")
        assert atom.species == "some_id"

    def test_promoted_subunit_key_can_be_constructed(self):
        atom = pd2af.predicates.promotedSubunitKey(species="some_id")
        assert atom.species == "some_id"

    def test_activity_wraps_kept_species_key(self):
        atom = pd2af.predicates.activity(
            key=pd2af.predicates.keptSpeciesKey(species="x")
        )
        assert isinstance(atom.key, pd2af.predicates.keptSpeciesKey)
        assert atom.key.species == "x"

    def test_activity_wraps_promoted_subunit_key(self):
        atom = pd2af.predicates.activity(
            key=pd2af.predicates.promotedSubunitKey(species="x")
        )
        assert isinstance(atom.key, pd2af.predicates.promotedSubunitKey)
        assert atom.key.species == "x"

    def test_positively_influences_with_key_wrappers(self):
        atom = pd2af.predicates.positivelyInfluences(
            source=pd2af.predicates.keptSpeciesKey(species="a"),
            target=pd2af.predicates.promotedSubunitKey(species="x"),
        )
        assert isinstance(atom.source, pd2af.predicates.keptSpeciesKey)
        assert isinstance(atom.target, pd2af.predicates.promotedSubunitKey)

    def test_negatively_influences_with_key_wrappers(self):
        atom = pd2af.predicates.negativelyInfluences(
            source=pd2af.predicates.keptSpeciesKey(species="a"),
            target=pd2af.predicates.keptSpeciesKey(species="b"),
        )
        assert isinstance(atom.source, pd2af.predicates.keptSpeciesKey)
        assert isinstance(atom.target, pd2af.predicates.keptSpeciesKey)

    def test_new_wraps_activity(self):
        atom = pd2af.predicates.new(
            object_=pd2af.predicates.activity(
                key=pd2af.predicates.keptSpeciesKey(species="x")
            )
        )
        assert isinstance(atom.object_, pd2af.predicates.activity)


class TestTypedInfluencePredicates:
    """The eight typed influence predicates emitted into ``new(...)``."""

    PREDICATE_NAMES = (
        "positivelyInfluences",
        "negativelyInfluences",
        "modulates",
        "triggers",
        "unknownPositivelyInfluences",
        "unknownNegativelyInfluences",
        "unknownModulates",
        "unknownTriggers",
    )

    def test_all_are_clorm_predicates(self):
        for name in self.PREDICATE_NAMES:
            predicate = getattr(pd2af.predicates, name)
            assert issubclass(predicate, clorm.Predicate), name

    def test_can_be_constructed_with_key_wrappers(self):
        for name in self.PREDICATE_NAMES:
            predicate = getattr(pd2af.predicates, name)
            atom = predicate(
                source=pd2af.predicates.keptSpeciesKey(species="a"),
                target=pd2af.predicates.keptSpeciesKey(species="b"),
            )
            assert isinstance(atom.source, pd2af.predicates.keptSpeciesKey)
            assert isinstance(atom.target, pd2af.predicates.keptSpeciesKey)

    def test_new_wraps_each_typed_influence(self):
        for name in self.PREDICATE_NAMES:
            predicate = getattr(pd2af.predicates, name)
            atom = pd2af.predicates.new(
                object_=predicate(
                    source=pd2af.predicates.keptSpeciesKey(species="a"),
                    target=pd2af.predicates.keptSpeciesKey(species="b"),
                )
            )
            assert isinstance(atom.object_, predicate), name


class TestMappingDicts:
    def test_predicate_to_model_element_class(self):
        mapping = pd2af.predicates.predicate_to_model_element_class
        celldesigner = momapy.celldesigner
        assert mapping == {
            pd2af.predicates.positivelyInfluences: celldesigner.PositiveInfluence,
            pd2af.predicates.negativelyInfluences: celldesigner.NegativeInfluence,
            pd2af.predicates.modulates: celldesigner.Modulation,
            pd2af.predicates.triggers: celldesigner.Triggering,
            pd2af.predicates.unknownPositivelyInfluences: celldesigner.UnknownPositiveInfluence,
            pd2af.predicates.unknownNegativelyInfluences: celldesigner.UnknownNegativeInfluence,
            pd2af.predicates.unknownModulates: celldesigner.UnknownModulation,
            pd2af.predicates.unknownTriggers: celldesigner.UnknownTriggering,
        }
