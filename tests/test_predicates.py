import clorm

import momapy.celldesigner

import pd2af.predicates


class TestPredicateClass:
    def test_kept_species_is_clorm_predicate(self):
        assert issubclass(pd2af.predicates.kept_species, clorm.Predicate)

    def test_derived_proteoform_class_is_clorm_predicate(self):
        assert issubclass(
            pd2af.predicates.derived_proteoform_class, clorm.Predicate
        )

    def test_activity_is_clorm_predicate(self):
        assert issubclass(pd2af.predicates.activity, clorm.Predicate)

    def test_positively_influences_is_clorm_predicate(self):
        assert issubclass(pd2af.predicates.positivelyInfluences, clorm.Predicate)

    def test_negatively_influences_is_clorm_predicate(self):
        assert issubclass(pd2af.predicates.negativelyInfluences, clorm.Predicate)

    def test_new_is_clorm_predicate(self):
        assert issubclass(pd2af.predicates.new, clorm.Predicate)

    def test_kept_species_can_be_constructed(self):
        atom = pd2af.predicates.kept_species(species="some_id")
        assert atom.species == "some_id"

    def test_derived_proteoform_class_can_be_constructed(self):
        atom = pd2af.predicates.derived_proteoform_class(
            template="t", compartment="c"
        )
        assert atom.template == "t"
        assert atom.compartment == "c"

    def test_activity_wraps_kept_species(self):
        atom = pd2af.predicates.activity(
            key=pd2af.predicates.kept_species(species="x")
        )
        assert isinstance(atom.key, pd2af.predicates.kept_species)
        assert atom.key.species == "x"

    def test_activity_wraps_derived_proteoform_class(self):
        atom = pd2af.predicates.activity(
            key=pd2af.predicates.derived_proteoform_class(
                template="t", compartment="c"
            )
        )
        assert isinstance(atom.key, pd2af.predicates.derived_proteoform_class)
        assert atom.key.template == "t"
        assert atom.key.compartment == "c"

    def test_positively_influences_with_key_wrappers(self):
        atom = pd2af.predicates.positivelyInfluences(
            source=pd2af.predicates.kept_species(species="a"),
            target=pd2af.predicates.derived_proteoform_class(
                template="t", compartment="c"
            ),
        )
        assert isinstance(atom.source, pd2af.predicates.kept_species)
        assert isinstance(atom.target, pd2af.predicates.derived_proteoform_class)

    def test_negatively_influences_with_key_wrappers(self):
        atom = pd2af.predicates.negativelyInfluences(
            source=pd2af.predicates.kept_species(species="a"),
            target=pd2af.predicates.kept_species(species="b"),
        )
        assert isinstance(atom.source, pd2af.predicates.kept_species)
        assert isinstance(atom.target, pd2af.predicates.kept_species)

    def test_new_wraps_activity(self):
        atom = pd2af.predicates.new(
            object_=pd2af.predicates.activity(
                key=pd2af.predicates.kept_species(species="x")
            )
        )
        assert isinstance(atom.object_, pd2af.predicates.activity)


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
