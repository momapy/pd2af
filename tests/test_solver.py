import clorm
import pytest

import momapy.celldesigner

import pd2af.predicates
import pd2af.solver


@pytest.fixture(scope="module")
def solved_keep_species(example_cd_map):
    return pd2af.solver.solve(example_cd_map, mode="keep-species")


@pytest.fixture(scope="module")
def new_model_keep_species(solved_keep_species):
    clingo_model, id_to = solved_keep_species
    return clingo_model, id_to, pd2af.solver.make_new_cd_model(clingo_model, id_to)


class TestSolve:
    def test_solve_keep_species_returns_factbase_and_id_map(
        self, example_cd_map, solved_keep_species
    ):
        clingo_model, id_to_model_element = solved_keep_species
        assert isinstance(clingo_model, clorm.FactBase)
        assert isinstance(id_to_model_element, dict)
        assert len(id_to_model_element) > 0
        mapped_species_ids = {
            obj.id_
            for obj in id_to_model_element.values()
            if isinstance(obj, momapy.celldesigner.Species)
        }
        for species in example_cd_map.model.species:
            assert species.id_ in mapped_species_ids

    def test_solve_keep_species_finds_five_activity_atoms(
        self, solved_keep_species
    ):
        # B, D, E, F, G — active subunit C of complex D is subsumed.
        clingo_model, _ = solved_keep_species
        atoms = pd2af.solver._get_activity_atoms(clingo_model)
        assert len(atoms) == 5

    def test_solve_keep_species_finds_four_influence_atoms(
        self, solved_keep_species
    ):
        clingo_model, _ = solved_keep_species
        atoms = pd2af.solver._get_influence_atoms(clingo_model)
        assert len(atoms) == 4

    def test_solve_keep_species_no_complex_excludes_complex_with_active_subunit(
        self, example_cd_map
    ):
        clingo_model, id_to_model_element = pd2af.solver.solve(
            example_cd_map, mode="keep-species-no-complex"
        )
        atoms = pd2af.solver._get_activity_atoms(clingo_model)
        names = sorted(
            id_to_model_element[atom.key.species].name for atom in atoms
        )
        assert "D" not in names
        assert "C" in names

    def test_solve_unknown_mode_raises(self, example_cd_map):
        with pytest.raises(ValueError):
            pd2af.solver.solve(example_cd_map, mode="not-a-mode")


class TestMakeNewCdModel:
    def test_returns_celldesigner_model(self, new_model_keep_species):
        _, _, new_model = new_model_keep_species
        assert isinstance(new_model, momapy.celldesigner.CellDesignerModel)

    def test_species_count_matches_activity_atoms(self, new_model_keep_species):
        clingo_model, _, new_model = new_model_keep_species
        atoms = pd2af.solver._get_activity_atoms(clingo_model)
        assert len(new_model.species) == len(atoms)

    def test_modulations_have_expected_classes(self, new_model_keep_species):
        _, _, new_model = new_model_keep_species
        allowed = (
            momapy.celldesigner.PositiveInfluence,
            momapy.celldesigner.Inhibition,
        )
        assert len(new_model.modulations) > 0
        for mod in new_model.modulations:
            assert isinstance(mod, allowed)

    def test_modulation_endpoints_are_species_in_new_model(self, new_model_keep_species):
        _, _, new_model = new_model_keep_species
        species_ids = {s.id_ for s in new_model.species}
        for mod in new_model.modulations:
            assert mod.source.id_ in species_ids
            assert mod.target.id_ in species_ids

    def test_compartments_collected_from_kept_species(self, new_model_keep_species):
        _, _, new_model = new_model_keep_species
        expected = {
            s.compartment for s in new_model.species if s.compartment is not None
        }
        assert expected <= set(new_model.compartments)


def test_supported_modes():
    assert pd2af.solver._VALID_MODES == frozenset(
        {
            "normal",
            "no-complex",
            "keep-species",
            "keep-species-no-complex",
            "casq",
        }
    )


@pytest.fixture(scope="module")
def solved_casq(example_cd_map):
    return pd2af.solver.solve(example_cd_map, mode="casq")


class TestSolveCasq:
    def test_returns_factbase_and_id_map(self, solved_casq):
        clingo_model, id_to_model_element = solved_casq
        assert isinstance(clingo_model, clorm.FactBase)
        assert isinstance(id_to_model_element, dict)

    def test_activity_keys_are_kept_species(self, solved_casq):
        clingo_model, _ = solved_casq
        atoms = pd2af.solver._get_activity_atoms(clingo_model)
        assert len(atoms) > 0
        for atom in atoms:
            assert isinstance(atom.key, pd2af.predicates.kept_species)

    def test_emits_some_influences(self, solved_casq):
        clingo_model, _ = solved_casq
        atoms = pd2af.solver._get_influence_atoms(clingo_model)
        assert len(atoms) > 0

    def test_make_new_cd_model(self, solved_casq):
        clingo_model, id_to = solved_casq
        new_model = pd2af.solver.make_new_cd_model(clingo_model, id_to)
        assert isinstance(new_model, momapy.celldesigner.CellDesignerModel)
        species_ids = {s.id_ for s in new_model.species}
        for mod in new_model.modulations:
            assert mod.source.id_ in species_ids
            assert mod.target.id_ in species_ids
