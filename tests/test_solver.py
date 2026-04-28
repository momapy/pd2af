import clorm
import pytest

import momapy.celldesigner

import pd2af.predicates
import pd2af.solver


@pytest.fixture(scope="module")
def solved_default(example_cd_map):
    return pd2af.solver.solve(example_cd_map, mode="normal")


@pytest.fixture(scope="module")
def new_model_default(solved_default):
    clingo_model, id_to = solved_default
    return clingo_model, id_to, pd2af.solver.make_new_cd_model(clingo_model, id_to)


class TestSolve:
    def test_solve_default_mode_returns_factbase_and_id_map(
        self, example_cd_map, solved_default
    ):
        clingo_model, id_to_model_element = solved_default
        assert isinstance(clingo_model, clorm.FactBase)
        assert isinstance(id_to_model_element, dict)
        assert len(id_to_model_element) > 0
        # The id map should include every species from the source model.
        mapped_species_ids = {
            obj.id_
            for obj in id_to_model_element.values()
            if isinstance(obj, momapy.celldesigner.Species)
        }
        for species in example_cd_map.model.species:
            assert species.id_ in mapped_species_ids

    def test_solve_default_mode_finds_six_activity_atoms(self, solved_default):
        clingo_model, _ = solved_default
        atoms = pd2af.solver._get_activity_atoms(clingo_model)
        assert len(atoms) == 6

    def test_solve_default_mode_finds_four_influence_atoms(self, solved_default):
        clingo_model, _ = solved_default
        atoms = pd2af.solver._get_influence_atoms(clingo_model)
        assert len(atoms) == 4

    def test_solve_no_complex_excludes_complex_with_active_subunit(
        self, example_cd_map
    ):
        clingo_model, id_to_model_element = pd2af.solver.solve(
            example_cd_map, mode="no-complex"
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
    def test_returns_celldesigner_model(self, new_model_default):
        _, _, new_model = new_model_default
        assert isinstance(new_model, momapy.celldesigner.CellDesignerModel)

    def test_species_count_matches_activity_atoms(self, new_model_default):
        clingo_model, _, new_model = new_model_default
        atoms = pd2af.solver._get_activity_atoms(clingo_model)
        assert len(new_model.species) == len(atoms)

    def test_modulations_have_expected_classes(self, new_model_default):
        _, _, new_model = new_model_default
        allowed = (
            momapy.celldesigner.PositiveInfluence,
            momapy.celldesigner.Inhibition,
        )
        assert len(new_model.modulations) > 0
        for mod in new_model.modulations:
            assert isinstance(mod, allowed)

    def test_modulation_endpoints_are_species_in_new_model(self, new_model_default):
        _, _, new_model = new_model_default
        species_ids = {s.id_ for s in new_model.species}
        for mod in new_model.modulations:
            assert mod.source.id_ in species_ids
            assert mod.target.id_ in species_ids

    def test_compartments_collected_from_kept_species(self, new_model_default):
        _, _, new_model = new_model_default
        expected = {
            s.compartment for s in new_model.species if s.compartment is not None
        }
        assert expected <= set(new_model.compartments)


def test_supported_modes():
    assert set(pd2af.solver._MODE_TO_PROFILE) == {
        "normal",
        "no-complex",
        "pure-af",
    }
