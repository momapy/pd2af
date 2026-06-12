import clorm
import pytest

import momapy.celldesigner

import pd2af.core
import pd2af.predicates
import pd2af.solver


_INFLUENCE_PREDICATES = (
    pd2af.predicates.positivelyInfluences,
    pd2af.predicates.negativelyInfluences,
    pd2af.predicates.modulates,
    pd2af.predicates.triggers,
    pd2af.predicates.unknownPositivelyInfluences,
    pd2af.predicates.unknownNegativelyInfluences,
    pd2af.predicates.unknownModulates,
    pd2af.predicates.unknownTriggers,
)


def _activity_atoms(clingo_model):
    return [
        fact.object_
        for fact in clingo_model
        if isinstance(fact, pd2af.predicates.new)
        and isinstance(fact.object_, pd2af.predicates.activity)
    ]


def _influence_atoms(clingo_model):
    return [
        fact.object_
        for fact in clingo_model
        if isinstance(fact, pd2af.predicates.new)
        and isinstance(fact.object_, _INFLUENCE_PREDICATES)
    ]


@pytest.fixture(scope="module")
def solved_keep_species(example_cd_map):
    return pd2af.solver.solve(example_cd_map, mode="keep-species")


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
        # B, D, E, F, G. Active subunit C of complex D no longer gets its own
        # activity -- it resolves to D's top-level `keptSpeciesKey(D)`, the same
        # key D itself carries, so it adds no distinct activity atom.
        clingo_model, _ = solved_keep_species
        atoms = _activity_atoms(clingo_model)
        assert len(atoms) == 5

    def test_solve_keep_species_finds_four_influence_atoms(
        self, solved_keep_species
    ):
        clingo_model, _ = solved_keep_species
        atoms = _influence_atoms(clingo_model)
        assert len(atoms) == 4

    def test_solve_keep_species_no_complex_excludes_complex_with_active_subunit(
        self, example_cd_map
    ):
        clingo_model, id_to_model_element = pd2af.solver.solve(
            example_cd_map, mode="keep-species-no-complex"
        )
        atoms = _activity_atoms(clingo_model)
        names = sorted(
            id_to_model_element[atom.key.species].name for atom in atoms
        )
        assert "D" not in names
        assert "C" in names

    def test_solve_unknown_mode_raises(self, example_cd_map):
        with pytest.raises(ValueError):
            pd2af.solver.solve(example_cd_map, mode="not-a-mode")


def test_supported_modes():
    assert pd2af.core._TRANSFORMATION_MODES == frozenset(
        {
            "normal",
            "normal-no-complex",
            "keep-species",
            "keep-species-no-complex",
            "casq",
        }
    )


class TestSolveSbgnPdMergedModes:
    """Regression guard for the SBGN-PD carrier bug: `normal`/`normal-no-complex`
    over SBGN-PD input must emit a non-empty influence set (was zero before
    the entity-pool carrier was added to these profiles)."""

    @pytest.mark.parametrize("mode", ("normal", "normal-no-complex"))
    def test_sbgn_pd_merged_mode_emits_influences(self, sbgn_example_map, mode):
        clingo_model, _ = pd2af.solver.solve(sbgn_example_map, mode=mode)
        assert len(_activity_atoms(clingo_model)) > 0
        assert len(_influence_atoms(clingo_model)) > 0


@pytest.fixture(scope="module")
def solved_casq(example_cd_map):
    return pd2af.solver.solve(example_cd_map, mode="casq")


class TestSolveCasq:
    def test_returns_factbase_and_id_map(self, solved_casq):
        clingo_model, id_to_model_element = solved_casq
        assert isinstance(clingo_model, clorm.FactBase)
        assert isinstance(id_to_model_element, dict)

    def test_activity_keys_are_kept_species_key(self, solved_casq):
        clingo_model, _ = solved_casq
        atoms = _activity_atoms(clingo_model)
        assert len(atoms) > 0
        for atom in atoms:
            assert isinstance(atom.key, pd2af.predicates.keptSpeciesKey)

    def test_emits_some_influences(self, solved_casq):
        clingo_model, _ = solved_casq
        atoms = _influence_atoms(clingo_model)
        assert len(atoms) > 0
