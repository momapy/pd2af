import clorm
import pytest

import momapy.celldesigner

import pd2af
import pd2af.modes
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

    def test_solve_keep_species_finds_five_activity_atoms(self, solved_keep_species):
        # B, D, E, F, G. Active subunit C of complex D no longer gets its own
        # activity -- it resolves to D's top-level `keptSpeciesKey(D)`, the same
        # key D itself carries, so it adds no distinct activity atom.
        clingo_model, _ = solved_keep_species
        atoms = _activity_atoms(clingo_model)
        assert len(atoms) == 5

    def test_solve_keep_species_finds_four_influence_atoms(self, solved_keep_species):
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
        names = sorted(id_to_model_element[atom.key.species].name for atom in atoms)
        assert "D" not in names
        assert "C" in names

    def test_solve_unknown_mode_raises(self, example_cd_map):
        with pytest.raises(ValueError):
            pd2af.solver.solve(example_cd_map, mode="not-a-mode")


def test_supported_modes():
    # The entry point is the only way a sixth mode can appear: pd2af itself
    # declares exactly these five, and every mode it offers is one of them
    # unless something contributed it.
    builtin_names = {mode.name for mode in pd2af.modes._BUILTIN_TRANSFORMATION_MODES}
    assert builtin_names == {
        "normal",
        "normal-no-complex",
        "keep-species",
        "keep-species-no-complex",
        "keep-reactions",
    }
    assert set(pd2af.get_transformation_modes()) >= builtin_names


class TestSolveSbgnPdMergedModes:
    """`normal`/`normal-no-complex` over SBGN-PD input must emit a non-empty
    influence set, which takes the entity-pool carrier these modes use."""

    @pytest.mark.parametrize("mode", ("normal", "normal-no-complex"))
    def test_sbgn_pd_merged_mode_emits_influences(self, sbgn_example_map, mode):
        clingo_model, _ = pd2af.solver.solve(sbgn_example_map, mode=mode)
        assert len(_activity_atoms(clingo_model)) > 0
        assert len(_influence_atoms(clingo_model)) > 0


def _named_influence_edges(clingo_model, id_to_model_element):
    """Return the set of (influence-class-name, source-name, target-name)
    tuples for every emitted influence atom, resolving the activity keys
    back to species names."""
    edges = set()
    for atom in _influence_atoms(clingo_model):
        source = id_to_model_element[atom.source.species].name
        target = id_to_model_element[atom.target.species].name
        edges.add((type(atom).__name__, source, target))
    return edges


def _species(name, active=False):
    return momapy.celldesigner.Species(name=name, active=active)


def _reaction(reactant, product, modifiers=()):
    return momapy.celldesigner.Reaction(
        reversible=False,
        reactants=frozenset({momapy.celldesigner.Reactant(referred_element=reactant)}),
        products=frozenset({momapy.celldesigner.Product(referred_element=product)}),
        modifiers=frozenset(modifiers),
    )


def _map_from(species, reactions):
    model = momapy.celldesigner.CellDesignerModel(
        species=frozenset(species), reactions=frozenset(reactions)
    )
    return momapy.celldesigner.CellDesignerMap(model=model)


class TestCycleAwareInfluences:
    """A source feeding a production cycle must not leak influence back
    around the loop onto members it directly depletes; transitivity across
    non-cycle edges is untouched. See plans/cycle-aware-influence-paths.md."""

    def test_cyclic_production_blocks_leaked_influence(self):
        # A catalyses B->C; C->B closes the production cycle; B and C active.
        species_a = _species("A")
        species_b = _species("B", active=True)
        species_c = _species("C", active=True)
        catalyzer = momapy.celldesigner.Catalyzer(referred_element=species_a)
        cyclic_map = _map_from(
            (species_a, species_b, species_c),
            (
                _reaction(species_b, species_c, modifiers=(catalyzer,)),
                _reaction(species_c, species_b),
            ),
        )
        clingo_model, id_to_model_element = pd2af.solver.solve(
            cyclic_map, mode="normal"
        )
        edges = _named_influence_edges(clingo_model, id_to_model_element)
        # Consumption survives: A depletes its catalysed reactant B.
        assert ("negativelyInfluences", "A", "B") in edges
        # Base path survives: A activates the product C.
        assert ("positivelyInfluences", "A", "C") in edges
        # Cyclic leak removed: the A->C->B hop crosses a within-cycle edge.
        assert ("positivelyInfluences", "A", "B") not in edges

    def test_acyclic_production_keeps_transitive_influence(self):
        # A catalyses B->C, then C->D (no cycle); transitivity must reach D.
        species_a = _species("A")
        species_b = _species("B", active=True)
        species_c = _species("C", active=True)
        species_d = _species("D", active=True)
        catalyzer = momapy.celldesigner.Catalyzer(referred_element=species_a)
        acyclic_map = _map_from(
            (species_a, species_b, species_c, species_d),
            (
                _reaction(species_b, species_c, modifiers=(catalyzer,)),
                _reaction(species_c, species_d),
            ),
        )
        clingo_model, id_to_model_element = pd2af.solver.solve(
            acyclic_map, mode="normal"
        )
        edges = _named_influence_edges(clingo_model, id_to_model_element)
        # Transitivity across the non-cycle C->D edge is preserved.
        assert ("positivelyInfluences", "A", "D") in edges


class TestSetActive:
    """User-declared active elements: `set_active` injects
    `hasActivity(..., isInputParameter)` for a species that carries no
    structural activity signal of its own."""

    def test_set_active_adds_activity_for_non_active_species(
        self, example_cd_map, solved_keep_species
    ):
        # Species A (id `s1`) is a bare, non-active top-level species: it is
        # absent from the baseline activities (B, D, E, F, G).
        baseline_model, _ = solved_keep_species
        baseline_count = len(_activity_atoms(baseline_model))
        clingo_model, id_to_model_element = pd2af.solver.solve(
            example_cd_map, mode="keep-species", set_active=["s1"]
        )
        atoms = _activity_atoms(clingo_model)
        names = {id_to_model_element[atom.key.species].name for atom in atoms}
        assert "A" in names
        assert len(atoms) == baseline_count + 1

    def test_unknown_set_active_id_raises(self, example_cd_map):
        with pytest.raises(ValueError):
            pd2af.solver.solve(
                example_cd_map, mode="keep-species", set_active=["not-an-id"]
            )


class TestSetInactive:
    """User-declared inactive elements: `set_inactive` suppresses *any*
    `hasActivity` for a species that would otherwise be active by default."""

    def test_set_inactive_suppresses_default_active_species(
        self, example_cd_map, solved_keep_species
    ):
        # Species B (id `s2`) is active by default (it is one of the baseline
        # activities B, D, E, F, G). Marking it inactive drops it entirely.
        baseline_model, _ = solved_keep_species
        baseline_count = len(_activity_atoms(baseline_model))
        clingo_model, id_to_model_element = pd2af.solver.solve(
            example_cd_map, mode="keep-species", set_inactive=["s2"]
        )
        atoms = _activity_atoms(clingo_model)
        names = {id_to_model_element[atom.key.species].name for atom in atoms}
        assert "B" not in names
        assert len(atoms) == baseline_count - 1

    def test_set_inactive_conflicting_with_set_active_raises(self, example_cd_map):
        with pytest.raises(ValueError):
            pd2af.solver.solve(
                example_cd_map,
                mode="keep-species",
                set_active=["s1"],
                set_inactive=["s1"],
            )

    def test_unknown_set_inactive_id_raises(self, example_cd_map):
        with pytest.raises(ValueError):
            pd2af.solver.solve(
                example_cd_map,
                mode="keep-species",
                set_inactive=["not-an-id"],
            )


class TestSetAllActive:
    """`set_all_active` turns every top-level species into an activity, while
    per-id `set_inactive` still carves out exceptions and both global toggles
    are mutually exclusive."""

    def test_set_all_active_activates_every_top_level_species(self, example_cd_map):
        clingo_model, id_to_model_element = pd2af.solver.solve(
            example_cd_map, mode="keep-species", set_all_active=True
        )
        atoms = _activity_atoms(clingo_model)
        names = {id_to_model_element[atom.key.species].name for atom in atoms}
        # A is bare (never active by default); the toggle surfaces it.
        assert "A" in names

    def test_set_inactive_overrides_set_all_active(self, example_cd_map):
        # Per-id > global: B is globally activated but explicitly suppressed.
        clingo_model, id_to_model_element = pd2af.solver.solve(
            example_cd_map,
            mode="keep-species",
            set_all_active=True,
            set_inactive=["s2"],
        )
        atoms = _activity_atoms(clingo_model)
        names = {id_to_model_element[atom.key.species].name for atom in atoms}
        assert "B" not in names

    def test_both_global_toggles_raises(self, example_cd_map):
        with pytest.raises(ValueError):
            pd2af.solver.solve(
                example_cd_map,
                mode="keep-species",
                set_all_active=True,
                set_all_inactive=True,
            )


class TestSetAllInactive:
    """`set_all_inactive` suppresses every activity, while per-id `set_active`
    still forces the named ids back on."""

    def test_set_all_inactive_suppresses_every_activity(self, example_cd_map):
        clingo_model, _ = pd2af.solver.solve(
            example_cd_map, mode="keep-species", set_all_inactive=True
        )
        assert not _activity_atoms(clingo_model)

    def test_set_active_overrides_set_all_inactive(self, example_cd_map):
        # Per-id > global: everything is suppressed except the forced id.
        clingo_model, id_to_model_element = pd2af.solver.solve(
            example_cd_map,
            mode="keep-species",
            set_all_inactive=True,
            set_active=["s1"],
        )
        atoms = _activity_atoms(clingo_model)
        names = {id_to_model_element[atom.key.species].name for atom in atoms}
        assert names == {"A"}
