import os

import pytest

import momapy.celldesigner

import pd2af

from tests._helpers import (
    MAPS_DIR,
    has_dot_binary,
    modulation_set,
    read_cd_map,
    species_names,
)


def _assert_recursively_stripped(species_iterable):
    """Every species (and subunit, recursively) carries no PTM decoration and
    keeps a real, reader-resolvable id (no synthesized prefix)."""
    for species in species_iterable:
        assert not species.id_.startswith("new_species_from_template")
        assert getattr(species, "homomultimer", 1) == 1
        assert not getattr(species, "structural_states", frozenset())
        assert not getattr(species, "modifications", frozenset())
        template = getattr(species, "template", None)
        if template is not None:
            assert not getattr(template, "modification_residues", frozenset())
            assert not getattr(template, "regions", frozenset())
        _assert_recursively_stripped(getattr(species, "subunits", ()) or ())


@pytest.fixture(scope="module")
def out_keep_species(example_cd_map):
    return pd2af.transform(
        example_cd_map, mode="keep-species", layout_mode="plain"
    )


@pytest.fixture(scope="module")
def out_keep_species_no_complex(example_cd_map):
    return pd2af.transform(
        example_cd_map, mode="keep-species-no-complex", layout_mode="plain"
    )


@pytest.fixture(scope="module")
def out_no_complex(example_cd_map):
    if not has_dot_binary():
        pytest.skip("graphviz `dot` binary not on PATH")
    return pd2af.transform(
        example_cd_map, mode="no-complex", layout_mode="auto"
    )


@pytest.fixture(scope="module")
def out_normal(example_cd_map):
    if not has_dot_binary():
        pytest.skip("graphviz `dot` binary not on PATH")
    return pd2af.transform(example_cd_map, mode="normal", layout_mode="auto")


class TestTransformExampleKeepSpeciesMode:
    """Golden test: example.xml under keep-species mode and plain layout."""

    def test_returns_celldesigner_map(self, out_keep_species):
        assert isinstance(out_keep_species, momapy.celldesigner.CellDesignerMap)

    def test_preserves_celldesigner_model_type(self, out_keep_species):
        assert isinstance(
            out_keep_species.model, momapy.celldesigner.CellDesignerModel
        )

    def test_expected_active_species(self, out_keep_species):
        # Active subunit C is subsumed into its containing complex D under
        # keep-species (which keeps complexes), so it does not appear as a
        # top-level activity.
        assert species_names(out_keep_species.model) == [
            "B", "D", "E", "F", "G",
        ]

    def test_expected_modulations(self, out_keep_species):
        assert modulation_set(out_keep_species.model) == {
            ("PositiveInfluence", "B", "D"),
            ("PositiveInfluence", "D", "F"),
            ("NegativeInfluence", "B", "E"),
            ("NegativeInfluence", "G", "B"),
        }

    def test_layout_present_with_plain_mode(self, out_keep_species):
        assert out_keep_species.layout is not None
        assert len(out_keep_species.layout.layout_elements) > 0


class TestTransformExampleKeepSpeciesNoComplexMode:
    """Golden test: example.xml under keep-species-no-complex mode."""

    def test_complex_subunits_replace_complex(self, out_keep_species_no_complex):
        # Complex D contains subunits A and C; under keep-species-no-complex,
        # D drops out and influences route through its active subunits
        # instead. Only C has activity, so D->F becomes C->F and B->D becomes
        # B->C.
        names = species_names(out_keep_species_no_complex.model)
        assert "D" not in names
        assert "C" in names

    def test_expected_modulations(self, out_keep_species_no_complex):
        assert modulation_set(out_keep_species_no_complex.model) == {
            ("PositiveInfluence", "B", "C"),
            ("PositiveInfluence", "C", "F"),
            ("NegativeInfluence", "B", "E"),
            ("NegativeInfluence", "G", "B"),
        }


class TestTransformExampleNoComplexMode:
    """Golden test: example.xml under no-complex mode (merge proteoforms,
    drop complexes)."""

    def test_returns_celldesigner_map(self, out_no_complex):
        assert isinstance(out_no_complex, momapy.celldesigner.CellDesignerMap)

    def test_complex_drops_out(self, out_no_complex):
        # Complex D has an active subunit; no-complex drops D in favour of
        # the active subunit C.
        names = species_names(out_no_complex.model)
        assert "D" not in names

    def test_synthetic_species_have_clean_template(self, out_no_complex):
        for species in out_no_complex.model.species:
            template = getattr(species, "template", None)
            if template is None:
                continue
            if not template.id_.startswith("merged_template__"):
                continue
            if hasattr(template, "modification_residues"):
                assert template.modification_residues == frozenset()
            if hasattr(template, "regions"):
                assert template.regions == frozenset()

    def test_synthetic_species_have_no_states_or_modifications(
        self, out_no_complex
    ):
        for species in out_no_complex.model.species:
            if not species.id_.startswith("merged__"):
                continue
            assert species.homomultimer == 1
            structural_states = getattr(species, "structural_states", None)
            if structural_states is not None:
                assert structural_states == frozenset()
            modifications = getattr(species, "modifications", None)
            if modifications is not None:
                assert modifications == frozenset()


class TestTransformExampleNormalMode:
    """Golden test: example.xml under normal mode (merge proteoforms, keep
    complexes)."""

    def test_returns_celldesigner_map(self, out_normal):
        assert isinstance(out_normal, momapy.celldesigner.CellDesignerMap)

    def test_complex_is_kept(self, out_normal):
        # Complex D is kept as its own activity. Its active subunit C is no
        # longer a separate top-level activity -- it is a structural component
        # of D, so any influence it carries routes to D (same as keep-species).
        names = species_names(out_normal.model)
        assert "D" in names
        assert "C" not in names

    def test_complex_routes_through_itself(self, out_normal):
        # Influences involving the complex go through the complex (kept_species)
        # rather than through subunits.
        mods = modulation_set(out_normal.model)
        assert ("PositiveInfluence", "B", "D") in mods
        assert ("PositiveInfluence", "D", "F") in mods

    def test_complex_kept_with_kept_species_id(self, out_normal):
        # The complex survives with its original id (not a synthesized one).
        complex_species = [
            s
            for s in out_normal.model.species
            if isinstance(s, momapy.celldesigner.Complex)
        ]
        assert complex_species
        for s in complex_species:
            assert not s.id_.startswith("new_species_from_template")

    def test_merged_species_carry_no_decorations(self, out_normal):
        # Merged modes strip every PTM decoration (recursively, incl. subunits)
        # and keep each species' real id -- there is no synthesized id prefix.
        _assert_recursively_stripped(out_normal.model.species)


class TestMergedModeStripping:
    """The merged modes strip all PTM decorations (recursively); keep-species
    keeps them. Exercised on a committed map that carries real modifications and
    structural states -- example.xml is decoration-free, so it cannot prove the
    strip on its own."""

    @pytest.fixture(scope="class")
    def rich_map(self):
        return read_cd_map(os.path.join(MAPS_DIR, "Apoptosis_pathway.xml"))

    def test_normal_strips_all_decorations(self, rich_map):
        out = pd2af.transform(rich_map, mode="normal", layout_mode=None)
        _assert_recursively_stripped(out.model.species)

    def test_no_complex_strips_all_decorations(self, rich_map):
        out = pd2af.transform(rich_map, mode="no-complex", layout_mode=None)
        _assert_recursively_stripped(out.model.species)

    def test_keep_species_retains_decorations(self, rich_map):
        out = pd2af.transform(rich_map, mode="keep-species", layout_mode=None)
        decorated = [
            s
            for s in out.model.species
            if getattr(s, "modifications", None)
            or getattr(s, "structural_states", None)
        ]
        assert decorated  # keep-species preserves PTM decorations


class TestTransformLayoutModes:
    def test_no_layout_mode_yields_map_without_layout(self, example_cd_map):
        out = pd2af.transform(
            example_cd_map, mode="keep-species", layout_mode=None
        )
        assert isinstance(out, momapy.celldesigner.CellDesignerMap)
        assert out.layout is None
        assert len(out.model.species) > 0

    def test_plain_mode_attaches_layout(self, example_cd_map):
        out = pd2af.transform(
            example_cd_map, mode="keep-species", layout_mode="plain"
        )
        assert out.layout is not None
        assert len(out.layout.layout_elements) > 0

    @pytest.mark.skipif(
        not has_dot_binary(), reason="graphviz `dot` binary not on PATH"
    )
    def test_auto_mode_runs_with_dot(self, example_cd_map):
        out = pd2af.transform(
            example_cd_map, mode="keep-species", layout_mode="auto"
        )
        assert out.layout is not None
        assert len(out.layout.layout_elements) > 0


class TestTransformErrors:
    def test_unknown_mode_raises_value_error(self, example_cd_map):
        with pytest.raises(ValueError):
            pd2af.transform(
                example_cd_map, mode="not-a-mode", layout_mode="plain"
            )

    @pytest.mark.parametrize("mode", ["normal", "no-complex"])
    @pytest.mark.parametrize("layout_mode", ["plain", "overlay"])
    def test_merged_modes_reject_input_derived_layout(
        self, example_cd_map, mode, layout_mode
    ):
        with pytest.raises(ValueError):
            pd2af.transform(example_cd_map, mode=mode, layout_mode=layout_mode)

    @pytest.mark.parametrize(
        "mode", ["keep-species", "keep-species-no-complex"]
    )
    @pytest.mark.parametrize("layout_mode", ["plain", "overlay"])
    def test_keep_species_modes_accept_non_auto_layout(
        self, example_cd_map, mode, layout_mode
    ):
        out = pd2af.transform(
            example_cd_map, mode=mode, layout_mode=layout_mode
        )
        assert out.layout is not None


class TestTransformIsPure:
    """transform should not mutate the input map's model."""

    def test_original_model_unchanged(self, example_map_path):
        cd_map = read_cd_map(example_map_path)
        before_species = sorted(s.id_ for s in cd_map.model.species)
        before_reactions = sorted(r.id_ for r in cd_map.model.reactions)
        pd2af.transform(cd_map, mode="keep-species", layout_mode=None)
        after_species = sorted(s.id_ for s in cd_map.model.species)
        after_reactions = sorted(r.id_ for r in cd_map.model.reactions)
        assert before_species == after_species
        assert before_reactions == after_reactions

    def test_no_complex_does_not_mutate_template(self, example_map_path):
        if not has_dot_binary():
            pytest.skip("graphviz `dot` binary not on PATH")
        cd_map = read_cd_map(example_map_path)
        before_template_ids = sorted(
            t.id_ for t in cd_map.model.species_templates
        )
        before_template_residues = {
            t.id_: getattr(t, "modification_residues", None)
            for t in cd_map.model.species_templates
        }
        pd2af.transform(cd_map, mode="no-complex", layout_mode="auto")
        after_template_ids = sorted(
            t.id_ for t in cd_map.model.species_templates
        )
        after_template_residues = {
            t.id_: getattr(t, "modification_residues", None)
            for t in cd_map.model.species_templates
        }
        assert before_template_ids == after_template_ids
        assert before_template_residues == after_template_residues
