import pytest

import momapy.celldesigner

import pd2af

from tests._helpers import (
    has_dot_binary,
    modulation_set,
    read_cd_map,
    species_names,
)


@pytest.fixture(scope="module")
def out_default(example_cd_map):
    return pd2af.transform(example_cd_map, mode="normal", layout_mode="plain")


@pytest.fixture(scope="module")
def out_no_complex(example_cd_map):
    return pd2af.transform(
        example_cd_map, mode="no-complex", layout_mode="plain"
    )


@pytest.fixture(scope="module")
def out_pure_af(example_cd_map):
    if not has_dot_binary():
        pytest.skip("graphviz `dot` binary not on PATH")
    return pd2af.transform(
        example_cd_map, mode="pure-af", layout_mode="auto"
    )


class TestTransformExampleDefaultMode:
    """Golden test: example.xml under default mode and plain layout."""

    def test_returns_celldesigner_map(self, out_default):
        assert isinstance(out_default, momapy.celldesigner.CellDesignerMap)

    def test_preserves_celldesigner_model_type(self, out_default):
        assert isinstance(out_default.model, momapy.celldesigner.CellDesignerModel)

    def test_expected_active_species(self, out_default):
        assert species_names(out_default.model) == ["B", "C", "D", "E", "F", "G"]

    def test_expected_modulations(self, out_default):
        assert modulation_set(out_default.model) == {
            ("PositiveInfluence", "B", "D"),
            ("PositiveInfluence", "D", "F"),
            ("Inhibition", "B", "E"),
            ("Inhibition", "G", "B"),
        }

    def test_layout_present_with_plain_mode(self, out_default):
        assert out_default.layout is not None
        # one layout element per species + one arc per modulation
        assert len(out_default.layout.layout_elements) == (
            len(out_default.model.species) + len(out_default.model.modulations)
        )


class TestTransformExampleNoComplexMode:
    """Golden test: example.xml under no-complex mode."""

    def test_complex_subunits_replace_complex(self, out_no_complex):
        # Complex D contains subunits A and C; under no-complex, D drops out
        # and influences route through its active subunits instead. Of A and
        # C only C has activity (active structural state), so D->F becomes
        # C->F and B->D becomes B->C.
        names = species_names(out_no_complex.model)
        assert "D" not in names
        assert "C" in names

    def test_expected_modulations(self, out_no_complex):
        assert modulation_set(out_no_complex.model) == {
            ("PositiveInfluence", "B", "C"),
            ("PositiveInfluence", "C", "F"),
            ("Inhibition", "B", "E"),
            ("Inhibition", "G", "B"),
        }


class TestTransformExamplePureAfMode:
    """Golden test: example.xml under pure-af mode."""

    def test_returns_celldesigner_map(self, out_pure_af):
        assert isinstance(out_pure_af, momapy.celldesigner.CellDesignerMap)

    def test_complex_drops_out_via_no_complex_inheritance(self, out_pure_af):
        # Complex D has an active subunit; pure-af inherits the no-complex
        # behaviour and drops D in favour of the active subunit C.
        names = species_names(out_pure_af.model)
        assert "D" not in names

    def test_synthetic_species_have_clean_template(self, out_pure_af):
        for species in out_pure_af.model.species:
            template = getattr(species, "template", None)
            if template is None:
                continue
            # synthesized species' template id is prefixed; original
            # template ids in the input map are short ("pr1", ...).
            if not template.id_.startswith("pure_af_template__"):
                continue
            if hasattr(template, "modification_residues"):
                assert template.modification_residues == frozenset()
            if hasattr(template, "regions"):
                assert template.regions == frozenset()

    def test_synthetic_species_have_no_states_or_modifications(self, out_pure_af):
        for species in out_pure_af.model.species:
            if not species.id_.startswith("pure_af__"):
                continue
            assert species.homomultimer == 1
            structural_states = getattr(species, "structural_states", None)
            if structural_states is not None:
                assert structural_states == frozenset()
            modifications = getattr(species, "modifications", None)
            if modifications is not None:
                assert modifications == frozenset()


class TestTransformLayoutModes:
    def test_no_layout_mode_yields_map_without_layout(self, example_cd_map):
        out = pd2af.transform(example_cd_map, mode="normal", layout_mode=None)
        assert isinstance(out, momapy.celldesigner.CellDesignerMap)
        assert out.layout is None
        assert len(out.model.species) > 0

    def test_plain_mode_attaches_layout(self, example_cd_map):
        out = pd2af.transform(example_cd_map, mode="normal", layout_mode="plain")
        assert out.layout is not None
        assert len(out.layout.layout_elements) > 0

    @pytest.mark.skipif(
        not has_dot_binary(), reason="graphviz `dot` binary not on PATH"
    )
    def test_auto_mode_runs_with_dot(self, example_cd_map):
        out = pd2af.transform(example_cd_map, mode="normal", layout_mode="auto")
        assert out.layout is not None
        assert len(out.layout.layout_elements) > 0


class TestTransformErrors:
    def test_unknown_mode_raises_value_error(self, example_cd_map):
        with pytest.raises(ValueError):
            pd2af.transform(example_cd_map, mode="not-a-mode", layout_mode="plain")

    def test_unknown_layout_mode_raises_value_error(self, example_cd_map):
        with pytest.raises(ValueError):
            pd2af.transform(example_cd_map, mode="normal", layout_mode="not-a-layout")

    @pytest.mark.parametrize("layout_mode", ["plain", "overlay", None])
    def test_pure_af_requires_auto_layout(self, example_cd_map, layout_mode):
        with pytest.raises(ValueError):
            pd2af.transform(
                example_cd_map, mode="pure-af", layout_mode=layout_mode
            )


class TestTransformIsPure:
    """transform should not mutate the input map's model."""

    def test_original_model_unchanged(self, example_map_path):
        cd_map = read_cd_map(example_map_path)
        before_species = sorted(s.id_ for s in cd_map.model.species)
        before_reactions = sorted(r.id_ for r in cd_map.model.reactions)
        pd2af.transform(cd_map, mode="normal", layout_mode=None)
        after_species = sorted(s.id_ for s in cd_map.model.species)
        after_reactions = sorted(r.id_ for r in cd_map.model.reactions)
        assert before_species == after_species
        assert before_reactions == after_reactions

    def test_pure_af_does_not_mutate_template(self, example_map_path):
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
        pd2af.transform(cd_map, mode="pure-af", layout_mode="auto")
        after_template_ids = sorted(
            t.id_ for t in cd_map.model.species_templates
        )
        after_template_residues = {
            t.id_: getattr(t, "modification_residues", None)
            for t in cd_map.model.species_templates
        }
        assert before_template_ids == after_template_ids
        assert before_template_residues == after_template_residues
