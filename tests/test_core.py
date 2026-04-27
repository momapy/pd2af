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
    return pd2af.transform(example_cd_map, mode="pd2af", layout_mode="plain")


@pytest.fixture(scope="module")
def out_no_complex(example_cd_map):
    return pd2af.transform(
        example_cd_map, mode="pd2af-no-complex", layout_mode="plain"
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
    """Golden test: example.xml under pd2af-no-complex mode."""

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


class TestTransformLayoutModes:
    def test_no_layout_mode_yields_map_without_layout(self, example_cd_map):
        out = pd2af.transform(example_cd_map, mode="pd2af", layout_mode=None)
        assert isinstance(out, momapy.celldesigner.CellDesignerMap)
        assert out.layout is None
        assert len(out.model.species) > 0

    def test_plain_mode_attaches_layout(self, example_cd_map):
        out = pd2af.transform(example_cd_map, mode="pd2af", layout_mode="plain")
        assert out.layout is not None
        assert len(out.layout.layout_elements) > 0

    @pytest.mark.skipif(
        not has_dot_binary(), reason="graphviz `dot` binary not on PATH"
    )
    def test_auto_mode_runs_with_dot(self, example_cd_map):
        out = pd2af.transform(example_cd_map, mode="pd2af", layout_mode="auto")
        assert out.layout is not None
        assert len(out.layout.layout_elements) > 0


class TestTransformErrors:
    def test_unknown_mode_raises_value_error(self, example_cd_map):
        with pytest.raises(ValueError):
            pd2af.transform(example_cd_map, mode="not-a-mode", layout_mode="plain")

    def test_unknown_layout_mode_raises_value_error(self, example_cd_map):
        with pytest.raises(ValueError):
            pd2af.transform(example_cd_map, mode="pd2af", layout_mode="not-a-layout")


class TestTransformIsPure:
    """transform should not mutate the input map's model."""

    def test_original_model_unchanged(self, example_map_path):
        cd_map = read_cd_map(example_map_path)
        before_species = sorted(s.id_ for s in cd_map.model.species)
        before_reactions = sorted(r.id_ for r in cd_map.model.reactions)
        pd2af.transform(cd_map, mode="pd2af", layout_mode=None)
        after_species = sorted(s.id_ for s in cd_map.model.species)
        after_reactions = sorted(r.id_ for r in cd_map.model.reactions)
        assert before_species == after_species
        assert before_reactions == after_reactions
