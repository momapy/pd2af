import pytest

import momapy.celldesigner

import pd2af.layouts
import pd2af.solver

from tests._helpers import has_dot_binary


def _solve_and_make_model(cd_map, mode="normal"):
    clingo_model, id_to = pd2af.solver.solve(cd_map, mode=mode)
    return pd2af.solver.make_new_cd_model(clingo_model, id_to)


@pytest.fixture(scope="module")
def new_model(example_cd_map):
    return _solve_and_make_model(example_cd_map)


class TestBuildMapDispatch:
    def test_none_mode_returns_no_layout_map(self, example_cd_map, new_model):
        out = pd2af.layouts.build_map(example_cd_map, new_model, None)
        assert isinstance(out, momapy.celldesigner.CellDesignerMap)
        assert out.model is new_model
        assert out.layout is None

    def test_unsupported_mode_raises(self, example_cd_map, new_model):
        with pytest.raises(ValueError):
            pd2af.layouts.build_map(example_cd_map, new_model, "not-a-mode")


class TestMakePlain:
    @pytest.fixture(scope="class")
    def out(self, example_cd_map, new_model):
        return pd2af.layouts.make_plain(example_cd_map, new_model)

    def test_returns_celldesigner_map(self, out):
        assert isinstance(out, momapy.celldesigner.CellDesignerMap)

    def test_layout_present(self, out):
        assert out.layout is not None

    def test_layout_has_species_and_arc_elements(self, out, new_model):
        # Plain mode emits one node per kept species and one arc per modulation.
        expected = len(new_model.species) + len(new_model.modulations)
        assert len(out.layout.layout_elements) == expected

    def test_mapping_resolves_each_species(self, out, new_model):
        for species in new_model.species:
            mapped = out.layout_model_mapping.get_mapping(species)
            assert mapped is not None, f"species {species.name!r} not mapped"


class TestMakeNoLayout:
    def test_returns_map_without_layout(self, new_model):
        out = pd2af.layouts.make_no_layout(new_model)
        assert isinstance(out, momapy.celldesigner.CellDesignerMap)
        assert out.model is new_model
        assert out.layout is None


class TestMakeAuto:
    @pytest.mark.skipif(
        not has_dot_binary(), reason="graphviz `dot` binary not on PATH"
    )
    def test_auto_layout_produces_positioned_layout(
        self, example_cd_map, new_model
    ):
        out = pd2af.layouts.make_auto(example_cd_map, new_model)
        assert out.layout is not None
        assert len(out.layout.layout_elements) > 0
