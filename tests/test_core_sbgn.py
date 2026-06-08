"""SBGN-PD -> SBGN-AF transform: `normal` and `no-complex` merged modes.

Two layers of coverage:

* Programmatic fixtures (built with momapy builders) pin the merged-mode
  semantics exactly: proteoform collapse, complex handling, and the
  subunit-promotion parity between `normal` and `no-complex`.
* Real committed SBGN-PD maps exercise the full pipeline including
  graphviz `auto` layout and `.sbgn` read-back (the strongest regression
  guard for the carrier bug, which produced zero influences).
"""

import os
import tempfile

import pytest

import momapy.io.core
import momapy.sbgn.af
import momapy.sbgn.pd

import pd2af

from tests._helpers import (
    SBGN_MAPS_DIR,
    SBGN_WITH_COMPARTMENTS_MAP_PATH,
    has_dot_binary,
    read_sbgn_map,
)


def _activity_labels(model):
    return sorted(activity.label for activity in model.activities)


def _has_complex_unit_of_information(model):
    return any(
        isinstance(unit, momapy.sbgn.af.ComplexUnitOfInformation)
        for activity in model.activities
        for unit in getattr(activity, "units_of_information", frozenset())
    )


@pytest.fixture
def proteoform_map():
    """Two proteoforms of one macromolecule (same name, both active, different
    state) -- they must collapse to a single merged activity in the merged
    modes and stay distinct under keep-species."""
    state_active = momapy.sbgn.pd.StateVariable(
        variable="r0", value="active", order=0
    )
    state_phosphorylated = momapy.sbgn.pd.StateVariable(
        variable="r1", value="P", order=1
    )
    state_unphosphorylated = momapy.sbgn.pd.StateVariable(
        variable="r1", value=None, order=1
    )
    proteoform_one = momapy.sbgn.pd.Macromolecule(
        label="AKT",
        state_variables=frozenset([state_active, state_phosphorylated]),
    )
    proteoform_two = momapy.sbgn.pd.Macromolecule(
        label="AKT",
        state_variables=frozenset([state_active, state_unphosphorylated]),
    )
    model = momapy.sbgn.pd.SBGNPDModel(
        entity_pools=frozenset([proteoform_one, proteoform_two])
    )
    return momapy.sbgn.pd.SBGNPDMap(model=model)


@pytest.fixture
def active_subunit_complex_map():
    """A complex whose single macromolecule subunit is active. The complex is
    therefore suppressed in no-complex (subunit promoted) and kept in normal
    (complex activity + promoted subunit, in parity with CellDesigner)."""
    state_active = momapy.sbgn.pd.StateVariable(
        variable="r0", value="active", order=0
    )
    subunit = momapy.sbgn.pd.MacromoleculeSubunit(
        label="RAF", state_variables=frozenset([state_active])
    )
    complex_ = momapy.sbgn.pd.Complex(label=None, subunits=frozenset([subunit]))
    model = momapy.sbgn.pd.SBGNPDModel(entity_pools=frozenset([complex_]))
    return momapy.sbgn.pd.SBGNPDMap(model=model)


class TestProteoformMerging:
    def test_keep_species_keeps_proteoforms_distinct(self, proteoform_map):
        out = pd2af.transform(proteoform_map, mode="keep-species", layout_mode=None)
        assert len(out.model.activities) == 2

    @pytest.mark.parametrize("mode", ("normal", "no-complex"))
    def test_merged_modes_collapse_proteoforms(self, proteoform_map, mode):
        out = pd2af.transform(proteoform_map, mode=mode, layout_mode=None)
        assert _activity_labels(out.model) == ["AKT"]


class TestComplexHandling:
    def test_normal_keeps_complex_and_promotes_subunit(
        self, active_subunit_complex_map
    ):
        out = pd2af.transform(
            active_subunit_complex_map, mode="normal", layout_mode=None
        )
        # The complex becomes a ComplexUnitOfInformation activity, and its
        # active mergeable subunit is promoted to its own activity (parity).
        assert _has_complex_unit_of_information(out.model)
        assert any(
            isinstance(unit, momapy.sbgn.af.MacromoleculeUnitOfInformation)
            for activity in out.model.activities
            for unit in activity.units_of_information
        )

    def test_no_complex_breaks_complex_into_subunit(
        self, active_subunit_complex_map
    ):
        out = pd2af.transform(
            active_subunit_complex_map, mode="no-complex", layout_mode=None
        )
        # The suppressed complex contributes no ComplexUnitOfInformation; only
        # the promoted subunit survives.
        assert not _has_complex_unit_of_information(out.model)
        assert _activity_labels(out.model) == ["RAF"]


_SBGN_MAP_NAMES = (
    "mapk_cascade",
    "insulin-like_growth_factor_signaling",
    "neuronal_muscle_signalling",
)


class TestRealMapIntegration:
    """Full pipeline over real maps: influences emitted (carrier-bug guard),
    auto layout built, and `.sbgn` round-trips through the sbgnml reader."""

    @pytest.fixture(scope="class", params=_SBGN_MAP_NAMES)
    def sbgn_map(self, request):
        return read_sbgn_map(
            os.path.join(SBGN_MAPS_DIR, f"{request.param}.sbgn")
        )

    @pytest.mark.parametrize("mode", ("normal", "no-complex"))
    def test_merged_mode_emits_influences(self, sbgn_map, mode):
        out = pd2af.transform(sbgn_map, mode=mode, layout_mode=None)
        # mapk_cascade and the others all carry modulation arcs, so a correct
        # carrier must yield at least one influence (was zero before the fix).
        assert len(out.model.influences) > 0

    @pytest.mark.parametrize("mode", ("normal", "no-complex"))
    def test_merged_mode_auto_layout_round_trips(self, sbgn_map, mode):
        if not has_dot_binary():
            pytest.skip("graphviz `dot` binary not on PATH")
        out = pd2af.transform(sbgn_map, mode=mode, layout_mode="auto")
        assert out.layout is not None
        path = os.path.join(tempfile.gettempdir(), "pd2af_test_sbgn.sbgn")
        momapy.io.core.write(out, path, writer="sbgnml")
        # Read-back must not raise (the integration invariant).
        momapy.io.core.read(path, reader="sbgnml")


def _compartment_layouts(layout):
    return [
        element
        for element in layout.layout_elements
        if isinstance(element, momapy.sbgn.af.CompartmentLayout)
    ]


class TestCompartments:
    """Compartments must survive into the AF output (model + layout) when the
    input carries `compartmentRef`. The model pass already handled this; these
    guard the output/layout path: the plain-mode input-compartment lookup and
    the auto-layout `outside`-free clustering (SBGN has no outside compartment).
    """

    @pytest.fixture(scope="class")
    def map_with_compartments(self):
        return read_sbgn_map(SBGN_WITH_COMPARTMENTS_MAP_PATH)

    @pytest.mark.parametrize(
        "mode",
        ("keep-species", "keep-species-no-complex", "normal", "no-complex"),
    )
    def test_model_carries_compartments(self, map_with_compartments, mode):
        out = pd2af.transform(map_with_compartments, mode=mode, layout_mode=None)
        assert len(out.model.compartments) == 1
        assert all(
            activity.compartment is not None for activity in out.model.activities
        )

    def test_plain_layout_renders_compartment(self, map_with_compartments):
        # Regression: the plain-mode input-compartment lookup used to return
        # None, so zero CompartmentLayouts were emitted.
        out = pd2af.transform(
            map_with_compartments, mode="keep-species", layout_mode="plain"
        )
        assert len(_compartment_layouts(out.layout)) == 1

    @pytest.mark.parametrize("mode", ("keep-species", "normal"))
    def test_auto_layout_renders_compartment(self, map_with_compartments, mode):
        # Regression: auto-layout used to crash on `compartment.outside`
        # (a CellDesigner-only relation; SBGN has no outside compartment).
        if not has_dot_binary():
            pytest.skip("graphviz `dot` binary not on PATH")
        out = pd2af.transform(map_with_compartments, mode=mode, layout_mode="auto")
        assert len(_compartment_layouts(out.layout)) == 1

    @pytest.mark.parametrize(
        "mode,layout_mode",
        (("keep-species", "plain"), ("normal", "auto")),
    )
    def test_compartments_round_trip(
        self, map_with_compartments, mode, layout_mode
    ):
        if layout_mode == "auto" and not has_dot_binary():
            pytest.skip("graphviz `dot` binary not on PATH")
        out = pd2af.transform(
            map_with_compartments, mode=mode, layout_mode=layout_mode
        )
        path = os.path.join(tempfile.gettempdir(), "pd2af_test_compartments.sbgn")
        momapy.io.core.write(out, path, writer="sbgnml")
        back = momapy.io.core.read(path, reader="sbgnml").obj
        assert len(back.model.compartments) == len(out.model.compartments) == 1
