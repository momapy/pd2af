"""SBGN-PD -> SBGN-AF transform: `normal` and `normal-no-complex` merged modes.

Two layers of coverage:

* Programmatic fixtures (built with momapy builders) pin the merged-mode
  semantics exactly: proteoform collapse, complex handling, and the
  subunit-promotion parity between `normal` and `normal-no-complex`.
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
    therefore suppressed in normal-no-complex (subunit promoted) and kept in normal
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


@pytest.fixture
def stateful_complex_map():
    """A complex carrying its *own* state variable (``tense``) plus a stateful
    subunit -- mirrors the actin:myosin case. Merged modes strip both the
    complex's and the subunit's state (recursively); keep-species keeps them."""
    complex_state = momapy.sbgn.pd.StateVariable(
        variable="r0", value="tense", order=0
    )
    subunit_state = momapy.sbgn.pd.StateVariable(
        variable="r0", value="active", order=0
    )
    subunit = momapy.sbgn.pd.MacromoleculeSubunit(
        label="RAF", state_variables=frozenset([subunit_state])
    )
    complex_ = momapy.sbgn.pd.Complex(
        label=None,
        state_variables=frozenset([complex_state]),
        subunits=frozenset([subunit]),
    )
    model = momapy.sbgn.pd.SBGNPDModel(entity_pools=frozenset([complex_]))
    return momapy.sbgn.pd.SBGNPDMap(model=model)


class TestProteoformMerging:
    def test_keep_species_keeps_proteoforms_distinct(self, proteoform_map):
        out = pd2af.transform(
            proteoform_map, mode="keep-species", layout_mode=None
        ).obj
        assert len(out.model.activities) == 2

    @pytest.mark.parametrize("mode", ("normal", "normal-no-complex"))
    def test_merged_modes_collapse_proteoforms(self, proteoform_map, mode):
        out = pd2af.transform(proteoform_map, mode=mode, layout_mode=None).obj
        assert _activity_labels(out.model) == ["AKT"]


class TestProvenance:
    """The TransformerResult.provenance maps each input element to the output
    AF elements derived from it; `.inverse` recovers the input elements behind
    a given output element, collapsing the many-to-one merge."""

    def test_merged_proteoforms_trace_back_to_both_inputs(self, proteoform_map):
        result = pd2af.transform(proteoform_map, mode="normal", layout_mode=None)
        # normal mode merges the two proteoforms into a single AKT activity.
        merged_activity = next(iter(result.obj.model.activities))
        input_proteoforms = frozenset(proteoform_map.model.entity_pools)
        # Forward: each input proteoform traces to the one merged activity.
        for proteoform in input_proteoforms:
            assert result.provenance[proteoform] == frozenset([merged_activity])
        # Inverse (many-to-one): the merged activity traces back to both inputs.
        assert (
            result.provenance.inverse[id(merged_activity)] == input_proteoforms
        )

    def test_keep_species_keeps_provenance_one_to_one(self, proteoform_map):
        result = pd2af.transform(
            proteoform_map, mode="keep-species", layout_mode=None
        )
        # keep-species keeps the proteoforms distinct: each input maps to its
        # own output activity, and every provenance value is in the model.
        model_activities = set(result.obj.model.activities)
        for output_elements in result.provenance.values():
            assert len(output_elements) == 1
            assert output_elements <= model_activities


class TestComplexHandling:
    def test_normal_keeps_complex_without_promoting_subunit(
        self, active_subunit_complex_map
    ):
        out = pd2af.transform(
            active_subunit_complex_map, mode="normal", layout_mode=None
        ).obj
        # The complex becomes a single ComplexUnitOfInformation activity. Its
        # active subunit is NOT promoted to its own activity -- it is a
        # structural component of the complex (carried in the composed label).
        assert _has_complex_unit_of_information(out.model)
        assert not any(
            isinstance(unit, momapy.sbgn.af.MacromoleculeUnitOfInformation)
            for activity in out.model.activities
            for unit in activity.units_of_information
        )
        assert len(out.model.activities) == 1

    def test_normal_no_complex_breaks_complex_into_subunit(
        self, active_subunit_complex_map
    ):
        out = pd2af.transform(
            active_subunit_complex_map, mode="normal-no-complex", layout_mode=None
        ).obj
        # The suppressed complex contributes no ComplexUnitOfInformation; only
        # the promoted subunit survives.
        assert not _has_complex_unit_of_information(out.model)
        assert _activity_labels(out.model) == ["RAF"]

    @pytest.mark.parametrize("mode", ("normal", "normal-no-complex"))
    def test_merged_modes_strip_complex_and_subunit_state(
        self, stateful_complex_map, mode
    ):
        out = pd2af.transform(
            stateful_complex_map, mode=mode, layout_mode=None
        ).obj
        labels = _activity_labels(out.model)
        # No state-variable bracket survives -- neither the complex's own
        # `tense` nor the subunit's `active`.
        assert all("tense" not in label for label in labels)
        assert all("active" not in label for label in labels)

    def test_keep_species_retains_complex_state(self, stateful_complex_map):
        out = pd2af.transform(
            stateful_complex_map, mode="keep-species", layout_mode=None
        ).obj
        labels = _activity_labels(out.model)
        assert any("tense" in label for label in labels)


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

    @pytest.mark.parametrize("mode", ("normal", "normal-no-complex"))
    def test_merged_mode_emits_influences(self, sbgn_map, mode):
        out = pd2af.transform(sbgn_map, mode=mode, layout_mode=None).obj
        # mapk_cascade and the others all carry modulation arcs, so a correct
        # carrier must yield at least one influence (was zero before the fix).
        assert len(out.model.influences) > 0

    @pytest.mark.parametrize("mode", ("normal", "normal-no-complex"))
    def test_merged_mode_auto_layout_round_trips(self, sbgn_map, mode):
        if not has_dot_binary():
            pytest.skip("graphviz `dot` binary not on PATH")
        out = pd2af.transform(sbgn_map, mode=mode, layout_mode="auto").obj
        assert out.layout is not None
        path = os.path.join(tempfile.gettempdir(), "pd2af_test_sbgn.sbgn")
        momapy.io.core.write(out, path, writer="sbgnml")
        # Read-back must not raise (the integration invariant).
        momapy.io.core.read(path, reader="sbgnml")


class TestSbgnPhenotypeActivity:
    """An SBGN-PD phenotype is a `Process`, not an `EntityPool`. It must still
    surface as an activity in every mode: the `keep-species`/`normal` key comes
    from `resolvesToTopLevel`, which needs a phenotype self-rule (the entity-pool
    self-rule cannot key a process); the `*-no-complex` modes already key it via
    `not isSubunit`/`not delete`."""

    @pytest.fixture(scope="class")
    def phenotype_map(self):
        # insulin-like_growth_factor_signaling carries one phenotype,
        # labelled "gene\ntranscription".
        return read_sbgn_map(
            os.path.join(
                SBGN_MAPS_DIR, "insulin-like_growth_factor_signaling.sbgn"
            )
        )

    @pytest.mark.parametrize(
        "mode",
        ("keep-species", "normal", "keep-species-no-complex", "normal-no-complex"),
    )
    def test_phenotype_is_an_activity(self, phenotype_map, mode):
        out = pd2af.transform(phenotype_map, mode=mode, layout_mode=None).obj
        assert any(
            label and "transcription" in label
            for label in _activity_labels(out.model)
        )


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
        ("keep-species", "keep-species-no-complex", "normal", "normal-no-complex"),
    )
    def test_model_carries_compartments(self, map_with_compartments, mode):
        out = pd2af.transform(
            map_with_compartments, mode=mode, layout_mode=None
        ).obj
        assert len(out.model.compartments) == 1
        assert all(
            activity.compartment is not None for activity in out.model.activities
        )

    def test_plain_layout_renders_compartment(self, map_with_compartments):
        # Regression: the plain-mode input-compartment lookup used to return
        # None, so zero CompartmentLayouts were emitted.
        out = pd2af.transform(
            map_with_compartments, mode="keep-species", layout_mode="plain"
        ).obj
        assert len(_compartment_layouts(out.layout)) == 1

    @pytest.mark.parametrize("mode", ("keep-species", "normal"))
    def test_auto_layout_renders_compartment(self, map_with_compartments, mode):
        # Regression: auto-layout used to crash on `compartment.outside`
        # (a CellDesigner-only relation; SBGN has no outside compartment).
        if not has_dot_binary():
            pytest.skip("graphviz `dot` binary not on PATH")
        out = pd2af.transform(
            map_with_compartments, mode=mode, layout_mode="auto"
        ).obj
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
        ).obj
        path = os.path.join(tempfile.gettempdir(), "pd2af_test_compartments.sbgn")
        momapy.io.core.write(out, path, writer="sbgnml")
        back = momapy.io.core.read(path, reader="sbgnml").obj
        assert len(back.model.compartments) == len(out.model.compartments) == 1


class TestTransformModelInput:
    """A bare SBGN-PD model in `normal` mode returns a bare SBGN-AF model with
    `layout_mode` forced to None -- so no graphviz `dot` is required."""

    def test_returns_sbgn_af_model(self, proteoform_map):
        result = pd2af.transform(proteoform_map.model, mode="normal")
        assert isinstance(result.obj, momapy.sbgn.af.SBGNAFModel)

    def test_model_output_matches_map_output(self, proteoform_map):
        from_model = pd2af.transform(proteoform_map.model, mode="normal").obj
        from_map = pd2af.transform(
            proteoform_map, mode="normal", layout_mode=None
        ).obj
        assert _activity_labels(from_model) == _activity_labels(from_map.model)

    @pytest.mark.parametrize("layout_mode", ["dot", "plain"])
    def test_rejects_explicit_layout_mode(self, proteoform_map, layout_mode):
        with pytest.raises(ValueError):
            pd2af.transform(
                proteoform_map.model, mode="normal", layout_mode=layout_mode
            )
