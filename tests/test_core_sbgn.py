"""SBGN-PD -> SBGN-AF transform: the `normal` and `no-complex` modes.

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

import momapy.core.layout
import momapy.io.core
import momapy.sbgn.af
import momapy.sbgn.pd
import pytest
from momapy.sbml.model import BQBiol, RDFAnnotation

import pd2af
from tests._helpers import (
    BINDING_ACTIVATION_MAPS_DIR,
    MAPS_DIR,
    SBGN_EXAMPLE_MAP_PATH,
    SBGN_MAPS_DIR,
    SBGN_WITH_COMPARTMENTS_MAP_PATH,
    has_dot_binary,
    read_cd_map,
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
    reading and stay distinct under `keep_species`."""
    state_active = momapy.sbgn.pd.StateVariable(variable="r0", value="active", order=0)
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
    state_active = momapy.sbgn.pd.StateVariable(variable="r0", value="active", order=0)
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
    complex's and the subunit's state (recursively); `keep_species` keeps them."""
    complex_state = momapy.sbgn.pd.StateVariable(variable="r0", value="tense", order=0)
    subunit_state = momapy.sbgn.pd.StateVariable(variable="r0", value="active", order=0)
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
            proteoform_map, mode="normal", keep_species=True, layout_mode=None
        ).obj
        assert len(out.model.activities) == 2

    @pytest.mark.parametrize("mode", ("normal", "no-complex"))
    def test_merged_modes_collapse_proteoforms(self, proteoform_map, mode):
        out = pd2af.transform(proteoform_map, mode=mode, layout_mode=None).obj
        assert _activity_labels(out.model) == ["AKT"]


class TestProvenance:
    """The TransformerResult.output_element_to_input_elements maps each output AF element to the
    input elements it derives from, collapsing the many-to-one merge; its
    `.inverse` recovers the output element behind a given input."""

    def test_merged_proteoforms_trace_back_to_both_inputs(self, proteoform_map):
        result = pd2af.transform(proteoform_map, mode="normal", layout_mode=None)
        # normal mode merges the two proteoforms into a single AKT activity.
        merged_activity = next(iter(result.obj.model.activities))
        input_proteoforms = frozenset(proteoform_map.model.entity_pools)
        # Forward (many-to-one): the merged activity traces back to both inputs.
        assert (
            result.output_element_to_input_elements[merged_activity]
            == input_proteoforms
        )
        # Inverse: each input proteoform recovers the one merged activity.
        for proteoform in input_proteoforms:
            assert result.output_element_to_input_elements.inverse[
                id(proteoform)
            ] == frozenset([merged_activity])

    def test_keep_species_keeps_provenance_one_to_one(self, proteoform_map):
        result = pd2af.transform(
            proteoform_map, mode="normal", keep_species=True, layout_mode=None
        )
        # `keep_species` keeps the forms distinct: each output activity has
        # a single input source, and every provenance key is an output element.
        model_activities = set(result.obj.model.activities)
        activity_sources = [
            input_elements
            for output_element, input_elements in result.output_element_to_input_elements.items()
            if output_element in model_activities
        ]
        assert activity_sources
        for input_elements in activity_sources:
            assert len(input_elements) == 1


class TestAnnotationCarry:
    """The merged modes union the annotations of every collapsed proteoform
    onto their single merged activity -- the case where content-match alone
    fails (the merged output's content differs from every input)."""

    def test_merged_mode_unions_proteoform_annotations(self, proteoform_map):
        proteoforms = list(proteoform_map.model.entity_pools)
        annotation_one = RDFAnnotation(
            qualifier=BQBiol.IS, resources=frozenset(["urn:one"])
        )
        annotation_two = RDFAnnotation(
            qualifier=BQBiol.IS, resources=frozenset(["urn:two"])
        )
        result = pd2af.transform(
            proteoform_map,
            mode="normal",
            layout_mode=None,
            element_to_annotations={
                proteoforms[0]: frozenset([annotation_one]),
                proteoforms[1]: frozenset([annotation_two]),
            },
            element_to_notes={},
        )
        merged_activity = next(iter(result.obj.model.activities))
        assert result.element_to_annotations[merged_activity] == frozenset(
            [annotation_one, annotation_two]
        )

    def test_map_level_annotation_survives_sbgnml_round_trip(self, tmp_path):
        reader_result = momapy.io.core.read(SBGN_EXAMPLE_MAP_PATH, reader="sbgnml")
        map_annotation = RDFAnnotation(
            qualifier=BQBiol.IS, resources=frozenset(["urn:map:level"])
        )
        element_to_annotations = dict(reader_result.element_to_annotations or {})
        element_to_annotations[reader_result.obj] = frozenset([map_annotation])
        result = pd2af.transform(
            reader_result.obj,
            mode="normal",
            keep_species=True,
            layout_mode="plain",
            element_to_annotations=element_to_annotations,
            element_to_notes=reader_result.element_to_notes,
        )
        # re-keyed onto the OUTPUT map object
        assert result.element_to_annotations.get(result.obj) == frozenset(
            [map_annotation]
        )
        path = str(tmp_path / "out.sbgn")
        momapy.io.core.write(
            result.obj,
            path,
            writer="sbgnml",
            element_to_annotations=result.element_to_annotations,
            element_to_notes=result.element_to_notes,
        )
        reread = momapy.io.core.read(path, reader="sbgnml")  # must not raise
        assert map_annotation in (reread.element_to_annotations or {}).get(
            reread.obj, frozenset()
        )


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
            active_subunit_complex_map, mode="no-complex", layout_mode=None
        ).obj
        # The suppressed complex contributes no ComplexUnitOfInformation; only
        # the promoted subunit survives.
        assert not _has_complex_unit_of_information(out.model)
        assert _activity_labels(out.model) == ["RAF"]

    @pytest.mark.parametrize("mode", ("normal", "no-complex"))
    def test_merged_modes_strip_complex_and_subunit_state(
        self, stateful_complex_map, mode
    ):
        out = pd2af.transform(stateful_complex_map, mode=mode, layout_mode=None).obj
        labels = _activity_labels(out.model)
        # No state-variable bracket survives -- neither the complex's own
        # `tense` nor the subunit's `active`.
        assert all("tense" not in label for label in labels)
        assert all("active" not in label for label in labels)

    def test_keep_species_retains_complex_state(self, stateful_complex_map):
        out = pd2af.transform(
            stateful_complex_map, mode="normal", keep_species=True, layout_mode=None
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
        return read_sbgn_map(os.path.join(SBGN_MAPS_DIR, f"{request.param}.sbgn"))

    @pytest.mark.parametrize("mode", ("normal", "no-complex"))
    def test_merged_mode_emits_influences(self, sbgn_map, mode):
        out = pd2af.transform(sbgn_map, mode=mode, layout_mode=None).obj
        # mapk_cascade and the others all carry modulation arcs, so a correct
        # carrier must yield at least one influence (was zero before the fix).
        assert len(out.model.influences) > 0

    @pytest.mark.parametrize("mode", ("normal", "no-complex"))
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
    surface as an activity in every mode: the complex-keeping key comes
    from `resolvesToTopLevel`, which needs a phenotype self-rule (the entity-pool
    self-rule cannot key a process); the `no-complex` mode already keys it via
    `not isSubunit`/`not delete`."""

    @pytest.fixture(scope="class")
    def phenotype_map(self):
        # insulin-like_growth_factor_signaling carries one phenotype,
        # labelled "gene\ntranscription".
        return read_sbgn_map(
            os.path.join(SBGN_MAPS_DIR, "insulin-like_growth_factor_signaling.sbgn")
        )

    @pytest.mark.parametrize("mode", ("normal", "no-complex"))
    @pytest.mark.parametrize("keep_species", (True, False))
    def test_phenotype_is_an_activity(self, phenotype_map, mode, keep_species):
        out = pd2af.transform(
            phenotype_map, mode=mode, keep_species=keep_species, layout_mode=None
        ).obj
        assert any(
            label and "transcription" in label for label in _activity_labels(out.model)
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

    @pytest.mark.parametrize("mode", ("normal", "no-complex"))
    @pytest.mark.parametrize("keep_species", (True, False))
    def test_model_carries_compartments(
        self, map_with_compartments, mode, keep_species
    ):
        out = pd2af.transform(
            map_with_compartments,
            mode=mode,
            keep_species=keep_species,
            layout_mode=None,
        ).obj
        assert len(out.model.compartments) == 1
        assert all(
            activity.compartment is not None for activity in out.model.activities
        )

    def test_plain_layout_renders_compartment(self, map_with_compartments):
        # Regression: the plain-mode input-compartment lookup used to return
        # None, so zero CompartmentLayouts were emitted.
        out = pd2af.transform(
            map_with_compartments, mode="normal", keep_species=True, layout_mode="plain"
        ).obj
        assert len(_compartment_layouts(out.layout)) == 1

    @pytest.mark.parametrize("keep_species", (True, False))
    def test_auto_layout_renders_compartment(self, map_with_compartments, keep_species):
        # Regression: auto-layout used to crash on `compartment.outside`
        # (a CellDesigner-only relation; SBGN has no outside compartment).
        if not has_dot_binary():
            pytest.skip("graphviz `dot` binary not on PATH")
        out = pd2af.transform(
            map_with_compartments,
            mode="normal",
            keep_species=keep_species,
            layout_mode="auto",
        ).obj
        assert len(_compartment_layouts(out.layout)) == 1

    @pytest.mark.parametrize("keep_species", (True, False))
    def test_auto_layout_compartment_encloses_all_members(
        self, map_with_compartments, keep_species
    ):
        # `set_all_active` surfaces S as an activity; with the consumption
        # group excluded nothing induces an influence on it, so it is the
        # member-with-no-influences case (this test is about layout, not
        # influence inference).
        if not has_dot_binary():
            pytest.skip("graphviz `dot` binary not on PATH")
        out = pd2af.transform(
            map_with_compartments,
            mode="normal",
            keep_species=keep_species,
            layout_mode="auto",
            set_all_active=True,
            exclude_groups=("influences:consumption",),
        ).obj
        mapping = out.layout_model_mapping
        (compartment_layout,) = _compartment_layouts(out.layout)
        compartment = mapping.get_mapping(compartment_layout)
        members = [
            element
            for element in out.layout.layout_elements
            if isinstance(element, momapy.core.layout.Node)
            and not isinstance(element, momapy.sbgn.af.CompartmentLayout)
            and getattr(mapping.get_mapping(element), "compartment", None)
            is compartment
        ]
        assert len(members) == 3  # A, B and the influence-free S
        influenced = {
            endpoint
            for influence in out.model.influences
            for endpoint in (influence.source, influence.target)
        }
        assert any(mapping.get_mapping(member) not in influenced for member in members)
        for member in members:
            assert (
                compartment_layout.position.x - compartment_layout.width / 2
                <= member.position.x
                <= compartment_layout.position.x + compartment_layout.width / 2
            )
            assert (
                compartment_layout.position.y - compartment_layout.height / 2
                <= member.position.y
                <= compartment_layout.position.y + compartment_layout.height / 2
            )

    @pytest.mark.parametrize(
        "keep_species,layout_mode",
        ((True, "plain"), (False, "auto"), (True, "auto")),
    )
    def test_compartments_round_trip(
        self, map_with_compartments, keep_species, layout_mode
    ):
        if layout_mode == "auto" and not has_dot_binary():
            pytest.skip("graphviz `dot` binary not on PATH")
        out = pd2af.transform(
            map_with_compartments,
            mode="normal",
            keep_species=keep_species,
            layout_mode=layout_mode,
        ).obj
        path = os.path.join(tempfile.gettempdir(), "pd2af_test_compartments.sbgn")
        momapy.io.core.write(out, path, writer="sbgnml")
        back = momapy.io.core.read(path, reader="sbgnml").obj
        assert len(back.model.compartments) == len(out.model.compartments) == 1


class TestDropCompartments:
    """`drop_compartments` leaves an SBGN-AF map with no compartment at all.

    SBGN-AF has no default compartment to fall back on, so every activity ends
    up with `compartment=None` and the model carries no compartment.
    """

    MODES = ("normal", "no-complex")

    @pytest.fixture(scope="class")
    def map_with_compartments(self):
        return read_sbgn_map(SBGN_WITH_COMPARTMENTS_MAP_PATH)

    @pytest.mark.parametrize("mode", MODES)
    def test_model_carries_no_compartment_at_all(self, map_with_compartments, mode):
        out = pd2af.transform(
            map_with_compartments, mode=mode, layout_mode=None, drop_compartments=True
        ).obj
        assert not out.model.compartments
        assert out.model.activities
        assert all(activity.compartment is None for activity in out.model.activities)

    @pytest.mark.parametrize("mode", MODES)
    def test_layout_renders_no_compartment_at_all(self, map_with_compartments, mode):
        if not has_dot_binary():
            pytest.skip("graphviz `dot` binary not on PATH")
        out = pd2af.transform(
            map_with_compartments, mode=mode, layout_mode="auto", drop_compartments=True
        ).obj
        assert _compartment_layouts(out.layout) == []

    @pytest.mark.parametrize("mode", MODES)
    def test_output_round_trips(self, map_with_compartments, mode):
        if not has_dot_binary():
            pytest.skip("graphviz `dot` binary not on PATH")
        out = pd2af.transform(
            map_with_compartments, mode=mode, layout_mode="auto", drop_compartments=True
        ).obj
        path = os.path.join(tempfile.gettempdir(), "pd2af_test_drop_compartments.sbgn")
        momapy.io.core.write(out, path, writer="sbgnml")
        back = momapy.io.core.read(path, reader="sbgnml").obj
        assert not back.model.compartments
        assert len(back.model.activities) == len(out.model.activities)


class TestNestedOperators:
    """SBGN-PD logical operators may feed other logical operators; the nested
    logic reaches the output instead of collapsing into an input-less
    operator. The committed nested_operators map is the stat1 map with an OR
    inserted between the AND's inputs and the AND."""

    @pytest.fixture(scope="class")
    def nested_map(self):
        return read_sbgn_map(os.path.join(SBGN_MAPS_DIR, "nested_operators.sbgn"))

    @pytest.mark.parametrize("keep_species", (True, False))
    def test_nested_operators_reach_model_and_layout(self, nested_map, keep_species):
        if not has_dot_binary():
            pytest.skip("graphviz `dot` binary not on PATH")
        out = pd2af.transform(
            nested_map,
            mode="normal",
            keep_species=keep_species,
            layout_mode="auto",
        ).obj
        operators = {type(o).__name__ for o in out.model.logical_operators}
        assert operators == {"AndOperator", "OrOperator"}
        and_operator = next(
            o
            for o in out.model.logical_operators
            if isinstance(o, momapy.sbgn.af.AndOperator)
        )
        assert any(
            isinstance(i.referred_element, momapy.sbgn.af.LogicalOperator)
            for i in and_operator.inputs
        )
        operator_layouts = [
            element
            for element in out.layout.layout_elements
            if isinstance(
                element,
                (
                    momapy.sbgn.af.AndOperatorLayout,
                    momapy.sbgn.af.OrOperatorLayout,
                    momapy.sbgn.af.NotOperatorLayout,
                ),
            )
        ]
        assert len(operator_layouts) == 2


def _consumption_map(modulation_class, reactant_active=True):
    """A -> B with `modulation_class` modulation onto the process; the source
    is a separate active entity pool E."""
    a = momapy.sbgn.pd.Macromolecule(
        id_="a",
        label="A",
        state_variables=frozenset([momapy.sbgn.pd.StateVariable(value="active")])
        if reactant_active
        else frozenset(),
    )
    b = momapy.sbgn.pd.Macromolecule(id_="b", label="B")
    source = momapy.sbgn.pd.Macromolecule(id_="e", label="E")
    process = momapy.sbgn.pd.GenericProcess(
        id_="proc",
        reactants=frozenset([momapy.sbgn.pd.Reactant(id_="r_a", referred_element=a)]),
        products=frozenset([momapy.sbgn.pd.Product(id_="pr_b", referred_element=b)]),
    )
    model = momapy.sbgn.pd.SBGNPDModel(
        id_="m",
        entity_pools=frozenset([a, b, source]),
        processes=frozenset([process]),
        modulations=frozenset(
            [modulation_class(id_="mod", source=source, target=process)]
        ),
    )
    return momapy.sbgn.pd.SBGNPDMap(id_="map", model=model)


def _phosphorylation_map(product_active):
    """A -> A-P under a stimulation from E. A is active; A-P is the same
    entity, active or not as asked."""
    product_state_variables = [momapy.sbgn.pd.StateVariable(value="P")]
    if product_active:
        product_state_variables.append(momapy.sbgn.pd.StateVariable(value="active"))
    a = momapy.sbgn.pd.Macromolecule(
        id_="a",
        label="A",
        state_variables=frozenset([momapy.sbgn.pd.StateVariable(value="active")]),
    )
    phosphorylated_a = momapy.sbgn.pd.Macromolecule(
        id_="a_p",
        label="A",
        state_variables=frozenset(product_state_variables),
    )
    source = momapy.sbgn.pd.Macromolecule(id_="e", label="E")
    process = momapy.sbgn.pd.GenericProcess(
        id_="proc",
        reactants=frozenset([momapy.sbgn.pd.Reactant(id_="r_a", referred_element=a)]),
        products=frozenset(
            [momapy.sbgn.pd.Product(id_="pr_a_p", referred_element=phosphorylated_a)]
        ),
    )
    model = momapy.sbgn.pd.SBGNPDModel(
        id_="m",
        entity_pools=frozenset([a, phosphorylated_a, source]),
        processes=frozenset([process]),
        modulations=frozenset(
            [momapy.sbgn.pd.Stimulation(id_="mod", source=source, target=process)]
        ),
    )
    return momapy.sbgn.pd.SBGNPDMap(id_="map", model=model)


def _influence_arrows(out):
    return {
        (type(influence).__name__, influence.source.label, influence.target.label)
        for influence in out.model.influences
    }


class TestConsumptionInfluences:
    """The SBGN consumption/sparing rules walk a modulation arc onto a process
    down to each reactant that is itself an activity: a stimulation (catalysis
    and necessary stimulation included) consumes -- negative -- and an
    inhibition spares -- positive. Bare modulations draw no edge."""

    @pytest.mark.parametrize(
        "modulation_class,expected_name",
        (
            (momapy.sbgn.pd.Stimulation, "NegativeInfluence"),
            (momapy.sbgn.pd.Catalysis, "NegativeInfluence"),
            (momapy.sbgn.pd.NecessaryStimulation, "NegativeInfluence"),
            (momapy.sbgn.pd.Inhibition, "PositiveInfluence"),
        ),
    )
    def test_modulation_kind_drives_the_edge_sign(
        self, modulation_class, expected_name
    ):
        out = pd2af.transform(
            _consumption_map(modulation_class),
            mode="normal",
            keep_species=True,
            layout_mode=None,
        ).obj
        assert _influence_arrows(out) == {(expected_name, "E", "[active]A")}

    def test_bare_modulation_draws_no_consumption_edge(self):
        out = pd2af.transform(
            _consumption_map(momapy.sbgn.pd.Modulation),
            mode="normal",
            keep_species=True,
            layout_mode=None,
        ).obj
        assert _influence_arrows(out) == set()

    def test_inactive_reactant_draws_no_consumption_edge(self):
        out = pd2af.transform(
            _consumption_map(momapy.sbgn.pd.Stimulation, reactant_active=False),
            mode="normal",
            keep_species=True,
            layout_mode=None,
        ).obj
        assert _influence_arrows(out) == set()

    @pytest.mark.parametrize("keep_species", (True, False))
    def test_excluding_the_group_removes_only_these_edges(self, keep_species):
        with_consumption = pd2af.transform(
            _consumption_map(momapy.sbgn.pd.Stimulation),
            mode="normal",
            keep_species=keep_species,
            layout_mode=None,
        )
        without = pd2af.transform(
            _consumption_map(momapy.sbgn.pd.Stimulation),
            mode="normal",
            keep_species=keep_species,
            layout_mode=None,
            exclude_groups=("influences:consumption",),
        )
        expected_label = "[active]A" if keep_species else "A"
        assert _influence_arrows(with_consumption.obj) == {
            ("NegativeInfluence", "E", expected_label)
        }
        assert _influence_arrows(without.obj) == set()

    @pytest.mark.parametrize("keep_species", (True, False))
    def test_reactant_handed_back_as_an_activity_draws_no_edge(self, keep_species):
        out = pd2af.transform(
            _phosphorylation_map(product_active=True),
            mode="normal",
            keep_species=keep_species,
            layout_mode=None,
        ).obj
        assert not any(
            name == "NegativeInfluence" for name, _, _ in _influence_arrows(out)
        )

    def test_reactant_handed_back_inactive_keeps_its_edge(self):
        out = pd2af.transform(
            _phosphorylation_map(product_active=False),
            mode="normal",
            keep_species=True,
            layout_mode=None,
        ).obj
        assert _influence_arrows(out) == {("NegativeInfluence", "E", "[active]A")}

    def test_excluding_the_guard_rule_restores_the_edge(self):
        out = pd2af.transform(
            _phosphorylation_map(product_active=True),
            mode="normal",
            layout_mode=None,
            exclude_rules=("influences:consumption:returns_reactant_as_activity",),
        ).obj
        assert _influence_arrows(out) == {
            ("PositiveInfluence", "E", "A"),
            ("NegativeInfluence", "E", "A"),
        }

    def test_self_consumption_is_a_self_edge(self):
        # The modulation source IS the reactant: E consumes itself.
        a = momapy.sbgn.pd.Macromolecule(
            id_="a",
            label="A",
            state_variables=frozenset([momapy.sbgn.pd.StateVariable(value="active")]),
        )
        b = momapy.sbgn.pd.Macromolecule(id_="b", label="B")
        process = momapy.sbgn.pd.GenericProcess(
            id_="proc",
            reactants=frozenset(
                [momapy.sbgn.pd.Reactant(id_="r_a", referred_element=a)]
            ),
            products=frozenset(
                [momapy.sbgn.pd.Product(id_="pr_b", referred_element=b)]
            ),
        )
        model = momapy.sbgn.pd.SBGNPDModel(
            id_="m",
            entity_pools=frozenset([a, b]),
            processes=frozenset([process]),
            modulations=frozenset(
                [momapy.sbgn.pd.Stimulation(id_="mod", source=a, target=process)]
            ),
        )
        out = pd2af.transform(
            momapy.sbgn.pd.SBGNPDMap(id_="map", model=model),
            mode="normal",
            keep_species=True,
            layout_mode=None,
        ).obj
        assert _influence_arrows(out) == {
            ("NegativeInfluence", "[active]A", "[active]A")
        }

    def test_deleted_complex_reactant_is_not_rerouted(self):
        """Parity pin with CellDesigner: the direct consumption rules do not
        pass through `paths:complex_traversal`, so a deleted complex reactant
        yields no consumption edge to a promoted subunit in the no-complex
        modes."""
        subunit = momapy.sbgn.pd.MacromoleculeSubunit(id_="su", label="S")
        complex_pool = momapy.sbgn.pd.Complex(
            id_="c", label="C", subunits=frozenset([subunit])
        )
        active_subunit = momapy.sbgn.pd.Macromolecule(
            id_="s",
            label="S",
            state_variables=frozenset([momapy.sbgn.pd.StateVariable(value="active")]),
        )
        source = momapy.sbgn.pd.Macromolecule(
            id_="e",
            label="E",
            state_variables=frozenset([momapy.sbgn.pd.StateVariable(value="active")]),
        )
        product = momapy.sbgn.pd.Macromolecule(id_="p", label="P")
        process = momapy.sbgn.pd.GenericProcess(
            id_="proc",
            reactants=frozenset(
                [momapy.sbgn.pd.Reactant(id_="r_c", referred_element=complex_pool)]
            ),
            products=frozenset(
                [momapy.sbgn.pd.Product(id_="pr_p", referred_element=product)]
            ),
        )
        model = momapy.sbgn.pd.SBGNPDModel(
            id_="m",
            entity_pools=frozenset([complex_pool, source, product, active_subunit]),
            processes=frozenset([process]),
            modulations=frozenset(
                [momapy.sbgn.pd.Stimulation(id_="mod", source=source, target=process)]
            ),
        )
        sbgn_map = momapy.sbgn.pd.SBGNPDMap(id_="map", model=model)
        for mode in ("normal", "no-complex"):
            out = pd2af.transform(
                sbgn_map, mode=mode, keep_species=True, layout_mode=None
            ).obj
            assert _influence_arrows(out) == set()

    def test_cell_designer_consumes_reactants_the_same_way(self):
        # Parity: the same biology in CellDesigner (catalyzer modifier on a
        # reaction with an active reactant) also yields a negative influence
        # on the reactant.
        creb_map = read_cd_map(os.path.join(MAPS_DIR, "CREB_activity.xml"))
        out = pd2af.transform(
            creb_map, mode="normal", keep_species=True, layout_mode=None
        ).obj
        negative = [
            modulation
            for modulation in out.model.modulations
            if isinstance(modulation, momapy.celldesigner.NegativeInfluence)
        ]
        assert negative


class TestTransformModelInput:
    """A bare SBGN-PD model in `normal` mode returns a bare SBGN-AF model with
    `layout_mode` forced to None -- so no graphviz `dot` is required."""

    def test_returns_sbgn_af_model(self, proteoform_map):
        result = pd2af.transform(proteoform_map.model, mode="normal")
        assert isinstance(result.obj, momapy.sbgn.af.SBGNAFModel)

    def test_model_output_matches_map_output(self, proteoform_map):
        from_model = pd2af.transform(proteoform_map.model, mode="normal").obj
        from_map = pd2af.transform(proteoform_map, mode="normal", layout_mode=None).obj
        assert _activity_labels(from_model) == _activity_labels(from_map.model)

    @pytest.mark.parametrize("layout_mode", ["dot", "plain"])
    def test_rejects_explicit_layout_mode(self, proteoform_map, layout_mode):
        with pytest.raises(ValueError):
            pd2af.transform(
                proteoform_map.model, mode="normal", layout_mode=layout_mode
            )


class TestLanguageCompatibility:
    """A mode declares the input languages it accepts, and `transform`
    enforces it: `keep-reactions` is CellDesigner-only."""

    def test_celldesigner_only_mode_rejects_sbgn_pd_input(self, sbgn_example_map):
        with pytest.raises(ValueError) as excinfo:
            pd2af.transform(sbgn_example_map, mode="keep-reactions", layout_mode=None)
        message = str(excinfo.value)
        assert "keep-reactions" in message
        assert "sbgn_pd" in message


class TestMergingIsNotProteinOnly:
    """The merging applies to every entity class, not only macromolecules: two
    nucleic acid features differing by a state variable merge just as two
    macromolecule forms do."""

    @pytest.fixture
    def nucleic_acid_feature_map(self):
        state_active = momapy.sbgn.pd.StateVariable(
            variable="r0", value="active", order=0
        )
        state_methylated = momapy.sbgn.pd.StateVariable(
            variable="r1", value="Me", order=1
        )
        state_unmethylated = momapy.sbgn.pd.StateVariable(
            variable="r1", value=None, order=1
        )
        feature_one = momapy.sbgn.pd.NucleicAcidFeature(
            label="MYC", state_variables=frozenset([state_active, state_methylated])
        )
        feature_two = momapy.sbgn.pd.NucleicAcidFeature(
            label="MYC", state_variables=frozenset([state_active, state_unmethylated])
        )
        model = momapy.sbgn.pd.SBGNPDModel(
            entity_pools=frozenset([feature_one, feature_two])
        )
        return momapy.sbgn.pd.SBGNPDMap(model=model)

    def test_merged_by_default(self, nucleic_acid_feature_map):
        out = pd2af.transform(
            nucleic_acid_feature_map, mode="normal", layout_mode=None
        ).obj
        assert _activity_labels(out.model) == ["MYC"]

    def test_kept_apart_under_keep_species(self, nucleic_acid_feature_map):
        out = pd2af.transform(
            nucleic_acid_feature_map,
            mode="normal",
            keep_species=True,
            layout_mode=None,
        ).obj
        assert len(out.model.activities) == 2


def _binding_activation_map(name):
    return read_sbgn_map(os.path.join(BINDING_ACTIVATION_MAPS_DIR, f"{name}.sbgn"))


class TestBindingActivation:
    """The SBGN-PD side of the binding activation cases: the same small maps as
    the CellDesigner tests, with the two forms of an entity matched by entity
    kind and label instead of by template."""

    def _transform(self, name, mode, **kwargs):
        out = pd2af.transform(
            _binding_activation_map(name), mode=mode, layout_mode=None, **kwargs
        ).obj
        return _activity_labels(out.model), _influence_arrows(out)

    # M -> L, L + R -> L:R (R active only in the complex), L:R + X -> L:R:X.
    def test_normal_targets_the_complex(self):
        labels, influences = self._transform("ligand_receptor", "normal")
        assert labels == ["L", "L:R", "L:R:X", "M"]
        assert ("PositiveInfluence", "L", "L:R") in influences

    def test_no_complex_targets_the_promoted_subunit(self):
        labels, influences = self._transform("ligand_receptor", "no-complex")
        assert labels == ["L", "L:R:X", "M", "R"]
        assert ("PositiveInfluence", "L", "R") in influences

    @pytest.mark.parametrize("mode", ("normal", "no-complex"))
    def test_excluding_the_group_restores_the_plain_reading(self, mode):
        labels, influences = self._transform(
            "ligand_receptor", mode, exclude_groups=("influences:binding_activation",)
        )
        assert "L" not in labels
        assert not {one for one in influences if one[1] == "L"}

    @pytest.mark.parametrize("mode", ("normal", "no-complex"))
    def test_upstream_path_reaches_the_activator_but_the_edge_does_not_chain(
        self, mode
    ):
        _, influences = self._transform("ligand_receptor", mode)
        assert ("PositiveInfluence", "M", "L") in influences
        assert ("PositiveInfluence", "M", "L:R:X") in influences
        assert ("PositiveInfluence", "L", "L:R:X") not in influences

    def test_explicitly_active_pool_is_not_newly_activated(self):
        labels, influences = self._transform(
            "ligand_receptor", "no-complex", set_active=["r_model"]
        )
        assert "L" not in labels
        assert not {one for one in influences if one[1] == "L"}

    # Ras:GTP (Ras active) + Raf -> Ras:GTP:Raf (Ras and Raf active).
    def test_active_recruiter_is_a_source(self):
        _, influences = self._transform("active_recruiter", "normal")
        assert influences == {("PositiveInfluence", "Ras:GTP", "Ras:GTP:Raf")}
        labels, influences = self._transform("active_recruiter", "no-complex")
        assert labels == ["Raf", "Ras"]
        assert influences == {("PositiveInfluence", "Ras", "Raf")}

    # A + B -> A:B, both active in the complex.
    @pytest.mark.parametrize("mode", ("normal", "no-complex"))
    def test_two_newly_active_reactants_draw_no_edge(self, mode):
        labels, influences = self._transform("mutual", mode)
        assert labels == (["A:B"] if mode == "normal" else ["A", "B"])
        assert influences == set()

    # L + R:S (R inactive) -> L:(R:S) (R active, nested one level down).
    def test_nested_subunits(self):
        _, influences = self._transform("nested", "normal")
        assert influences == {("PositiveInfluence", "L", "L:(R:S)")}

    # X the simple chemical and X the macromolecule share a label, not a kind.
    @pytest.mark.parametrize("mode", ("normal", "no-complex"))
    def test_equal_labels_of_different_kinds_do_not_match(self, mode):
        labels, influences = self._transform("different_kinds", mode)
        assert "P" not in labels
        assert influences == set()

    # A macromolecule pool and a macromolecule subunit are the same kind.
    def test_pool_and_subunit_of_the_same_kind_match(self):
        _, influences = self._transform("reversible", "no-complex")
        assert influences == {("PositiveInfluence", "L", "R")}

    # R -> R:R, both subunits active.
    @pytest.mark.parametrize("mode", ("normal", "no-complex"))
    def test_same_entity_is_never_its_own_activator(self, mode):
        _, influences = self._transform("homodimer", mode)
        assert influences == set()

    # L (membrane) + R (cytosol) -> L:R (R active); R also stimulates Y -> Z.
    @pytest.mark.parametrize("mode", ("normal", "no-complex"))
    @pytest.mark.parametrize("keep_species", (False, True))
    @pytest.mark.parametrize("drop_compartments", (False, True))
    def test_no_self_influence(self, mode, keep_species, drop_compartments):
        _, influences = self._transform(
            "merge",
            mode,
            keep_species=keep_species,
            drop_compartments=drop_compartments,
        )
        # Under `keep_species` a label keeps its state prefix ("[active]R").
        target = "L:R" if mode == "normal" else "R"
        assert any(
            kind == "PositiveInfluence" and source == "L" and label.endswith(target)
            for kind, source, label in influences
        )
        assert not {one for one in influences if one[1] == one[2]}

    @pytest.mark.parametrize("mode", ("normal", "no-complex"))
    @pytest.mark.parametrize("name", ("ligand_receptor", "merge"))
    def test_output_round_trips(self, tmp_path, name, mode):
        out = pd2af.transform(
            _binding_activation_map(name),
            mode=mode,
            keep_species=True,
            layout_mode="plain",
        ).obj
        path = os.path.join(tmp_path, f"{name}_{mode}.sbgn")
        momapy.io.core.write(out, path, writer="sbgnml")
        momapy.io.core.read(path, reader="sbgnml")
