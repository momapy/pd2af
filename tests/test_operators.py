"""Authored logical operators: CellDesigner ``BooleanLogicGate`` and SBGN-PD
``LogicalOperator`` -> AF operator nodes, input edges and operator-sourced
influences.

Coverage mirrors the plan's verification list:

* predicate level -- each typed influence accepts a ``logicalOperatorKey``
  source; ``new`` wraps the operator node and input predicates;
* rule level -- the ``_GATES`` group is present in every path-inference mode,
  the gate-input activation clause is emitted per language, and the NOT token
  dodges the reserved ``not`` keyword;
* builder/layout level -- the operator-type and operator-class maps are total
  over the gate/operator classes they dispatch on;
* integration -- CREB (Shape A: gate as reaction modifier), the SBGN-PD
  operator map (Shape B: operator as modulation source), and a read-back sweep
  over every committed CellDesigner gate map. An operator whose target is not an
  activity yields no influence and is dropped (SRR), exactly as a species source
  to a non-activity target is.
"""

import glob
import os
import tempfile

import clorm
import pytest

import momapy.celldesigner
import momapy.io.core
import momapy.sbgn.af

import pd2af
import pd2af.celldesigner.building_layout
import pd2af.celldesigner.building_model
import pd2af.languages
import pd2af.predicates
import pd2af.rules
import pd2af.sbgn.building_layout
import pd2af.sbgn.building_model

from tests._helpers import (
    MAPS_DIR,
    SBGN_MAPS_DIR,
    has_dot_binary,
    read_cd_map,
    read_sbgn_map,
)


_PATH_INFERENCE_MODES = (
    "normal",
    "normal-no-complex",
    "keep-species",
    "keep-species-no-complex",
)

_TYPED_INFLUENCE_NAMES = (
    "positivelyInfluences",
    "negativelyInfluences",
    "modulates",
    "triggers",
    "unknownPositivelyInfluences",
    "unknownNegativelyInfluences",
    "unknownModulates",
    "unknownTriggers",
)

_CREB_MAP_PATH = os.path.join(MAPS_DIR, "CREB_activity.xml")
_SRR_MAP_PATH = os.path.join(MAPS_DIR, "SRR_signaling.xml")
_SBGN_OPERATOR_MAP_PATH = os.path.join(
    SBGN_MAPS_DIR, "activated_stat1alpha_induction_of_the_irf1_gene.sbgn"
)


def _points_close(first, second, tolerance=1e-6):
    return abs(first.x - second.x) < tolerance and abs(first.y - second.y) < tolerance


def _gate_map_paths():
    paths = []
    for path in sorted(glob.glob(os.path.join(MAPS_DIR, "*.xml"))):
        with open(path) as handle:
            content = handle.read()
        if "BOOLEAN_LOGIC_GATE" in content or "GateMember" in content:
            paths.append(path)
    return paths


class TestOperatorPredicates:
    def test_logical_operator_key_is_clorm_predicate(self):
        assert issubclass(pd2af.predicates.logicalOperatorKey, clorm.Predicate)

    def test_logical_operator_is_clorm_predicate(self):
        assert issubclass(pd2af.predicates.logicalOperator, clorm.Predicate)

    def test_logical_operator_input_is_clorm_predicate(self):
        assert issubclass(pd2af.predicates.logicalOperatorInput, clorm.Predicate)

    def test_logical_operator_key_can_be_constructed(self):
        atom = pd2af.predicates.logicalOperatorKey(gate="some_gate")
        assert atom.gate == "some_gate"

    def test_logical_operator_carries_type_token(self):
        atom = pd2af.predicates.logicalOperator(
            key=pd2af.predicates.logicalOperatorKey(gate="g"), type_="and"
        )
        assert atom.type_ == "and"
        assert isinstance(atom.key, pd2af.predicates.logicalOperatorKey)

    def test_each_typed_influence_accepts_operator_source(self):
        for name in _TYPED_INFLUENCE_NAMES:
            predicate = getattr(pd2af.predicates, name)
            atom = predicate(
                source=pd2af.predicates.logicalOperatorKey(gate="g"),
                target=pd2af.predicates.keptSpeciesKey(species="t"),
            )
            assert isinstance(atom.source, pd2af.predicates.logicalOperatorKey), name
            assert isinstance(atom.target, pd2af.predicates.keptSpeciesKey), name

    def test_new_wraps_operator_node_and_input(self):
        key = pd2af.predicates.logicalOperatorKey(gate="g")
        node = pd2af.predicates.new(
            object_=pd2af.predicates.logicalOperator(key=key, type_="or")
        )
        edge = pd2af.predicates.new(
            object_=pd2af.predicates.logicalOperatorInput(
                operator=key, input=pd2af.predicates.keptSpeciesKey(species="i")
            )
        )
        assert isinstance(node.object_, pd2af.predicates.logicalOperator)
        assert isinstance(edge.object_, pd2af.predicates.logicalOperatorInput)


class TestGateRules:
    @pytest.mark.parametrize("mode", _PATH_INFERENCE_MODES)
    @pytest.mark.parametrize("language", tuple(pd2af.languages.LANGUAGES))
    def test_gates_group_present_in_path_inference_modes(self, mode, language):
        program = pd2af.rules.build_program(mode, language=language)
        assert "new(logicalOperator(logicalOperatorKey(OPERATOR)," in program
        assert "new(logicalOperatorInput(logicalOperatorKey(OPERATOR)," in program
        assert (
            "influences(logicalOperatorKey(OPERATOR), TARGET_KEY, INFLUENCE_KIND)"
            in program
        )

    def test_celldesigner_activates_gate_inputs(self):
        program = pd2af.rules.build_program("keep-species", language="celldesigner")
        assert "hasActivityCandidate(ELEMENT, isGateInput)" in program
        assert "booleanLogicGateInput(INPUT)" in program

    def test_sbgn_pd_operator_input_activation_is_entity_pool_guarded(self):
        program = pd2af.rules.build_program("keep-species", language="sbgn_pd")
        assert "hasActivityCandidate(ELEMENT, isGateInput)" in program
        assert "logicalOperatorInput(INPUT)" in program
        assert "entityPool(ELEMENT)" in program

    def test_not_token_dodges_reserved_keyword(self):
        # bare `not` is a reserved clingo keyword, so the NOT token is `not_`.
        program = pd2af.rules.build_program("keep-species", language="celldesigner")
        assert (
            "logicalOperator(logicalOperatorKey(OPERATOR), not_)) :- notGate(OPERATOR)."
            in program
        )

    def test_gates_influence_rule_guards_on_umbrella(self):
        # Without the booleanLogicGate/logicalOperator umbrella guard, every
        # path/3 source (species included) would be read as an operator key.
        cd_program = pd2af.rules.build_program("keep-species", language="celldesigner")
        assert "booleanLogicGate(OPERATOR)" in cd_program
        sbgn_program = pd2af.rules.build_program("keep-species", language="sbgn_pd")
        assert "logicalOperator(OPERATOR)" in sbgn_program


class TestOperatorClassMapsAreTotal:
    def test_celldesigner_gate_type_map_covers_every_gate_class(self):
        gate_classes = set(
            pd2af.celldesigner.building_model._OPERATOR_TYPE_TO_GATE_CLASS.values()
        )
        assert gate_classes == {
            momapy.celldesigner.AndGate,
            momapy.celldesigner.OrGate,
            momapy.celldesigner.NotGate,
            momapy.celldesigner.UnknownGate,
        }

    def test_celldesigner_gate_layout_map_covers_every_gate_class(self):
        layout_map = pd2af.celldesigner.building_layout._GATE_CLASS_TO_LAYOUT_CLASS
        assert set(layout_map) == set(
            pd2af.celldesigner.building_model._OPERATOR_TYPE_TO_GATE_CLASS.values()
        )

    def test_sbgn_operator_type_map_covers_authored_operators(self):
        operator_classes = set(
            pd2af.sbgn.building_model._OPERATOR_TYPE_TO_OPERATOR_CLASS.values()
        )
        assert operator_classes == {
            momapy.sbgn.af.AndOperator,
            momapy.sbgn.af.OrOperator,
            momapy.sbgn.af.NotOperator,
        }

    def test_sbgn_operator_layout_map_covers_its_operator_classes(self):
        layout_map = pd2af.sbgn.building_layout._OPERATOR_CLASS_TO_LAYOUT_CLASS
        model_classes = set(
            pd2af.sbgn.building_model._OPERATOR_TYPE_TO_OPERATOR_CLASS.values()
        )
        assert model_classes.issubset(set(layout_map))

    def test_not_token_is_consistent_across_maps(self):
        assert "not_" in pd2af.celldesigner.building_model._OPERATOR_TYPE_TO_GATE_CLASS
        assert "not_" in pd2af.sbgn.building_model._OPERATOR_TYPE_TO_OPERATOR_CLASS


class TestCelldesignerGatesShapeA:
    """CREB carries one AND gate whose inputs feed a catalyzed reaction
    (Shape A); the gate must surface as a ``BooleanLogicGate`` sourcing one
    positive influence."""

    @pytest.fixture(scope="class")
    def creb_map(self):
        return read_cd_map(_CREB_MAP_PATH)

    def test_emits_one_and_gate_with_two_inputs(self, creb_map):
        out = pd2af.transform(creb_map, mode="keep-species", layout_mode=None).obj
        gates = list(out.model.boolean_logic_gates)
        assert len(gates) == 1
        assert isinstance(gates[0], momapy.celldesigner.AndGate)
        assert len(gates[0].inputs) == 2

    def test_gate_sources_exactly_one_modulation(self, creb_map):
        out = pd2af.transform(creb_map, mode="keep-species", layout_mode=None).obj
        gate_modulations = [
            modulation
            for modulation in out.model.modulations
            if isinstance(modulation.source, momapy.celldesigner.BooleanLogicGate)
        ]
        assert len(gate_modulations) == 1
        assert isinstance(gate_modulations[0], momapy.celldesigner.PositiveInfluence)

    def test_gate_inputs_are_activities_in_the_model(self, creb_map):
        out = pd2af.transform(creb_map, mode="keep-species", layout_mode=None).obj
        gate = next(iter(out.model.boolean_logic_gates))
        model_species = set(out.model.species)
        for gate_input in gate.inputs:
            assert gate_input.referred_element in model_species

    @pytest.mark.parametrize(
        "mode,layout_mode",
        (
            ("keep-species", "plain"),
            ("keep-species", "overlay"),
            ("keep-species", "auto"),
            ("normal", "auto"),
            ("normal-no-complex", "auto"),
        ),
    )
    def test_round_trips(self, creb_map, mode, layout_mode):
        if layout_mode == "auto" and not has_dot_binary():
            pytest.skip("graphviz `dot` binary not on PATH")
        out = pd2af.transform(creb_map, mode=mode, layout_mode=layout_mode).obj
        assert len(out.model.boolean_logic_gates) == 1
        path = os.path.join(tempfile.gettempdir(), "pd2af_test_gate.xml")
        momapy.io.core.write(out, path, writer="celldesigner")
        back = momapy.io.core.read(path).obj
        assert len(back.model.boolean_logic_gates) == 1


class TestCelldesignerGatesAreAdditive:
    """An operator whose target is not an activity yields no influence and is
    dropped -- exactly as a species source to a non-activity target is. SRR's
    gates all modulate pure-sink RNAs, so SRR emits zero gates: gates are
    purely additive, never spuriously emitted."""

    def test_srr_emits_no_surviving_gates(self):
        srr_map = read_cd_map(_SRR_MAP_PATH)
        out = pd2af.transform(srr_map, mode="keep-species", layout_mode=None).obj
        assert len(out.model.boolean_logic_gates) == 0


@pytest.mark.parametrize("path", _gate_map_paths())
def test_every_gate_map_round_trips(path):
    """Read-back invariant over every committed CellDesigner gate map: the
    gate count survives a write -> read cycle (plain layout, no graphviz)."""
    cd_map = read_cd_map(path)
    out = pd2af.transform(cd_map, mode="keep-species", layout_mode="plain").obj
    written = os.path.join(tempfile.gettempdir(), "pd2af_test_gate_map.xml")
    momapy.io.core.write(out, written, writer="celldesigner")
    back = momapy.io.core.read(written).obj
    assert len(back.model.boolean_logic_gates) == len(out.model.boolean_logic_gates)


class TestSbgnOperatorsShapeB:
    """The stat1alpha SBGN-PD map carries one AND operator feeding a necessary
    stimulation (Shape B: operator as modulation source) -> one AF
    ``AndOperator`` sourcing one ``NecessaryStimulation``."""

    @pytest.fixture(scope="class")
    def operator_map(self):
        return read_sbgn_map(_SBGN_OPERATOR_MAP_PATH)

    def test_emits_one_and_operator_with_two_inputs(self, operator_map):
        out = pd2af.transform(operator_map, mode="keep-species", layout_mode=None).obj
        operators = list(out.model.logical_operators)
        assert len(operators) == 1
        assert isinstance(operators[0], momapy.sbgn.af.AndOperator)
        assert len(operators[0].inputs) == 2

    def test_operator_sources_one_influence(self, operator_map):
        out = pd2af.transform(operator_map, mode="keep-species", layout_mode=None).obj
        operator_influences = [
            influence
            for influence in out.model.influences
            if isinstance(influence.source, momapy.sbgn.af.LogicalOperator)
        ]
        assert len(operator_influences) == 1
        assert isinstance(operator_influences[0], momapy.sbgn.af.NecessaryStimulation)

    @pytest.mark.parametrize(
        "mode,layout_mode",
        (
            ("keep-species", "plain"),
            ("keep-species", "auto"),
            ("normal", "auto"),
            ("normal-no-complex", "auto"),
        ),
    )
    def test_round_trips(self, operator_map, mode, layout_mode):
        if layout_mode == "auto" and not has_dot_binary():
            pytest.skip("graphviz `dot` binary not on PATH")
        out = pd2af.transform(operator_map, mode=mode, layout_mode=layout_mode).obj
        assert len(out.model.logical_operators) == 1
        path = os.path.join(tempfile.gettempdir(), "pd2af_test_operator.sbgn")
        momapy.io.core.write(out, path, writer="sbgnml")
        back = momapy.io.core.read(path, reader="sbgnml").obj
        assert len(back.model.logical_operators) == 1

    @pytest.mark.parametrize("layout_mode", ("plain", "auto"))
    def test_arcs_attach_to_operator_connectors(self, operator_map, layout_mode):
        """Logic arcs must meet the operator's input connector tip and the
        operator-sourced influence arc must leave its output connector tip --
        not the circle border. The two connectors are distinct points."""
        if layout_mode == "auto" and not has_dot_binary():
            pytest.skip("graphviz `dot` binary not on PATH")
        out = pd2af.transform(
            operator_map, mode="keep-species", layout_mode=layout_mode
        ).obj
        operator_layout = next(
            element
            for element in out.layout.layout_elements
            if isinstance(element, momapy.sbgn.af.AndOperatorLayout)
        )
        if operator_layout.left_to_right:
            input_tip = operator_layout.left_connector_tip()
            output_tip = operator_layout.right_connector_tip()
        else:
            input_tip = operator_layout.right_connector_tip()
            output_tip = operator_layout.left_connector_tip()
        assert input_tip != output_tip
        logic_arc_starts = []
        influence_arc_starts = []
        for element in out.layout.layout_elements:
            if isinstance(element, momapy.sbgn.af.LogicArcLayout):
                logic_arc_starts.append(element.points()[0])
            elif getattr(element, "source", None) is operator_layout:
                influence_arc_starts.append(element.points()[0])
        assert logic_arc_starts
        assert influence_arc_starts
        for start in logic_arc_starts:
            assert _points_close(start, input_tip)
        for start in influence_arc_starts:
            assert _points_close(start, output_tip)
