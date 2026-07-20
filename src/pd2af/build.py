"""Coordinator: drive the two-phase BuilderContext pipeline.

Pass 1 (:func:`pd2af.celldesigner.building_model.make_and_add_model`) walks clingo
activity / influence atoms and populates ``context.model`` with canonical,
content-deduped compartments, templates, species and modulations.
References between elements are wired to canonical instances at
construction time.

Pass 2 (:func:`pd2af.celldesigner.building_layout.make_and_add_layout`) -- skipped
entirely when ``layout_mode is None`` -- populates ``context.layout`` and
``context.layout_model_mapping``, branching on ``layout_mode``.

``build_map`` creates the ``BuilderContext``, runs the two passes, then
assembles the final CellDesignerMap from the three context slots and
returns the frozen map.
"""

import dataclasses

import momapy.builder
import momapy.celldesigner
import momapy.sbgn.af
import momapy.utils

import pd2af.celldesigner.building_layout
import pd2af.celldesigner.building_model
import pd2af.languages
import pd2af.sbgn.building_layout
import pd2af.sbgn.building_model
import pd2af.utils


@dataclasses.dataclass
class BuilderContext:
    # --- inputs ---
    input_map: object
    layout_mode: str | None
    clingo_id_to_model_element: dict
    influence_pairing: str = "cross"
    mode: str = "normal"

    # --- outputs being built ---
    model: object = None
    layout: object = None
    layout_model_mapping: object = None

    # --- Pass-1 -> Pass-2 handoff ---
    species_emissions: list = dataclasses.field(default_factory=list)
    input_model_element_to_canonical_model_element: dict = dataclasses.field(
        default_factory=dict
    )

    # --- Pass-1 scratch ---
    cache: dict = dataclasses.field(default_factory=dict)
    subunit_to_top_level: dict = None
    activity_atoms_by_key_class: dict = dataclasses.field(default_factory=dict)
    influence_atoms: list = dataclasses.field(default_factory=list)
    key_to_species: dict = dataclasses.field(default_factory=dict)

    # --- Logical-operator scratch (shared collection; per-language outputs) ---
    operator_atoms: list = dataclasses.field(default_factory=list)
    operator_input_atoms: list = dataclasses.field(default_factory=list)
    # CellDesigner: operator key -> BooleanLogicGate; emitted (gate, input_gate).
    key_to_gate: dict = dataclasses.field(default_factory=dict)
    gate_emissions: list = dataclasses.field(default_factory=list)
    # SBGN-AF: operator key -> LogicalOperator; emitted (operator, input_operator).
    key_to_operator: dict = dataclasses.field(default_factory=dict)
    operator_emissions: list = dataclasses.field(default_factory=list)

    # --- Pass-2 scratch ---
    model_element_to_layout_elements: dict = dataclasses.field(default_factory=dict)
    object_to_builder: dict = dataclasses.field(default_factory=dict)
    synthetic_index: int = 0

    # --- SBGN-AF pass scratch ---
    activity_atoms: list = dataclasses.field(default_factory=list)
    key_to_activity: dict = dataclasses.field(default_factory=dict)
    activity_emissions: list = dataclasses.field(default_factory=list)
    input_compartment_to_af_compartment: dict = dataclasses.field(
        default_factory=dict
    )
    af_compartment_to_input_compartment: dict = dataclasses.field(
        default_factory=dict
    )
    subunit_id_to_parent_compartment: dict = dataclasses.field(
        default_factory=dict
    )


def build_map(
    map_,
    layout_mode,
    clingo_model,
    clingo_id_to_model_element,
    influence_pairing="cross",
    mode="normal",
):
    context = BuilderContext(
        input_map=map_,
        layout_mode=layout_mode,
        clingo_id_to_model_element=clingo_id_to_model_element,
        influence_pairing=influence_pairing,
        mode=mode,
    )
    language = pd2af.languages.language_from_map(map_)
    if language == pd2af.languages.SBGN_PD:
        pd2af.sbgn.building_model.make_and_add_model(context, clingo_model)
        if layout_mode is not None:
            pd2af.sbgn.building_layout.make_and_add_layout(context)
        map_builder_class = momapy.builder.get_or_make_builder_cls(
            momapy.sbgn.af.SBGNAFMap
        )
    else:
        pd2af.celldesigner.building_model.make_and_add_model(context, clingo_model)
        if layout_mode is not None:
            pd2af.celldesigner.building_layout.make_and_add_layout(context)
        map_builder_class = momapy.builder.get_or_make_builder_cls(
            momapy.celldesigner.CellDesignerMap
        )

    map_builder = map_builder_class(
        model=context.model,
        layout=context.layout,
        layout_model_mapping=context.layout_model_mapping,
    )
    new_map = momapy.builder.object_from_builder(map_builder)

    # the "dot" mode (graphviz) repositions an already-built layout. The
    # compartment-layout classes differ per language (see
    # pd2af.utils.make_auto_layout).
    if layout_mode == "dot":
        if language == pd2af.languages.CELLDESIGNER:
            new_map = pd2af.utils.make_auto_layout(new_map)
        elif language == pd2af.languages.SBGN_PD:
            # SBGN-AF logical operators have input/output connectors: rank their
            # logic-arc inputs upstream (reversed_arc_classes) and re-attach the
            # operator arcs to the connector tips after graphviz repositions
            # (operator_arc_resolver). Both hooks are no-ops on operator-free
            # maps, so non-operator SBGN output is unchanged.
            new_map = pd2af.utils.make_auto_layout(
                new_map,
                compartment_layout_classes=(momapy.sbgn.af.CompartmentLayout,),
                reversed_arc_classes=(momapy.sbgn.af.LogicArcLayout,),
                operator_arc_resolver=(
                    pd2af.sbgn.building_layout.resolve_operator_arc_segments
                ),
            )

    # Imported here (not at module top) to break the core -> build import cycle.
    from pd2af.core import TransformerResult

    return TransformerResult(
        obj=new_map,
        provenance=make_provenance_from_context(context, language),
    )


def make_provenance_from_context(context, language):
    """Build the input-element -> output-element provenance mapping.

    Walks the ``key_to_*`` dicts the model pass populates for *every* activity
    and operator atom (the ``*_emissions`` lists only record the first element
    per identity and would miss the many-to-one dedup pairings). Each activity
    key carries the input element's clingo id in ``key.species`` and each
    operator key in ``key.gate``; both resolve back to the input model element
    through ``context.clingo_id_to_model_element``.

    Returns a :class:`momapy.utils.FrozenIdentityMultiDict` mapping each input
    model element to the ``frozenset`` of output elements derived from it.
    """
    if language == pd2af.languages.SBGN_PD:
        key_to_activity = context.key_to_activity
        key_to_operator = context.key_to_operator
    else:
        key_to_activity = context.key_to_species
        key_to_operator = context.key_to_gate

    input_element_to_output_elements = {}

    def record_provenance(clingo_id, output_element):
        if output_element is None:
            return
        input_element = context.clingo_id_to_model_element[clingo_id]
        input_element_to_output_elements.setdefault(input_element, set()).add(
            output_element
        )

    for activity_key, output_element in key_to_activity.items():
        record_provenance(activity_key.species, output_element)
    for operator_key, output_element in key_to_operator.items():
        record_provenance(operator_key.gate, output_element)

    return momapy.utils.FrozenIdentityMultiDict(
        {
            input_element: frozenset(output_elements)
            for input_element, output_elements in (
                input_element_to_output_elements.items()
            )
        }
    )
