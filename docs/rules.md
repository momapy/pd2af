# Rules reference

The pd2af transformation is driven by a set of rules (encoded in [ASP](https://en.wikipedia.org/wiki/Answer_set_programming)). Every rule is language-agnostic (a *base* rule) or specific to an input language (a *variant* rule, for CellDesigner or SBGN PD). Rules are organized in groups, each representing a coherent functional unit — the whole a `--exclude-group` can drop. Every group id reads `<family>:<member>`, and the family is the coarse concern the group addresses.

This page is generated from the live registry in `pd2af.rules`. It presents a [by-mode](#by-mode) view of which rules each mode uses, a [canonical reference](#rule-groups) documenting each rule once, and an alphabetical [index](#rule-index) of every rule.

## By mode {#by-mode}

Within a mode a group is *excludable* (nothing else the mode includes depends on it, so `--exclude-group` can drop it cleanly) or *mandatory* (another included group depends on it, directly or by filling a slot it requires). The same group can be excludable in one mode and mandatory in another.

### `normal`

**`activity:core`** (mandatory)

Mandatory activity machinery: the candidate->activity bridge and the two global toggles.

| Rule | Variant |
| --- | --- |
| [`activity:core:from_candidate`](#rule-activity-core-from_candidate) |  |
| [`activity:core:from_global_suppress`](#rule-activity-core-from_global_suppress) |  |
| [`activity:core:celldesigner:from_global_activate`](#rule-activity-core-celldesigner-from_global_activate) | CellDesigner |
| [`activity:core:sbgn_pd:from_global_activate`](#rule-activity-core-sbgn_pd-from_global_activate) | SBGN PD |

**`activity:phenotype`** (excludable)

Excludable: a phenotype is an activity candidate (reason `isPhenotype`).

| Rule | Variant |
| --- | --- |
| [`activity:phenotype:from_phenotype`](#rule-activity-phenotype-from_phenotype) |  |

**`activity:active_marker`** (excludable)

Excludable: an explicit active marker makes a species/entity pool an activity candidate (CellDesigner: `hasActive` flag or active structural state; SBGN-PD: active state variable, including on subunits).

| Rule | Variant |
| --- | --- |
| [`activity:active_marker:celldesigner:from_active_flag`](#rule-activity-active_marker-celldesigner-from_active_flag) | CellDesigner |
| [`activity:active_marker:celldesigner:from_active_structural_state`](#rule-activity-active_marker-celldesigner-from_active_structural_state) | CellDesigner |
| [`activity:active_marker:sbgn_pd:from_active_state_variable`](#rule-activity-active_marker-sbgn_pd-from_active_state_variable) | SBGN PD |
| [`activity:active_marker:sbgn_pd:from_active_subunit_state_variable`](#rule-activity-active_marker-sbgn_pd-from_active_subunit_state_variable) | SBGN PD |

**`activity:modulation_source`** (excludable)

Excludable: the source of a modulation arc (CellDesigner: a modulation arc or a reaction modifier, which in SBGN-PD folds into the same arc-to-process case) is an activity candidate.

| Rule | Variant |
| --- | --- |
| [`activity:modulation_source:celldesigner:from_modulation_source`](#rule-activity-modulation_source-celldesigner-from_modulation_source) | CellDesigner |
| [`activity:modulation_source:celldesigner:from_reaction_modulator`](#rule-activity-modulation_source-celldesigner-from_reaction_modulator) | CellDesigner |
| [`activity:modulation_source:sbgn_pd:from_modulation_source`](#rule-activity-modulation_source-sbgn_pd-from_modulation_source) | SBGN PD |

**`activity:gate_input`** (excludable)

Excludable: an input feeding a boolean logic gate / logical operator is an activity candidate (reason `isGateInput`).

| Rule | Variant |
| --- | --- |
| [`activity:gate_input:celldesigner:from_gate_input`](#rule-activity-gate_input-celldesigner-from_gate_input) | CellDesigner |
| [`activity:gate_input:sbgn_pd:from_operator_input`](#rule-activity-gate_input-sbgn_pd-from_operator_input) | SBGN PD |

**`topology:core`** (mandatory)

Mode-agnostic structural helpers shared by every mode: `isSubunit`, `hasActiveDescendantSubunit`.

| Rule | Variant |
| --- | --- |
| [`topology:core:is_subunit`](#rule-topology-core-is_subunit) |  |
| [`topology:core:has_active_descendant_direct`](#rule-topology-core-has_active_descendant_direct) |  |
| [`topology:core:has_active_descendant_transitive`](#rule-topology-core-has_active_descendant_transitive) |  |

**`topology:top_level`** (mandatory)

Resolves every species to its outermost top-level entity: a non-subunit resolves to itself; a subunit -- at any nesting depth -- resolves to the outermost complex that contains it.

| Rule | Variant |
| --- | --- |
| [`topology:top_level:recursive`](#rule-topology-top_level-recursive) |  |
| [`topology:top_level:celldesigner:self`](#rule-topology-top_level-celldesigner-self) | CellDesigner |
| [`topology:top_level:sbgn_pd:self`](#rule-topology-top_level-sbgn_pd-self) | SBGN PD |
| [`topology:top_level:sbgn_pd:phenotype_self`](#rule-topology-top_level-sbgn_pd-phenotype_self) | SBGN PD |

**`preparation:complex`** (mandatory)

The complex-keeping modes (`normal`, `keep-species`, `keep-reactions`) key a species with activity by the `keptSpeciesKey` of its top-level entity (the `topology:top_level` group): itself when top-level, its outermost complex when a subunit.

| Rule | Variant |
| --- | --- |
| [`preparation:complex:key`](#rule-preparation-complex-key) |  |
| [`preparation:complex:celldesigner:carrier`](#rule-preparation-complex-celldesigner-carrier) | CellDesigner |
| [`preparation:complex:sbgn_pd:carrier_entity_pool`](#rule-preparation-complex-sbgn_pd-carrier_entity_pool) | SBGN PD |
| [`preparation:complex:sbgn_pd:carrier_phenotype`](#rule-preparation-complex-sbgn_pd-carrier_phenotype) | SBGN PD |

**`paths:core`** (mandatory)

Builds the direct (single-hop) kinded `propagatesInfluence(SOURCE_SPECIES, TARGET_SPECIES, INFLUENCE_KIND)` relation from PD reactions and modulation arcs, plus the production-cycle relations that gate transitive extension.

| Rule | Variant |
| --- | --- |
| [`paths:core:is_transformed_to`](#rule-paths-core-is_transformed_to) |  |
| [`paths:core:is_cyclically_transformed_to`](#rule-paths-core-is_cyclically_transformed_to) |  |
| [`paths:core:composes_to`](#rule-paths-core-composes_to) |  |
| [`paths:core:celldesigner:catalyzer_to_product`](#rule-paths-core-celldesigner-catalyzer_to_product) | CellDesigner |
| [`paths:core:celldesigner:physical_stimulator_to_product`](#rule-paths-core-celldesigner-physical_stimulator_to_product) | CellDesigner |
| [`paths:core:celldesigner:trigger_to_product`](#rule-paths-core-celldesigner-trigger_to_product) | CellDesigner |
| [`paths:core:celldesigner:modulation_arc_influence`](#rule-paths-core-celldesigner-modulation_arc_influence) | CellDesigner |
| [`paths:core:celldesigner:inhibitor_to_product`](#rule-paths-core-celldesigner-inhibitor_to_product) | CellDesigner |
| [`paths:core:celldesigner:modulator_to_product`](#rule-paths-core-celldesigner-modulator_to_product) | CellDesigner |
| [`paths:core:celldesigner:unknown_catalyzer_to_product`](#rule-paths-core-celldesigner-unknown_catalyzer_to_product) | CellDesigner |
| [`paths:core:celldesigner:unknown_inhibitor_to_product`](#rule-paths-core-celldesigner-unknown_inhibitor_to_product) | CellDesigner |
| [`paths:core:celldesigner:is_directly_transformed_to`](#rule-paths-core-celldesigner-is_directly_transformed_to) | CellDesigner |
| [`paths:core:sbgn_pd:modulation_to_product`](#rule-paths-core-sbgn_pd-modulation_to_product) | SBGN PD |
| [`paths:core:sbgn_pd:modulation_to_phenotype`](#rule-paths-core-sbgn_pd-modulation_to_phenotype) | SBGN PD |
| [`paths:core:sbgn_pd:is_directly_transformed_to`](#rule-paths-core-sbgn_pd-is_directly_transformed_to) | SBGN PD |

**`paths:chaining`** (excludable)

Excludable multi-hop transitivity: extends a `propagatesInfluence` path through one more reactant->product hop (a reaction in CellDesigner, a process in SBGN-PD), carrying the kind via `composesTo` (triggering degrades to positivelyInfluences) and refusing to extend across a hop inside a production cycle (`not isCyclicallyTransformedTo`).

| Rule | Variant |
| --- | --- |
| [`paths:chaining:celldesigner:transitive_through_reaction`](#rule-paths-chaining-celldesigner-transitive_through_reaction) | CellDesigner |
| [`paths:chaining:sbgn_pd:transitive_through_process`](#rule-paths-chaining-sbgn_pd-transitive_through_process) | SBGN PD |

**`influences:kind`** (mandatory)

Maps each modulation arc to its influence kind via `hasInfluenceKind(MODULATION, INFLUENCE_KIND)`, the single place the arc-type->kind knowledge lives.

| Rule | Variant |
| --- | --- |
| [`influences:kind:celldesigner:catalysis`](#rule-influences-kind-celldesigner-catalysis) | CellDesigner |
| [`influences:kind:celldesigner:physical_stimulation`](#rule-influences-kind-celldesigner-physical_stimulation) | CellDesigner |
| [`influences:kind:celldesigner:positive_influence`](#rule-influences-kind-celldesigner-positive_influence) | CellDesigner |
| [`influences:kind:celldesigner:triggering`](#rule-influences-kind-celldesigner-triggering) | CellDesigner |
| [`influences:kind:celldesigner:inhibition`](#rule-influences-kind-celldesigner-inhibition) | CellDesigner |
| [`influences:kind:celldesigner:negative_influence`](#rule-influences-kind-celldesigner-negative_influence) | CellDesigner |
| [`influences:kind:celldesigner:modulation`](#rule-influences-kind-celldesigner-modulation) | CellDesigner |
| [`influences:kind:celldesigner:unknown_catalysis`](#rule-influences-kind-celldesigner-unknown_catalysis) | CellDesigner |
| [`influences:kind:celldesigner:unknown_physical_stimulation`](#rule-influences-kind-celldesigner-unknown_physical_stimulation) | CellDesigner |
| [`influences:kind:celldesigner:unknown_positive_influence`](#rule-influences-kind-celldesigner-unknown_positive_influence) | CellDesigner |
| [`influences:kind:celldesigner:unknown_triggering`](#rule-influences-kind-celldesigner-unknown_triggering) | CellDesigner |
| [`influences:kind:celldesigner:unknown_inhibition`](#rule-influences-kind-celldesigner-unknown_inhibition) | CellDesigner |
| [`influences:kind:celldesigner:unknown_negative_influence`](#rule-influences-kind-celldesigner-unknown_negative_influence) | CellDesigner |
| [`influences:kind:celldesigner:unknown_modulation`](#rule-influences-kind-celldesigner-unknown_modulation) | CellDesigner |
| [`influences:kind:sbgn_pd:necessary_stimulation`](#rule-influences-kind-sbgn_pd-necessary_stimulation) | SBGN PD |
| [`influences:kind:sbgn_pd:stimulation`](#rule-influences-kind-sbgn_pd-stimulation) | SBGN PD |
| [`influences:kind:sbgn_pd:inhibition`](#rule-influences-kind-sbgn_pd-inhibition) | SBGN PD |
| [`influences:kind:sbgn_pd:modulation`](#rule-influences-kind-sbgn_pd-modulation) | SBGN PD |

**`influences:core`** (excludable)

Derivation: emits `new(activity(KEY))` for every activity key, and lifts every kinded path into the internal `influences(SOURCE_KEY, TARGET_KEY, INFLUENCE_KIND)` relation by resolving both endpoints through their activity carrier and key.

| Rule | Variant |
| --- | --- |
| [`influences:core:activity`](#rule-influences-core-activity) |  |
| [`influences:core:path`](#rule-influences-core-path) |  |

**`influences:consumption`** (excludable)

Excludable consumption/sparing reasoning: a reaction depletes its reactants, so a modifier that drives the reaction also acts on every reactant that is itself an activity -- catalyzer/physicalStimulator/trigger negatively influence each consumed reactant, inhibitor positively influences each spared reactant, and the unknown modifiers contribute the unknown twins.

| Rule | Variant |
| --- | --- |
| [`influences:consumption:celldesigner:catalyzer_consumes_reactant`](#rule-influences-consumption-celldesigner-catalyzer_consumes_reactant) | CellDesigner |
| [`influences:consumption:celldesigner:physical_stimulator_consumes_reactant`](#rule-influences-consumption-celldesigner-physical_stimulator_consumes_reactant) | CellDesigner |
| [`influences:consumption:celldesigner:trigger_consumes_reactant`](#rule-influences-consumption-celldesigner-trigger_consumes_reactant) | CellDesigner |
| [`influences:consumption:celldesigner:inhibitor_spares_reactant`](#rule-influences-consumption-celldesigner-inhibitor_spares_reactant) | CellDesigner |
| [`influences:consumption:celldesigner:unknown_catalyzer_consumes_reactant`](#rule-influences-consumption-celldesigner-unknown_catalyzer_consumes_reactant) | CellDesigner |
| [`influences:consumption:celldesigner:unknown_inhibitor_spares_reactant`](#rule-influences-consumption-celldesigner-unknown_inhibitor_spares_reactant) | CellDesigner |

**`influences:output`** (excludable)

Shared fan-out from the internal `influences(SOURCE, TARGET, INFLUENCE_KIND)` relation to the typed `new(...)` influence heads — one rule per kind.

| Rule | Variant |
| --- | --- |
| [`influences:output:positive`](#rule-influences-output-positive) |  |
| [`influences:output:negative`](#rule-influences-output-negative) |  |
| [`influences:output:modulation`](#rule-influences-output-modulation) |  |
| [`influences:output:triggering`](#rule-influences-output-triggering) |  |
| [`influences:output:unknown_positive`](#rule-influences-output-unknown_positive) |  |
| [`influences:output:unknown_negative`](#rule-influences-output-unknown_negative) |  |
| [`influences:output:unknown_modulation`](#rule-influences-output-unknown_modulation) |  |
| [`influences:output:unknown_triggering`](#rule-influences-output-unknown_triggering) |  |

**`gates:core`** (excludable)

Authored logical operators (CellDesigner `BooleanLogicGate`, SBGN-PD `LogicalOperator`): emits the operator node (`logicalOperator/2`, token-typed), its input edges (`logicalOperatorInput/2`, each input resolved through carrier/key), and the operator-sourced influence written straight into the internal `influences/3` relation by reusing the existing `propagatesInfluence(OPERATOR, TARGET, INFLUENCE_KIND)` closure.

| Rule | Variant |
| --- | --- |
| [`gates:core:celldesigner:node_and`](#rule-gates-core-celldesigner-node_and) | CellDesigner |
| [`gates:core:celldesigner:node_or`](#rule-gates-core-celldesigner-node_or) | CellDesigner |
| [`gates:core:celldesigner:node_not`](#rule-gates-core-celldesigner-node_not) | CellDesigner |
| [`gates:core:celldesigner:node_unknown`](#rule-gates-core-celldesigner-node_unknown) | CellDesigner |
| [`gates:core:celldesigner:input_edge`](#rule-gates-core-celldesigner-input_edge) | CellDesigner |
| [`gates:core:celldesigner:influence`](#rule-gates-core-celldesigner-influence) | CellDesigner |
| [`gates:core:sbgn_pd:node_and`](#rule-gates-core-sbgn_pd-node_and) | SBGN PD |
| [`gates:core:sbgn_pd:node_or`](#rule-gates-core-sbgn_pd-node_or) | SBGN PD |
| [`gates:core:sbgn_pd:node_not`](#rule-gates-core-sbgn_pd-node_not) | SBGN PD |
| [`gates:core:sbgn_pd:input_edge`](#rule-gates-core-sbgn_pd-input_edge) | SBGN PD |
| [`gates:core:sbgn_pd:influence`](#rule-gates-core-sbgn_pd-influence) | SBGN PD |

### `normal-no-complex`

**`activity:core`** (mandatory)

Mandatory activity machinery: the candidate->activity bridge and the two global toggles.

| Rule | Variant |
| --- | --- |
| [`activity:core:from_candidate`](#rule-activity-core-from_candidate) |  |
| [`activity:core:from_global_suppress`](#rule-activity-core-from_global_suppress) |  |
| [`activity:core:celldesigner:from_global_activate`](#rule-activity-core-celldesigner-from_global_activate) | CellDesigner |
| [`activity:core:sbgn_pd:from_global_activate`](#rule-activity-core-sbgn_pd-from_global_activate) | SBGN PD |

**`activity:phenotype`** (excludable)

Excludable: a phenotype is an activity candidate (reason `isPhenotype`).

| Rule | Variant |
| --- | --- |
| [`activity:phenotype:from_phenotype`](#rule-activity-phenotype-from_phenotype) |  |

**`activity:active_marker`** (excludable)

Excludable: an explicit active marker makes a species/entity pool an activity candidate (CellDesigner: `hasActive` flag or active structural state; SBGN-PD: active state variable, including on subunits).

| Rule | Variant |
| --- | --- |
| [`activity:active_marker:celldesigner:from_active_flag`](#rule-activity-active_marker-celldesigner-from_active_flag) | CellDesigner |
| [`activity:active_marker:celldesigner:from_active_structural_state`](#rule-activity-active_marker-celldesigner-from_active_structural_state) | CellDesigner |
| [`activity:active_marker:sbgn_pd:from_active_state_variable`](#rule-activity-active_marker-sbgn_pd-from_active_state_variable) | SBGN PD |
| [`activity:active_marker:sbgn_pd:from_active_subunit_state_variable`](#rule-activity-active_marker-sbgn_pd-from_active_subunit_state_variable) | SBGN PD |

**`activity:modulation_source`** (excludable)

Excludable: the source of a modulation arc (CellDesigner: a modulation arc or a reaction modifier, which in SBGN-PD folds into the same arc-to-process case) is an activity candidate.

| Rule | Variant |
| --- | --- |
| [`activity:modulation_source:celldesigner:from_modulation_source`](#rule-activity-modulation_source-celldesigner-from_modulation_source) | CellDesigner |
| [`activity:modulation_source:celldesigner:from_reaction_modulator`](#rule-activity-modulation_source-celldesigner-from_reaction_modulator) | CellDesigner |
| [`activity:modulation_source:sbgn_pd:from_modulation_source`](#rule-activity-modulation_source-sbgn_pd-from_modulation_source) | SBGN PD |

**`activity:gate_input`** (excludable)

Excludable: an input feeding a boolean logic gate / logical operator is an activity candidate (reason `isGateInput`).

| Rule | Variant |
| --- | --- |
| [`activity:gate_input:celldesigner:from_gate_input`](#rule-activity-gate_input-celldesigner-from_gate_input) | CellDesigner |
| [`activity:gate_input:sbgn_pd:from_operator_input`](#rule-activity-gate_input-sbgn_pd-from_operator_input) | SBGN PD |

**`topology:core`** (mandatory)

Mode-agnostic structural helpers shared by every mode: `isSubunit`, `hasActiveDescendantSubunit`.

| Rule | Variant |
| --- | --- |
| [`topology:core:is_subunit`](#rule-topology-core-is_subunit) |  |
| [`topology:core:has_active_descendant_direct`](#rule-topology-core-has_active_descendant_direct) |  |
| [`topology:core:has_active_descendant_transitive`](#rule-topology-core-has_active_descendant_transitive) |  |

**`preparation:no_complex`** (mandatory)

The complex-dissolving modes (`normal-no-complex`, `keep-species-no-complex`): a complex with any (transitive) active descendant is deleted; a non-deleted top-level species with activity is keyed by `keptSpeciesKey(SELF)`; a subunit of a deleted complex is promoted to top level, keyed by `promotedSubunitKey(SELF)` (paths reach it via `paths:complex_traversal`).

| Rule | Variant |
| --- | --- |
| [`preparation:no_complex:deleted`](#rule-preparation-no_complex-deleted) |  |
| [`preparation:no_complex:subunit_of_deleted_direct`](#rule-preparation-no_complex-subunit_of_deleted_direct) |  |
| [`preparation:no_complex:subunit_of_deleted_transitive`](#rule-preparation-no_complex-subunit_of_deleted_transitive) |  |
| [`preparation:no_complex:key_top_level`](#rule-preparation-no_complex-key_top_level) |  |
| [`preparation:no_complex:key_promoted_subunit`](#rule-preparation-no_complex-key_promoted_subunit) |  |
| [`preparation:no_complex:celldesigner:carrier`](#rule-preparation-no_complex-celldesigner-carrier) | CellDesigner |
| [`preparation:no_complex:sbgn_pd:carrier_entity_pool`](#rule-preparation-no_complex-sbgn_pd-carrier_entity_pool) | SBGN PD |
| [`preparation:no_complex:sbgn_pd:carrier_phenotype`](#rule-preparation-no_complex-sbgn_pd-carrier_phenotype) | SBGN PD |
| [`preparation:no_complex:sbgn_pd:carrier_subunit`](#rule-preparation-no_complex-sbgn_pd-carrier_subunit) | SBGN PD |

**`paths:core`** (mandatory)

Builds the direct (single-hop) kinded `propagatesInfluence(SOURCE_SPECIES, TARGET_SPECIES, INFLUENCE_KIND)` relation from PD reactions and modulation arcs, plus the production-cycle relations that gate transitive extension.

| Rule | Variant |
| --- | --- |
| [`paths:core:is_transformed_to`](#rule-paths-core-is_transformed_to) |  |
| [`paths:core:is_cyclically_transformed_to`](#rule-paths-core-is_cyclically_transformed_to) |  |
| [`paths:core:composes_to`](#rule-paths-core-composes_to) |  |
| [`paths:core:celldesigner:catalyzer_to_product`](#rule-paths-core-celldesigner-catalyzer_to_product) | CellDesigner |
| [`paths:core:celldesigner:physical_stimulator_to_product`](#rule-paths-core-celldesigner-physical_stimulator_to_product) | CellDesigner |
| [`paths:core:celldesigner:trigger_to_product`](#rule-paths-core-celldesigner-trigger_to_product) | CellDesigner |
| [`paths:core:celldesigner:modulation_arc_influence`](#rule-paths-core-celldesigner-modulation_arc_influence) | CellDesigner |
| [`paths:core:celldesigner:inhibitor_to_product`](#rule-paths-core-celldesigner-inhibitor_to_product) | CellDesigner |
| [`paths:core:celldesigner:modulator_to_product`](#rule-paths-core-celldesigner-modulator_to_product) | CellDesigner |
| [`paths:core:celldesigner:unknown_catalyzer_to_product`](#rule-paths-core-celldesigner-unknown_catalyzer_to_product) | CellDesigner |
| [`paths:core:celldesigner:unknown_inhibitor_to_product`](#rule-paths-core-celldesigner-unknown_inhibitor_to_product) | CellDesigner |
| [`paths:core:celldesigner:is_directly_transformed_to`](#rule-paths-core-celldesigner-is_directly_transformed_to) | CellDesigner |
| [`paths:core:sbgn_pd:modulation_to_product`](#rule-paths-core-sbgn_pd-modulation_to_product) | SBGN PD |
| [`paths:core:sbgn_pd:modulation_to_phenotype`](#rule-paths-core-sbgn_pd-modulation_to_phenotype) | SBGN PD |
| [`paths:core:sbgn_pd:is_directly_transformed_to`](#rule-paths-core-sbgn_pd-is_directly_transformed_to) | SBGN PD |

**`paths:chaining`** (excludable)

Excludable multi-hop transitivity: extends a `propagatesInfluence` path through one more reactant->product hop (a reaction in CellDesigner, a process in SBGN-PD), carrying the kind via `composesTo` (triggering degrades to positivelyInfluences) and refusing to extend across a hop inside a production cycle (`not isCyclicallyTransformedTo`).

| Rule | Variant |
| --- | --- |
| [`paths:chaining:celldesigner:transitive_through_reaction`](#rule-paths-chaining-celldesigner-transitive_through_reaction) | CellDesigner |
| [`paths:chaining:sbgn_pd:transitive_through_process`](#rule-paths-chaining-sbgn_pd-transitive_through_process) | SBGN PD |

**`paths:complex_traversal`** (excludable)

Extends paths through complex containment for the ``*-no-complex`` modes: a path touching a complex is propagated to/from each of its subunits so influences reach the surviving subunit activities.

| Rule | Variant |
| --- | --- |
| [`paths:complex_traversal:into_subunits`](#rule-paths-complex_traversal-into_subunits) |  |
| [`paths:complex_traversal:from_subunits`](#rule-paths-complex_traversal-from_subunits) |  |

**`influences:kind`** (mandatory)

Maps each modulation arc to its influence kind via `hasInfluenceKind(MODULATION, INFLUENCE_KIND)`, the single place the arc-type->kind knowledge lives.

| Rule | Variant |
| --- | --- |
| [`influences:kind:celldesigner:catalysis`](#rule-influences-kind-celldesigner-catalysis) | CellDesigner |
| [`influences:kind:celldesigner:physical_stimulation`](#rule-influences-kind-celldesigner-physical_stimulation) | CellDesigner |
| [`influences:kind:celldesigner:positive_influence`](#rule-influences-kind-celldesigner-positive_influence) | CellDesigner |
| [`influences:kind:celldesigner:triggering`](#rule-influences-kind-celldesigner-triggering) | CellDesigner |
| [`influences:kind:celldesigner:inhibition`](#rule-influences-kind-celldesigner-inhibition) | CellDesigner |
| [`influences:kind:celldesigner:negative_influence`](#rule-influences-kind-celldesigner-negative_influence) | CellDesigner |
| [`influences:kind:celldesigner:modulation`](#rule-influences-kind-celldesigner-modulation) | CellDesigner |
| [`influences:kind:celldesigner:unknown_catalysis`](#rule-influences-kind-celldesigner-unknown_catalysis) | CellDesigner |
| [`influences:kind:celldesigner:unknown_physical_stimulation`](#rule-influences-kind-celldesigner-unknown_physical_stimulation) | CellDesigner |
| [`influences:kind:celldesigner:unknown_positive_influence`](#rule-influences-kind-celldesigner-unknown_positive_influence) | CellDesigner |
| [`influences:kind:celldesigner:unknown_triggering`](#rule-influences-kind-celldesigner-unknown_triggering) | CellDesigner |
| [`influences:kind:celldesigner:unknown_inhibition`](#rule-influences-kind-celldesigner-unknown_inhibition) | CellDesigner |
| [`influences:kind:celldesigner:unknown_negative_influence`](#rule-influences-kind-celldesigner-unknown_negative_influence) | CellDesigner |
| [`influences:kind:celldesigner:unknown_modulation`](#rule-influences-kind-celldesigner-unknown_modulation) | CellDesigner |
| [`influences:kind:sbgn_pd:necessary_stimulation`](#rule-influences-kind-sbgn_pd-necessary_stimulation) | SBGN PD |
| [`influences:kind:sbgn_pd:stimulation`](#rule-influences-kind-sbgn_pd-stimulation) | SBGN PD |
| [`influences:kind:sbgn_pd:inhibition`](#rule-influences-kind-sbgn_pd-inhibition) | SBGN PD |
| [`influences:kind:sbgn_pd:modulation`](#rule-influences-kind-sbgn_pd-modulation) | SBGN PD |

**`influences:core`** (excludable)

Derivation: emits `new(activity(KEY))` for every activity key, and lifts every kinded path into the internal `influences(SOURCE_KEY, TARGET_KEY, INFLUENCE_KIND)` relation by resolving both endpoints through their activity carrier and key.

| Rule | Variant |
| --- | --- |
| [`influences:core:activity`](#rule-influences-core-activity) |  |
| [`influences:core:path`](#rule-influences-core-path) |  |

**`influences:consumption`** (excludable)

Excludable consumption/sparing reasoning: a reaction depletes its reactants, so a modifier that drives the reaction also acts on every reactant that is itself an activity -- catalyzer/physicalStimulator/trigger negatively influence each consumed reactant, inhibitor positively influences each spared reactant, and the unknown modifiers contribute the unknown twins.

| Rule | Variant |
| --- | --- |
| [`influences:consumption:celldesigner:catalyzer_consumes_reactant`](#rule-influences-consumption-celldesigner-catalyzer_consumes_reactant) | CellDesigner |
| [`influences:consumption:celldesigner:physical_stimulator_consumes_reactant`](#rule-influences-consumption-celldesigner-physical_stimulator_consumes_reactant) | CellDesigner |
| [`influences:consumption:celldesigner:trigger_consumes_reactant`](#rule-influences-consumption-celldesigner-trigger_consumes_reactant) | CellDesigner |
| [`influences:consumption:celldesigner:inhibitor_spares_reactant`](#rule-influences-consumption-celldesigner-inhibitor_spares_reactant) | CellDesigner |
| [`influences:consumption:celldesigner:unknown_catalyzer_consumes_reactant`](#rule-influences-consumption-celldesigner-unknown_catalyzer_consumes_reactant) | CellDesigner |
| [`influences:consumption:celldesigner:unknown_inhibitor_spares_reactant`](#rule-influences-consumption-celldesigner-unknown_inhibitor_spares_reactant) | CellDesigner |

**`influences:output`** (excludable)

Shared fan-out from the internal `influences(SOURCE, TARGET, INFLUENCE_KIND)` relation to the typed `new(...)` influence heads — one rule per kind.

| Rule | Variant |
| --- | --- |
| [`influences:output:positive`](#rule-influences-output-positive) |  |
| [`influences:output:negative`](#rule-influences-output-negative) |  |
| [`influences:output:modulation`](#rule-influences-output-modulation) |  |
| [`influences:output:triggering`](#rule-influences-output-triggering) |  |
| [`influences:output:unknown_positive`](#rule-influences-output-unknown_positive) |  |
| [`influences:output:unknown_negative`](#rule-influences-output-unknown_negative) |  |
| [`influences:output:unknown_modulation`](#rule-influences-output-unknown_modulation) |  |
| [`influences:output:unknown_triggering`](#rule-influences-output-unknown_triggering) |  |

**`gates:core`** (excludable)

Authored logical operators (CellDesigner `BooleanLogicGate`, SBGN-PD `LogicalOperator`): emits the operator node (`logicalOperator/2`, token-typed), its input edges (`logicalOperatorInput/2`, each input resolved through carrier/key), and the operator-sourced influence written straight into the internal `influences/3` relation by reusing the existing `propagatesInfluence(OPERATOR, TARGET, INFLUENCE_KIND)` closure.

| Rule | Variant |
| --- | --- |
| [`gates:core:celldesigner:node_and`](#rule-gates-core-celldesigner-node_and) | CellDesigner |
| [`gates:core:celldesigner:node_or`](#rule-gates-core-celldesigner-node_or) | CellDesigner |
| [`gates:core:celldesigner:node_not`](#rule-gates-core-celldesigner-node_not) | CellDesigner |
| [`gates:core:celldesigner:node_unknown`](#rule-gates-core-celldesigner-node_unknown) | CellDesigner |
| [`gates:core:celldesigner:input_edge`](#rule-gates-core-celldesigner-input_edge) | CellDesigner |
| [`gates:core:celldesigner:influence`](#rule-gates-core-celldesigner-influence) | CellDesigner |
| [`gates:core:sbgn_pd:node_and`](#rule-gates-core-sbgn_pd-node_and) | SBGN PD |
| [`gates:core:sbgn_pd:node_or`](#rule-gates-core-sbgn_pd-node_or) | SBGN PD |
| [`gates:core:sbgn_pd:node_not`](#rule-gates-core-sbgn_pd-node_not) | SBGN PD |
| [`gates:core:sbgn_pd:input_edge`](#rule-gates-core-sbgn_pd-input_edge) | SBGN PD |
| [`gates:core:sbgn_pd:influence`](#rule-gates-core-sbgn_pd-influence) | SBGN PD |

### `keep-species`

**`activity:core`** (mandatory)

Mandatory activity machinery: the candidate->activity bridge and the two global toggles.

| Rule | Variant |
| --- | --- |
| [`activity:core:from_candidate`](#rule-activity-core-from_candidate) |  |
| [`activity:core:from_global_suppress`](#rule-activity-core-from_global_suppress) |  |
| [`activity:core:celldesigner:from_global_activate`](#rule-activity-core-celldesigner-from_global_activate) | CellDesigner |
| [`activity:core:sbgn_pd:from_global_activate`](#rule-activity-core-sbgn_pd-from_global_activate) | SBGN PD |

**`activity:phenotype`** (excludable)

Excludable: a phenotype is an activity candidate (reason `isPhenotype`).

| Rule | Variant |
| --- | --- |
| [`activity:phenotype:from_phenotype`](#rule-activity-phenotype-from_phenotype) |  |

**`activity:active_marker`** (excludable)

Excludable: an explicit active marker makes a species/entity pool an activity candidate (CellDesigner: `hasActive` flag or active structural state; SBGN-PD: active state variable, including on subunits).

| Rule | Variant |
| --- | --- |
| [`activity:active_marker:celldesigner:from_active_flag`](#rule-activity-active_marker-celldesigner-from_active_flag) | CellDesigner |
| [`activity:active_marker:celldesigner:from_active_structural_state`](#rule-activity-active_marker-celldesigner-from_active_structural_state) | CellDesigner |
| [`activity:active_marker:sbgn_pd:from_active_state_variable`](#rule-activity-active_marker-sbgn_pd-from_active_state_variable) | SBGN PD |
| [`activity:active_marker:sbgn_pd:from_active_subunit_state_variable`](#rule-activity-active_marker-sbgn_pd-from_active_subunit_state_variable) | SBGN PD |

**`activity:modulation_source`** (excludable)

Excludable: the source of a modulation arc (CellDesigner: a modulation arc or a reaction modifier, which in SBGN-PD folds into the same arc-to-process case) is an activity candidate.

| Rule | Variant |
| --- | --- |
| [`activity:modulation_source:celldesigner:from_modulation_source`](#rule-activity-modulation_source-celldesigner-from_modulation_source) | CellDesigner |
| [`activity:modulation_source:celldesigner:from_reaction_modulator`](#rule-activity-modulation_source-celldesigner-from_reaction_modulator) | CellDesigner |
| [`activity:modulation_source:sbgn_pd:from_modulation_source`](#rule-activity-modulation_source-sbgn_pd-from_modulation_source) | SBGN PD |

**`activity:gate_input`** (excludable)

Excludable: an input feeding a boolean logic gate / logical operator is an activity candidate (reason `isGateInput`).

| Rule | Variant |
| --- | --- |
| [`activity:gate_input:celldesigner:from_gate_input`](#rule-activity-gate_input-celldesigner-from_gate_input) | CellDesigner |
| [`activity:gate_input:sbgn_pd:from_operator_input`](#rule-activity-gate_input-sbgn_pd-from_operator_input) | SBGN PD |

**`topology:core`** (mandatory)

Mode-agnostic structural helpers shared by every mode: `isSubunit`, `hasActiveDescendantSubunit`.

| Rule | Variant |
| --- | --- |
| [`topology:core:is_subunit`](#rule-topology-core-is_subunit) |  |
| [`topology:core:has_active_descendant_direct`](#rule-topology-core-has_active_descendant_direct) |  |
| [`topology:core:has_active_descendant_transitive`](#rule-topology-core-has_active_descendant_transitive) |  |

**`topology:top_level`** (mandatory)

Resolves every species to its outermost top-level entity: a non-subunit resolves to itself; a subunit -- at any nesting depth -- resolves to the outermost complex that contains it.

| Rule | Variant |
| --- | --- |
| [`topology:top_level:recursive`](#rule-topology-top_level-recursive) |  |
| [`topology:top_level:celldesigner:self`](#rule-topology-top_level-celldesigner-self) | CellDesigner |
| [`topology:top_level:sbgn_pd:self`](#rule-topology-top_level-sbgn_pd-self) | SBGN PD |
| [`topology:top_level:sbgn_pd:phenotype_self`](#rule-topology-top_level-sbgn_pd-phenotype_self) | SBGN PD |

**`preparation:complex`** (mandatory)

The complex-keeping modes (`normal`, `keep-species`, `keep-reactions`) key a species with activity by the `keptSpeciesKey` of its top-level entity (the `topology:top_level` group): itself when top-level, its outermost complex when a subunit.

| Rule | Variant |
| --- | --- |
| [`preparation:complex:key`](#rule-preparation-complex-key) |  |
| [`preparation:complex:celldesigner:carrier`](#rule-preparation-complex-celldesigner-carrier) | CellDesigner |
| [`preparation:complex:sbgn_pd:carrier_entity_pool`](#rule-preparation-complex-sbgn_pd-carrier_entity_pool) | SBGN PD |
| [`preparation:complex:sbgn_pd:carrier_phenotype`](#rule-preparation-complex-sbgn_pd-carrier_phenotype) | SBGN PD |

**`paths:core`** (mandatory)

Builds the direct (single-hop) kinded `propagatesInfluence(SOURCE_SPECIES, TARGET_SPECIES, INFLUENCE_KIND)` relation from PD reactions and modulation arcs, plus the production-cycle relations that gate transitive extension.

| Rule | Variant |
| --- | --- |
| [`paths:core:is_transformed_to`](#rule-paths-core-is_transformed_to) |  |
| [`paths:core:is_cyclically_transformed_to`](#rule-paths-core-is_cyclically_transformed_to) |  |
| [`paths:core:composes_to`](#rule-paths-core-composes_to) |  |
| [`paths:core:celldesigner:catalyzer_to_product`](#rule-paths-core-celldesigner-catalyzer_to_product) | CellDesigner |
| [`paths:core:celldesigner:physical_stimulator_to_product`](#rule-paths-core-celldesigner-physical_stimulator_to_product) | CellDesigner |
| [`paths:core:celldesigner:trigger_to_product`](#rule-paths-core-celldesigner-trigger_to_product) | CellDesigner |
| [`paths:core:celldesigner:modulation_arc_influence`](#rule-paths-core-celldesigner-modulation_arc_influence) | CellDesigner |
| [`paths:core:celldesigner:inhibitor_to_product`](#rule-paths-core-celldesigner-inhibitor_to_product) | CellDesigner |
| [`paths:core:celldesigner:modulator_to_product`](#rule-paths-core-celldesigner-modulator_to_product) | CellDesigner |
| [`paths:core:celldesigner:unknown_catalyzer_to_product`](#rule-paths-core-celldesigner-unknown_catalyzer_to_product) | CellDesigner |
| [`paths:core:celldesigner:unknown_inhibitor_to_product`](#rule-paths-core-celldesigner-unknown_inhibitor_to_product) | CellDesigner |
| [`paths:core:celldesigner:is_directly_transformed_to`](#rule-paths-core-celldesigner-is_directly_transformed_to) | CellDesigner |
| [`paths:core:sbgn_pd:modulation_to_product`](#rule-paths-core-sbgn_pd-modulation_to_product) | SBGN PD |
| [`paths:core:sbgn_pd:modulation_to_phenotype`](#rule-paths-core-sbgn_pd-modulation_to_phenotype) | SBGN PD |
| [`paths:core:sbgn_pd:is_directly_transformed_to`](#rule-paths-core-sbgn_pd-is_directly_transformed_to) | SBGN PD |

**`paths:chaining`** (excludable)

Excludable multi-hop transitivity: extends a `propagatesInfluence` path through one more reactant->product hop (a reaction in CellDesigner, a process in SBGN-PD), carrying the kind via `composesTo` (triggering degrades to positivelyInfluences) and refusing to extend across a hop inside a production cycle (`not isCyclicallyTransformedTo`).

| Rule | Variant |
| --- | --- |
| [`paths:chaining:celldesigner:transitive_through_reaction`](#rule-paths-chaining-celldesigner-transitive_through_reaction) | CellDesigner |
| [`paths:chaining:sbgn_pd:transitive_through_process`](#rule-paths-chaining-sbgn_pd-transitive_through_process) | SBGN PD |

**`influences:kind`** (mandatory)

Maps each modulation arc to its influence kind via `hasInfluenceKind(MODULATION, INFLUENCE_KIND)`, the single place the arc-type->kind knowledge lives.

| Rule | Variant |
| --- | --- |
| [`influences:kind:celldesigner:catalysis`](#rule-influences-kind-celldesigner-catalysis) | CellDesigner |
| [`influences:kind:celldesigner:physical_stimulation`](#rule-influences-kind-celldesigner-physical_stimulation) | CellDesigner |
| [`influences:kind:celldesigner:positive_influence`](#rule-influences-kind-celldesigner-positive_influence) | CellDesigner |
| [`influences:kind:celldesigner:triggering`](#rule-influences-kind-celldesigner-triggering) | CellDesigner |
| [`influences:kind:celldesigner:inhibition`](#rule-influences-kind-celldesigner-inhibition) | CellDesigner |
| [`influences:kind:celldesigner:negative_influence`](#rule-influences-kind-celldesigner-negative_influence) | CellDesigner |
| [`influences:kind:celldesigner:modulation`](#rule-influences-kind-celldesigner-modulation) | CellDesigner |
| [`influences:kind:celldesigner:unknown_catalysis`](#rule-influences-kind-celldesigner-unknown_catalysis) | CellDesigner |
| [`influences:kind:celldesigner:unknown_physical_stimulation`](#rule-influences-kind-celldesigner-unknown_physical_stimulation) | CellDesigner |
| [`influences:kind:celldesigner:unknown_positive_influence`](#rule-influences-kind-celldesigner-unknown_positive_influence) | CellDesigner |
| [`influences:kind:celldesigner:unknown_triggering`](#rule-influences-kind-celldesigner-unknown_triggering) | CellDesigner |
| [`influences:kind:celldesigner:unknown_inhibition`](#rule-influences-kind-celldesigner-unknown_inhibition) | CellDesigner |
| [`influences:kind:celldesigner:unknown_negative_influence`](#rule-influences-kind-celldesigner-unknown_negative_influence) | CellDesigner |
| [`influences:kind:celldesigner:unknown_modulation`](#rule-influences-kind-celldesigner-unknown_modulation) | CellDesigner |
| [`influences:kind:sbgn_pd:necessary_stimulation`](#rule-influences-kind-sbgn_pd-necessary_stimulation) | SBGN PD |
| [`influences:kind:sbgn_pd:stimulation`](#rule-influences-kind-sbgn_pd-stimulation) | SBGN PD |
| [`influences:kind:sbgn_pd:inhibition`](#rule-influences-kind-sbgn_pd-inhibition) | SBGN PD |
| [`influences:kind:sbgn_pd:modulation`](#rule-influences-kind-sbgn_pd-modulation) | SBGN PD |

**`influences:core`** (excludable)

Derivation: emits `new(activity(KEY))` for every activity key, and lifts every kinded path into the internal `influences(SOURCE_KEY, TARGET_KEY, INFLUENCE_KIND)` relation by resolving both endpoints through their activity carrier and key.

| Rule | Variant |
| --- | --- |
| [`influences:core:activity`](#rule-influences-core-activity) |  |
| [`influences:core:path`](#rule-influences-core-path) |  |

**`influences:consumption`** (excludable)

Excludable consumption/sparing reasoning: a reaction depletes its reactants, so a modifier that drives the reaction also acts on every reactant that is itself an activity -- catalyzer/physicalStimulator/trigger negatively influence each consumed reactant, inhibitor positively influences each spared reactant, and the unknown modifiers contribute the unknown twins.

| Rule | Variant |
| --- | --- |
| [`influences:consumption:celldesigner:catalyzer_consumes_reactant`](#rule-influences-consumption-celldesigner-catalyzer_consumes_reactant) | CellDesigner |
| [`influences:consumption:celldesigner:physical_stimulator_consumes_reactant`](#rule-influences-consumption-celldesigner-physical_stimulator_consumes_reactant) | CellDesigner |
| [`influences:consumption:celldesigner:trigger_consumes_reactant`](#rule-influences-consumption-celldesigner-trigger_consumes_reactant) | CellDesigner |
| [`influences:consumption:celldesigner:inhibitor_spares_reactant`](#rule-influences-consumption-celldesigner-inhibitor_spares_reactant) | CellDesigner |
| [`influences:consumption:celldesigner:unknown_catalyzer_consumes_reactant`](#rule-influences-consumption-celldesigner-unknown_catalyzer_consumes_reactant) | CellDesigner |
| [`influences:consumption:celldesigner:unknown_inhibitor_spares_reactant`](#rule-influences-consumption-celldesigner-unknown_inhibitor_spares_reactant) | CellDesigner |

**`influences:output`** (excludable)

Shared fan-out from the internal `influences(SOURCE, TARGET, INFLUENCE_KIND)` relation to the typed `new(...)` influence heads — one rule per kind.

| Rule | Variant |
| --- | --- |
| [`influences:output:positive`](#rule-influences-output-positive) |  |
| [`influences:output:negative`](#rule-influences-output-negative) |  |
| [`influences:output:modulation`](#rule-influences-output-modulation) |  |
| [`influences:output:triggering`](#rule-influences-output-triggering) |  |
| [`influences:output:unknown_positive`](#rule-influences-output-unknown_positive) |  |
| [`influences:output:unknown_negative`](#rule-influences-output-unknown_negative) |  |
| [`influences:output:unknown_modulation`](#rule-influences-output-unknown_modulation) |  |
| [`influences:output:unknown_triggering`](#rule-influences-output-unknown_triggering) |  |

**`gates:core`** (excludable)

Authored logical operators (CellDesigner `BooleanLogicGate`, SBGN-PD `LogicalOperator`): emits the operator node (`logicalOperator/2`, token-typed), its input edges (`logicalOperatorInput/2`, each input resolved through carrier/key), and the operator-sourced influence written straight into the internal `influences/3` relation by reusing the existing `propagatesInfluence(OPERATOR, TARGET, INFLUENCE_KIND)` closure.

| Rule | Variant |
| --- | --- |
| [`gates:core:celldesigner:node_and`](#rule-gates-core-celldesigner-node_and) | CellDesigner |
| [`gates:core:celldesigner:node_or`](#rule-gates-core-celldesigner-node_or) | CellDesigner |
| [`gates:core:celldesigner:node_not`](#rule-gates-core-celldesigner-node_not) | CellDesigner |
| [`gates:core:celldesigner:node_unknown`](#rule-gates-core-celldesigner-node_unknown) | CellDesigner |
| [`gates:core:celldesigner:input_edge`](#rule-gates-core-celldesigner-input_edge) | CellDesigner |
| [`gates:core:celldesigner:influence`](#rule-gates-core-celldesigner-influence) | CellDesigner |
| [`gates:core:sbgn_pd:node_and`](#rule-gates-core-sbgn_pd-node_and) | SBGN PD |
| [`gates:core:sbgn_pd:node_or`](#rule-gates-core-sbgn_pd-node_or) | SBGN PD |
| [`gates:core:sbgn_pd:node_not`](#rule-gates-core-sbgn_pd-node_not) | SBGN PD |
| [`gates:core:sbgn_pd:input_edge`](#rule-gates-core-sbgn_pd-input_edge) | SBGN PD |
| [`gates:core:sbgn_pd:influence`](#rule-gates-core-sbgn_pd-influence) | SBGN PD |

### `keep-species-no-complex`

**`activity:core`** (mandatory)

Mandatory activity machinery: the candidate->activity bridge and the two global toggles.

| Rule | Variant |
| --- | --- |
| [`activity:core:from_candidate`](#rule-activity-core-from_candidate) |  |
| [`activity:core:from_global_suppress`](#rule-activity-core-from_global_suppress) |  |
| [`activity:core:celldesigner:from_global_activate`](#rule-activity-core-celldesigner-from_global_activate) | CellDesigner |
| [`activity:core:sbgn_pd:from_global_activate`](#rule-activity-core-sbgn_pd-from_global_activate) | SBGN PD |

**`activity:phenotype`** (excludable)

Excludable: a phenotype is an activity candidate (reason `isPhenotype`).

| Rule | Variant |
| --- | --- |
| [`activity:phenotype:from_phenotype`](#rule-activity-phenotype-from_phenotype) |  |

**`activity:active_marker`** (excludable)

Excludable: an explicit active marker makes a species/entity pool an activity candidate (CellDesigner: `hasActive` flag or active structural state; SBGN-PD: active state variable, including on subunits).

| Rule | Variant |
| --- | --- |
| [`activity:active_marker:celldesigner:from_active_flag`](#rule-activity-active_marker-celldesigner-from_active_flag) | CellDesigner |
| [`activity:active_marker:celldesigner:from_active_structural_state`](#rule-activity-active_marker-celldesigner-from_active_structural_state) | CellDesigner |
| [`activity:active_marker:sbgn_pd:from_active_state_variable`](#rule-activity-active_marker-sbgn_pd-from_active_state_variable) | SBGN PD |
| [`activity:active_marker:sbgn_pd:from_active_subunit_state_variable`](#rule-activity-active_marker-sbgn_pd-from_active_subunit_state_variable) | SBGN PD |

**`activity:modulation_source`** (excludable)

Excludable: the source of a modulation arc (CellDesigner: a modulation arc or a reaction modifier, which in SBGN-PD folds into the same arc-to-process case) is an activity candidate.

| Rule | Variant |
| --- | --- |
| [`activity:modulation_source:celldesigner:from_modulation_source`](#rule-activity-modulation_source-celldesigner-from_modulation_source) | CellDesigner |
| [`activity:modulation_source:celldesigner:from_reaction_modulator`](#rule-activity-modulation_source-celldesigner-from_reaction_modulator) | CellDesigner |
| [`activity:modulation_source:sbgn_pd:from_modulation_source`](#rule-activity-modulation_source-sbgn_pd-from_modulation_source) | SBGN PD |

**`activity:gate_input`** (excludable)

Excludable: an input feeding a boolean logic gate / logical operator is an activity candidate (reason `isGateInput`).

| Rule | Variant |
| --- | --- |
| [`activity:gate_input:celldesigner:from_gate_input`](#rule-activity-gate_input-celldesigner-from_gate_input) | CellDesigner |
| [`activity:gate_input:sbgn_pd:from_operator_input`](#rule-activity-gate_input-sbgn_pd-from_operator_input) | SBGN PD |

**`topology:core`** (mandatory)

Mode-agnostic structural helpers shared by every mode: `isSubunit`, `hasActiveDescendantSubunit`.

| Rule | Variant |
| --- | --- |
| [`topology:core:is_subunit`](#rule-topology-core-is_subunit) |  |
| [`topology:core:has_active_descendant_direct`](#rule-topology-core-has_active_descendant_direct) |  |
| [`topology:core:has_active_descendant_transitive`](#rule-topology-core-has_active_descendant_transitive) |  |

**`preparation:no_complex`** (mandatory)

The complex-dissolving modes (`normal-no-complex`, `keep-species-no-complex`): a complex with any (transitive) active descendant is deleted; a non-deleted top-level species with activity is keyed by `keptSpeciesKey(SELF)`; a subunit of a deleted complex is promoted to top level, keyed by `promotedSubunitKey(SELF)` (paths reach it via `paths:complex_traversal`).

| Rule | Variant |
| --- | --- |
| [`preparation:no_complex:deleted`](#rule-preparation-no_complex-deleted) |  |
| [`preparation:no_complex:subunit_of_deleted_direct`](#rule-preparation-no_complex-subunit_of_deleted_direct) |  |
| [`preparation:no_complex:subunit_of_deleted_transitive`](#rule-preparation-no_complex-subunit_of_deleted_transitive) |  |
| [`preparation:no_complex:key_top_level`](#rule-preparation-no_complex-key_top_level) |  |
| [`preparation:no_complex:key_promoted_subunit`](#rule-preparation-no_complex-key_promoted_subunit) |  |
| [`preparation:no_complex:celldesigner:carrier`](#rule-preparation-no_complex-celldesigner-carrier) | CellDesigner |
| [`preparation:no_complex:sbgn_pd:carrier_entity_pool`](#rule-preparation-no_complex-sbgn_pd-carrier_entity_pool) | SBGN PD |
| [`preparation:no_complex:sbgn_pd:carrier_phenotype`](#rule-preparation-no_complex-sbgn_pd-carrier_phenotype) | SBGN PD |
| [`preparation:no_complex:sbgn_pd:carrier_subunit`](#rule-preparation-no_complex-sbgn_pd-carrier_subunit) | SBGN PD |

**`paths:core`** (mandatory)

Builds the direct (single-hop) kinded `propagatesInfluence(SOURCE_SPECIES, TARGET_SPECIES, INFLUENCE_KIND)` relation from PD reactions and modulation arcs, plus the production-cycle relations that gate transitive extension.

| Rule | Variant |
| --- | --- |
| [`paths:core:is_transformed_to`](#rule-paths-core-is_transformed_to) |  |
| [`paths:core:is_cyclically_transformed_to`](#rule-paths-core-is_cyclically_transformed_to) |  |
| [`paths:core:composes_to`](#rule-paths-core-composes_to) |  |
| [`paths:core:celldesigner:catalyzer_to_product`](#rule-paths-core-celldesigner-catalyzer_to_product) | CellDesigner |
| [`paths:core:celldesigner:physical_stimulator_to_product`](#rule-paths-core-celldesigner-physical_stimulator_to_product) | CellDesigner |
| [`paths:core:celldesigner:trigger_to_product`](#rule-paths-core-celldesigner-trigger_to_product) | CellDesigner |
| [`paths:core:celldesigner:modulation_arc_influence`](#rule-paths-core-celldesigner-modulation_arc_influence) | CellDesigner |
| [`paths:core:celldesigner:inhibitor_to_product`](#rule-paths-core-celldesigner-inhibitor_to_product) | CellDesigner |
| [`paths:core:celldesigner:modulator_to_product`](#rule-paths-core-celldesigner-modulator_to_product) | CellDesigner |
| [`paths:core:celldesigner:unknown_catalyzer_to_product`](#rule-paths-core-celldesigner-unknown_catalyzer_to_product) | CellDesigner |
| [`paths:core:celldesigner:unknown_inhibitor_to_product`](#rule-paths-core-celldesigner-unknown_inhibitor_to_product) | CellDesigner |
| [`paths:core:celldesigner:is_directly_transformed_to`](#rule-paths-core-celldesigner-is_directly_transformed_to) | CellDesigner |
| [`paths:core:sbgn_pd:modulation_to_product`](#rule-paths-core-sbgn_pd-modulation_to_product) | SBGN PD |
| [`paths:core:sbgn_pd:modulation_to_phenotype`](#rule-paths-core-sbgn_pd-modulation_to_phenotype) | SBGN PD |
| [`paths:core:sbgn_pd:is_directly_transformed_to`](#rule-paths-core-sbgn_pd-is_directly_transformed_to) | SBGN PD |

**`paths:chaining`** (excludable)

Excludable multi-hop transitivity: extends a `propagatesInfluence` path through one more reactant->product hop (a reaction in CellDesigner, a process in SBGN-PD), carrying the kind via `composesTo` (triggering degrades to positivelyInfluences) and refusing to extend across a hop inside a production cycle (`not isCyclicallyTransformedTo`).

| Rule | Variant |
| --- | --- |
| [`paths:chaining:celldesigner:transitive_through_reaction`](#rule-paths-chaining-celldesigner-transitive_through_reaction) | CellDesigner |
| [`paths:chaining:sbgn_pd:transitive_through_process`](#rule-paths-chaining-sbgn_pd-transitive_through_process) | SBGN PD |

**`paths:complex_traversal`** (excludable)

Extends paths through complex containment for the ``*-no-complex`` modes: a path touching a complex is propagated to/from each of its subunits so influences reach the surviving subunit activities.

| Rule | Variant |
| --- | --- |
| [`paths:complex_traversal:into_subunits`](#rule-paths-complex_traversal-into_subunits) |  |
| [`paths:complex_traversal:from_subunits`](#rule-paths-complex_traversal-from_subunits) |  |

**`influences:kind`** (mandatory)

Maps each modulation arc to its influence kind via `hasInfluenceKind(MODULATION, INFLUENCE_KIND)`, the single place the arc-type->kind knowledge lives.

| Rule | Variant |
| --- | --- |
| [`influences:kind:celldesigner:catalysis`](#rule-influences-kind-celldesigner-catalysis) | CellDesigner |
| [`influences:kind:celldesigner:physical_stimulation`](#rule-influences-kind-celldesigner-physical_stimulation) | CellDesigner |
| [`influences:kind:celldesigner:positive_influence`](#rule-influences-kind-celldesigner-positive_influence) | CellDesigner |
| [`influences:kind:celldesigner:triggering`](#rule-influences-kind-celldesigner-triggering) | CellDesigner |
| [`influences:kind:celldesigner:inhibition`](#rule-influences-kind-celldesigner-inhibition) | CellDesigner |
| [`influences:kind:celldesigner:negative_influence`](#rule-influences-kind-celldesigner-negative_influence) | CellDesigner |
| [`influences:kind:celldesigner:modulation`](#rule-influences-kind-celldesigner-modulation) | CellDesigner |
| [`influences:kind:celldesigner:unknown_catalysis`](#rule-influences-kind-celldesigner-unknown_catalysis) | CellDesigner |
| [`influences:kind:celldesigner:unknown_physical_stimulation`](#rule-influences-kind-celldesigner-unknown_physical_stimulation) | CellDesigner |
| [`influences:kind:celldesigner:unknown_positive_influence`](#rule-influences-kind-celldesigner-unknown_positive_influence) | CellDesigner |
| [`influences:kind:celldesigner:unknown_triggering`](#rule-influences-kind-celldesigner-unknown_triggering) | CellDesigner |
| [`influences:kind:celldesigner:unknown_inhibition`](#rule-influences-kind-celldesigner-unknown_inhibition) | CellDesigner |
| [`influences:kind:celldesigner:unknown_negative_influence`](#rule-influences-kind-celldesigner-unknown_negative_influence) | CellDesigner |
| [`influences:kind:celldesigner:unknown_modulation`](#rule-influences-kind-celldesigner-unknown_modulation) | CellDesigner |
| [`influences:kind:sbgn_pd:necessary_stimulation`](#rule-influences-kind-sbgn_pd-necessary_stimulation) | SBGN PD |
| [`influences:kind:sbgn_pd:stimulation`](#rule-influences-kind-sbgn_pd-stimulation) | SBGN PD |
| [`influences:kind:sbgn_pd:inhibition`](#rule-influences-kind-sbgn_pd-inhibition) | SBGN PD |
| [`influences:kind:sbgn_pd:modulation`](#rule-influences-kind-sbgn_pd-modulation) | SBGN PD |

**`influences:core`** (excludable)

Derivation: emits `new(activity(KEY))` for every activity key, and lifts every kinded path into the internal `influences(SOURCE_KEY, TARGET_KEY, INFLUENCE_KIND)` relation by resolving both endpoints through their activity carrier and key.

| Rule | Variant |
| --- | --- |
| [`influences:core:activity`](#rule-influences-core-activity) |  |
| [`influences:core:path`](#rule-influences-core-path) |  |

**`influences:consumption`** (excludable)

Excludable consumption/sparing reasoning: a reaction depletes its reactants, so a modifier that drives the reaction also acts on every reactant that is itself an activity -- catalyzer/physicalStimulator/trigger negatively influence each consumed reactant, inhibitor positively influences each spared reactant, and the unknown modifiers contribute the unknown twins.

| Rule | Variant |
| --- | --- |
| [`influences:consumption:celldesigner:catalyzer_consumes_reactant`](#rule-influences-consumption-celldesigner-catalyzer_consumes_reactant) | CellDesigner |
| [`influences:consumption:celldesigner:physical_stimulator_consumes_reactant`](#rule-influences-consumption-celldesigner-physical_stimulator_consumes_reactant) | CellDesigner |
| [`influences:consumption:celldesigner:trigger_consumes_reactant`](#rule-influences-consumption-celldesigner-trigger_consumes_reactant) | CellDesigner |
| [`influences:consumption:celldesigner:inhibitor_spares_reactant`](#rule-influences-consumption-celldesigner-inhibitor_spares_reactant) | CellDesigner |
| [`influences:consumption:celldesigner:unknown_catalyzer_consumes_reactant`](#rule-influences-consumption-celldesigner-unknown_catalyzer_consumes_reactant) | CellDesigner |
| [`influences:consumption:celldesigner:unknown_inhibitor_spares_reactant`](#rule-influences-consumption-celldesigner-unknown_inhibitor_spares_reactant) | CellDesigner |

**`influences:output`** (excludable)

Shared fan-out from the internal `influences(SOURCE, TARGET, INFLUENCE_KIND)` relation to the typed `new(...)` influence heads — one rule per kind.

| Rule | Variant |
| --- | --- |
| [`influences:output:positive`](#rule-influences-output-positive) |  |
| [`influences:output:negative`](#rule-influences-output-negative) |  |
| [`influences:output:modulation`](#rule-influences-output-modulation) |  |
| [`influences:output:triggering`](#rule-influences-output-triggering) |  |
| [`influences:output:unknown_positive`](#rule-influences-output-unknown_positive) |  |
| [`influences:output:unknown_negative`](#rule-influences-output-unknown_negative) |  |
| [`influences:output:unknown_modulation`](#rule-influences-output-unknown_modulation) |  |
| [`influences:output:unknown_triggering`](#rule-influences-output-unknown_triggering) |  |

**`gates:core`** (excludable)

Authored logical operators (CellDesigner `BooleanLogicGate`, SBGN-PD `LogicalOperator`): emits the operator node (`logicalOperator/2`, token-typed), its input edges (`logicalOperatorInput/2`, each input resolved through carrier/key), and the operator-sourced influence written straight into the internal `influences/3` relation by reusing the existing `propagatesInfluence(OPERATOR, TARGET, INFLUENCE_KIND)` closure.

| Rule | Variant |
| --- | --- |
| [`gates:core:celldesigner:node_and`](#rule-gates-core-celldesigner-node_and) | CellDesigner |
| [`gates:core:celldesigner:node_or`](#rule-gates-core-celldesigner-node_or) | CellDesigner |
| [`gates:core:celldesigner:node_not`](#rule-gates-core-celldesigner-node_not) | CellDesigner |
| [`gates:core:celldesigner:node_unknown`](#rule-gates-core-celldesigner-node_unknown) | CellDesigner |
| [`gates:core:celldesigner:input_edge`](#rule-gates-core-celldesigner-input_edge) | CellDesigner |
| [`gates:core:celldesigner:influence`](#rule-gates-core-celldesigner-influence) | CellDesigner |
| [`gates:core:sbgn_pd:node_and`](#rule-gates-core-sbgn_pd-node_and) | SBGN PD |
| [`gates:core:sbgn_pd:node_or`](#rule-gates-core-sbgn_pd-node_or) | SBGN PD |
| [`gates:core:sbgn_pd:node_not`](#rule-gates-core-sbgn_pd-node_not) | SBGN PD |
| [`gates:core:sbgn_pd:input_edge`](#rule-gates-core-sbgn_pd-input_edge) | SBGN PD |
| [`gates:core:sbgn_pd:influence`](#rule-gates-core-sbgn_pd-influence) | SBGN PD |

### `keep-reactions`

**`activity:core`** (mandatory)

Mandatory activity machinery: the candidate->activity bridge and the two global toggles.

| Rule | Variant |
| --- | --- |
| [`activity:core:from_candidate`](#rule-activity-core-from_candidate) |  |
| [`activity:core:from_global_suppress`](#rule-activity-core-from_global_suppress) |  |
| [`activity:core:celldesigner:from_global_activate`](#rule-activity-core-celldesigner-from_global_activate) | CellDesigner |
| [`activity:core:sbgn_pd:from_global_activate`](#rule-activity-core-sbgn_pd-from_global_activate) | SBGN PD |

**`topology:core`** (mandatory)

Mode-agnostic structural helpers shared by every mode: `isSubunit`, `hasActiveDescendantSubunit`.

| Rule | Variant |
| --- | --- |
| [`topology:core:is_subunit`](#rule-topology-core-is_subunit) |  |
| [`topology:core:has_active_descendant_direct`](#rule-topology-core-has_active_descendant_direct) |  |
| [`topology:core:has_active_descendant_transitive`](#rule-topology-core-has_active_descendant_transitive) |  |

**`topology:top_level`** (mandatory)

Resolves every species to its outermost top-level entity: a non-subunit resolves to itself; a subunit -- at any nesting depth -- resolves to the outermost complex that contains it.

| Rule | Variant |
| --- | --- |
| [`topology:top_level:recursive`](#rule-topology-top_level-recursive) |  |
| [`topology:top_level:celldesigner:self`](#rule-topology-top_level-celldesigner-self) | CellDesigner |
| [`topology:top_level:sbgn_pd:self`](#rule-topology-top_level-sbgn_pd-self) | SBGN PD |
| [`topology:top_level:sbgn_pd:phenotype_self`](#rule-topology-top_level-sbgn_pd-phenotype_self) | SBGN PD |

**`preparation:complex`** (mandatory)

The complex-keeping modes (`normal`, `keep-species`, `keep-reactions`) key a species with activity by the `keptSpeciesKey` of its top-level entity (the `topology:top_level` group): itself when top-level, its outermost complex when a subunit.

| Rule | Variant |
| --- | --- |
| [`preparation:complex:key`](#rule-preparation-complex-key) |  |
| [`preparation:complex:celldesigner:carrier`](#rule-preparation-complex-celldesigner-carrier) | CellDesigner |
| [`preparation:complex:sbgn_pd:carrier_entity_pool`](#rule-preparation-complex-sbgn_pd-carrier_entity_pool) | SBGN PD |
| [`preparation:complex:sbgn_pd:carrier_phenotype`](#rule-preparation-complex-sbgn_pd-carrier_phenotype) | SBGN PD |

**`paths:core`** (mandatory)

Builds the direct (single-hop) kinded `propagatesInfluence(SOURCE_SPECIES, TARGET_SPECIES, INFLUENCE_KIND)` relation from PD reactions and modulation arcs, plus the production-cycle relations that gate transitive extension.

| Rule | Variant |
| --- | --- |
| [`paths:core:is_transformed_to`](#rule-paths-core-is_transformed_to) |  |
| [`paths:core:is_cyclically_transformed_to`](#rule-paths-core-is_cyclically_transformed_to) |  |
| [`paths:core:composes_to`](#rule-paths-core-composes_to) |  |
| [`paths:core:celldesigner:catalyzer_to_product`](#rule-paths-core-celldesigner-catalyzer_to_product) | CellDesigner |
| [`paths:core:celldesigner:physical_stimulator_to_product`](#rule-paths-core-celldesigner-physical_stimulator_to_product) | CellDesigner |
| [`paths:core:celldesigner:trigger_to_product`](#rule-paths-core-celldesigner-trigger_to_product) | CellDesigner |
| [`paths:core:celldesigner:modulation_arc_influence`](#rule-paths-core-celldesigner-modulation_arc_influence) | CellDesigner |
| [`paths:core:celldesigner:inhibitor_to_product`](#rule-paths-core-celldesigner-inhibitor_to_product) | CellDesigner |
| [`paths:core:celldesigner:modulator_to_product`](#rule-paths-core-celldesigner-modulator_to_product) | CellDesigner |
| [`paths:core:celldesigner:unknown_catalyzer_to_product`](#rule-paths-core-celldesigner-unknown_catalyzer_to_product) | CellDesigner |
| [`paths:core:celldesigner:unknown_inhibitor_to_product`](#rule-paths-core-celldesigner-unknown_inhibitor_to_product) | CellDesigner |
| [`paths:core:celldesigner:is_directly_transformed_to`](#rule-paths-core-celldesigner-is_directly_transformed_to) | CellDesigner |
| [`paths:core:sbgn_pd:modulation_to_product`](#rule-paths-core-sbgn_pd-modulation_to_product) | SBGN PD |
| [`paths:core:sbgn_pd:modulation_to_phenotype`](#rule-paths-core-sbgn_pd-modulation_to_phenotype) | SBGN PD |
| [`paths:core:sbgn_pd:is_directly_transformed_to`](#rule-paths-core-sbgn_pd-is_directly_transformed_to) | SBGN PD |

**`influences:kind`** (mandatory)

Maps each modulation arc to its influence kind via `hasInfluenceKind(MODULATION, INFLUENCE_KIND)`, the single place the arc-type->kind knowledge lives.

| Rule | Variant |
| --- | --- |
| [`influences:kind:celldesigner:catalysis`](#rule-influences-kind-celldesigner-catalysis) | CellDesigner |
| [`influences:kind:celldesigner:physical_stimulation`](#rule-influences-kind-celldesigner-physical_stimulation) | CellDesigner |
| [`influences:kind:celldesigner:positive_influence`](#rule-influences-kind-celldesigner-positive_influence) | CellDesigner |
| [`influences:kind:celldesigner:triggering`](#rule-influences-kind-celldesigner-triggering) | CellDesigner |
| [`influences:kind:celldesigner:inhibition`](#rule-influences-kind-celldesigner-inhibition) | CellDesigner |
| [`influences:kind:celldesigner:negative_influence`](#rule-influences-kind-celldesigner-negative_influence) | CellDesigner |
| [`influences:kind:celldesigner:modulation`](#rule-influences-kind-celldesigner-modulation) | CellDesigner |
| [`influences:kind:celldesigner:unknown_catalysis`](#rule-influences-kind-celldesigner-unknown_catalysis) | CellDesigner |
| [`influences:kind:celldesigner:unknown_physical_stimulation`](#rule-influences-kind-celldesigner-unknown_physical_stimulation) | CellDesigner |
| [`influences:kind:celldesigner:unknown_positive_influence`](#rule-influences-kind-celldesigner-unknown_positive_influence) | CellDesigner |
| [`influences:kind:celldesigner:unknown_triggering`](#rule-influences-kind-celldesigner-unknown_triggering) | CellDesigner |
| [`influences:kind:celldesigner:unknown_inhibition`](#rule-influences-kind-celldesigner-unknown_inhibition) | CellDesigner |
| [`influences:kind:celldesigner:unknown_negative_influence`](#rule-influences-kind-celldesigner-unknown_negative_influence) | CellDesigner |
| [`influences:kind:celldesigner:unknown_modulation`](#rule-influences-kind-celldesigner-unknown_modulation) | CellDesigner |
| [`influences:kind:sbgn_pd:necessary_stimulation`](#rule-influences-kind-sbgn_pd-necessary_stimulation) | SBGN PD |
| [`influences:kind:sbgn_pd:stimulation`](#rule-influences-kind-sbgn_pd-stimulation) | SBGN PD |
| [`influences:kind:sbgn_pd:inhibition`](#rule-influences-kind-sbgn_pd-inhibition) | SBGN PD |
| [`influences:kind:sbgn_pd:modulation`](#rule-influences-kind-sbgn_pd-modulation) | SBGN PD |

**`influences:core`** (mandatory)

Derivation: emits `new(activity(KEY))` for every activity key, and lifts every kinded path into the internal `influences(SOURCE_KEY, TARGET_KEY, INFLUENCE_KIND)` relation by resolving both endpoints through their activity carrier and key.

| Rule | Variant |
| --- | --- |
| [`influences:core:activity`](#rule-influences-core-activity) |  |
| [`influences:core:path`](#rule-influences-core-path) |  |

**`influences:output`** (excludable)

Shared fan-out from the internal `influences(SOURCE, TARGET, INFLUENCE_KIND)` relation to the typed `new(...)` influence heads — one rule per kind.

| Rule | Variant |
| --- | --- |
| [`influences:output:positive`](#rule-influences-output-positive) |  |
| [`influences:output:negative`](#rule-influences-output-negative) |  |
| [`influences:output:modulation`](#rule-influences-output-modulation) |  |
| [`influences:output:triggering`](#rule-influences-output-triggering) |  |
| [`influences:output:unknown_positive`](#rule-influences-output-unknown_positive) |  |
| [`influences:output:unknown_negative`](#rule-influences-output-unknown_negative) |  |
| [`influences:output:unknown_modulation`](#rule-influences-output-unknown_modulation) |  |
| [`influences:output:unknown_triggering`](#rule-influences-output-unknown_triggering) |  |

**`gates:core`** (excludable)

Authored logical operators (CellDesigner `BooleanLogicGate`, SBGN-PD `LogicalOperator`): emits the operator node (`logicalOperator/2`, token-typed), its input edges (`logicalOperatorInput/2`, each input resolved through carrier/key), and the operator-sourced influence written straight into the internal `influences/3` relation by reusing the existing `propagatesInfluence(OPERATOR, TARGET, INFLUENCE_KIND)` closure.

| Rule | Variant |
| --- | --- |
| [`gates:core:celldesigner:node_and`](#rule-gates-core-celldesigner-node_and) | CellDesigner |
| [`gates:core:celldesigner:node_or`](#rule-gates-core-celldesigner-node_or) | CellDesigner |
| [`gates:core:celldesigner:node_not`](#rule-gates-core-celldesigner-node_not) | CellDesigner |
| [`gates:core:celldesigner:node_unknown`](#rule-gates-core-celldesigner-node_unknown) | CellDesigner |
| [`gates:core:celldesigner:input_edge`](#rule-gates-core-celldesigner-input_edge) | CellDesigner |
| [`gates:core:celldesigner:influence`](#rule-gates-core-celldesigner-influence) | CellDesigner |
| [`gates:core:sbgn_pd:node_and`](#rule-gates-core-sbgn_pd-node_and) | SBGN PD |
| [`gates:core:sbgn_pd:node_or`](#rule-gates-core-sbgn_pd-node_or) | SBGN PD |
| [`gates:core:sbgn_pd:node_not`](#rule-gates-core-sbgn_pd-node_not) | SBGN PD |
| [`gates:core:sbgn_pd:input_edge`](#rule-gates-core-sbgn_pd-input_edge) | SBGN PD |
| [`gates:core:sbgn_pd:influence`](#rule-gates-core-sbgn_pd-influence) | SBGN PD |

**`keep_reactions:activity`** (excludable)

The `keep-reactions` premise that every species is an activity, replacing the structural-reason `activity:*` feature-groups (each of which is inert once every species is a candidate, so the mode omits them).

| Rule | Variant |
| --- | --- |
| [`keep_reactions:activity:from_species`](#rule-keep_reactions-activity-from_species) |  |

**`keep_reactions:influences`** (excludable)

The `keep-reactions` premise that every reaction is kept: each (reactant, product) pair of a reaction becomes one positive influence.

| Rule | Variant |
| --- | --- |
| [`keep_reactions:influences:reactant_to_product`](#rule-keep_reactions-influences-reactant_to_product) |  |

## Rule groups {#rule-groups}

### `activity` {#concern-activity}

#### `activity:core`

Mandatory activity machinery: the candidate->activity bridge and the two global toggles. A single bridging rule promotes a `hasActivityCandidate(ELEMENT, REASON)` to `hasActivity` unless the element is `suppressActivity` (the `--set-inactive` veto). `globalSuppress` (`--set-all-inactive`) suppresses every candidate except those the solver marked `forceActive` (the per-id `--set-active` override), realising the precedence per-id > global > rules; `globalActivate` (`--set-all-active`) turns every top-level species/entity pool into a candidate. Every activity feature-group and every downstream group depends on this, so the dependency graph forbids excluding it -- the global toggles stay wired to their CLI flags. The structural-reason rules (phenotype, active marker, modulation source, gate input) live in the excludable `activity:*` feature-groups.

| Rule | Variant | Documentation |
| --- | --- | --- |
| <a id="rule-activity-core-from_candidate"></a>`activity:core:from_candidate` |  | An activity candidate becomes an actual activity unless it has been vetoed. This is the single interception point for the `--set-inactive` veto (the `suppressActivity` marker; blanket suppression, so it also blocks the `--set-active` `isInputParameter` candidate, which is injected as a candidate). All downstream body references key on `hasActivity`, so they automatically respect suppression. |
| <a id="rule-activity-core-from_global_suppress"></a>`activity:core:from_global_suppress` |  | The set-all-inactive toggle suppresses every activity candidate, except those pinned active per id. This realises the precedence per-id > global > rules for the inactive toggle: `globalSuppress` (`--set-all-inactive`) blankets everything, and `forceActive` (the per-id `--set-active` override) carves out the exceptions. It keys on `hasActivityCandidate`, so it also silences subunit-derived candidates (subunits are suppressed under `--set-all-inactive`). |
| <a id="rule-activity-core-celldesigner-from_global_activate"></a>`activity:core:celldesigner:from_global_activate` | CellDesigner | Under the set-all-active toggle, every top-level species is an activity candidate. It fires with reason `isGlobalActive` when `globalActivate` (`--set-all-active`) is set. The `not hasSubunit(_, SPECIES)` guard excludes subunits (in CellDesigner a subunit is a species): a subunit is a structural component of its complex, never a top-level activity, so the toggle activates the complex, not its parts. This is the subunit asymmetry -- subunits are *not* activated under `--set-all-active`, though they *are* suppressed under `--set-all-inactive`. |
| <a id="rule-activity-core-sbgn_pd-from_global_activate"></a>`activity:core:sbgn_pd:from_global_activate` | SBGN PD | Under the set-all-active toggle, every entity pool is an activity candidate. It fires with reason `isGlobalActive` when `globalActivate` (`--set-all-active`) is set -- the SBGN PD parallel of CellDesigner's global-activate rule. No `not hasSubunit` guard is needed: in SBGN PD a subunit is an `SBGNAuxiliaryUnit`, not an `entityPool`, so the `entityPool` guard already excludes subunits. This preserves the subunit asymmetry -- subunits are not activated under `--set-all-active`. |

#### `activity:phenotype`

Excludable: a phenotype is an activity candidate (reason `isPhenotype`). Exclude with `--exclude-group activity:phenotype` to stop treating phenotypes as activities.

| Rule | Variant | Documentation |
| --- | --- | --- |
| <a id="rule-activity-phenotype-from_phenotype"></a>`activity:phenotype:from_phenotype` |  | A phenotype is always an activity candidate. Both languages emit the `phenotype` signal, so this rule is shared; it fires with reason `isPhenotype`. |

#### `activity:active_marker`

Excludable: an explicit active marker makes a species/entity pool an activity candidate (CellDesigner: `hasActive` flag or active structural state; SBGN-PD: active state variable, including on subunits). Exclude with `--exclude-group activity:active_marker`.

| Rule | Variant | Documentation |
| --- | --- | --- |
| <a id="rule-activity-active_marker-celldesigner-from_active_flag"></a>`activity:active_marker:celldesigner:from_active_flag` | CellDesigner | A species explicitly marked as active is an activity candidate. This fires when its `hasActive` flag is 1 (reason `isActive`). |
| <a id="rule-activity-active_marker-celldesigner-from_active_structural_state"></a>`activity:active_marker:celldesigner:from_active_structural_state` | CellDesigner | A species in an active structural state is an activity candidate. This fires when it carries a structural state whose value is "active" (reason `hasActiveStructuralState`). |
| <a id="rule-activity-active_marker-sbgn_pd-from_active_state_variable"></a>`activity:active_marker:sbgn_pd:from_active_state_variable` | SBGN PD | An entity pool in an active state is an activity candidate. This fires when it carries a state variable whose value is "active" (reason `hasActiveStateVariable`) -- the SBGN PD parallel of CellDesigner's active structural state. |
| <a id="rule-activity-active_marker-sbgn_pd-from_active_subunit_state_variable"></a>`activity:active_marker:sbgn_pd:from_active_subunit_state_variable` | SBGN PD | A subunit in an active state is an activity candidate. This fires when the subunit carries a state variable whose value is "active". Its complex therefore inherits activity (keep-species), and in the ``*-no-complex`` modes the subunit can be promoted. It is the parallel of CellDesigner, where subunits are species and so the active-structural-state rule already covers them. |

#### `activity:modulation_source`

Excludable: the source of a modulation arc (CellDesigner: a modulation arc or a reaction modifier, which in SBGN-PD folds into the same arc-to-process case) is an activity candidate. Exclude with `--exclude-group activity:modulation_source`.

| Rule | Variant | Documentation |
| --- | --- | --- |
| <a id="rule-activity-modulation_source-celldesigner-from_modulation_source"></a>`activity:modulation_source:celldesigner:from_modulation_source` | CellDesigner | A species that is the source of a modulation arc is an activity candidate. This covers both known and unknown modulations that have some target (reason `isModulationSource`); keying on `knownOrUnknownModulation` rather than `modulation` is what lets the source of an unknown modulation become an activity node. |
| <a id="rule-activity-modulation_source-celldesigner-from_reaction_modulator"></a>`activity:modulation_source:celldesigner:from_reaction_modulator` | CellDesigner | A species that modulates a reaction is an activity candidate. This fires when the species is referred to by a reaction modulator, known or unknown, that modifies some target reaction (reason `isReactionModulator`); keying on `knownOrUnknownModulator` rather than `modulator` is what lets an unknown catalyzer/inhibitor become an activity node. |
| <a id="rule-activity-modulation_source-sbgn_pd-from_modulation_source"></a>`activity:modulation_source:sbgn_pd:from_modulation_source` | SBGN PD | An entity pool that is the source of a modulation arc is an activity candidate. This fires when the arc's target is a process (reason `isModulationSource`) and is the SBGN PD parallel of both the CellDesigner modulation-arc and reaction-modifier rules. |

#### `activity:gate_input`

Excludable: an input feeding a boolean logic gate / logical operator is an activity candidate (reason `isGateInput`). Exclude with `--exclude-group activity:gate_input`.

| Rule | Variant | Documentation |
| --- | --- | --- |
| <a id="rule-activity-gate_input-celldesigner-from_gate_input"></a>`activity:gate_input:celldesigner:from_gate_input` | CellDesigner | A species feeding a boolean logic gate is an activity candidate. It fires with reason `isGateInput`. A gate is structurally always an influence/modulation source or reaction modifier, so each of its inputs is an active driver of the downstream target -- a semantic guarantee, not a fallback. The element is activated regardless of kind (species, complex, ion, ...); the carrier rules route each kind. `hasReferredElement` in CellDesigner only relates a gate input to its species, so no extra guard is needed. |
| <a id="rule-activity-gate_input-sbgn_pd-from_operator_input"></a>`activity:gate_input:sbgn_pd:from_operator_input` | SBGN PD | An entity pool feeding a logical operator is an activity candidate. It fires with reason `isGateInput` -- the SBGN PD parallel of CellDesigner's gate-input activation. The `entityPool` guard excludes the deferred nested-operator case: an operator feeding another operator has no activity carrier, so such an input simply dangles. In SBGN PD `hasReferredElement` relates many roles (reactant, product, modulation participant) to entities, so both the `logicalOperatorInput` and `entityPool` guards are needed. |

### `topology` {#concern-topology}

#### `topology:core`

Mode-agnostic structural helpers shared by every mode: `isSubunit`, `hasActiveDescendantSubunit`.

| Rule | Variant | Documentation |
| --- | --- | --- |
| <a id="rule-topology-core-is_subunit"></a>`topology:core:is_subunit` |  | A species is a subunit if it appears on the right side of any `hasSubunit` relation. |
| <a id="rule-topology-core-has_active_descendant_direct"></a>`topology:core:has_active_descendant_direct` |  | A complex has an active descendant if any direct subunit has activity. |
| <a id="rule-topology-core-has_active_descendant_transitive"></a>`topology:core:has_active_descendant_transitive` |  | The `hasActiveDescendantSubunit` relation is transitive through complex containment. |

#### `topology:top_level`

Resolves every species to its outermost top-level entity: a non-subunit resolves to itself; a subunit -- at any nesting depth -- resolves to the outermost complex that contains it. The complex-keeping modes (`keep-species`, `normal`, `keep-reactions`) key a species by its top-level entity, so a subunit is never its own activity and its influences attach to its top-level complex (a subunit is a structural component, not an independent influencer).

| Rule | Variant | Documentation |
| --- | --- | --- |
| <a id="rule-topology-top_level-recursive"></a>`topology:top_level:recursive` |  | A subunit resolves to the same top-level entity as its parent complex, recursively through nested complexes -- so a subunit at any depth resolves to its outermost complex. |
| <a id="rule-topology-top_level-celldesigner-self"></a>`topology:top_level:celldesigner:self` | CellDesigner | CellDesigner: a species that is not a subunit of any complex is its own top-level entity. |
| <a id="rule-topology-top_level-sbgn_pd-self"></a>`topology:top_level:sbgn_pd:self` | SBGN PD | SBGN-PD: an entity pool that is not a subunit of any complex is its own top-level entity. |
| <a id="rule-topology-top_level-sbgn_pd-phenotype_self"></a>`topology:top_level:sbgn_pd:phenotype_self` | SBGN PD | SBGN-PD: a phenotype is a process, not an entity pool, so the entity-pool self-rule never keys it; a phenotype is never a subunit, so it is always its own top-level entity. Without this an SBGN phenotype gets `hasActivity` but no activity key and is silently dropped from `keep-species`/`normal` output. The `*-no-complex` modes key via `not isSubunit`/`not delete` and already include phenotypes; CellDesigner phenotypes are species (covered by the species self-rule). |

### `preparation` {#concern-preparation}

#### `preparation:complex`

The complex-keeping modes (`normal`, `keep-species`, `keep-reactions`) key a species with activity by the `keptSpeciesKey` of its top-level entity (the `topology:top_level` group): itself when top-level, its outermost complex when a subunit. A subunit therefore contributes no activity of its own -- it is a structural component of its complex, and any influence it carries attaches to the top-level complex (the active descendant directly keys the complex, so the assembly is represented without a separate inherit-activity rule). Carriers are identity; the rerouting lives entirely in the key. `normal` and `keep-species` share these keys and differ only at the build stage, where `normal` strips PTM decorations and merges content-equal results while `keep-species` keeps the decorations.

| Rule | Variant | Documentation |
| --- | --- | --- |
| <a id="rule-preparation-complex-key"></a>`preparation:complex:key` |  | A species with activity is keyed by the `keptSpeciesKey` of its top-level entity (itself when top-level; its outermost complex when a subunit). |
| <a id="rule-preparation-complex-celldesigner-carrier"></a>`preparation:complex:celldesigner:carrier` | CellDesigner | Each species is its own activity carrier (no rerouting; a subunit endpoint is rerouted to its top-level complex by the `topology:top_level` key, not by the carrier). |
| <a id="rule-preparation-complex-sbgn_pd-carrier_entity_pool"></a>`preparation:complex:sbgn_pd:carrier_entity_pool` | SBGN PD | SBGN-PD: each entity pool is its own activity carrier. |
| <a id="rule-preparation-complex-sbgn_pd-carrier_phenotype"></a>`preparation:complex:sbgn_pd:carrier_phenotype` | SBGN PD | SBGN-PD: each phenotype process is its own activity carrier (phenotypes are activities, not entity pools). |

#### `preparation:no_complex`

The complex-dissolving modes (`normal-no-complex`, `keep-species-no-complex`): a complex with any (transitive) active descendant is deleted; a non-deleted top-level species with activity is keyed by `keptSpeciesKey(SELF)`; a subunit of a deleted complex is promoted to top level, keyed by `promotedSubunitKey(SELF)` (paths reach it via `paths:complex_traversal`). `normal-no-complex` and `keep-species-no-complex` share these keys and differ only at the build stage, where `normal-no-complex` strips PTM decorations and merges content-equal results while `keep-species-no-complex` keeps the decorations.

| Rule | Variant | Documentation |
| --- | --- | --- |
| <a id="rule-preparation-no_complex-deleted"></a>`preparation:no_complex:deleted` |  | A complex with any (transitive) active descendant is deleted. |
| <a id="rule-preparation-no_complex-subunit_of_deleted_direct"></a>`preparation:no_complex:subunit_of_deleted_direct` |  | A direct subunit of a deleted complex is itself subunit-of-deleted. |
| <a id="rule-preparation-no_complex-subunit_of_deleted_transitive"></a>`preparation:no_complex:subunit_of_deleted_transitive` |  | Subunit-of-deleted is transitive through nested complexes. |
| <a id="rule-preparation-no_complex-key_top_level"></a>`preparation:no_complex:key_top_level` |  | A non-deleted top-level species with activity is keyed by `keptSpeciesKey(SELF)`. |
| <a id="rule-preparation-no_complex-key_promoted_subunit"></a>`preparation:no_complex:key_promoted_subunit` |  | A subunit of a deleted complex with activity is promoted to top level, keyed by `promotedSubunitKey(SELF)`. |
| <a id="rule-preparation-no_complex-celldesigner-carrier"></a>`preparation:no_complex:celldesigner:carrier` | CellDesigner | Each species is its own activity carrier. |
| <a id="rule-preparation-no_complex-sbgn_pd-carrier_entity_pool"></a>`preparation:no_complex:sbgn_pd:carrier_entity_pool` | SBGN PD | SBGN-PD: each entity pool is its own activity carrier. |
| <a id="rule-preparation-no_complex-sbgn_pd-carrier_phenotype"></a>`preparation:no_complex:sbgn_pd:carrier_phenotype` | SBGN PD | SBGN-PD: each phenotype process is its own activity carrier. |
| <a id="rule-preparation-no_complex-sbgn_pd-carrier_subunit"></a>`preparation:no_complex:sbgn_pd:carrier_subunit` | SBGN PD | SBGN-PD: each subunit is its own activity carrier, so a promoted subunit can be an influence endpoint (paths reach it via `paths:complex_traversal`). |

### `paths` {#concern-paths}

#### `paths:core`

Builds the direct (single-hop) kinded `propagatesInfluence(SOURCE_SPECIES, TARGET_SPECIES, INFLUENCE_KIND)` relation from PD reactions and modulation arcs, plus the production-cycle relations that gate transitive extension. INFLUENCE_KIND is one of `positive`, `negative`, `triggering`, `modulation` and their `unknown_*` twins. Reaction modifiers and the matching species→species modulation arcs map to the *same* kind (e.g. a trigger modifier and a triggering arc both give `triggering`; catalysis and physical stimulation both give `positive`; inhibition gives `negative`). Each variant defines the single reactant->product hop `isDirectlyTransformedTo/2`, the shared `isTransformedTo/2` is its transitive closure, and `isCyclicallyTransformedTo/2` marks the hops that lie inside a cycle so transitivity never propagates influence back around a production loop. The multi-hop reactant-chained transitivity itself is the excludable `paths:chaining` group.

| Rule | Variant | Documentation |
| --- | --- | --- |
| <a id="rule-paths-core-is_transformed_to"></a>`paths:core:is_transformed_to` |  | Transitive closure of the single-hop `isDirectlyTransformedTo/2`: `isTransformedTo(UPSTREAM, DOWNSTREAM)` holds when DOWNSTREAM is reachable from UPSTREAM through one or more reactant->product hops (the direct hop included, as the base case). Language-agnostic; `isDirectlyTransformedTo/2` is supplied per language variant. Depends only on `isDirectlyTransformedTo`, never on `propagatesInfluence`, so the negation that gates transitivity stays stratified. |
| <a id="rule-paths-core-is_cyclically_transformed_to"></a>`paths:core:is_cyclically_transformed_to` |  | UPSTREAM and DOWNSTREAM are mutually reachable through the production graph -- both lie on a common cycle. At every transitivity guard site the hop in question is already a direct reactant->product edge, so requiring mutual reachability there is exactly equivalent to 'this hop is on a cycle'. The transitivity rules forbid extending a path across such a hop, blocking influence from leaking around production loops. |
| <a id="rule-paths-core-composes_to"></a>`paths:core:composes_to` |  | `composesTo(INCOMING_INFLUENCE_KIND, OUTGOING_INFLUENCE_KIND)`: how an influence kind transforms when a path is extended by one reactant→product hop (a reaction in CellDesigner, a process in SBGN-PD). Identity for every kind except `triggering` (a necessary-stimulation relationship is a property of the direct edge; composed through a reaction it weakens to a plain `positive` influence) and its unknown twin `unknown_triggering` → `unknown_positive`. Modulation composes like the signed kinds (stays `modulation`). |
| <a id="rule-paths-core-celldesigner-catalyzer_to_product"></a>`paths:core:celldesigner:catalyzer_to_product` | CellDesigner | If a species is referred to by a catalyzer of a reaction and another species is referred to by a product of that reaction, then there is a positive path from the first to the second. |
| <a id="rule-paths-core-celldesigner-physical_stimulator_to_product"></a>`paths:core:celldesigner:physical_stimulator_to_product` | CellDesigner | If a species is referred to by a physical stimulator of a reaction and another species is referred to by a product of that reaction, then there is a positive path from the first to the second. |
| <a id="rule-paths-core-celldesigner-trigger_to_product"></a>`paths:core:celldesigner:trigger_to_product` | CellDesigner | If a species is referred to by a trigger of a reaction and another species is referred to by a product of that reaction, then there is a triggering path from the first to the second (a trigger→product edge is direct, so it keeps the `triggering` kind). |
| <a id="rule-paths-core-celldesigner-modulation_arc_influence"></a>`paths:core:celldesigner:modulation_arc_influence` | CellDesigner | A modulation arc influences its target with the arc's kind (`hasInfluenceKind`, from the `influences:kind` group). This covers every arc type -- catalysis/physical stimulation/positive influence give `positive`, triggering gives `triggering`, inhibition/negative influence give `negative`, a bare modulation gives `modulation`, and the `unknown*` twins give the `unknown_*` kinds. The reaction modifier->product rules and the transitive/cycle rules remain separate. |
| <a id="rule-paths-core-celldesigner-inhibitor_to_product"></a>`paths:core:celldesigner:inhibitor_to_product` | CellDesigner | If a species is referred to by an inhibitor of a reaction and another species is referred to by a product of that reaction, then there is a negative path from the first to the second. |
| <a id="rule-paths-core-celldesigner-modulator_to_product"></a>`paths:core:celldesigner:modulator_to_product` | CellDesigner | If a species is referred to by a *bare* modulator of a reaction (a generic MODULATION modifier — not a physical stimulator, inhibitor or trigger; catalyzers are physical stimulators) and another species is referred to by a product of that reaction, then there is a modulation path from the first to the second. |
| <a id="rule-paths-core-celldesigner-unknown_catalyzer_to_product"></a>`paths:core:celldesigner:unknown_catalyzer_to_product` | CellDesigner | If a species is referred to by an unknown catalyzer of a reaction and another species is referred to by a product of that reaction, then there is an unknown-positive path from the first to the second. |
| <a id="rule-paths-core-celldesigner-unknown_inhibitor_to_product"></a>`paths:core:celldesigner:unknown_inhibitor_to_product` | CellDesigner | If a species is referred to by an unknown inhibitor of a reaction and another species is referred to by a product of that reaction, then there is an unknown-negative path from the first to the second. |
| <a id="rule-paths-core-celldesigner-is_directly_transformed_to"></a>`paths:core:celldesigner:is_directly_transformed_to` | CellDesigner | The single reactant->product hop in the production graph: the upstream species is referred to by a reactant and the downstream species by a product of the same reaction. Passive voice (the reaction does the transforming, not the species) keeps it language-neutral. Feeds the shared `isTransformedTo`/`isCyclicallyTransformedTo` cycle relations that gate transitive path extension; carries no influence kind itself. |
| <a id="rule-paths-core-sbgn_pd-modulation_to_product"></a>`paths:core:sbgn_pd:modulation_to_product` | SBGN PD | A modulation arc's entity-pool source influences each product of its target process, carrying the arc's kind. |
| <a id="rule-paths-core-sbgn_pd-modulation_to_phenotype"></a>`paths:core:sbgn_pd:modulation_to_phenotype` | SBGN PD | A modulation arc whose target is a phenotype influences the phenotype itself (a phenotype process has no products; it is the activity). |
| <a id="rule-paths-core-sbgn_pd-is_directly_transformed_to"></a>`paths:core:sbgn_pd:is_directly_transformed_to` | SBGN PD | The single reactant->product hop in the production graph: the upstream entity pool is an element of a reactant and the downstream pool an element of a product of the same process. Passive voice (the process does the transforming, not the pool) keeps it language-neutral. Feeds the shared `isTransformedTo`/`isCyclicallyTransformedTo` cycle relations that gate transitive path extension; carries no influence kind itself. |

#### `paths:chaining`

Excludable multi-hop transitivity: extends a `propagatesInfluence` path through one more reactant->product hop (a reaction in CellDesigner, a process in SBGN-PD), carrying the kind via `composesTo` (triggering degrades to positivelyInfluences) and refusing to extend across a hop inside a production cycle (`not isCyclicallyTransformedTo`). Exclude with `--exclude-group paths:chaining` to keep only direct single-hop influences. Depends on `paths:core`, which supplies both `composesTo` and the cycle relations.

| Rule | Variant | Documentation |
| --- | --- | --- |
| <a id="rule-paths-chaining-celldesigner-transitive_through_reaction"></a>`paths:chaining:celldesigner:transitive_through_reaction` | CellDesigner | If there is a path from a species to an intermediate species, and the intermediate species is referred to by a reactant of a reaction whose product refers to another species, then there is a path from the first species to the second. The kind is carried through `composesTo`, which degrades `triggering`→`positive` (and `unknown_triggering`→`unknown_positive`) at the reaction hop while leaving every other kind unchanged. Extension is suppressed across a reactant->product hop that lies inside a cycle (`isCyclicallyTransformedTo`), so a source feeding a production cycle does not leak influence back around the loop onto members it directly depletes. |
| <a id="rule-paths-chaining-sbgn_pd-transitive_through_process"></a>`paths:chaining:sbgn_pd:transitive_through_process` | SBGN PD | Extends a path through a process reactant->product hop, carrying the kind via `composesTo` (triggering degrades to positivelyInfluences). Extension is suppressed across a reactant->product hop that lies inside a cycle (`isCyclicallyTransformedTo`), so a source feeding a production cycle does not leak influence back around the loop onto members it directly depletes. |

#### `paths:complex_traversal`

Extends paths through complex containment for the ``*-no-complex`` modes: a path touching a complex is propagated to/from each of its subunits so influences reach the surviving subunit activities.

| Rule | Variant | Documentation |
| --- | --- | --- |
| <a id="rule-paths-complex_traversal-into_subunits"></a>`paths:complex_traversal:into_subunits` |  | If an influence of a given kind propagates from a source to a complex, then it also propagates from the source to each subunit of the complex. |
| <a id="rule-paths-complex_traversal-from_subunits"></a>`paths:complex_traversal:from_subunits` |  | If an influence of a given kind propagates from a complex to a target, then it also propagates from each subunit of the complex to the target. |

### `influences` {#concern-influences}

#### `influences:kind`

Maps each modulation arc to its influence kind via `hasInfluenceKind(MODULATION, INFLUENCE_KIND)`, the single place the arc-type->kind knowledge lives. Every group that emits a modulation-arc influence reads this relation rather than re-encoding the mapping. The mapping is per-language: CellDesigner arcs carry the sign in the arc type (catalysis, inhibition, ...), while an SBGN-PD arc's kind comes from its stimulation/inhibition/necessary-stimulation classification.

| Rule | Variant | Documentation |
| --- | --- | --- |
| <a id="rule-influences-kind-celldesigner-catalysis"></a>`influences:kind:celldesigner:catalysis` | CellDesigner | A catalysis arc contributes a `positive` kind. |
| <a id="rule-influences-kind-celldesigner-physical_stimulation"></a>`influences:kind:celldesigner:physical_stimulation` | CellDesigner | A physical stimulation arc contributes a `positive` kind. |
| <a id="rule-influences-kind-celldesigner-positive_influence"></a>`influences:kind:celldesigner:positive_influence` | CellDesigner | A positive influence arc contributes a `positive` kind. |
| <a id="rule-influences-kind-celldesigner-triggering"></a>`influences:kind:celldesigner:triggering` | CellDesigner | A triggering arc contributes a `triggering` kind (kept on the direct edge; it degrades to `positive` when composed through a reaction). |
| <a id="rule-influences-kind-celldesigner-inhibition"></a>`influences:kind:celldesigner:inhibition` | CellDesigner | An inhibition arc contributes a `negative` kind. |
| <a id="rule-influences-kind-celldesigner-negative_influence"></a>`influences:kind:celldesigner:negative_influence` | CellDesigner | A negative influence arc contributes a `negative` kind. |
| <a id="rule-influences-kind-celldesigner-modulation"></a>`influences:kind:celldesigner:modulation` | CellDesigner | A bare modulation arc (the generic MODULATION arc, not one of the signed/typed subtypes) contributes an unknown-sign `modulation` kind. The negations exclude the subtypes, which `modulation` is the umbrella over. |
| <a id="rule-influences-kind-celldesigner-unknown_catalysis"></a>`influences:kind:celldesigner:unknown_catalysis` | CellDesigner | An unknown catalysis arc contributes an `unknown_positive` kind. |
| <a id="rule-influences-kind-celldesigner-unknown_physical_stimulation"></a>`influences:kind:celldesigner:unknown_physical_stimulation` | CellDesigner | An unknown physical stimulation arc contributes an `unknown_positive` kind. |
| <a id="rule-influences-kind-celldesigner-unknown_positive_influence"></a>`influences:kind:celldesigner:unknown_positive_influence` | CellDesigner | An unknown positive influence arc contributes an `unknown_positive` kind. |
| <a id="rule-influences-kind-celldesigner-unknown_triggering"></a>`influences:kind:celldesigner:unknown_triggering` | CellDesigner | An unknown triggering arc contributes an `unknown_triggering` kind (kept on the direct edge; it degrades to `unknown_positive` when composed through a reaction). |
| <a id="rule-influences-kind-celldesigner-unknown_inhibition"></a>`influences:kind:celldesigner:unknown_inhibition` | CellDesigner | An unknown inhibition arc contributes an `unknown_negative` kind. |
| <a id="rule-influences-kind-celldesigner-unknown_negative_influence"></a>`influences:kind:celldesigner:unknown_negative_influence` | CellDesigner | An unknown negative influence arc contributes an `unknown_negative` kind. |
| <a id="rule-influences-kind-celldesigner-unknown_modulation"></a>`influences:kind:celldesigner:unknown_modulation` | CellDesigner | A bare unknown modulation arc (not one of the unknown subtypes) contributes an `unknown_modulation` kind. The negations exclude the subtypes, which `unknownModulation` is the umbrella over. |
| <a id="rule-influences-kind-sbgn_pd-necessary_stimulation"></a>`influences:kind:sbgn_pd:necessary_stimulation` | SBGN PD | An SBGN-PD necessary stimulation contributes a `triggering` kind. |
| <a id="rule-influences-kind-sbgn_pd-stimulation"></a>`influences:kind:sbgn_pd:stimulation` | SBGN PD | A stimulation (catalysis included) that is not a necessary stimulation contributes a `positive` kind. |
| <a id="rule-influences-kind-sbgn_pd-inhibition"></a>`influences:kind:sbgn_pd:inhibition` | SBGN PD | An SBGN-PD inhibition contributes a `negative` kind. |
| <a id="rule-influences-kind-sbgn_pd-modulation"></a>`influences:kind:sbgn_pd:modulation` | SBGN PD | A bare modulation (neither stimulation nor inhibition) contributes an unknown-sign `modulation` kind. |

#### `influences:core`

Derivation: emits `new(activity(KEY))` for every activity key, and lifts every kinded path into the internal `influences(SOURCE_KEY, TARGET_KEY, INFLUENCE_KIND)` relation by resolving both endpoints through their activity carrier and key. Both rules are mode- and language-agnostic; the consumption/sparing reasoning lives in the separate `influences:consumption` group. The internal `influences/3` relation is fanned out to the typed `new(...)` heads by the shared `influences:output` group.

| Rule | Variant | Documentation |
| --- | --- | --- |
| <a id="rule-influences-core-activity"></a>`influences:core:activity` |  | Every species with an activity key emits an activity node with that key. |
| <a id="rule-influences-core-path"></a>`influences:core:path` |  | A propagated influence of any kind between two raw nodes yields an influence of that same kind between their activity-key images (via their carriers). |

#### `influences:consumption`

Excludable consumption/sparing reasoning: a reaction depletes its reactants, so a modifier that drives the reaction also acts on every reactant that is itself an activity -- catalyzer/physicalStimulator/trigger negatively influence each consumed reactant, inhibitor positively influences each spared reactant, and the unknown modifiers contribute the unknown twins. This is inference beyond what the map draws, so it is a group of its own: exclude with `--exclude-group influences:consumption` to keep only the influences the map states. The `keep-reactions` mode omits it, since that mode renders each reaction directly instead of reasoning about it.

| Rule | Variant | Documentation |
| --- | --- | --- |
| <a id="rule-influences-consumption-celldesigner-catalyzer_consumes_reactant"></a>`influences:consumption:celldesigner:catalyzer_consumes_reactant` | CellDesigner | A catalyzer of a reaction negatively influences each reactant that is itself an activity (consumption depletes the reactant — a negative influence regardless of the modifier's positive role on the product). |
| <a id="rule-influences-consumption-celldesigner-physical_stimulator_consumes_reactant"></a>`influences:consumption:celldesigner:physical_stimulator_consumes_reactant` | CellDesigner | A physical stimulator of a reaction negatively influences each reactant that is itself an activity (consumption). |
| <a id="rule-influences-consumption-celldesigner-trigger_consumes_reactant"></a>`influences:consumption:celldesigner:trigger_consumes_reactant` | CellDesigner | A trigger of a reaction negatively influences each reactant that is itself an activity (consumption is depletion, hence negative — not triggers, which is only the trigger→product relationship). |
| <a id="rule-influences-consumption-celldesigner-inhibitor_spares_reactant"></a>`influences:consumption:celldesigner:inhibitor_spares_reactant` | CellDesigner | An inhibitor of a reaction positively influences each reactant that is itself an activity (sparing). |
| <a id="rule-influences-consumption-celldesigner-unknown_catalyzer_consumes_reactant"></a>`influences:consumption:celldesigner:unknown_catalyzer_consumes_reactant` | CellDesigner | An unknown catalyzer of a reaction unknown-negatively influences each reactant that is itself an activity (consumption, uncertain). |
| <a id="rule-influences-consumption-celldesigner-unknown_inhibitor_spares_reactant"></a>`influences:consumption:celldesigner:unknown_inhibitor_spares_reactant` | CellDesigner | An unknown inhibitor of a reaction unknown-positively influences each reactant that is itself an activity (sparing, uncertain). |

#### `influences:output`

Shared fan-out from the internal `influences(SOURCE, TARGET, INFLUENCE_KIND)` relation to the typed `new(...)` influence heads — one rule per kind. Every pipeline converges on `influences/3`; this group is the single place that turns a kind into its output predicate.

| Rule | Variant | Documentation |
| --- | --- | --- |
| <a id="rule-influences-output-positive"></a>`influences:output:positive` |  | A `positive` influence emits a `positivelyInfluences` edge (PositiveInfluence). |
| <a id="rule-influences-output-negative"></a>`influences:output:negative` |  | A `negative` influence emits a `negativelyInfluences` edge (NegativeInfluence). |
| <a id="rule-influences-output-modulation"></a>`influences:output:modulation` |  | A `modulation` influence emits a `modulates` edge (Modulation). |
| <a id="rule-influences-output-triggering"></a>`influences:output:triggering` |  | A `triggering` influence emits a `triggers` edge (Triggering). |
| <a id="rule-influences-output-unknown_positive"></a>`influences:output:unknown_positive` |  | An `unknown_positive` influence emits an `unknownPositivelyInfluences` edge (UnknownPositiveInfluence). |
| <a id="rule-influences-output-unknown_negative"></a>`influences:output:unknown_negative` |  | An `unknown_negative` influence emits an `unknownNegativelyInfluences` edge (UnknownNegativeInfluence). |
| <a id="rule-influences-output-unknown_modulation"></a>`influences:output:unknown_modulation` |  | An `unknown_modulation` influence emits an `unknownModulates` edge (UnknownModulation). |
| <a id="rule-influences-output-unknown_triggering"></a>`influences:output:unknown_triggering` |  | An `unknown_triggering` influence emits an `unknownTriggers` edge (UnknownTriggering). |

### `gates` {#concern-gates}

#### `gates:core`

Authored logical operators (CellDesigner `BooleanLogicGate`, SBGN-PD `LogicalOperator`): emits the operator node (`logicalOperator/2`, token-typed), its input edges (`logicalOperatorInput/2`, each input resolved through carrier/key), and the operator-sourced influence written straight into the internal `influences/3` relation by reusing the existing `propagatesInfluence(OPERATOR, TARGET, INFLUENCE_KIND)` closure. The widened influence `source` union (`predicates._INFLUENCE_SOURCE`) lets `influences:output` fan these out with no change. Provenance-agnostic.

| Rule | Variant | Documentation |
| --- | --- | --- |
| <a id="rule-gates-core-celldesigner-node_and"></a>`gates:core:celldesigner:node_and` | CellDesigner | An `andGate` emits an AND logical-operator node. |
| <a id="rule-gates-core-celldesigner-node_or"></a>`gates:core:celldesigner:node_or` | CellDesigner | An `orGate` emits an OR logical-operator node. |
| <a id="rule-gates-core-celldesigner-node_not"></a>`gates:core:celldesigner:node_not` | CellDesigner | A `notGate` emits a NOT logical-operator node. The token is `not_` because bare `not` is a reserved clingo keyword. |
| <a id="rule-gates-core-celldesigner-node_unknown"></a>`gates:core:celldesigner:node_unknown` | CellDesigner | An `unknownGate` emits an unknown-type logical-operator node. |
| <a id="rule-gates-core-celldesigner-input_edge"></a>`gates:core:celldesigner:input_edge` | CellDesigner | Each gate input is resolved through its activity carrier and key (mirroring `influences:core:path`), so an input that is a subunit resolves to its top-level complex's key. The `booleanLogicGate` umbrella matches every gate type via the ontology's isa rules. |
| <a id="rule-gates-core-celldesigner-influence"></a>`gates:core:celldesigner:influence` | CellDesigner | One rule covering both Shape A (gate is a reaction modifier) and Shape B (gate is a modulation source): the `paths:core` rules already bind a `propagatesInfluence/3` whose source is the gate, so the gate only resolves its TARGET through carrier/key and writes the influence keyed by the operator. Inherits the transitive closure (Decision D1). The `booleanLogicGate` guard is essential -- without it every species `propagatesInfluence/3` source would be read as an operator key. |
| <a id="rule-gates-core-sbgn_pd-node_and"></a>`gates:core:sbgn_pd:node_and` | SBGN PD | An `andOperator` emits an AND logical-operator node. |
| <a id="rule-gates-core-sbgn_pd-node_or"></a>`gates:core:sbgn_pd:node_or` | SBGN PD | An `orOperator` emits an OR logical-operator node. |
| <a id="rule-gates-core-sbgn_pd-node_not"></a>`gates:core:sbgn_pd:node_not` | SBGN PD | A `notOperator` emits a NOT logical-operator node. The token is `not_` because bare `not` is a reserved clingo keyword. (SBGN-PD has no unknown-operator type.) |
| <a id="rule-gates-core-sbgn_pd-input_edge"></a>`gates:core:sbgn_pd:input_edge` | SBGN PD | SBGN-PD parallel of the CellDesigner input-edge rule, guarded on the `logicalOperator` umbrella (derived from the per-type operators via the ontology isa rules). |
| <a id="rule-gates-core-sbgn_pd-influence"></a>`gates:core:sbgn_pd:influence` | SBGN PD | SBGN-PD parallel of the CellDesigner operator-sourced influence rule (Shape B: the operator is a modulation source). Reuses the `propagatesInfluence(OPERATOR, TARGET, INFLUENCE_KIND)` closure emitted by `paths:core`, resolving only the target through carrier/key. |

### `keep_reactions` {#concern-keep_reactions}

#### `keep_reactions:activity`

The `keep-reactions` premise that every species is an activity, replacing the structural-reason `activity:*` feature-groups (each of which is inert once every species is a candidate, so the mode omits them). Unlike the `--set-all-active` toggle this rule carries no `not hasSubunit` guard: a subunit is a candidate too, which does not make it an activity of its own -- `preparation:complex:key` keys it by the `keptSpeciesKey` of its outermost complex -- but it does mean a reaction, modulation arc or gate input touching a subunit routes to the containing complex instead of being dropped for want of an activity key.

| Rule | Variant | Documentation |
| --- | --- | --- |
| <a id="rule-keep_reactions-activity-from_species"></a>`keep_reactions:activity:from_species` |  | Every species is an activity candidate, for the sole reason that it is a species (reason `isSpecies`). It still passes through the `activity:core` bridge, so `--set-inactive` and `--set-all-inactive` veto candidates in this mode exactly as in any other. |

#### `keep_reactions:influences`

The `keep-reactions` premise that every reaction is kept: each (reactant, product) pair of a reaction becomes one positive influence. `positivelyInfluences` rather than `triggers` -- consuming a reactant to make a product is a contribution to it, not the necessary-stimulation relationship `triggers` (CellDesigner `Triggering`, SBGN-AF `NecessaryStimulation`) asserts. The rule writes `propagatesInfluence/3` rather than `influences/3` so the shared `influences:core:path` bridge resolves both endpoints through their carrier and key, routing a subunit endpoint to its top-level complex. Reactant and product may resolve to the same activity (a reaction whose participants share a complex, or a state transition drawn on a single species), and the resulting self-influence is emitted like any other.

| Rule | Variant | Documentation |
| --- | --- | --- |
| <a id="rule-keep_reactions-influences-reactant_to_product"></a>`keep_reactions:influences:reactant_to_product` |  | Each reactant of a reaction positively influences each product of that reaction. This is the one rule that makes a reaction itself an influence: in the other modes the reactant->product hop only feeds the cycle relations, and a bare reactant is never an influence source. |

## Rule index {#rule-index}

All rules in alphabetical order (identifiers are group- and variant-prefixed).

| Rule | Group | Variant | Modes | Summary |
| --- | --- | --- | --- | --- |
| [`activity:active_marker:celldesigner:from_active_flag`](#rule-activity-active_marker-celldesigner-from_active_flag) | `activity:active_marker` | CellDesigner | `normal`, `normal-no-complex`, `keep-species`, `keep-species-no-complex` | A species explicitly marked as active is an activity candidate. |
| [`activity:active_marker:celldesigner:from_active_structural_state`](#rule-activity-active_marker-celldesigner-from_active_structural_state) | `activity:active_marker` | CellDesigner | `normal`, `normal-no-complex`, `keep-species`, `keep-species-no-complex` | A species in an active structural state is an activity candidate. |
| [`activity:active_marker:sbgn_pd:from_active_state_variable`](#rule-activity-active_marker-sbgn_pd-from_active_state_variable) | `activity:active_marker` | SBGN PD | `normal`, `normal-no-complex`, `keep-species`, `keep-species-no-complex` | An entity pool in an active state is an activity candidate. |
| [`activity:active_marker:sbgn_pd:from_active_subunit_state_variable`](#rule-activity-active_marker-sbgn_pd-from_active_subunit_state_variable) | `activity:active_marker` | SBGN PD | `normal`, `normal-no-complex`, `keep-species`, `keep-species-no-complex` | A subunit in an active state is an activity candidate. |
| [`activity:core:celldesigner:from_global_activate`](#rule-activity-core-celldesigner-from_global_activate) | `activity:core` | CellDesigner | all | Under the set-all-active toggle, every top-level species is an activity candidate. |
| [`activity:core:from_candidate`](#rule-activity-core-from_candidate) | `activity:core` |  | all | An activity candidate becomes an actual activity unless it has been vetoed. |
| [`activity:core:from_global_suppress`](#rule-activity-core-from_global_suppress) | `activity:core` |  | all | The set-all-inactive toggle suppresses every activity candidate, except those pinned active per id. |
| [`activity:core:sbgn_pd:from_global_activate`](#rule-activity-core-sbgn_pd-from_global_activate) | `activity:core` | SBGN PD | all | Under the set-all-active toggle, every entity pool is an activity candidate. |
| [`activity:gate_input:celldesigner:from_gate_input`](#rule-activity-gate_input-celldesigner-from_gate_input) | `activity:gate_input` | CellDesigner | `normal`, `normal-no-complex`, `keep-species`, `keep-species-no-complex` | A species feeding a boolean logic gate is an activity candidate. |
| [`activity:gate_input:sbgn_pd:from_operator_input`](#rule-activity-gate_input-sbgn_pd-from_operator_input) | `activity:gate_input` | SBGN PD | `normal`, `normal-no-complex`, `keep-species`, `keep-species-no-complex` | An entity pool feeding a logical operator is an activity candidate. |
| [`activity:modulation_source:celldesigner:from_modulation_source`](#rule-activity-modulation_source-celldesigner-from_modulation_source) | `activity:modulation_source` | CellDesigner | `normal`, `normal-no-complex`, `keep-species`, `keep-species-no-complex` | A species that is the source of a modulation arc is an activity candidate. |
| [`activity:modulation_source:celldesigner:from_reaction_modulator`](#rule-activity-modulation_source-celldesigner-from_reaction_modulator) | `activity:modulation_source` | CellDesigner | `normal`, `normal-no-complex`, `keep-species`, `keep-species-no-complex` | A species that modulates a reaction is an activity candidate. |
| [`activity:modulation_source:sbgn_pd:from_modulation_source`](#rule-activity-modulation_source-sbgn_pd-from_modulation_source) | `activity:modulation_source` | SBGN PD | `normal`, `normal-no-complex`, `keep-species`, `keep-species-no-complex` | An entity pool that is the source of a modulation arc is an activity candidate. |
| [`activity:phenotype:from_phenotype`](#rule-activity-phenotype-from_phenotype) | `activity:phenotype` |  | `normal`, `normal-no-complex`, `keep-species`, `keep-species-no-complex` | A phenotype is always an activity candidate. |
| [`gates:core:celldesigner:influence`](#rule-gates-core-celldesigner-influence) | `gates:core` | CellDesigner | all | One rule covering both Shape A (gate is a reaction modifier) and Shape B (gate is a modulation source): the `paths:core` rules already bind a `propagatesInfluence/3` whose source is the gate, so the gate only resolves its TARGET through carrier/key and writes the influence keyed by the operator. |
| [`gates:core:celldesigner:input_edge`](#rule-gates-core-celldesigner-input_edge) | `gates:core` | CellDesigner | all | Each gate input is resolved through its activity carrier and key (mirroring `influences:core:path`), so an input that is a subunit resolves to its top-level complex's key. |
| [`gates:core:celldesigner:node_and`](#rule-gates-core-celldesigner-node_and) | `gates:core` | CellDesigner | all | An `andGate` emits an AND logical-operator node. |
| [`gates:core:celldesigner:node_not`](#rule-gates-core-celldesigner-node_not) | `gates:core` | CellDesigner | all | A `notGate` emits a NOT logical-operator node. |
| [`gates:core:celldesigner:node_or`](#rule-gates-core-celldesigner-node_or) | `gates:core` | CellDesigner | all | An `orGate` emits an OR logical-operator node. |
| [`gates:core:celldesigner:node_unknown`](#rule-gates-core-celldesigner-node_unknown) | `gates:core` | CellDesigner | all | An `unknownGate` emits an unknown-type logical-operator node. |
| [`gates:core:sbgn_pd:influence`](#rule-gates-core-sbgn_pd-influence) | `gates:core` | SBGN PD | all | SBGN-PD parallel of the CellDesigner operator-sourced influence rule (Shape B: the operator is a modulation source). |
| [`gates:core:sbgn_pd:input_edge`](#rule-gates-core-sbgn_pd-input_edge) | `gates:core` | SBGN PD | all | SBGN-PD parallel of the CellDesigner input-edge rule, guarded on the `logicalOperator` umbrella (derived from the per-type operators via the ontology isa rules). |
| [`gates:core:sbgn_pd:node_and`](#rule-gates-core-sbgn_pd-node_and) | `gates:core` | SBGN PD | all | An `andOperator` emits an AND logical-operator node. |
| [`gates:core:sbgn_pd:node_not`](#rule-gates-core-sbgn_pd-node_not) | `gates:core` | SBGN PD | all | A `notOperator` emits a NOT logical-operator node. |
| [`gates:core:sbgn_pd:node_or`](#rule-gates-core-sbgn_pd-node_or) | `gates:core` | SBGN PD | all | An `orOperator` emits an OR logical-operator node. |
| [`influences:consumption:celldesigner:catalyzer_consumes_reactant`](#rule-influences-consumption-celldesigner-catalyzer_consumes_reactant) | `influences:consumption` | CellDesigner | `normal`, `normal-no-complex`, `keep-species`, `keep-species-no-complex` | A catalyzer of a reaction negatively influences each reactant that is itself an activity (consumption depletes the reactant — a negative influence regardless of the modifier's positive role on the product). |
| [`influences:consumption:celldesigner:inhibitor_spares_reactant`](#rule-influences-consumption-celldesigner-inhibitor_spares_reactant) | `influences:consumption` | CellDesigner | `normal`, `normal-no-complex`, `keep-species`, `keep-species-no-complex` | An inhibitor of a reaction positively influences each reactant that is itself an activity (sparing). |
| [`influences:consumption:celldesigner:physical_stimulator_consumes_reactant`](#rule-influences-consumption-celldesigner-physical_stimulator_consumes_reactant) | `influences:consumption` | CellDesigner | `normal`, `normal-no-complex`, `keep-species`, `keep-species-no-complex` | A physical stimulator of a reaction negatively influences each reactant that is itself an activity (consumption). |
| [`influences:consumption:celldesigner:trigger_consumes_reactant`](#rule-influences-consumption-celldesigner-trigger_consumes_reactant) | `influences:consumption` | CellDesigner | `normal`, `normal-no-complex`, `keep-species`, `keep-species-no-complex` | A trigger of a reaction negatively influences each reactant that is itself an activity (consumption is depletion, hence negative — not triggers, which is only the trigger→product relationship). |
| [`influences:consumption:celldesigner:unknown_catalyzer_consumes_reactant`](#rule-influences-consumption-celldesigner-unknown_catalyzer_consumes_reactant) | `influences:consumption` | CellDesigner | `normal`, `normal-no-complex`, `keep-species`, `keep-species-no-complex` | An unknown catalyzer of a reaction unknown-negatively influences each reactant that is itself an activity (consumption, uncertain). |
| [`influences:consumption:celldesigner:unknown_inhibitor_spares_reactant`](#rule-influences-consumption-celldesigner-unknown_inhibitor_spares_reactant) | `influences:consumption` | CellDesigner | `normal`, `normal-no-complex`, `keep-species`, `keep-species-no-complex` | An unknown inhibitor of a reaction unknown-positively influences each reactant that is itself an activity (sparing, uncertain). |
| [`influences:core:activity`](#rule-influences-core-activity) | `influences:core` |  | all | Every species with an activity key emits an activity node with that key. |
| [`influences:core:path`](#rule-influences-core-path) | `influences:core` |  | all | A propagated influence of any kind between two raw nodes yields an influence of that same kind between their activity-key images (via their carriers). |
| [`influences:kind:celldesigner:catalysis`](#rule-influences-kind-celldesigner-catalysis) | `influences:kind` | CellDesigner | all | A catalysis arc contributes a `positive` kind. |
| [`influences:kind:celldesigner:inhibition`](#rule-influences-kind-celldesigner-inhibition) | `influences:kind` | CellDesigner | all | An inhibition arc contributes a `negative` kind. |
| [`influences:kind:celldesigner:modulation`](#rule-influences-kind-celldesigner-modulation) | `influences:kind` | CellDesigner | all | A bare modulation arc (the generic MODULATION arc, not one of the signed/typed subtypes) contributes an unknown-sign `modulation` kind. |
| [`influences:kind:celldesigner:negative_influence`](#rule-influences-kind-celldesigner-negative_influence) | `influences:kind` | CellDesigner | all | A negative influence arc contributes a `negative` kind. |
| [`influences:kind:celldesigner:physical_stimulation`](#rule-influences-kind-celldesigner-physical_stimulation) | `influences:kind` | CellDesigner | all | A physical stimulation arc contributes a `positive` kind. |
| [`influences:kind:celldesigner:positive_influence`](#rule-influences-kind-celldesigner-positive_influence) | `influences:kind` | CellDesigner | all | A positive influence arc contributes a `positive` kind. |
| [`influences:kind:celldesigner:triggering`](#rule-influences-kind-celldesigner-triggering) | `influences:kind` | CellDesigner | all | A triggering arc contributes a `triggering` kind (kept on the direct edge; it degrades to `positive` when composed through a reaction). |
| [`influences:kind:celldesigner:unknown_catalysis`](#rule-influences-kind-celldesigner-unknown_catalysis) | `influences:kind` | CellDesigner | all | An unknown catalysis arc contributes an `unknown_positive` kind. |
| [`influences:kind:celldesigner:unknown_inhibition`](#rule-influences-kind-celldesigner-unknown_inhibition) | `influences:kind` | CellDesigner | all | An unknown inhibition arc contributes an `unknown_negative` kind. |
| [`influences:kind:celldesigner:unknown_modulation`](#rule-influences-kind-celldesigner-unknown_modulation) | `influences:kind` | CellDesigner | all | A bare unknown modulation arc (not one of the unknown subtypes) contributes an `unknown_modulation` kind. |
| [`influences:kind:celldesigner:unknown_negative_influence`](#rule-influences-kind-celldesigner-unknown_negative_influence) | `influences:kind` | CellDesigner | all | An unknown negative influence arc contributes an `unknown_negative` kind. |
| [`influences:kind:celldesigner:unknown_physical_stimulation`](#rule-influences-kind-celldesigner-unknown_physical_stimulation) | `influences:kind` | CellDesigner | all | An unknown physical stimulation arc contributes an `unknown_positive` kind. |
| [`influences:kind:celldesigner:unknown_positive_influence`](#rule-influences-kind-celldesigner-unknown_positive_influence) | `influences:kind` | CellDesigner | all | An unknown positive influence arc contributes an `unknown_positive` kind. |
| [`influences:kind:celldesigner:unknown_triggering`](#rule-influences-kind-celldesigner-unknown_triggering) | `influences:kind` | CellDesigner | all | An unknown triggering arc contributes an `unknown_triggering` kind (kept on the direct edge; it degrades to `unknown_positive` when composed through a reaction). |
| [`influences:kind:sbgn_pd:inhibition`](#rule-influences-kind-sbgn_pd-inhibition) | `influences:kind` | SBGN PD | all | An SBGN-PD inhibition contributes a `negative` kind. |
| [`influences:kind:sbgn_pd:modulation`](#rule-influences-kind-sbgn_pd-modulation) | `influences:kind` | SBGN PD | all | A bare modulation (neither stimulation nor inhibition) contributes an unknown-sign `modulation` kind. |
| [`influences:kind:sbgn_pd:necessary_stimulation`](#rule-influences-kind-sbgn_pd-necessary_stimulation) | `influences:kind` | SBGN PD | all | An SBGN-PD necessary stimulation contributes a `triggering` kind. |
| [`influences:kind:sbgn_pd:stimulation`](#rule-influences-kind-sbgn_pd-stimulation) | `influences:kind` | SBGN PD | all | A stimulation (catalysis included) that is not a necessary stimulation contributes a `positive` kind. |
| [`influences:output:modulation`](#rule-influences-output-modulation) | `influences:output` |  | all | A `modulation` influence emits a `modulates` edge (Modulation). |
| [`influences:output:negative`](#rule-influences-output-negative) | `influences:output` |  | all | A `negative` influence emits a `negativelyInfluences` edge (NegativeInfluence). |
| [`influences:output:positive`](#rule-influences-output-positive) | `influences:output` |  | all | A `positive` influence emits a `positivelyInfluences` edge (PositiveInfluence). |
| [`influences:output:triggering`](#rule-influences-output-triggering) | `influences:output` |  | all | A `triggering` influence emits a `triggers` edge (Triggering). |
| [`influences:output:unknown_modulation`](#rule-influences-output-unknown_modulation) | `influences:output` |  | all | An `unknown_modulation` influence emits an `unknownModulates` edge (UnknownModulation). |
| [`influences:output:unknown_negative`](#rule-influences-output-unknown_negative) | `influences:output` |  | all | An `unknown_negative` influence emits an `unknownNegativelyInfluences` edge (UnknownNegativeInfluence). |
| [`influences:output:unknown_positive`](#rule-influences-output-unknown_positive) | `influences:output` |  | all | An `unknown_positive` influence emits an `unknownPositivelyInfluences` edge (UnknownPositiveInfluence). |
| [`influences:output:unknown_triggering`](#rule-influences-output-unknown_triggering) | `influences:output` |  | all | An `unknown_triggering` influence emits an `unknownTriggers` edge (UnknownTriggering). |
| [`keep_reactions:activity:from_species`](#rule-keep_reactions-activity-from_species) | `keep_reactions:activity` |  | `keep-reactions` | Every species is an activity candidate, for the sole reason that it is a species (reason `isSpecies`). |
| [`keep_reactions:influences:reactant_to_product`](#rule-keep_reactions-influences-reactant_to_product) | `keep_reactions:influences` |  | `keep-reactions` | Each reactant of a reaction positively influences each product of that reaction. |
| [`paths:chaining:celldesigner:transitive_through_reaction`](#rule-paths-chaining-celldesigner-transitive_through_reaction) | `paths:chaining` | CellDesigner | `normal`, `normal-no-complex`, `keep-species`, `keep-species-no-complex` | If there is a path from a species to an intermediate species, and the intermediate species is referred to by a reactant of a reaction whose product refers to another species, then there is a path from the first species to the second. |
| [`paths:chaining:sbgn_pd:transitive_through_process`](#rule-paths-chaining-sbgn_pd-transitive_through_process) | `paths:chaining` | SBGN PD | `normal`, `normal-no-complex`, `keep-species`, `keep-species-no-complex` | Extends a path through a process reactant->product hop, carrying the kind via `composesTo` (triggering degrades to positivelyInfluences). |
| [`paths:complex_traversal:from_subunits`](#rule-paths-complex_traversal-from_subunits) | `paths:complex_traversal` |  | `normal-no-complex`, `keep-species-no-complex` | If an influence of a given kind propagates from a complex to a target, then it also propagates from each subunit of the complex to the target. |
| [`paths:complex_traversal:into_subunits`](#rule-paths-complex_traversal-into_subunits) | `paths:complex_traversal` |  | `normal-no-complex`, `keep-species-no-complex` | If an influence of a given kind propagates from a source to a complex, then it also propagates from the source to each subunit of the complex. |
| [`paths:core:celldesigner:catalyzer_to_product`](#rule-paths-core-celldesigner-catalyzer_to_product) | `paths:core` | CellDesigner | all | If a species is referred to by a catalyzer of a reaction and another species is referred to by a product of that reaction, then there is a positive path from the first to the second. |
| [`paths:core:celldesigner:inhibitor_to_product`](#rule-paths-core-celldesigner-inhibitor_to_product) | `paths:core` | CellDesigner | all | If a species is referred to by an inhibitor of a reaction and another species is referred to by a product of that reaction, then there is a negative path from the first to the second. |
| [`paths:core:celldesigner:is_directly_transformed_to`](#rule-paths-core-celldesigner-is_directly_transformed_to) | `paths:core` | CellDesigner | all | The single reactant->product hop in the production graph: the upstream species is referred to by a reactant and the downstream species by a product of the same reaction. |
| [`paths:core:celldesigner:modulation_arc_influence`](#rule-paths-core-celldesigner-modulation_arc_influence) | `paths:core` | CellDesigner | all | A modulation arc influences its target with the arc's kind (`hasInfluenceKind`, from the `influences:kind` group). |
| [`paths:core:celldesigner:modulator_to_product`](#rule-paths-core-celldesigner-modulator_to_product) | `paths:core` | CellDesigner | all | If a species is referred to by a *bare* modulator of a reaction (a generic MODULATION modifier — not a physical stimulator, inhibitor or trigger; catalyzers are physical stimulators) and another species is referred to by a product of that reaction, then there is a modulation path from the first to the second. |
| [`paths:core:celldesigner:physical_stimulator_to_product`](#rule-paths-core-celldesigner-physical_stimulator_to_product) | `paths:core` | CellDesigner | all | If a species is referred to by a physical stimulator of a reaction and another species is referred to by a product of that reaction, then there is a positive path from the first to the second. |
| [`paths:core:celldesigner:trigger_to_product`](#rule-paths-core-celldesigner-trigger_to_product) | `paths:core` | CellDesigner | all | If a species is referred to by a trigger of a reaction and another species is referred to by a product of that reaction, then there is a triggering path from the first to the second (a trigger→product edge is direct, so it keeps the `triggering` kind). |
| [`paths:core:celldesigner:unknown_catalyzer_to_product`](#rule-paths-core-celldesigner-unknown_catalyzer_to_product) | `paths:core` | CellDesigner | all | If a species is referred to by an unknown catalyzer of a reaction and another species is referred to by a product of that reaction, then there is an unknown-positive path from the first to the second. |
| [`paths:core:celldesigner:unknown_inhibitor_to_product`](#rule-paths-core-celldesigner-unknown_inhibitor_to_product) | `paths:core` | CellDesigner | all | If a species is referred to by an unknown inhibitor of a reaction and another species is referred to by a product of that reaction, then there is an unknown-negative path from the first to the second. |
| [`paths:core:composes_to`](#rule-paths-core-composes_to) | `paths:core` |  | all | `composesTo(INCOMING_INFLUENCE_KIND, OUTGOING_INFLUENCE_KIND)`: how an influence kind transforms when a path is extended by one reactant→product hop (a reaction in CellDesigner, a process in SBGN-PD). |
| [`paths:core:is_cyclically_transformed_to`](#rule-paths-core-is_cyclically_transformed_to) | `paths:core` |  | all | UPSTREAM and DOWNSTREAM are mutually reachable through the production graph -- both lie on a common cycle. |
| [`paths:core:is_transformed_to`](#rule-paths-core-is_transformed_to) | `paths:core` |  | all | Transitive closure of the single-hop `isDirectlyTransformedTo/2`: `isTransformedTo(UPSTREAM, DOWNSTREAM)` holds when DOWNSTREAM is reachable from UPSTREAM through one or more reactant->product hops (the direct hop included, as the base case). |
| [`paths:core:sbgn_pd:is_directly_transformed_to`](#rule-paths-core-sbgn_pd-is_directly_transformed_to) | `paths:core` | SBGN PD | all | The single reactant->product hop in the production graph: the upstream entity pool is an element of a reactant and the downstream pool an element of a product of the same process. |
| [`paths:core:sbgn_pd:modulation_to_phenotype`](#rule-paths-core-sbgn_pd-modulation_to_phenotype) | `paths:core` | SBGN PD | all | A modulation arc whose target is a phenotype influences the phenotype itself (a phenotype process has no products; it is the activity). |
| [`paths:core:sbgn_pd:modulation_to_product`](#rule-paths-core-sbgn_pd-modulation_to_product) | `paths:core` | SBGN PD | all | A modulation arc's entity-pool source influences each product of its target process, carrying the arc's kind. |
| [`preparation:complex:celldesigner:carrier`](#rule-preparation-complex-celldesigner-carrier) | `preparation:complex` | CellDesigner | `normal`, `keep-species`, `keep-reactions` | Each species is its own activity carrier (no rerouting; a subunit endpoint is rerouted to its top-level complex by the `topology:top_level` key, not by the carrier). |
| [`preparation:complex:key`](#rule-preparation-complex-key) | `preparation:complex` |  | `normal`, `keep-species`, `keep-reactions` | A species with activity is keyed by the `keptSpeciesKey` of its top-level entity (itself when top-level; its outermost complex when a subunit). |
| [`preparation:complex:sbgn_pd:carrier_entity_pool`](#rule-preparation-complex-sbgn_pd-carrier_entity_pool) | `preparation:complex` | SBGN PD | `normal`, `keep-species`, `keep-reactions` | SBGN-PD: each entity pool is its own activity carrier. |
| [`preparation:complex:sbgn_pd:carrier_phenotype`](#rule-preparation-complex-sbgn_pd-carrier_phenotype) | `preparation:complex` | SBGN PD | `normal`, `keep-species`, `keep-reactions` | SBGN-PD: each phenotype process is its own activity carrier (phenotypes are activities, not entity pools). |
| [`preparation:no_complex:celldesigner:carrier`](#rule-preparation-no_complex-celldesigner-carrier) | `preparation:no_complex` | CellDesigner | `normal-no-complex`, `keep-species-no-complex` | Each species is its own activity carrier. |
| [`preparation:no_complex:deleted`](#rule-preparation-no_complex-deleted) | `preparation:no_complex` |  | `normal-no-complex`, `keep-species-no-complex` | A complex with any (transitive) active descendant is deleted. |
| [`preparation:no_complex:key_promoted_subunit`](#rule-preparation-no_complex-key_promoted_subunit) | `preparation:no_complex` |  | `normal-no-complex`, `keep-species-no-complex` | A subunit of a deleted complex with activity is promoted to top level, keyed by `promotedSubunitKey(SELF)`. |
| [`preparation:no_complex:key_top_level`](#rule-preparation-no_complex-key_top_level) | `preparation:no_complex` |  | `normal-no-complex`, `keep-species-no-complex` | A non-deleted top-level species with activity is keyed by `keptSpeciesKey(SELF)`. |
| [`preparation:no_complex:sbgn_pd:carrier_entity_pool`](#rule-preparation-no_complex-sbgn_pd-carrier_entity_pool) | `preparation:no_complex` | SBGN PD | `normal-no-complex`, `keep-species-no-complex` | SBGN-PD: each entity pool is its own activity carrier. |
| [`preparation:no_complex:sbgn_pd:carrier_phenotype`](#rule-preparation-no_complex-sbgn_pd-carrier_phenotype) | `preparation:no_complex` | SBGN PD | `normal-no-complex`, `keep-species-no-complex` | SBGN-PD: each phenotype process is its own activity carrier. |
| [`preparation:no_complex:sbgn_pd:carrier_subunit`](#rule-preparation-no_complex-sbgn_pd-carrier_subunit) | `preparation:no_complex` | SBGN PD | `normal-no-complex`, `keep-species-no-complex` | SBGN-PD: each subunit is its own activity carrier, so a promoted subunit can be an influence endpoint (paths reach it via `paths:complex_traversal`). |
| [`preparation:no_complex:subunit_of_deleted_direct`](#rule-preparation-no_complex-subunit_of_deleted_direct) | `preparation:no_complex` |  | `normal-no-complex`, `keep-species-no-complex` | A direct subunit of a deleted complex is itself subunit-of-deleted. |
| [`preparation:no_complex:subunit_of_deleted_transitive`](#rule-preparation-no_complex-subunit_of_deleted_transitive) | `preparation:no_complex` |  | `normal-no-complex`, `keep-species-no-complex` | Subunit-of-deleted is transitive through nested complexes. |
| [`topology:core:has_active_descendant_direct`](#rule-topology-core-has_active_descendant_direct) | `topology:core` |  | all | A complex has an active descendant if any direct subunit has activity. |
| [`topology:core:has_active_descendant_transitive`](#rule-topology-core-has_active_descendant_transitive) | `topology:core` |  | all | The `hasActiveDescendantSubunit` relation is transitive through complex containment. |
| [`topology:core:is_subunit`](#rule-topology-core-is_subunit) | `topology:core` |  | all | A species is a subunit if it appears on the right side of any `hasSubunit` relation. |
| [`topology:top_level:celldesigner:self`](#rule-topology-top_level-celldesigner-self) | `topology:top_level` | CellDesigner | `normal`, `keep-species`, `keep-reactions` | CellDesigner: a species that is not a subunit of any complex is its own top-level entity. |
| [`topology:top_level:recursive`](#rule-topology-top_level-recursive) | `topology:top_level` |  | `normal`, `keep-species`, `keep-reactions` | A subunit resolves to the same top-level entity as its parent complex, recursively through nested complexes -- so a subunit at any depth resolves to its outermost complex. |
| [`topology:top_level:sbgn_pd:phenotype_self`](#rule-topology-top_level-sbgn_pd-phenotype_self) | `topology:top_level` | SBGN PD | `normal`, `keep-species`, `keep-reactions` | SBGN-PD: a phenotype is a process, not an entity pool, so the entity-pool self-rule never keys it; a phenotype is never a subunit, so it is always its own top-level entity. |
| [`topology:top_level:sbgn_pd:self`](#rule-topology-top_level-sbgn_pd-self) | `topology:top_level` | SBGN PD | `normal`, `keep-species`, `keep-reactions` | SBGN-PD: an entity pool that is not a subunit of any complex is its own top-level entity. |

