"""ASP rule composition for the pd2af transformation modes.

Which groups a mode is made of is declared by the mode itself
(:mod:`pd2af.modes`); this module owns the groups and composes the program
text for a named mode. The rules are organised in layers:

* **topology** (shared) — structural helpers: ``isSubunit``,
  ``hasActiveDescendantSubunit``.
* **top_level** (the complex-keeping modes ``keep-species`` and
  ``normal``) — ``resolvesToTopLevel(SPECIES, TOPLEVEL)`` resolves every species to its
  outermost top-level entity, so a subunit is keyed by (and its
  influences routed to) its top-level complex rather than itself.
* **preparation** (a ``slot``: exactly one filler group per mode, and
  every consumer of the activity-key predicates depends on the slot
  rather than on a particular filler) — emits
  ``hasActivityCarrier(RAW_SPECIES, ACTIVITY_BEARER)`` and
  ``hasActivityKey(ACTIVITY_BEARER, KEY)`` where ``KEY`` is one of two
  per-species wrappers: ``keptSpeciesKey/1`` (top-level species, and the
  top-level complex a subunit resolves to) or ``promotedSubunitKey/1`` (a
  subunit promoted to top level when its complex is dissolved in the
  ``*-no-complex`` modes). Proteoform/PTM stripping for the merged modes
  (``normal``/``normal-no-complex``) happens at the build stage, not in the key.
* **derivation** (shared) — emits
  ``new(activity(KEY))`` and ``new(positivelyInfluences(...))`` /
  ``new(negativelyInfluences(...))`` from ``hasActivityCarrier`` /
  ``hasActivityKey``. The inference layers on top of it — multi-hop
  ``paths:chaining`` and ``influences_consumption`` — are carried by the four
  path-inference modes only.

The ``keep-reactions`` mode reuses that whole scaffolding and swaps
discovery and inference for two rules: every species is an activity
(``keep_reactions:activity``) and every (reactant, product) pair of a reaction
is a positive influence (``keep_reactions:influences``). It keeps the single-hop
modulation influences of ``paths:core`` and drops the inference layers.
"""

from textwrap import dedent

from aspcompose import (
    CollectionPlan,
    PlanInvalidError,
    Rule,
    RuleGroup,
    RuleRegistry,
)

import pd2af.modes

# Activity discovery drives the path-inference modes: `keep-reactions` keeps the
# mandatory `activity:core` bridge but replaces the structural-reason
# feature-groups with a single rule making every species a candidate
# (`keep_reactions:activity`), which renders each feature-group inert.
#
# Activity discovery is parallel across languages: only the `phenotype` rule is
# language-agnostic (both languages emit the `phenotype` functor); every other
# signal is a per-language variant. CellDesigner and SBGN-PD share the same
# signals -- an "active" marker, a modulation source, and a phenotype -- except
# that SBGN-PD has no explicit `hasActive` flag, CellDesigner's active
# *structural state* becomes SBGN-PD's active *state variable*, and a reaction
# modifier (CellDesigner-only) is, in SBGN-PD, just a modulation arc whose
# target is a process (so it folds into the modulation-source rule). A bare
# reactant/product is never an activity in either language.
#
# The mandatory `activity:core` group holds the candidate->activity bridge and
# the two global toggles (`--set-all-active`/`--set-all-inactive`); everything
# depends on it, so the dependency graph makes it non-excludable. Each
# structural reason (phenotype, active marker, modulation source, gate input) is
# its own excludable feature-group depending on `activity:core`, so
# `--exclude-group activity:phenotype` drops exactly the phenotype-as-activity
# behavior and leaves the rest of the program intact.


def _activity_feature(name, *, base=(), cd=(), sbgn=(), docs=""):
    """One activity feature-group: a structural reason that derives a
    `hasActivityCandidate`, excludable independently, depending on
    `activity:core` (its candidate is only meaningful through the bridge)."""
    variants = {}
    if cd:
        variants["celldesigner"] = cd
    if sbgn:
        variants["sbgn_pd"] = sbgn
    return RuleGroup(
        identifier=f"activity:{name}",
        depends_on=frozenset({"activity:core"}),
        rules=base,
        variants=variants,
        docs=docs,
    )


_ACTIVITY_CORE = RuleGroup(
    identifier="activity:core",
    docs="Mandatory activity machinery: the candidate->activity bridge and the two global toggles. A single bridging rule promotes a `hasActivityCandidate(ELEMENT, REASON)` to `hasActivity` unless the element is `suppressActivity` (the `--set-inactive` veto). `globalSuppress` (`--set-all-inactive`) suppresses every candidate except those the solver marked `forceActive` (the per-id `--set-active` override), realising the precedence per-id > global > rules; `globalActivate` (`--set-all-active`) turns every top-level species/entity pool into a candidate. Every activity feature-group and every downstream group depends on this, so the dependency graph forbids excluding it -- the global toggles stay wired to their CLI flags. The structural-reason rules (phenotype, active marker, modulation source, gate input) live in the excludable `activity:*` feature-groups.",
    rules=(
        Rule(
            identifier="activity:core:from_candidate",
            text=dedent("""\
                hasActivity(ELEMENT, REASON) :-
                    hasActivityCandidate(ELEMENT, REASON),
                    not suppressActivity(ELEMENT)."""),
            docs="An activity candidate becomes an actual activity unless it has been vetoed. This is the single interception point for the `--set-inactive` veto (the `suppressActivity` marker; blanket suppression, so it also blocks the `--set-active` `isInputParameter` candidate, which is injected as a candidate). All downstream body references key on `hasActivity`, so they automatically respect suppression.",
        ),
        Rule(
            identifier="activity:core:from_global_suppress",
            text=dedent("""\
                suppressActivity(ELEMENT) :-
                    hasActivityCandidate(ELEMENT, _),
                    globalSuppress,
                    not forceActive(ELEMENT)."""),
            docs="The set-all-inactive toggle suppresses every activity candidate, except those pinned active per id. This realises the precedence per-id > global > rules for the inactive toggle: `globalSuppress` (`--set-all-inactive`) blankets everything, and `forceActive` (the per-id `--set-active` override) carves out the exceptions. It keys on `hasActivityCandidate`, so it also silences subunit-derived candidates (subunits are suppressed under `--set-all-inactive`).",
        ),
    ),
    variants={
        "celldesigner": (
            Rule(
                identifier="activity:core:celldesigner:from_global_activate",
                text=dedent("""\
                    hasActivityCandidate(SPECIES, isGlobalActive) :-
                        species(SPECIES),
                        globalActivate,
                        not hasSubunit(_, SPECIES)."""),
                docs="Under the set-all-active toggle, every top-level species is an activity candidate. It fires with reason `isGlobalActive` when `globalActivate` (`--set-all-active`) is set. The `not hasSubunit(_, SPECIES)` guard excludes subunits (in CellDesigner a subunit is a species): a subunit is a structural component of its complex, never a top-level activity, so the toggle activates the complex, not its parts. This is the subunit asymmetry -- subunits are *not* activated under `--set-all-active`, though they *are* suppressed under `--set-all-inactive`.",
            ),
        ),
        "sbgn_pd": (
            Rule(
                identifier="activity:core:sbgn_pd:from_global_activate",
                text=dedent("""\
                    hasActivityCandidate(ENTITY_POOL, isGlobalActive) :-
                        entityPool(ENTITY_POOL),
                        globalActivate."""),
                docs="Under the set-all-active toggle, every entity pool is an activity candidate. It fires with reason `isGlobalActive` when `globalActivate` (`--set-all-active`) is set -- the SBGN PD parallel of CellDesigner's global-activate rule. No `not hasSubunit` guard is needed: in SBGN PD a subunit is an `SBGNAuxiliaryUnit`, not an `entityPool`, so the `entityPool` guard already excludes subunits. This preserves the subunit asymmetry -- subunits are not activated under `--set-all-active`.",
            ),
        ),
    },
)

_ACTIVITY_PHENOTYPE = _activity_feature(
    "phenotype",
    base=(
        Rule(
            identifier="activity:phenotype:from_phenotype",
            text="hasActivityCandidate(PHENOTYPE, isPhenotype) :- phenotype(PHENOTYPE).",
            docs="A phenotype is always an activity candidate. Both languages emit the `phenotype` signal, so this rule is shared; it fires with reason `isPhenotype`.",
        ),
    ),
    docs="Excludable: a phenotype is an activity candidate (reason `isPhenotype`). Exclude with `--exclude-group activity:phenotype` to stop treating phenotypes as activities.",
)

_ACTIVITY_ACTIVE_MARKER = _activity_feature(
    "active_marker",
    cd=(
        Rule(
            identifier="activity:active_marker:celldesigner:from_active_flag",
            text="hasActivityCandidate(SPECIES, isActive) :- species(SPECIES), hasActive(SPECIES, 1).",
            docs="A species explicitly marked as active is an activity candidate. This fires when its `hasActive` flag is 1 (reason `isActive`).",
        ),
        Rule(
            identifier="activity:active_marker:celldesigner:from_active_structural_state",
            text=dedent("""\
                hasActivityCandidate(SPECIES, hasActiveStructuralState) :-
                    species(SPECIES),
                    hasStructuralState(SPECIES, STRUCTURAL_STATE),
                    hasValue(STRUCTURAL_STATE, "active")."""),
            docs='A species in an active structural state is an activity candidate. This fires when it carries a structural state whose value is "active" (reason `hasActiveStructuralState`).',
        ),
    ),
    sbgn=(
        Rule(
            identifier="activity:active_marker:sbgn_pd:from_active_state_variable",
            text=dedent("""\
                hasActivityCandidate(ENTITY_POOL, hasActiveStateVariable) :-
                    entityPool(ENTITY_POOL),
                    hasStateVariable(ENTITY_POOL, STATE_VARIABLE),
                    hasValue(STATE_VARIABLE, "active")."""),
            docs='An entity pool in an active state is an activity candidate. This fires when it carries a state variable whose value is "active" (reason `hasActiveStateVariable`) -- the SBGN PD parallel of CellDesigner\'s active structural state.',
        ),
        Rule(
            identifier="activity:active_marker:sbgn_pd:from_active_subunit_state_variable",
            text=dedent("""\
                hasActivityCandidate(SUBUNIT, hasActiveStateVariable) :-
                    isSubunit(SUBUNIT),
                    hasStateVariable(SUBUNIT, STATE_VARIABLE),
                    hasValue(STATE_VARIABLE, "active")."""),
            docs='A subunit in an active state is an activity candidate. This fires when the subunit carries a state variable whose value is "active". Its complex therefore inherits activity (keep-species), and in the ``*-no-complex`` modes the subunit can be promoted. It is the parallel of CellDesigner, where subunits are species and so the active-structural-state rule already covers them.',
        ),
    ),
    docs="Excludable: an explicit active marker makes a species/entity pool an activity candidate (CellDesigner: `hasActive` flag or active structural state; SBGN-PD: active state variable, including on subunits). Exclude with `--exclude-group activity:active_marker`.",
)

_ACTIVITY_MODULATION_SOURCE = _activity_feature(
    "modulation_source",
    cd=(
        Rule(
            identifier="activity:modulation_source:celldesigner:from_modulation_source",
            text=dedent("""\
                hasActivityCandidate(SOURCE, isModulationSource) :-
                    species(SOURCE),
                    knownOrUnknownModulation(MODULATION),
                    hasSource(MODULATION, SOURCE),
                    hasTarget(MODULATION, _)."""),
            docs="A species that is the source of a modulation arc is an activity candidate. This covers both known and unknown modulations that have some target (reason `isModulationSource`); keying on `knownOrUnknownModulation` rather than `modulation` is what lets the source of an unknown modulation become an activity node.",
        ),
        Rule(
            identifier="activity:modulation_source:celldesigner:from_reaction_modulator",
            text=dedent("""\
                hasActivityCandidate(SOURCE, isReactionModulator) :-
                    species(SOURCE),
                    knownOrUnknownModulator(MODULATOR),
                    hasReferredElement(MODULATOR, SOURCE),
                    hasModifier(_, MODULATOR)."""),
            docs="A species that modulates a reaction is an activity candidate. This fires when the species is referred to by a reaction modulator, known or unknown, that modifies some target reaction (reason `isReactionModulator`); keying on `knownOrUnknownModulator` rather than `modulator` is what lets an unknown catalyzer/inhibitor become an activity node.",
        ),
    ),
    sbgn=(
        Rule(
            identifier="activity:modulation_source:sbgn_pd:from_modulation_source",
            text=dedent("""\
                hasActivityCandidate(SOURCE, isModulationSource) :-
                    entityPool(SOURCE),
                    modulation(MODULATION),
                    hasSource(MODULATION, SOURCE),
                    hasTarget(MODULATION, _)."""),
            docs="An entity pool that is the source of a modulation arc is an activity candidate. This fires when the arc's target is a process (reason `isModulationSource`) and is the SBGN PD parallel of both the CellDesigner modulation-arc and reaction-modifier rules.",
        ),
    ),
    docs="Excludable: the source of a modulation arc (CellDesigner: a modulation arc or a reaction modifier, which in SBGN-PD folds into the same arc-to-process case) is an activity candidate. Exclude with `--exclude-group activity:modulation_source`.",
)

_ACTIVITY_GATE_INPUT = _activity_feature(
    "gate_input",
    cd=(
        Rule(
            identifier="activity:gate_input:celldesigner:from_gate_input",
            text=dedent("""\
                hasActivityCandidate(ELEMENT, isGateInput) :-
                    booleanLogicGateInput(INPUT),
                    hasReferredElement(INPUT, ELEMENT)."""),
            docs="A species feeding a boolean logic gate is an activity candidate. It fires with reason `isGateInput`. A gate is structurally always an influence/modulation source or reaction modifier, so each of its inputs is an active driver of the downstream target -- a semantic guarantee, not a fallback. The element is activated regardless of kind (species, complex, ion, ...); the carrier rules route each kind. `hasReferredElement` in CellDesigner only relates a gate input to its species, so no extra guard is needed.",
        ),
    ),
    sbgn=(
        Rule(
            identifier="activity:gate_input:sbgn_pd:from_operator_input",
            text=dedent("""\
                hasActivityCandidate(ELEMENT, isGateInput) :-
                    logicalOperatorInput(INPUT),
                    hasReferredElement(INPUT, ELEMENT),
                    entityPool(ELEMENT)."""),
            docs="An entity pool feeding a logical operator is an activity candidate. It fires with reason `isGateInput` -- the SBGN PD parallel of CellDesigner's gate-input activation. The `entityPool` guard excludes the deferred nested-operator case: an operator feeding another operator has no activity carrier, so such an input simply dangles. In SBGN PD `hasReferredElement` relates many roles (reactant, product, modulation participant) to entities, so both the `logicalOperatorInput` and `entityPool` guards are needed.",
        ),
    ),
    docs="Excludable: an input feeding a boolean logic gate / logical operator is an activity candidate (reason `isGateInput`). Exclude with `--exclude-group activity:gate_input`.",
)

_TOPOLOGY = RuleGroup(
    identifier="topology",
    depends_on=frozenset({"activity:core"}),
    docs="Mode-agnostic structural helpers shared by every mode: `isSubunit`, `hasActiveDescendantSubunit`.",
    rules=(
        Rule(
            identifier="topology:is_subunit",
            text="isSubunit(SUBUNIT) :- hasSubunit(_, SUBUNIT).",
            docs="A species is a subunit if it appears on the right side of any `hasSubunit` relation.",
        ),
        Rule(
            identifier="topology:has_active_descendant_direct",
            text=dedent("""\
                hasActiveDescendantSubunit(COMPLEX) :-
                    complex(COMPLEX),
                    hasSubunit(COMPLEX, SUBUNIT),
                    hasActivity(SUBUNIT, _)."""),
            docs="A complex has an active descendant if any direct subunit has activity.",
        ),
        Rule(
            identifier="topology:has_active_descendant_transitive",
            text=dedent("""\
                hasActiveDescendantSubunit(COMPLEX) :-
                    complex(COMPLEX),
                    hasSubunit(COMPLEX, NESTED_COMPLEX),
                    hasActiveDescendantSubunit(NESTED_COMPLEX)."""),
            docs="The `hasActiveDescendantSubunit` relation is transitive through complex containment.",
        ),
    ),
)

_TOP_LEVEL = RuleGroup(
    identifier="top_level",
    depends_on=frozenset({"topology"}),
    docs="Resolves every species to its outermost top-level entity: a non-subunit resolves to itself; a subunit -- at any nesting depth -- resolves to the outermost complex that contains it. The complex-keeping modes (`keep-species`, `normal`, `keep-reactions`) key a species by its top-level entity, so a subunit is never its own activity and its influences attach to its top-level complex (a subunit is a structural component, not an independent influencer).",
    rules=(
        Rule(
            identifier="top_level:recursive",
            text=dedent("""\
                resolvesToTopLevel(SUBUNIT, TOPLEVEL) :-
                    hasSubunit(PARENT_COMPLEX, SUBUNIT),
                    resolvesToTopLevel(PARENT_COMPLEX, TOPLEVEL)."""),
            docs="A subunit resolves to the same top-level entity as its parent complex, recursively through nested complexes -- so a subunit at any depth resolves to its outermost complex.",
        ),
    ),
    variants={
        "celldesigner": (
            Rule(
                identifier="top_level:celldesigner:self",
                text=dedent("""\
                    resolvesToTopLevel(SPECIES, SPECIES) :-
                        species(SPECIES),
                        not hasSubunit(_, SPECIES)."""),
                docs="CellDesigner: a species that is not a subunit of any complex is its own top-level entity.",
            ),
        ),
        "sbgn_pd": (
            Rule(
                identifier="top_level:sbgn_pd:self",
                text=dedent("""\
                    resolvesToTopLevel(ENTITY_POOL, ENTITY_POOL) :-
                        entityPool(ENTITY_POOL),
                        not hasSubunit(_, ENTITY_POOL)."""),
                docs="SBGN-PD: an entity pool that is not a subunit of any complex is its own top-level entity.",
            ),
            Rule(
                identifier="top_level:sbgn_pd:phenotype_self",
                text="resolvesToTopLevel(PHENOTYPE, PHENOTYPE) :- phenotype(PHENOTYPE).",
                docs="SBGN-PD: a phenotype is a process, not an entity pool, so the entity-pool self-rule never keys it; a phenotype is never a subunit, so it is always its own top-level entity. Without this an SBGN phenotype gets `hasActivity` but no activity key and is silently dropped from `keep-species`/`normal` output. The `*-no-complex` modes key via `not isSubunit`/`not delete` and already include phenotypes; CellDesigner phenotypes are species (covered by the species self-rule).",
            ),
        ),
    },
)

_PREPARATION_COMPLEX = RuleGroup(
    identifier="preparation:complex",
    slot="preparation",
    depends_on=frozenset({"activity:core", "topology", "top_level"}),
    docs="The complex-keeping modes (`normal`, `keep-species`, `keep-reactions`) key a species with activity by the `keptSpeciesKey` of its top-level entity (the `top_level` group): itself when top-level, its outermost complex when a subunit. A subunit therefore contributes no activity of its own -- it is a structural component of its complex, and any influence it carries attaches to the top-level complex (the active descendant directly keys the complex, so the assembly is represented without a separate inherit-activity rule). Carriers are identity; the rerouting lives entirely in the key. `normal` and `keep-species` share these keys and differ only at the build stage, where `normal` strips PTM decorations and merges content-equal results while `keep-species` keeps the decorations.",
    rules=(
        Rule(
            identifier="preparation:complex:key",
            text=dedent("""\
                hasActivityKey(SPECIES, keptSpeciesKey(TOPLEVEL)) :-
                    hasActivity(SPECIES, _),
                    resolvesToTopLevel(SPECIES, TOPLEVEL)."""),
            docs="A species with activity is keyed by the `keptSpeciesKey` of its top-level entity (itself when top-level; its outermost complex when a subunit).",
        ),
    ),
    # The activity carrier is a per-language variant: a CellDesigner species is
    # its own carrier; an SBGN-PD entity pool or phenotype is its own carrier
    # (phenotypes are processes, not entity pools, so they need their own rule).
    variants={
        "celldesigner": (
            Rule(
                identifier="preparation:complex:celldesigner:carrier",
                text="hasActivityCarrier(SPECIES, SPECIES) :- species(SPECIES).",
                docs="Each species is its own activity carrier (no rerouting; a subunit endpoint is rerouted to its top-level complex by the `top_level` key, not by the carrier).",
            ),
        ),
        "sbgn_pd": (
            Rule(
                identifier="preparation:complex:sbgn_pd:carrier_entity_pool",
                text="hasActivityCarrier(ENTITY_POOL, ENTITY_POOL) :- entityPool(ENTITY_POOL).",
                docs="SBGN-PD: each entity pool is its own activity carrier.",
            ),
            Rule(
                identifier="preparation:complex:sbgn_pd:carrier_phenotype",
                text="hasActivityCarrier(PHENOTYPE, PHENOTYPE) :- phenotype(PHENOTYPE).",
                docs="SBGN-PD: each phenotype process is its own activity carrier (phenotypes are activities, not entity pools).",
            ),
        ),
    },
)

_PREPARATION_NO_COMPLEX = RuleGroup(
    identifier="preparation:no_complex",
    slot="preparation",
    depends_on=frozenset({"activity:core", "topology"}),
    docs="The complex-dissolving modes (`normal-no-complex`, `keep-species-no-complex`): a complex with any (transitive) active descendant is deleted; a non-deleted top-level species with activity is keyed by `keptSpeciesKey(SELF)`; a subunit of a deleted complex is promoted to top level, keyed by `promotedSubunitKey(SELF)` (paths reach it via `paths_complex_traversal`). `normal-no-complex` and `keep-species-no-complex` share these keys and differ only at the build stage, where `normal-no-complex` strips PTM decorations and merges content-equal results while `keep-species-no-complex` keeps the decorations.",
    rules=(
        Rule(
            identifier="preparation:no_complex:deleted",
            text=dedent("""\
                delete(COMPLEX) :-
                    complex(COMPLEX),
                    hasActiveDescendantSubunit(COMPLEX)."""),
            docs="A complex with any (transitive) active descendant is deleted.",
        ),
        Rule(
            identifier="preparation:no_complex:subunit_of_deleted_direct",
            text=dedent("""\
                isDescendantSubunitOfDeleted(SUBUNIT) :-
                    hasSubunit(COMPLEX, SUBUNIT),
                    delete(COMPLEX)."""),
            docs="A direct subunit of a deleted complex is itself subunit-of-deleted.",
        ),
        Rule(
            identifier="preparation:no_complex:subunit_of_deleted_transitive",
            text=dedent("""\
                isDescendantSubunitOfDeleted(SUBUNIT) :-
                    hasSubunit(PARENT_COMPLEX, SUBUNIT),
                    isDescendantSubunitOfDeleted(PARENT_COMPLEX)."""),
            docs="Subunit-of-deleted is transitive through nested complexes.",
        ),
        Rule(
            identifier="preparation:no_complex:key_top_level",
            text=dedent("""\
                hasActivityKey(SPECIES, keptSpeciesKey(SPECIES)) :-
                    hasActivity(SPECIES, _),
                    not isSubunit(SPECIES),
                    not delete(SPECIES)."""),
            docs="A non-deleted top-level species with activity is keyed by `keptSpeciesKey(SELF)`.",
        ),
        Rule(
            identifier="preparation:no_complex:key_promoted_subunit",
            text=dedent("""\
                hasActivityKey(SPECIES, promotedSubunitKey(SPECIES)) :-
                    hasActivity(SPECIES, _),
                    isDescendantSubunitOfDeleted(SPECIES)."""),
            docs="A subunit of a deleted complex with activity is promoted to top level, keyed by `promotedSubunitKey(SELF)`.",
        ),
    ),
    # The activity carrier is a per-language variant. In SBGN-PD a promoted
    # subunit must also be its own carrier so a rerouted path (a complex's
    # influence propagated to its subunit via `paths_complex_traversal`)
    # becomes an influence on the promoted subunit.
    variants={
        "celldesigner": (
            Rule(
                identifier="preparation:no_complex:celldesigner:carrier",
                text="hasActivityCarrier(SPECIES, SPECIES) :- species(SPECIES).",
                docs="Each species is its own activity carrier.",
            ),
        ),
        "sbgn_pd": (
            Rule(
                identifier="preparation:no_complex:sbgn_pd:carrier_entity_pool",
                text="hasActivityCarrier(ENTITY_POOL, ENTITY_POOL) :- entityPool(ENTITY_POOL).",
                docs="SBGN-PD: each entity pool is its own activity carrier.",
            ),
            Rule(
                identifier="preparation:no_complex:sbgn_pd:carrier_phenotype",
                text="hasActivityCarrier(PHENOTYPE, PHENOTYPE) :- phenotype(PHENOTYPE).",
                docs="SBGN-PD: each phenotype process is its own activity carrier.",
            ),
            Rule(
                identifier="preparation:no_complex:sbgn_pd:carrier_subunit",
                text="hasActivityCarrier(SUBUNIT, SUBUNIT) :- isSubunit(SUBUNIT).",
                docs="SBGN-PD: each subunit is its own activity carrier, so a promoted subunit can be an influence endpoint (paths reach it via `paths_complex_traversal`).",
            ),
        ),
    },
)

_MODULATION_KIND = RuleGroup(
    identifier="modulation_kind",
    docs="Maps each modulation arc to its influence kind via `hasModulationKind(MODULATION, INFLUENCE_KIND)`, the single place the arc-type->kind knowledge lives. Every group that emits a modulation-arc influence reads this relation rather than re-encoding the mapping. The mapping is per-language: CellDesigner arcs carry the sign in the arc type (catalysis, inhibition, ...), while an SBGN-PD arc's kind comes from its stimulation/inhibition/necessary-stimulation classification.",
    rules=(),
    variants={
        "celldesigner": (
            Rule(
                identifier="modulation_kind:celldesigner:catalysis",
                text="hasModulationKind(MODULATION, positivelyInfluences) :- catalysis(MODULATION).",
                docs="A catalysis arc contributes a `positive` kind.",
            ),
            Rule(
                identifier="modulation_kind:celldesigner:physical_stimulation",
                text="hasModulationKind(MODULATION, positivelyInfluences) :- physicalStimulation(MODULATION).",
                docs="A physical stimulation arc contributes a `positive` kind.",
            ),
            Rule(
                identifier="modulation_kind:celldesigner:positive_influence",
                text="hasModulationKind(MODULATION, positivelyInfluences) :- positiveInfluence(MODULATION).",
                docs="A positive influence arc contributes a `positive` kind.",
            ),
            Rule(
                identifier="modulation_kind:celldesigner:triggering",
                text="hasModulationKind(MODULATION, triggers) :- triggering(MODULATION).",
                docs="A triggering arc contributes a `triggering` kind (kept on the direct edge; it degrades to `positive` when composed through a reaction).",
            ),
            Rule(
                identifier="modulation_kind:celldesigner:inhibition",
                text="hasModulationKind(MODULATION, negativelyInfluences) :- inhibition(MODULATION).",
                docs="An inhibition arc contributes a `negative` kind.",
            ),
            Rule(
                identifier="modulation_kind:celldesigner:negative_influence",
                text="hasModulationKind(MODULATION, negativelyInfluences) :- negativeInfluence(MODULATION).",
                docs="A negative influence arc contributes a `negative` kind.",
            ),
            Rule(
                identifier="modulation_kind:celldesigner:modulation",
                text=dedent("""\
                    hasModulationKind(MODULATION, modulates) :-
                        modulation(MODULATION),
                        not catalysis(MODULATION),
                        not physicalStimulation(MODULATION),
                        not inhibition(MODULATION),
                        not triggering(MODULATION),
                        not positiveInfluence(MODULATION),
                        not negativeInfluence(MODULATION)."""),
                docs="A bare modulation arc (the generic MODULATION arc, not one of the signed/typed subtypes) contributes an unknown-sign `modulation` kind. The negations exclude the subtypes, which `modulation` is the umbrella over.",
            ),
            Rule(
                identifier="modulation_kind:celldesigner:unknown_catalysis",
                text="hasModulationKind(MODULATION, unknownPositivelyInfluences) :- unknownCatalysis(MODULATION).",
                docs="An unknown catalysis arc contributes an `unknown_positive` kind.",
            ),
            Rule(
                identifier="modulation_kind:celldesigner:unknown_physical_stimulation",
                text="hasModulationKind(MODULATION, unknownPositivelyInfluences) :- unknownPhysicalStimulation(MODULATION).",
                docs="An unknown physical stimulation arc contributes an `unknown_positive` kind.",
            ),
            Rule(
                identifier="modulation_kind:celldesigner:unknown_positive_influence",
                text="hasModulationKind(MODULATION, unknownPositivelyInfluences) :- unknownPositiveInfluence(MODULATION).",
                docs="An unknown positive influence arc contributes an `unknown_positive` kind.",
            ),
            Rule(
                identifier="modulation_kind:celldesigner:unknown_triggering",
                text="hasModulationKind(MODULATION, unknownTriggers) :- unknownTriggering(MODULATION).",
                docs="An unknown triggering arc contributes an `unknown_triggering` kind (kept on the direct edge; it degrades to `unknown_positive` when composed through a reaction).",
            ),
            Rule(
                identifier="modulation_kind:celldesigner:unknown_inhibition",
                text="hasModulationKind(MODULATION, unknownNegativelyInfluences) :- unknownInhibition(MODULATION).",
                docs="An unknown inhibition arc contributes an `unknown_negative` kind.",
            ),
            Rule(
                identifier="modulation_kind:celldesigner:unknown_negative_influence",
                text="hasModulationKind(MODULATION, unknownNegativelyInfluences) :- unknownNegativeInfluence(MODULATION).",
                docs="An unknown negative influence arc contributes an `unknown_negative` kind.",
            ),
            Rule(
                identifier="modulation_kind:celldesigner:unknown_modulation",
                text=dedent("""\
                    hasModulationKind(MODULATION, unknownModulates) :-
                        unknownModulation(MODULATION),
                        not unknownCatalysis(MODULATION),
                        not unknownPhysicalStimulation(MODULATION),
                        not unknownInhibition(MODULATION),
                        not unknownTriggering(MODULATION),
                        not unknownPositiveInfluence(MODULATION),
                        not unknownNegativeInfluence(MODULATION)."""),
                docs="A bare unknown modulation arc (not one of the unknown subtypes) contributes an `unknown_modulation` kind. The negations exclude the subtypes, which `unknownModulation` is the umbrella over.",
            ),
        ),
        "sbgn_pd": (
            Rule(
                identifier="modulation_kind:sbgn_pd:necessary_stimulation",
                text="hasModulationKind(MODULATION, triggers) :- necessaryStimulation(MODULATION).",
                docs="An SBGN-PD necessary stimulation contributes a `triggering` kind.",
            ),
            Rule(
                identifier="modulation_kind:sbgn_pd:stimulation",
                text=dedent("""\
                    hasModulationKind(MODULATION, positivelyInfluences) :-
                        stimulation(MODULATION),
                        not necessaryStimulation(MODULATION)."""),
                docs="A stimulation (catalysis included) that is not a necessary stimulation contributes a `positive` kind.",
            ),
            Rule(
                identifier="modulation_kind:sbgn_pd:inhibition",
                text="hasModulationKind(MODULATION, negativelyInfluences) :- inhibition(MODULATION).",
                docs="An SBGN-PD inhibition contributes a `negative` kind.",
            ),
            Rule(
                identifier="modulation_kind:sbgn_pd:modulation",
                text=dedent("""\
                    hasModulationKind(MODULATION, modulates) :-
                        modulation(MODULATION),
                        not stimulation(MODULATION),
                        not inhibition(MODULATION)."""),
                docs="A bare modulation (neither stimulation nor inhibition) contributes an unknown-sign `modulation` kind.",
            ),
        ),
    },
)

# SBGN-PD path rules. Unlike CellDesigner (where a modulation arc connects two
# species), an SBGN-PD modulation arc's source is an entity pool and its target
# is a *process*; the influence therefore propagates to the process's products
# (or, when the target is a phenotype, to the phenotype itself). These live in
# the `sbgn_pd` variant of `paths:core` because the CellDesigner modulation-arc
# rules would misfire on PD facts (their `catalysis`/`inhibition`/`modulation`
# functors match PD modulations, but PD targets are processes, yielding spurious
# edges).
# The SBGN-PD `paths:core` variant rules: the single reactant->product hop and
# the two modulation-arc emitters. The multi-hop `transitive_through_process`
# rule is excludable and lives in the `paths:chaining` group's SBGN-PD variant.
_PATHS_CORE_SBGN_PD = (
    Rule(
        identifier="paths:core:sbgn_pd:modulation_to_product",
        text=dedent("""\
            propagatesInfluence(SOURCE_ENTITY_POOL, TARGET_ENTITY_POOL, INFLUENCE_KIND) :-
                hasModulationKind(MODULATION, INFLUENCE_KIND),
                hasSource(MODULATION, SOURCE_ENTITY_POOL),
                hasTarget(MODULATION, PROCESS),
                hasProduct(PROCESS, PRODUCT),
                hasReferredElement(PRODUCT, TARGET_ENTITY_POOL)."""),
        docs="A modulation arc's entity-pool source influences each product of its target process, carrying the arc's kind.",
    ),
    Rule(
        identifier="paths:core:sbgn_pd:modulation_to_phenotype",
        text=dedent("""\
            propagatesInfluence(SOURCE_ENTITY_POOL, TARGET_ENTITY_POOL, INFLUENCE_KIND) :-
                hasModulationKind(MODULATION, INFLUENCE_KIND),
                hasSource(MODULATION, SOURCE_ENTITY_POOL),
                hasTarget(MODULATION, TARGET_ENTITY_POOL),
                phenotype(TARGET_ENTITY_POOL)."""),
        docs="A modulation arc whose target is a phenotype influences the phenotype itself (a phenotype process has no products; it is the activity).",
    ),
    Rule(
        identifier="paths:core:sbgn_pd:is_directly_transformed_to",
        text=dedent("""\
            isDirectlyTransformedTo(UPSTREAM_ENTITY_POOL, DOWNSTREAM_ENTITY_POOL) :-
                hasReactant(PROCESS, REACTANT),
                hasReferredElement(REACTANT, UPSTREAM_ENTITY_POOL),
                hasProduct(PROCESS, PRODUCT),
                hasReferredElement(PRODUCT, DOWNSTREAM_ENTITY_POOL)."""),
        docs="The single reactant->product hop in the production graph: the upstream entity pool is an element of a reactant and the downstream pool an element of a product of the same process. Passive voice (the process does the transforming, not the pool) keeps it language-neutral. Feeds the shared `isTransformedTo`/`isCyclicallyTransformedTo` cycle relations that gate transitive path extension; carries no influence kind itself.",
    ),
)


_PATHS_CORE = RuleGroup(
    identifier="paths:core",
    depends_on=frozenset({"modulation_kind"}),
    docs="Builds the direct (single-hop) kinded `propagatesInfluence(SOURCE_SPECIES, TARGET_SPECIES, INFLUENCE_KIND)` relation from PD reactions and modulation arcs, plus the production-cycle relations that gate transitive extension. INFLUENCE_KIND is one of `positive`, `negative`, `triggering`, `modulation` and their `unknown_*` twins. Reaction modifiers and the matching species→species modulation arcs map to the *same* kind (e.g. a trigger modifier and a triggering arc both give `triggering`; catalysis and physical stimulation both give `positive`; inhibition gives `negative`). Each variant defines the single reactant->product hop `isDirectlyTransformedTo/2`, the shared `isTransformedTo/2` is its transitive closure, and `isCyclicallyTransformedTo/2` marks the hops that lie inside a cycle so transitivity never propagates influence back around a production loop. The multi-hop reactant-chained transitivity itself is the excludable `paths:chaining` group.",
    rules=(
        Rule(
            identifier="paths:core:is_transformed_to",
            text=dedent("""\
                isTransformedTo(UPSTREAM, DOWNSTREAM) :- isDirectlyTransformedTo(UPSTREAM, DOWNSTREAM).
                isTransformedTo(UPSTREAM, DOWNSTREAM) :-
                    isTransformedTo(UPSTREAM, INTERMEDIATE),
                    isDirectlyTransformedTo(INTERMEDIATE, DOWNSTREAM)."""),
            docs="Transitive closure of the single-hop `isDirectlyTransformedTo/2`: `isTransformedTo(UPSTREAM, DOWNSTREAM)` holds when DOWNSTREAM is reachable from UPSTREAM through one or more reactant->product hops (the direct hop included, as the base case). Language-agnostic; `isDirectlyTransformedTo/2` is supplied per language variant. Depends only on `isDirectlyTransformedTo`, never on `propagatesInfluence`, so the negation that gates transitivity stays stratified.",
        ),
        Rule(
            identifier="paths:core:is_cyclically_transformed_to",
            text=dedent("""\
                isCyclicallyTransformedTo(UPSTREAM, DOWNSTREAM) :-
                    isTransformedTo(UPSTREAM, DOWNSTREAM),
                    isTransformedTo(DOWNSTREAM, UPSTREAM)."""),
            docs="UPSTREAM and DOWNSTREAM are mutually reachable through the production graph -- both lie on a common cycle. At every transitivity guard site the hop in question is already a direct reactant->product edge, so requiring mutual reachability there is exactly equivalent to 'this hop is on a cycle'. The transitivity rules forbid extending a path across such a hop, blocking influence from leaking around production loops.",
        ),
        Rule(
            identifier="paths:core:composes_to",
            text=dedent("""\
                composesTo(positivelyInfluences, positivelyInfluences).
                composesTo(negativelyInfluences, negativelyInfluences).
                composesTo(modulates, modulates).
                composesTo(triggers, positivelyInfluences).
                composesTo(unknownPositivelyInfluences, unknownPositivelyInfluences).
                composesTo(unknownNegativelyInfluences, unknownNegativelyInfluences).
                composesTo(unknownModulates, unknownModulates).
                composesTo(unknownTriggers, unknownPositivelyInfluences)."""),
            docs="`composesTo(INCOMING_INFLUENCE_KIND, OUTGOING_INFLUENCE_KIND)`: how an influence kind transforms when a path is extended by one reactant→product hop (a reaction in CellDesigner, a process in SBGN-PD). Identity for every kind except `triggering` (a necessary-stimulation relationship is a property of the direct edge; composed through a reaction it weakens to a plain `positive` influence) and its unknown twin `unknown_triggering` → `unknown_positive`. Modulation composes like the signed kinds (stays `modulation`).",
        ),
    ),
    variants={
        "sbgn_pd": _PATHS_CORE_SBGN_PD,
        "celldesigner": (
        Rule(
            identifier="paths:core:celldesigner:catalyzer_to_product",
            text=dedent("""\
                propagatesInfluence(SOURCE_SPECIES, TARGET_SPECIES, positivelyInfluences) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER),
                    catalyzer(MODIFIER),
                    hasReferredElement(MODIFIER, SOURCE_SPECIES),
                    hasProduct(REACTION, PRODUCT),
                    hasReferredElement(PRODUCT, TARGET_SPECIES)."""),
            docs="If a species is referred to by a catalyzer of a reaction and another species is referred to by a product of that reaction, then there is a positive path from the first to the second.",
        ),
        Rule(
            identifier="paths:core:celldesigner:physical_stimulator_to_product",
            text=dedent("""\
                propagatesInfluence(SOURCE_SPECIES, TARGET_SPECIES, positivelyInfluences) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER),
                    physicalStimulator(MODIFIER),
                    hasReferredElement(MODIFIER, SOURCE_SPECIES),
                    hasProduct(REACTION, PRODUCT),
                    hasReferredElement(PRODUCT, TARGET_SPECIES)."""),
            docs="If a species is referred to by a physical stimulator of a reaction and another species is referred to by a product of that reaction, then there is a positive path from the first to the second.",
        ),
        Rule(
            identifier="paths:core:celldesigner:trigger_to_product",
            text=dedent("""\
                propagatesInfluence(SOURCE_SPECIES, TARGET_SPECIES, triggers) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER),
                    trigger(MODIFIER),
                    hasReferredElement(MODIFIER, SOURCE_SPECIES),
                    hasProduct(REACTION, PRODUCT),
                    hasReferredElement(PRODUCT, TARGET_SPECIES)."""),
            docs="If a species is referred to by a trigger of a reaction and another species is referred to by a product of that reaction, then there is a triggering path from the first to the second (a trigger→product edge is direct, so it keeps the `triggering` kind).",
        ),
        Rule(
            identifier="paths:core:celldesigner:modulation_arc_influence",
            text=dedent("""\
                propagatesInfluence(SOURCE_SPECIES, TARGET_SPECIES, INFLUENCE_KIND) :-
                    hasModulationKind(MODULATION, INFLUENCE_KIND),
                    hasSource(MODULATION, SOURCE_SPECIES),
                    hasTarget(MODULATION, TARGET_SPECIES)."""),
            docs="A modulation arc influences its target with the arc's kind (`hasModulationKind`, from the `modulation_kind` group). This covers every arc type -- catalysis/physical stimulation/positive influence give `positive`, triggering gives `triggering`, inhibition/negative influence give `negative`, a bare modulation gives `modulation`, and the `unknown*` twins give the `unknown_*` kinds. The reaction modifier->product rules and the transitive/cycle rules remain separate.",
        ),
        Rule(
            identifier="paths:core:celldesigner:inhibitor_to_product",
            text=dedent("""\
                propagatesInfluence(SOURCE_SPECIES, TARGET_SPECIES, negativelyInfluences) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER),
                    inhibitor(MODIFIER),
                    hasReferredElement(MODIFIER, SOURCE_SPECIES),
                    hasProduct(REACTION, PRODUCT),
                    hasReferredElement(PRODUCT, TARGET_SPECIES)."""),
            docs="If a species is referred to by an inhibitor of a reaction and another species is referred to by a product of that reaction, then there is a negative path from the first to the second.",
        ),
        Rule(
            identifier="paths:core:celldesigner:modulator_to_product",
            text=dedent("""\
                propagatesInfluence(SOURCE_SPECIES, TARGET_SPECIES, modulates) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER),
                    modulator(MODIFIER),
                    not physicalStimulator(MODIFIER),
                    not inhibitor(MODIFIER),
                    not trigger(MODIFIER),
                    hasReferredElement(MODIFIER, SOURCE_SPECIES),
                    hasProduct(REACTION, PRODUCT),
                    hasReferredElement(PRODUCT, TARGET_SPECIES)."""),
            docs="If a species is referred to by a *bare* modulator of a reaction (a generic MODULATION modifier — not a physical stimulator, inhibitor or trigger; catalyzers are physical stimulators) and another species is referred to by a product of that reaction, then there is a modulation path from the first to the second.",
        ),
        Rule(
            identifier="paths:core:celldesigner:unknown_catalyzer_to_product",
            text=dedent("""\
                propagatesInfluence(SOURCE_SPECIES, TARGET_SPECIES, unknownPositivelyInfluences) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER),
                    unknownCatalyzer(MODIFIER),
                    hasReferredElement(MODIFIER, SOURCE_SPECIES),
                    hasProduct(REACTION, PRODUCT),
                    hasReferredElement(PRODUCT, TARGET_SPECIES)."""),
            docs="If a species is referred to by an unknown catalyzer of a reaction and another species is referred to by a product of that reaction, then there is an unknown-positive path from the first to the second.",
        ),
        Rule(
            identifier="paths:core:celldesigner:unknown_inhibitor_to_product",
            text=dedent("""\
                propagatesInfluence(SOURCE_SPECIES, TARGET_SPECIES, unknownNegativelyInfluences) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER),
                    unknownInhibitor(MODIFIER),
                    hasReferredElement(MODIFIER, SOURCE_SPECIES),
                    hasProduct(REACTION, PRODUCT),
                    hasReferredElement(PRODUCT, TARGET_SPECIES)."""),
            docs="If a species is referred to by an unknown inhibitor of a reaction and another species is referred to by a product of that reaction, then there is an unknown-negative path from the first to the second.",
        ),
        Rule(
            identifier="paths:core:celldesigner:is_directly_transformed_to",
            text=dedent("""\
                isDirectlyTransformedTo(UPSTREAM_SPECIES, DOWNSTREAM_SPECIES) :-
                    reaction(REACTION),
                    hasReactant(REACTION, REACTANT),
                    hasReferredElement(REACTANT, UPSTREAM_SPECIES),
                    hasProduct(REACTION, PRODUCT),
                    hasReferredElement(PRODUCT, DOWNSTREAM_SPECIES)."""),
            docs="The single reactant->product hop in the production graph: the upstream species is referred to by a reactant and the downstream species by a product of the same reaction. Passive voice (the reaction does the transforming, not the species) keeps it language-neutral. Feeds the shared `isTransformedTo`/`isCyclicallyTransformedTo` cycle relations that gate transitive path extension; carries no influence kind itself.",
        ),
        ),
    },
)

# Multi-hop path chaining: extends a `propagatesInfluence` path across one more
# reactant->product hop. Excludable -- dropping it leaves only single-hop
# influences. It consumes `composesTo` and `not isCyclicallyTransformedTo`, both
# supplied by `paths:core`, so excluding it never breaks the core relations.
_PATHS_CHAINING = RuleGroup(
    identifier="paths:chaining",
    depends_on=frozenset({"paths:core"}),
    docs="Excludable multi-hop transitivity: extends a `propagatesInfluence` path through one more reactant->product hop (a reaction in CellDesigner, a process in SBGN-PD), carrying the kind via `composesTo` (triggering degrades to positivelyInfluences) and refusing to extend across a hop inside a production cycle (`not isCyclicallyTransformedTo`). Exclude with `--exclude-group paths:chaining` to keep only direct single-hop influences. Depends on `paths:core`, which supplies both `composesTo` and the cycle relations.",
    rules=(),
    variants={
        "celldesigner": (
            Rule(
                identifier="paths:chaining:celldesigner:transitive_through_reaction",
                text=dedent("""\
                    propagatesInfluence(SOURCE_SPECIES, TARGET_SPECIES, OUTGOING_INFLUENCE_KIND) :-
                        propagatesInfluence(SOURCE_SPECIES, INTERMEDIATE_SPECIES, INCOMING_INFLUENCE_KIND),
                        composesTo(INCOMING_INFLUENCE_KIND, OUTGOING_INFLUENCE_KIND),
                        reaction(REACTION),
                        hasReactant(REACTION, REACTANT),
                        hasReferredElement(REACTANT, INTERMEDIATE_SPECIES),
                        hasProduct(REACTION, PRODUCT),
                        hasReferredElement(PRODUCT, TARGET_SPECIES),
                        not isCyclicallyTransformedTo(INTERMEDIATE_SPECIES, TARGET_SPECIES)."""),
                docs="If there is a path from a species to an intermediate species, and the intermediate species is referred to by a reactant of a reaction whose product refers to another species, then there is a path from the first species to the second. The kind is carried through `composesTo`, which degrades `triggering`→`positive` (and `unknown_triggering`→`unknown_positive`) at the reaction hop while leaving every other kind unchanged. Extension is suppressed across a reactant->product hop that lies inside a cycle (`isCyclicallyTransformedTo`), so a source feeding a production cycle does not leak influence back around the loop onto members it directly depletes.",
            ),
        ),
        "sbgn_pd": (
            Rule(
                identifier="paths:chaining:sbgn_pd:transitive_through_process",
                text=dedent("""\
                    propagatesInfluence(SOURCE_ENTITY_POOL, TARGET_ENTITY_POOL, OUTGOING_INFLUENCE_KIND) :-
                        propagatesInfluence(SOURCE_ENTITY_POOL, INTERMEDIATE_ENTITY_POOL, INCOMING_INFLUENCE_KIND),
                        composesTo(INCOMING_INFLUENCE_KIND, OUTGOING_INFLUENCE_KIND),
                        hasReactant(PROCESS, REACTANT),
                        hasReferredElement(REACTANT, INTERMEDIATE_ENTITY_POOL),
                        hasProduct(PROCESS, PRODUCT),
                        hasReferredElement(PRODUCT, TARGET_ENTITY_POOL),
                        not isCyclicallyTransformedTo(INTERMEDIATE_ENTITY_POOL, TARGET_ENTITY_POOL)."""),
                docs="Extends a path through a process reactant->product hop, carrying the kind via `composesTo` (triggering degrades to positivelyInfluences). Extension is suppressed across a reactant->product hop that lies inside a cycle (`isCyclicallyTransformedTo`), so a source feeding a production cycle does not leak influence back around the loop onto members it directly depletes.",
            ),
        ),
    },
)

_PATHS_COMPLEX_TRAVERSAL = RuleGroup(
    identifier="paths_complex_traversal",
    depends_on=frozenset({"paths:core"}),
    docs="Extends paths through complex containment for the ``*-no-complex`` profiles: a path touching a complex is propagated to/from each of its subunits so influences reach the surviving subunit activities.",
    rules=(
        Rule(
            identifier="paths_complex_traversal:into_subunits",
            text=dedent("""\
                propagatesInfluence(SOURCE, SUBUNIT, INFLUENCE_KIND) :-
                    propagatesInfluence(SOURCE, COMPLEX, INFLUENCE_KIND),
                    complex(COMPLEX),
                    hasSubunit(COMPLEX, SUBUNIT)."""),
            docs="If an influence of a given kind propagates from a source to a complex, then it also propagates from the source to each subunit of the complex.",
        ),
        Rule(
            identifier="paths_complex_traversal:from_subunits",
            text=dedent("""\
                propagatesInfluence(SUBUNIT, TARGET, INFLUENCE_KIND) :-
                    propagatesInfluence(COMPLEX, TARGET, INFLUENCE_KIND),
                    complex(COMPLEX),
                    hasSubunit(COMPLEX, SUBUNIT)."""),
            docs="If an influence of a given kind propagates from a complex to a target, then it also propagates from each subunit of the complex to the target.",
        ),
    ),
)

_INFLUENCES_DERIVATION = RuleGroup(
    identifier="influences_derivation",
    depends_on=frozenset({"paths:core", "preparation"}),
    docs="Derivation: emits `new(activity(KEY))` for every activity key, and lifts every kinded path into the internal `influences(SOURCE_KEY, TARGET_KEY, INFLUENCE_KIND)` relation by resolving both endpoints through their activity carrier and key. Both rules are mode- and language-agnostic; the consumption/sparing reasoning lives in the separate `influences_consumption` group. The internal `influences/3` relation is fanned out to the typed `new(...)` heads by the shared `influence_output` group.",
    rules=(
        Rule(
            identifier="influences_derivation:activity",
            text="new(activity(KEY)) :- hasActivityKey(_, KEY).",
            docs="Every species with an activity key emits an activity node with that key.",
        ),
        Rule(
            identifier="influences_derivation:path",
            text=dedent("""\
                influences(SOURCE_KEY, TARGET_KEY, INFLUENCE_KIND) :-
                    propagatesInfluence(SOURCE, TARGET, INFLUENCE_KIND),
                    hasActivityCarrier(SOURCE, SOURCE_CARRIER),
                    hasActivityCarrier(TARGET, TARGET_CARRIER),
                    hasActivityKey(SOURCE_CARRIER, SOURCE_KEY),
                    hasActivityKey(TARGET_CARRIER, TARGET_KEY)."""),
            docs="A propagated influence of any kind between two raw nodes yields an influence of that same kind between their activity-key images (via their carriers).",
        ),
    ),
)

# The consumption/sparing rules are CellDesigner-only (they reason over reaction
# modifiers); SBGN-PD has no analog yet, so its variant is empty. This is a known
# gap to implement, not a deliberate design choice: the same biology in
# CellDesigner vs SBGN-PD currently yields different AF influences.
_INFLUENCES_CONSUMPTION = RuleGroup(
    identifier="influences_consumption",
    depends_on=frozenset({"activity:core", "preparation"}),
    docs="Excludable consumption/sparing reasoning: a reaction depletes its reactants, so a modifier that drives the reaction also acts on every reactant that is itself an activity -- catalyzer/physicalStimulator/trigger negatively influence each consumed reactant, inhibitor positively influences each spared reactant, and the unknown modifiers contribute the unknown twins. This is inference beyond what the map draws, so it is a group of its own: exclude with `--exclude-group influences_consumption` to keep only the influences the map states. The `keep-reactions` profile omits it, since that mode renders each reaction directly instead of reasoning about it.",
    rules=(),
    variants={
        "sbgn_pd": (),
        "celldesigner": (
        Rule(
            identifier="influences_consumption:celldesigner:catalyzer_consumes_reactant",
            text=dedent("""\
                influences(SOURCE_KEY, TARGET_KEY, negativelyInfluences) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), catalyzer(MODIFIER),
                    hasReferredElement(MODIFIER, SOURCE_SPECIES),
                    hasReactant(REACTION, REACTANT),
                    hasReferredElement(REACTANT, TARGET_SPECIES),
                    hasActivityCarrier(SOURCE_SPECIES, SOURCE_CARRIER),
                    hasActivityCarrier(TARGET_SPECIES, TARGET_CARRIER),
                    hasActivityKey(SOURCE_CARRIER, SOURCE_KEY),
                    hasActivityKey(TARGET_CARRIER, TARGET_KEY)."""),
            docs="A catalyzer of a reaction negatively influences each reactant that is itself an activity (consumption depletes the reactant — a negative influence regardless of the modifier's positive role on the product).",
        ),
        Rule(
            identifier="influences_consumption:celldesigner:physical_stimulator_consumes_reactant",
            text=dedent("""\
                influences(SOURCE_KEY, TARGET_KEY, negativelyInfluences) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), physicalStimulator(MODIFIER),
                    hasReferredElement(MODIFIER, SOURCE_SPECIES),
                    hasReactant(REACTION, REACTANT),
                    hasReferredElement(REACTANT, TARGET_SPECIES),
                    hasActivityCarrier(SOURCE_SPECIES, SOURCE_CARRIER),
                    hasActivityCarrier(TARGET_SPECIES, TARGET_CARRIER),
                    hasActivityKey(SOURCE_CARRIER, SOURCE_KEY),
                    hasActivityKey(TARGET_CARRIER, TARGET_KEY)."""),
            docs="A physical stimulator of a reaction negatively influences each reactant that is itself an activity (consumption).",
        ),
        Rule(
            identifier="influences_consumption:celldesigner:trigger_consumes_reactant",
            text=dedent("""\
                influences(SOURCE_KEY, TARGET_KEY, negativelyInfluences) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), trigger(MODIFIER),
                    hasReferredElement(MODIFIER, SOURCE_SPECIES),
                    hasReactant(REACTION, REACTANT),
                    hasReferredElement(REACTANT, TARGET_SPECIES),
                    hasActivityCarrier(SOURCE_SPECIES, SOURCE_CARRIER),
                    hasActivityCarrier(TARGET_SPECIES, TARGET_CARRIER),
                    hasActivityKey(SOURCE_CARRIER, SOURCE_KEY),
                    hasActivityKey(TARGET_CARRIER, TARGET_KEY)."""),
            docs="A trigger of a reaction negatively influences each reactant that is itself an activity (consumption is depletion, hence negative — not triggers, which is only the trigger→product relationship).",
        ),
        Rule(
            identifier="influences_consumption:celldesigner:inhibitor_spares_reactant",
            text=dedent("""\
                influences(SOURCE_KEY, TARGET_KEY, positivelyInfluences) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), inhibitor(MODIFIER),
                    hasReferredElement(MODIFIER, SOURCE_SPECIES),
                    hasReactant(REACTION, REACTANT),
                    hasReferredElement(REACTANT, TARGET_SPECIES),
                    hasActivityCarrier(SOURCE_SPECIES, SOURCE_CARRIER),
                    hasActivityCarrier(TARGET_SPECIES, TARGET_CARRIER),
                    hasActivityKey(SOURCE_CARRIER, SOURCE_KEY),
                    hasActivityKey(TARGET_CARRIER, TARGET_KEY)."""),
            docs="An inhibitor of a reaction positively influences each reactant that is itself an activity (sparing).",
        ),
        Rule(
            identifier="influences_consumption:celldesigner:unknown_catalyzer_consumes_reactant",
            text=dedent("""\
                influences(SOURCE_KEY, TARGET_KEY, unknownNegativelyInfluences) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), unknownCatalyzer(MODIFIER),
                    hasReferredElement(MODIFIER, SOURCE_SPECIES),
                    hasReactant(REACTION, REACTANT),
                    hasReferredElement(REACTANT, TARGET_SPECIES),
                    hasActivityCarrier(SOURCE_SPECIES, SOURCE_CARRIER),
                    hasActivityCarrier(TARGET_SPECIES, TARGET_CARRIER),
                    hasActivityKey(SOURCE_CARRIER, SOURCE_KEY),
                    hasActivityKey(TARGET_CARRIER, TARGET_KEY)."""),
            docs="An unknown catalyzer of a reaction unknown-negatively influences each reactant that is itself an activity (consumption, uncertain).",
        ),
        Rule(
            identifier="influences_consumption:celldesigner:unknown_inhibitor_spares_reactant",
            text=dedent("""\
                influences(SOURCE_KEY, TARGET_KEY, unknownPositivelyInfluences) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), unknownInhibitor(MODIFIER),
                    hasReferredElement(MODIFIER, SOURCE_SPECIES),
                    hasReactant(REACTION, REACTANT),
                    hasReferredElement(REACTANT, TARGET_SPECIES),
                    hasActivityCarrier(SOURCE_SPECIES, SOURCE_CARRIER),
                    hasActivityCarrier(TARGET_SPECIES, TARGET_CARRIER),
                    hasActivityKey(SOURCE_CARRIER, SOURCE_KEY),
                    hasActivityKey(TARGET_CARRIER, TARGET_KEY)."""),
            docs="An unknown inhibitor of a reaction unknown-positively influences each reactant that is itself an activity (sparing, uncertain).",
        ),
        ),
    },
)

_INFLUENCE_OUTPUT = RuleGroup(
    identifier="influence_output",
    docs="Shared fan-out from the internal `influences(SOURCE, TARGET, INFLUENCE_KIND)` relation to the typed `new(...)` influence heads — one rule per kind. Every pipeline converges on `influences/3`; this group is the single place that turns a kind into its output predicate.",
    rules=(
        Rule(
            identifier="influence_output:positive",
            text="new(positivelyInfluences(SOURCE, TARGET)) :- influences(SOURCE, TARGET, positivelyInfluences).",
            docs="A `positive` influence emits a `positivelyInfluences` edge (PositiveInfluence).",
        ),
        Rule(
            identifier="influence_output:negative",
            text="new(negativelyInfluences(SOURCE, TARGET)) :- influences(SOURCE, TARGET, negativelyInfluences).",
            docs="A `negative` influence emits a `negativelyInfluences` edge (NegativeInfluence).",
        ),
        Rule(
            identifier="influence_output:modulation",
            text="new(modulates(SOURCE, TARGET)) :- influences(SOURCE, TARGET, modulates).",
            docs="A `modulation` influence emits a `modulates` edge (Modulation).",
        ),
        Rule(
            identifier="influence_output:triggering",
            text="new(triggers(SOURCE, TARGET)) :- influences(SOURCE, TARGET, triggers).",
            docs="A `triggering` influence emits a `triggers` edge (Triggering).",
        ),
        Rule(
            identifier="influence_output:unknown_positive",
            text="new(unknownPositivelyInfluences(SOURCE, TARGET)) :- influences(SOURCE, TARGET, unknownPositivelyInfluences).",
            docs="An `unknown_positive` influence emits an `unknownPositivelyInfluences` edge (UnknownPositiveInfluence).",
        ),
        Rule(
            identifier="influence_output:unknown_negative",
            text="new(unknownNegativelyInfluences(SOURCE, TARGET)) :- influences(SOURCE, TARGET, unknownNegativelyInfluences).",
            docs="An `unknown_negative` influence emits an `unknownNegativelyInfluences` edge (UnknownNegativeInfluence).",
        ),
        Rule(
            identifier="influence_output:unknown_modulation",
            text="new(unknownModulates(SOURCE, TARGET)) :- influences(SOURCE, TARGET, unknownModulates).",
            docs="An `unknown_modulation` influence emits an `unknownModulates` edge (UnknownModulation).",
        ),
        Rule(
            identifier="influence_output:unknown_triggering",
            text="new(unknownTriggers(SOURCE, TARGET)) :- influences(SOURCE, TARGET, unknownTriggers).",
            docs="An `unknown_triggering` influence emits an `unknownTriggers` edge (UnknownTriggering).",
        ),
    ),
)


# Authored logical operators (CellDesigner `BooleanLogicGate`, SBGN-PD
# `LogicalOperator`). Today a gate yields `propagatesInfluence(OPERATOR, ...)`, but a gate id has no
# `hasActivityCarrier`, so `influences_derivation:path` never matches and the gate
# is silently dropped. `_GATES` carries the gate through three rule kinds:
#
#   (a) operator node -- one head per gate type, carrying a type token;
#   (b) input edges -- each gate input resolved through carrier/key (mirroring
#       `influences_derivation:path`);
#   (c) operator-sourced influence -- ONE rule that *reuses* `propagatesInfluence/3`: the
#       `paths:core` rules already emit `propagatesInfluence(OPERATOR, TARGET, INFLUENCE_KIND)` (the source
#       binds to whatever the modifier refers to / the modulation source -- gate
#       included), so `_GATES` only resolves the TARGET through carrier/key and
#       writes `influences(logicalOperatorKey(OPERATOR), TARGET_KEY, INFLUENCE_KIND)` directly,
#       bypassing the carrier-keyed `influences_derivation:path`. Because `propagatesInfluence/3`
#       is the transitive closure, an operator influences its direct target and
#       everything transitively downstream -- consistent with how a species
#       activity source already behaves (Decision D1).
#
# The `booleanLogicGate(OPERATOR)` / `logicalOperator(OPERATOR)` umbrella guard in (b) and (c)
# is essential -- it is derived from the per-type facts by the input ontology's
# isa rules, and without it (c) would treat every `propagatesInfluence/3` source (species
# included) as an operator key. Provenance-agnostic: a future derived-operator
# layer emits the same predicates and reuses this group's builder/layout/output
# path unchanged.
_GATES = RuleGroup(
    identifier="gates",
    depends_on=frozenset({"activity:core", "paths:core", "preparation"}),
    docs="Authored logical operators (CellDesigner `BooleanLogicGate`, SBGN-PD `LogicalOperator`): emits the operator node (`logicalOperator/2`, token-typed), its input edges (`logicalOperatorInput/2`, each input resolved through carrier/key), and the operator-sourced influence written straight into the internal `influences/3` relation by reusing the existing `propagatesInfluence(OPERATOR, TARGET, INFLUENCE_KIND)` closure. The widened influence `source` union (`predicates._INFLUENCE_SOURCE`) lets `influence_output` fan these out with no change. Provenance-agnostic.",
    rules=(),
    variants={
        "celldesigner": (
            Rule(
                identifier="gates:celldesigner:node_and",
                text="new(logicalOperator(logicalOperatorKey(OPERATOR), and)) :- andGate(OPERATOR).",
                docs="An `andGate` emits an AND logical-operator node.",
            ),
            Rule(
                identifier="gates:celldesigner:node_or",
                text="new(logicalOperator(logicalOperatorKey(OPERATOR), or)) :- orGate(OPERATOR).",
                docs="An `orGate` emits an OR logical-operator node.",
            ),
            Rule(
                identifier="gates:celldesigner:node_not",
                text="new(logicalOperator(logicalOperatorKey(OPERATOR), not_)) :- notGate(OPERATOR).",
                docs="A `notGate` emits a NOT logical-operator node. The token is `not_` because bare `not` is a reserved clingo keyword.",
            ),
            Rule(
                identifier="gates:celldesigner:node_unknown",
                text="new(logicalOperator(logicalOperatorKey(OPERATOR), unknown)) :- unknownGate(OPERATOR).",
                docs="An `unknownGate` emits an unknown-type logical-operator node.",
            ),
            Rule(
                identifier="gates:celldesigner:input_edge",
                text=dedent("""\
                    new(logicalOperatorInput(logicalOperatorKey(OPERATOR), INPUT_KEY)) :-
                        booleanLogicGate(OPERATOR),
                        hasInput(OPERATOR, INPUT),
                        hasReferredElement(INPUT, INPUT_SPECIES),
                        hasActivityCarrier(INPUT_SPECIES, INPUT_CARRIER),
                        hasActivityKey(INPUT_CARRIER, INPUT_KEY)."""),
                docs="Each gate input is resolved through its activity carrier and key (mirroring `influences_derivation:path`), so an input that is a subunit resolves to its top-level complex's key. The `booleanLogicGate` umbrella matches every gate type via the ontology's isa rules.",
            ),
            Rule(
                identifier="gates:celldesigner:influence",
                text=dedent("""\
                    influences(logicalOperatorKey(OPERATOR), TARGET_KEY, INFLUENCE_KIND) :-
                        booleanLogicGate(OPERATOR),
                        propagatesInfluence(OPERATOR, TARGET_SPECIES, INFLUENCE_KIND),
                        hasActivityCarrier(TARGET_SPECIES, TARGET_CARRIER),
                        hasActivityKey(TARGET_CARRIER, TARGET_KEY)."""),
                docs="One rule covering both Shape A (gate is a reaction modifier) and Shape B (gate is a modulation source): the `paths:core` rules already bind a `propagatesInfluence/3` whose source is the gate, so the gate only resolves its TARGET through carrier/key and writes the influence keyed by the operator. Inherits the transitive closure (Decision D1). The `booleanLogicGate` guard is essential -- without it every species `propagatesInfluence/3` source would be read as an operator key.",
            ),
        ),
        "sbgn_pd": (
            Rule(
                identifier="gates:sbgn_pd:node_and",
                text="new(logicalOperator(logicalOperatorKey(OPERATOR), and)) :- andOperator(OPERATOR).",
                docs="An `andOperator` emits an AND logical-operator node.",
            ),
            Rule(
                identifier="gates:sbgn_pd:node_or",
                text="new(logicalOperator(logicalOperatorKey(OPERATOR), or)) :- orOperator(OPERATOR).",
                docs="An `orOperator` emits an OR logical-operator node.",
            ),
            Rule(
                identifier="gates:sbgn_pd:node_not",
                text="new(logicalOperator(logicalOperatorKey(OPERATOR), not_)) :- notOperator(OPERATOR).",
                docs="A `notOperator` emits a NOT logical-operator node. The token is `not_` because bare `not` is a reserved clingo keyword. (SBGN-PD has no unknown-operator type.)",
            ),
            Rule(
                identifier="gates:sbgn_pd:input_edge",
                text=dedent("""\
                    new(logicalOperatorInput(logicalOperatorKey(OPERATOR), INPUT_KEY)) :-
                        logicalOperator(OPERATOR),
                        hasInput(OPERATOR, INPUT),
                        hasReferredElement(INPUT, INPUT_ENTITY_POOL),
                        hasActivityCarrier(INPUT_ENTITY_POOL, INPUT_CARRIER),
                        hasActivityKey(INPUT_CARRIER, INPUT_KEY)."""),
                docs="SBGN-PD parallel of the CellDesigner input-edge rule, guarded on the `logicalOperator` umbrella (derived from the per-type operators via the ontology isa rules).",
            ),
            Rule(
                identifier="gates:sbgn_pd:influence",
                text=dedent("""\
                    influences(logicalOperatorKey(OPERATOR), TARGET_KEY, INFLUENCE_KIND) :-
                        logicalOperator(OPERATOR),
                        propagatesInfluence(OPERATOR, TARGET_ENTITY_POOL, INFLUENCE_KIND),
                        hasActivityCarrier(TARGET_ENTITY_POOL, TARGET_CARRIER),
                        hasActivityKey(TARGET_CARRIER, TARGET_KEY)."""),
                docs="SBGN-PD parallel of the CellDesigner operator-sourced influence rule (Shape B: the operator is a modulation source). Reuses the `propagatesInfluence(OPERATOR, TARGET, INFLUENCE_KIND)` closure emitted by `paths:core`, resolving only the target through carrier/key.",
            ),
        ),
    },
)


# The `keep_reactions` profile: every species is an activity and every reaction
# is rendered as reactant->product positive influences, so the map's own topology
# *is* the influence network. It shares the whole scaffolding (activity
# bridge, top-level resolution, `keptSpeciesKey` keying, the single-hop
# modulation influences of `paths:core`, gates, `influence_output`) and adds only the two
# rules below, while omitting the inference layers (`paths:chaining`,
# `influences_consumption`) that the other modes carry. These rules are
# CellDesigner-only: the mode declares `celldesigner` as its sole compatible
# language.
_KEEP_REACTIONS_ACTIVITY = RuleGroup(
    identifier="keep_reactions:activity",
    depends_on=frozenset({"activity:core"}),
    docs="The `keep-reactions` premise that every species is an activity, replacing the structural-reason `activity:*` feature-groups (each of which is inert once every species is a candidate, so the profile omits them). Unlike the `--set-all-active` toggle this rule carries no `not hasSubunit` guard: a subunit is a candidate too, which does not make it an activity of its own -- `preparation:complex:key` keys it by the `keptSpeciesKey` of its outermost complex -- but it does mean a reaction, modulation arc or gate input touching a subunit routes to the containing complex instead of being dropped for want of an activity key.",
    rules=(
        Rule(
            identifier="keep_reactions:activity:from_species",
            text="hasActivityCandidate(SPECIES, isSpecies) :- species(SPECIES).",
            docs="Every species is an activity candidate, for the sole reason that it is a species (reason `isSpecies`). It still passes through the `activity:core` bridge, so `--set-inactive` and `--set-all-inactive` veto candidates in this mode exactly as in any other.",
        ),
    ),
)

_KEEP_REACTIONS_INFLUENCES = RuleGroup(
    identifier="keep_reactions:influences",
    depends_on=frozenset({"influences_derivation"}),
    docs="The `keep-reactions` premise that every reaction is kept: each (reactant, product) pair of a reaction becomes one positive influence. `positivelyInfluences` rather than `triggers` -- consuming a reactant to make a product is a contribution to it, not the necessary-stimulation relationship `triggers` (CellDesigner `Triggering`, SBGN-AF `NecessaryStimulation`) asserts. The rule writes `propagatesInfluence/3` rather than `influences/3` so the shared `influences_derivation:path` bridge resolves both endpoints through their carrier and key, routing a subunit endpoint to its top-level complex. Reactant and product may resolve to the same activity (a reaction whose participants share a complex, or a state transition drawn on a single species), and the resulting self-influence is emitted like any other.",
    rules=(
        Rule(
            identifier="keep_reactions:influences:reactant_to_product",
            text=dedent("""\
                propagatesInfluence(REACTANT_SPECIES, PRODUCT_SPECIES, positivelyInfluences) :-
                    reaction(REACTION),
                    hasReactant(REACTION, REACTANT),
                    hasReferredElement(REACTANT, REACTANT_SPECIES),
                    hasProduct(REACTION, PRODUCT),
                    hasReferredElement(PRODUCT, PRODUCT_SPECIES)."""),
            docs="Each reactant of a reaction positively influences each product of that reaction. This is the one rule that makes a reaction itself an influence: in the other modes the reactant->product hop only feeds the cycle relations, and a bare reactant is never an influence source.",
        ),
    ),
)


_BUILTIN_GROUPS = (
    _ACTIVITY_CORE,
    _ACTIVITY_PHENOTYPE,
    _ACTIVITY_ACTIVE_MARKER,
    _ACTIVITY_MODULATION_SOURCE,
    _ACTIVITY_GATE_INPUT,
    _TOPOLOGY,
    _TOP_LEVEL,
    _PREPARATION_COMPLEX,
    _PREPARATION_NO_COMPLEX,
    _MODULATION_KIND,
    _PATHS_CORE,
    _PATHS_CHAINING,
    _PATHS_COMPLEX_TRAVERSAL,
    _INFLUENCES_DERIVATION,
    _INFLUENCES_CONSUMPTION,
    _INFLUENCE_OUTPUT,
    _GATES,
    _KEEP_REACTIONS_ACTIVITY,
    _KEEP_REACTIONS_INFLUENCES,
)


def build_registry(modes=None) -> RuleRegistry:
    """Register the built-in rule groups plus every group the modes define.

    ``modes`` is an iterable of :class:`pd2af.modes.TransformationMode`,
    defaulting to every known mode. Passing a narrower set — the built-in modes
    alone, say — yields a registry that ignores contributed groups.
    """
    if modes is None:
        modes = pd2af.modes.get_transformation_modes().values()
    contributed_groups = []
    mode_by_contributed_group_id: dict[str, str] = {}
    for mode in modes:
        for group in mode.group_definitions:
            already_contributed_by = mode_by_contributed_group_id.get(
                group.identifier
            )
            if already_contributed_by is not None:
                raise RuntimeError(
                    f"transformation modes {already_contributed_by!r} and "
                    f"{mode.name!r} both define rule group "
                    f"{group.identifier!r}"
                )
            mode_by_contributed_group_id[group.identifier] = mode.name
            contributed_groups.append(group)
    registry = RuleRegistry()
    registry.register(list(_BUILTIN_GROUPS) + contributed_groups)
    issues = registry.validate()
    if issues:
        raise RuntimeError(
            "pd2af rule registry is invalid:\n"
            + "\n".join(f"  - {issue}" for issue in issues)
        )
    return registry


def get_excludable_groups(mode_name: str) -> tuple[frozenset[str], frozenset[str]]:
    """Return ``(excludable, mandatory)`` group ids for a transformation mode.

    A group is *mandatory* when another included group depends on it (the
    dependency graph forbids excluding it — ``resolve`` would raise
    ``excluded_dependency``); every other included group is an *excludable*
    leaf that ``--exclude-group`` can drop cleanly. A dependency on a *slot*
    resolves to the included group filling it, so a filler is mandatory
    exactly like a directly named dependency.
    """
    registry = build_registry()
    included = frozenset(pd2af.modes.get_transformation_mode(mode_name).group_ids)
    depended_on: set[str] = set()
    for group_id in included:
        for dependency in registry.groups[group_id].depends_on:
            if dependency in included:
                depended_on.add(dependency)
            else:
                depended_on.update(
                    registry.slots.get(dependency, frozenset()) & included
                )
    mandatory = frozenset(depended_on)
    excludable = included - mandatory
    return excludable, mandatory


def _friendly_exclusion_message(error: PlanInvalidError) -> str:
    """Turn an ``excluded_dependency`` failure into an actionable message."""
    required_by: dict[str, list[str]] = {}
    other_issues: list[str] = []
    for issue in error.issues:
        if issue.kind == "excluded_dependency":
            group_id, dependency = issue.groups
            required_by.setdefault(dependency, []).append(group_id)
        else:
            other_issues.append(issue.message)
    lines = []
    for dependency, dependents in sorted(required_by.items()):
        lines.append(
            f"cannot exclude '{dependency}' — required by "
            + ", ".join(sorted(dependents))
        )
    lines.extend(other_issues)
    return "\n".join(lines) or str(error)


def build_program(
    mode_name: str,
    language: str = "celldesigner",
    exclude_groups: tuple[str, ...] = (),
    exclude_rules: tuple[str, ...] = (),
) -> str:
    """Return the composed ASP program text for the named transformation mode
    and input ``language``.

    The mode names the groups its program is made of; the input language is an
    aspcompose *variant*. Language-agnostic rules live in each group's
    ``rules`` and are emitted for every language; language-specific rules live
    in ``variants={"celldesigner": ..., "sbgn_pd": ...}`` and are selected here
    by ``resolve(variant=language)``.

    ``exclude_groups`` drops whole rule groups (the primary toggle: each
    excludable group is a coherent functional unit); excluding a group that
    another included group depends on raises ``ValueError`` with a friendly
    message. ``exclude_rules`` is a finer scalpel that drops individual rules
    by identifier (used for the whole-with-scalpel table groups); an id that
    does not name a resolved rule raises ``ValueError``.
    """
    mode = pd2af.modes.get_transformation_mode(mode_name)
    plan = CollectionPlan(build_registry())
    for group_id in mode.group_ids:
        plan.add_group(group_id)
    for group_id in exclude_groups:
        plan.exclude_group(group_id)
    try:
        resolved = plan.resolve(variant=language)
    except PlanInvalidError as error:
        raise ValueError(_friendly_exclusion_message(error)) from error
    if exclude_rules:
        resolved_ids = {rule.identifier for rule in resolved}
        unknown = sorted(set(exclude_rules) - resolved_ids)
        if unknown:
            raise ValueError(
                "unknown rule id(s) for --disable-rule: " + ", ".join(unknown)
            )
        excluded = set(exclude_rules)
        resolved = [rule for rule in resolved if rule.identifier not in excluded]
    return "\n".join(rule.text for rule in resolved)
