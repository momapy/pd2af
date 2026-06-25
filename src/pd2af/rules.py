"""ASP rule composition for the pd2af transformation modes.

The non-CASQ rules are organised in layers:

* **topology** (shared across all non-CASQ profiles) — structural
  helpers: ``isSubunit``, ``hasActiveDescendantSubunit``.
* **top_level** (the complex-keeping modes ``keep-species`` and
  ``normal``) — ``resolvesToTopLevel(SPECIES, TOPLEVEL)`` resolves every species to its
  outermost top-level entity, so a subunit is keyed by (and its
  influences routed to) its top-level complex rather than itself.
* **preparation** (one rule group per non-CASQ mode) — emits
  ``hasActivityCarrier(RAW_SPECIES, ACTIVITY_BEARER)`` and
  ``hasActivityKey(ACTIVITY_BEARER, KEY)`` where ``KEY`` is one of two
  per-species wrappers: ``keptSpeciesKey/1`` (top-level species, and the
  top-level complex a subunit resolves to) or ``promotedSubunitKey/1`` (a
  subunit promoted to top level when its complex is dissolved in the
  ``*-no-complex`` modes). Proteoform/PTM stripping for the merged modes
  (``normal``/``normal-no-complex``) happens at the build stage, not in the key.
* **derivation** (shared across all non-CASQ profiles) — emits
  ``new(activity(KEY))`` and ``new(positivelyInfluences(...))`` /
  ``new(negativelyInfluences(...))`` from ``hasActivityCarrier`` /
  ``hasActivityKey``.

The ``casq`` profile keeps its own pipeline (deletion rules,
bridged-product rewiring, direct reactant→product / modifier→product
influence emission). It uses ``keptSpeciesKey/1`` for every surviving
species's activity key, sharing only the predicate type with the
non-CASQ layers.
"""

from textwrap import dedent

from aspcompose import CollectionPlan, Rule, RuleGroup, RuleRegistry

_NON_CASQ_PROFILES = frozenset(
    {"normal", "normal_no_complex", "keep_species", "keep_species_no_complex"}
)

_CASQ_PROFILES = frozenset({"casq"})

# Activity discovery is parallel across languages: only the `phenotype` rule is
# language-agnostic (both languages emit the `phenotype` functor); every other
# signal is a per-language variant. CellDesigner and SBGN-PD share the same
# signals -- an "active" marker, a modulation source, and a phenotype -- except
# that SBGN-PD has no explicit `hasActive` flag, CellDesigner's active
# *structural state* becomes SBGN-PD's active *state variable*, and a reaction
# modifier (CellDesigner-only) is, in SBGN-PD, just a modulation arc whose
# target is a process (so it folds into the modulation-source rule). A bare
# reactant/product is never an activity in either language.
_ACTIVITY_BASE = RuleGroup(
    identifier="activity_base",
    profiles=_NON_CASQ_PROFILES | _CASQ_PROFILES,
    docs="Base rules deriving `hasActivity(SPECIES, REASON)`. Only the phenotype rule is shared; the active-marker and modulation-source rules are per-language variants (CellDesigner: active flag, active structural state, modulation arc, reaction modifier; SBGN-PD: active state variable, modulation arc).",
    rules=(
        Rule(
            identifier="activity_base:from_phenotype",
            text="hasActivity(PHENOTYPE, isPhenotype) :- phenotype(PHENOTYPE).",
            docs="If a species/process is a phenotype, then it has activity, with reason `isPhenotype`. Shared: both languages emit the `phenotype` functor.",
        ),
    ),
    variants={
        "celldesigner": (
            Rule(
                identifier="activity_base:from_active_flag",
                text="hasActivity(SPECIES, isActive) :- species(SPECIES), hasActive(SPECIES, 1).",
                docs="If a species has its `hasActive` flag set to 1, then it has activity, with reason `isActive`.",
            ),
            Rule(
                identifier="activity_base:from_active_structural_state",
                text=dedent("""\
                    hasActivity(SPECIES, hasActiveStructuralState) :-
                        species(SPECIES),
                        hasStructuralState(SPECIES, STRUCTURAL_STATE),
                        hasValue(STRUCTURAL_STATE, "active")."""),
                docs='If a species carries a structural state whose value is "active", then it has activity, with reason `hasActiveStructuralState`.',
            ),
            Rule(
                identifier="activity_base:from_modulation_source",
                text=dedent("""\
                    hasActivity(SOURCE, isModulationSource) :-
                        species(SOURCE),
                        knownOrUnknownModulation(MODULATION),
                        hasSource(MODULATION, SOURCE),
                        hasTarget(MODULATION, _)."""),
                docs="If a species is the source of a modulation arc (known *or* unknown) with some target, then it has activity, with reason `isModulationSource`. Keying on `knownOrUnknownModulation` rather than `modulation` is what lets the source of an unknown modulation become an activity node.",
            ),
            Rule(
                identifier="activity_base:from_reaction_modulator",
                text=dedent("""\
                    hasActivity(SOURCE, isReactionModifier) :-
                        species(SOURCE),
                        knownOrUnknownModulator(MODULATOR),
                        hasReferredElement(MODULATOR, SOURCE),
                        hasModifier(_, MODULATOR)."""),
                docs="If a species is referred to by a reaction modulator (known *or* unknown) that modifies some target reaction, then it has activity, with reason `isReactionModifier`. Keying on `knownOrUnknownModulator` rather than `modulator` is what lets an unknown catalyzer/inhibitor become an activity node.",
            ),
            Rule(
                identifier="activity_base:from_gate_input",
                text=dedent("""\
                    hasActivity(ELEMENT, isGateInput) :-
                        booleanLogicGateInput(INPUT),
                        hasReferredElement(INPUT, ELEMENT)."""),
                docs="If a species feeds a boolean logic gate input, then it has activity, with reason `isGateInput`. A gate is structurally always an influence/modulation source or reaction modifier, so each of its inputs is an active driver of the downstream target -- a semantic guarantee, not a fallback. The element is activated regardless of kind (species, complex, ion, ...); the carrier rules route each kind. `hasReferredElement` in CellDesigner only relates a gate input to its species, so no extra guard is needed.",
            ),
        ),
        "sbgn_pd": (
            Rule(
                identifier="activity_base:sbgn_pd:from_active_state_variable",
                text=dedent("""\
                    hasActivity(ENTITY_POOL, hasActiveStateVariable) :-
                        entityPool(ENTITY_POOL),
                        hasStateVariable(ENTITY_POOL, STATE_VARIABLE),
                        hasValue(STATE_VARIABLE, "active")."""),
                docs='If an entity pool carries a state variable whose value is "active", then it has activity (the SBGN-PD parallel of CellDesigner\'s active structural state), with reason `hasActiveStateVariable`.',
            ),
            Rule(
                identifier="activity_base:sbgn_pd:from_active_subunit_state_variable",
                text=dedent("""\
                    hasActivity(SUBUNIT, hasActiveStateVariable) :-
                        isSubunit(SUBUNIT),
                        hasStateVariable(SUBUNIT, STATE_VARIABLE),
                        hasValue(STATE_VARIABLE, "active")."""),
                docs='If a subunit carries a state variable whose value is "active", then it has activity. Its complex therefore inherits activity (keep-species), and in the ``*-no-complex`` modes the subunit can be promoted. The parallel of CellDesigner, where subunits are species and so the active-structural-state rule already covers them.',
            ),
            Rule(
                identifier="activity_base:sbgn_pd:from_modulation_source",
                text=dedent("""\
                    hasActivity(SOURCE, isModulationSource) :-
                        entityPool(SOURCE),
                        modulation(MODULATION),
                        hasSource(MODULATION, SOURCE),
                        hasTarget(MODULATION, _)."""),
                docs="If an entity pool is the source of a modulation arc (whose target is a process), then it has activity, with reason `isModulationSource`. This is the SBGN-PD parallel of both the CellDesigner modulation-arc and reaction-modifier rules.",
            ),
            Rule(
                identifier="activity_base:sbgn_pd:from_operator_input",
                text=dedent("""\
                    hasActivity(ELEMENT, isGateInput) :-
                        logicalOperatorInput(INPUT),
                        hasReferredElement(INPUT, ELEMENT),
                        entityPool(ELEMENT)."""),
                docs="If an entity pool feeds a logical operator input, then it has activity, with reason `isGateInput` (the SBGN-PD parallel of CellDesigner's gate-input activation). The `entityPool` guard excludes the deferred nested-operator case: an operator feeding another operator has no activity carrier, so such an input simply dangles. In SBGN-PD `hasReferredElement` relates many roles (reactant, product, modulation participant) to entities, so both the `logicalOperatorInput` and `entityPool` guards are needed.",
            ),
        ),
    },
)

_TOPOLOGY = RuleGroup(
    identifier="topology",
    profiles=_NON_CASQ_PROFILES,
    docs="Mode-agnostic structural helpers shared by all non-CASQ profiles: `isSubunit`, `hasActiveDescendantSubunit`.",
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
    profiles=frozenset({"keep_species", "normal"}),
    depends_on=frozenset({"topology"}),
    docs="Resolves every species to its outermost top-level entity: a non-subunit resolves to itself; a subunit -- at any nesting depth -- resolves to the outermost complex that contains it. The complex-keeping modes (`keep-species`, `normal`) key a species by its top-level entity, so a subunit is never its own activity and its influences attach to its top-level complex (a subunit is a structural component, not an independent influencer). Mirrors the resolution `casq` performs in its own pipeline.",
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
        ),
    },
)

_PREPARATION_KEEP_SPECIES = RuleGroup(
    identifier="preparation:keep_species",
    profiles=frozenset({"keep_species"}),
    depends_on=frozenset({"activity_base", "topology", "top_level"}),
    docs="`keep-species` preparation: a species with activity is keyed by the `keptSpeciesKey` of its top-level entity (`top_level` group) -- itself when top-level, its outermost complex when a subunit. A subunit therefore contributes no activity of its own: it is a structural component of its complex, and any influence it carries attaches to the top-level complex (the active descendant directly keys the complex, so the assembly is represented without a separate inherit-activity rule). Carriers are identity; the rerouting lives entirely in the key.",
    rules=(
        Rule(
            identifier="preparation:keep_species:key",
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
                identifier="preparation:keep_species:carrier",
                text="hasActivityCarrier(SPECIES, SPECIES) :- species(SPECIES).",
                docs="Each species is its own activity carrier (no rerouting).",
            ),
        ),
        "sbgn_pd": (
            Rule(
                identifier="preparation:keep_species:sbgn_pd:carrier_entity_pool",
                text="hasActivityCarrier(ENTITY_POOL, ENTITY_POOL) :- entityPool(ENTITY_POOL).",
                docs="SBGN-PD: each entity pool is its own activity carrier.",
            ),
            Rule(
                identifier="preparation:keep_species:sbgn_pd:carrier_phenotype",
                text="hasActivityCarrier(PHENOTYPE, PHENOTYPE) :- phenotype(PHENOTYPE).",
                docs="SBGN-PD: each phenotype process is its own activity carrier (phenotypes are activities, not entity pools).",
            ),
        ),
    },
)

_PREPARATION_KEEP_SPECIES_NO_COMPLEX = RuleGroup(
    identifier="preparation:keep_species_no_complex",
    profiles=frozenset({"keep_species_no_complex"}),
    depends_on=frozenset({"activity_base", "topology"}),
    docs="`keep-species-no-complex` preparation: a complex with any active descendant is deleted; non-deleted top-level species with activity contribute their own `keptSpeciesKey(SELF)` activity; subunits of deleted complexes are promoted to top-level activities keyed by `promotedSubunitKey(SELF)`.",
    rules=(
        Rule(
            identifier="preparation:keep_species_no_complex:deleted",
            text=dedent("""\
                delete(COMPLEX) :-
                    complex(COMPLEX),
                    hasActiveDescendantSubunit(COMPLEX)."""),
            docs="A complex with any (transitive) active descendant is deleted.",
        ),
        Rule(
            identifier="preparation:keep_species_no_complex:subunit_of_deleted_direct",
            text=dedent("""\
                isDescendantSubunitOfDeleted(SUBUNIT) :-
                    hasSubunit(COMPLEX, SUBUNIT),
                    delete(COMPLEX)."""),
            docs="A direct subunit of a deleted complex is itself subunit-of-deleted.",
        ),
        Rule(
            identifier="preparation:keep_species_no_complex:subunit_of_deleted_transitive",
            text=dedent("""\
                isDescendantSubunitOfDeleted(SUBUNIT) :-
                    hasSubunit(PARENT_COMPLEX, SUBUNIT),
                    isDescendantSubunitOfDeleted(PARENT_COMPLEX)."""),
            docs="Subunit-of-deleted is transitive through nested complexes.",
        ),
        Rule(
            identifier="preparation:keep_species_no_complex:key_top_level",
            text=dedent("""\
                hasActivityKey(SPECIES, keptSpeciesKey(SPECIES)) :-
                    hasActivity(SPECIES, _),
                    not isSubunit(SPECIES),
                    not delete(SPECIES)."""),
            docs="A non-deleted top-level species with activity is keyed by `keptSpeciesKey(SELF)`.",
        ),
        Rule(
            identifier="preparation:keep_species_no_complex:key_promoted_subunit",
            text=dedent("""\
                hasActivityKey(SPECIES, promotedSubunitKey(SPECIES)) :-
                    hasActivity(SPECIES, _),
                    isDescendantSubunitOfDeleted(SPECIES)."""),
            docs="A subunit of a deleted complex with activity is keyed by `promotedSubunitKey(SELF)` and added at top level.",
        ),
    ),
    # The activity carrier is a per-language variant. In SBGN-PD a promoted
    # subunit must also be its own carrier so a rerouted path (a complex's
    # influence propagated to its subunit via `paths_complex_traversal`)
    # becomes an influence on the promoted subunit.
    variants={
        "celldesigner": (
            Rule(
                identifier="preparation:keep_species_no_complex:carrier",
                text="hasActivityCarrier(SPECIES, SPECIES) :- species(SPECIES).",
                docs="Each species is its own activity carrier.",
            ),
        ),
        "sbgn_pd": (
            Rule(
                identifier="preparation:keep_species_no_complex:sbgn_pd:carrier_entity_pool",
                text="hasActivityCarrier(ENTITY_POOL, ENTITY_POOL) :- entityPool(ENTITY_POOL).",
                docs="SBGN-PD: each entity pool is its own activity carrier.",
            ),
            Rule(
                identifier="preparation:keep_species_no_complex:sbgn_pd:carrier_phenotype",
                text="hasActivityCarrier(PHENOTYPE, PHENOTYPE) :- phenotype(PHENOTYPE).",
                docs="SBGN-PD: each phenotype process is its own activity carrier.",
            ),
            Rule(
                identifier="preparation:keep_species_no_complex:sbgn_pd:carrier_subunit",
                text="hasActivityCarrier(SUBUNIT, SUBUNIT) :- isSubunit(SUBUNIT).",
                docs="SBGN-PD: each subunit is its own activity carrier, so a promoted subunit can be an influence endpoint.",
            ),
        ),
    },
)

_PREPARATION_NO_COMPLEX = RuleGroup(
    identifier="preparation:normal_no_complex",
    profiles=frozenset({"normal_no_complex"}),
    depends_on=frozenset({"activity_base", "topology"}),
    docs="`normal-no-complex` preparation: a complex with any active descendant is deleted; non-deleted top-level species with activity are keyed by `keptSpeciesKey(SELF)`; subunits of deleted complexes are promoted to top level, keyed by `promotedSubunitKey(SELF)` (paths reach them via `paths_complex_traversal`). PTM stripping for the merged `normal-no-complex` mode happens at the build stage, not in the activity key.",
    rules=(
        Rule(
            identifier="preparation:normal_no_complex:deleted",
            text=dedent("""\
                delete(COMPLEX) :-
                    complex(COMPLEX),
                    hasActiveDescendantSubunit(COMPLEX)."""),
            docs="A complex with any (transitive) active descendant is deleted.",
        ),
        Rule(
            identifier="preparation:normal_no_complex:subunit_of_deleted_direct",
            text=dedent("""\
                isDescendantSubunitOfDeleted(SUBUNIT) :-
                    hasSubunit(COMPLEX, SUBUNIT),
                    delete(COMPLEX)."""),
            docs="A direct subunit of a deleted complex is itself subunit-of-deleted.",
        ),
        Rule(
            identifier="preparation:normal_no_complex:subunit_of_deleted_transitive",
            text=dedent("""\
                isDescendantSubunitOfDeleted(SUBUNIT) :-
                    hasSubunit(PARENT_COMPLEX, SUBUNIT),
                    isDescendantSubunitOfDeleted(PARENT_COMPLEX)."""),
            docs="Subunit-of-deleted is transitive through nested complexes.",
        ),
    ),
    variants={
        "celldesigner": (
            Rule(
                identifier="preparation:normal_no_complex:carrier",
                text="hasActivityCarrier(SPECIES, SPECIES) :- species(SPECIES).",
                docs="Each species is its own activity carrier.",
            ),
            Rule(
                identifier="preparation:normal_no_complex:key_top_level",
                text=dedent("""\
                    hasActivityKey(SPECIES, keptSpeciesKey(SPECIES)) :-
                        hasActivity(SPECIES, _),
                        not isSubunit(SPECIES),
                        not delete(SPECIES)."""),
                docs="A non-deleted top-level active species is keyed by `keptSpeciesKey(SELF)`.",
            ),
            Rule(
                identifier="preparation:normal_no_complex:key_promoted_subunit",
                text=dedent("""\
                    hasActivityKey(SPECIES, promotedSubunitKey(SPECIES)) :-
                        hasActivity(SPECIES, _),
                        isDescendantSubunitOfDeleted(SPECIES)."""),
                docs="A subunit of a deleted complex with activity is promoted to top level, keyed by `promotedSubunitKey(SELF)`.",
            ),
        ),
        "sbgn_pd": (
            Rule(
                identifier="preparation:normal_no_complex:sbgn_pd:carrier_entity_pool",
                text="hasActivityCarrier(ENTITY_POOL, ENTITY_POOL) :- entityPool(ENTITY_POOL).",
                docs="SBGN-PD: each entity pool is its own activity carrier.",
            ),
            Rule(
                identifier="preparation:normal_no_complex:sbgn_pd:carrier_phenotype",
                text="hasActivityCarrier(PHENOTYPE, PHENOTYPE) :- phenotype(PHENOTYPE).",
                docs="SBGN-PD: each phenotype process is its own activity carrier.",
            ),
            Rule(
                identifier="preparation:normal_no_complex:sbgn_pd:carrier_subunit",
                text="hasActivityCarrier(SUBUNIT, SUBUNIT) :- isSubunit(SUBUNIT).",
                docs="SBGN-PD: each subunit is its own activity carrier, so a promoted subunit can be an influence endpoint (paths reach it via `paths_complex_traversal`).",
            ),
            Rule(
                identifier="preparation:normal_no_complex:sbgn_pd:key_top_level",
                text=dedent("""\
                    hasActivityKey(SPECIES, keptSpeciesKey(SPECIES)) :-
                        hasActivity(SPECIES, _),
                        not isSubunit(SPECIES),
                        not delete(SPECIES)."""),
                docs="SBGN-PD: a non-deleted top-level active entity is keyed by `keptSpeciesKey(SELF)`.",
            ),
            Rule(
                identifier="preparation:normal_no_complex:sbgn_pd:key_promoted_subunit",
                text=dedent("""\
                    hasActivityKey(SPECIES, promotedSubunitKey(SPECIES)) :-
                        hasActivity(SPECIES, _),
                        isDescendantSubunitOfDeleted(SPECIES)."""),
                docs="SBGN-PD: a subunit of a deleted complex with activity is promoted to top level, keyed by `promotedSubunitKey(SELF)`.",
            ),
        ),
    },
)

_PREPARATION_NORMAL = RuleGroup(
    identifier="preparation:normal",
    profiles=frozenset({"normal"}),
    depends_on=frozenset({"activity_base", "topology", "top_level"}),
    docs="`normal` preparation: same activity keying as `keep-species` -- a species with activity is keyed by the `keptSpeciesKey` of its top-level entity (`top_level` group), so a subunit is never its own activity and its influences attach to its top-level complex. `normal` differs from `keep-species` only at the build stage (proteoform/PTM stripping in the merged modes), not in the activity keys. Carriers are identity.",
    rules=(
        Rule(
            identifier="preparation:normal:key",
            text=dedent("""\
                hasActivityKey(SPECIES, keptSpeciesKey(TOPLEVEL)) :-
                    hasActivity(SPECIES, _),
                    resolvesToTopLevel(SPECIES, TOPLEVEL)."""),
            docs="A species with activity is keyed by the `keptSpeciesKey` of its top-level entity (itself when top-level; its outermost complex when a subunit).",
        ),
    ),
    variants={
        "celldesigner": (
            Rule(
                identifier="preparation:normal:carrier",
                text="hasActivityCarrier(SPECIES, SPECIES) :- species(SPECIES).",
                docs="Each species is its own activity carrier (no rerouting; a subunit endpoint is rerouted to its top-level complex by the `top_level` key, not by the carrier).",
            ),
        ),
        "sbgn_pd": (
            Rule(
                identifier="preparation:normal:sbgn_pd:carrier_entity_pool",
                text="hasActivityCarrier(ENTITY_POOL, ENTITY_POOL) :- entityPool(ENTITY_POOL).",
                docs="SBGN-PD: each entity pool is its own activity carrier.",
            ),
            Rule(
                identifier="preparation:normal:sbgn_pd:carrier_phenotype",
                text="hasActivityCarrier(PHENOTYPE, PHENOTYPE) :- phenotype(PHENOTYPE).",
                docs="SBGN-PD: each phenotype process is its own activity carrier.",
            ),
        ),
    },
)

# SBGN-PD path rules. Unlike CellDesigner (where a modulation arc connects two
# species), an SBGN-PD modulation arc's source is an entity pool and its target
# is a *process*; the influence therefore propagates to the process's products
# (or, when the target is a phenotype, to the phenotype itself). These live in
# the `sbgn_pd` variant of `paths_base` because the CellDesigner modulation-arc
# rules would misfire on PD facts (their `catalysis`/`inhibition`/`modulation`
# functors match PD modulations, but PD targets are processes, yielding spurious
# edges).
_PATHS_BASE_SBGN_PD = (
    Rule(
        identifier="paths_base:sbgn_pd:modulation_kind_necessary_stimulation",
        text="hasModulationKind(MODULATION, triggers) :- necessaryStimulation(MODULATION).",
        docs="An SBGN-PD necessary stimulation contributes a `triggering` kind.",
    ),
    Rule(
        identifier="paths_base:sbgn_pd:modulation_kind_stimulation",
        text=dedent("""\
            hasModulationKind(MODULATION, positivelyInfluences) :-
                stimulation(MODULATION),
                not necessaryStimulation(MODULATION)."""),
        docs="A stimulation (catalysis included) that is not a necessary stimulation contributes a `positive` kind.",
    ),
    Rule(
        identifier="paths_base:sbgn_pd:modulation_kind_inhibition",
        text="hasModulationKind(MODULATION, negativelyInfluences) :- inhibition(MODULATION).",
        docs="An SBGN-PD inhibition contributes a `negative` kind.",
    ),
    Rule(
        identifier="paths_base:sbgn_pd:modulation_kind_modulation",
        text=dedent("""\
            hasModulationKind(MODULATION, modulates) :-
                modulation(MODULATION),
                not stimulation(MODULATION),
                not inhibition(MODULATION)."""),
        docs="A bare modulation (neither stimulation nor inhibition) contributes an unknown-sign `modulation` kind.",
    ),
    Rule(
        identifier="paths_base:sbgn_pd:modulation_to_product",
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
        identifier="paths_base:sbgn_pd:modulation_to_phenotype",
        text=dedent("""\
            propagatesInfluence(SOURCE_ENTITY_POOL, TARGET_ENTITY_POOL, INFLUENCE_KIND) :-
                hasModulationKind(MODULATION, INFLUENCE_KIND),
                hasSource(MODULATION, SOURCE_ENTITY_POOL),
                hasTarget(MODULATION, TARGET_ENTITY_POOL),
                phenotype(TARGET_ENTITY_POOL)."""),
        docs="A modulation arc whose target is a phenotype influences the phenotype itself (a phenotype process has no products; it is the activity).",
    ),
    Rule(
        identifier="paths_base:sbgn_pd:is_directly_transformed_to",
        text=dedent("""\
            isDirectlyTransformedTo(UPSTREAM_ENTITY_POOL, DOWNSTREAM_ENTITY_POOL) :-
                hasReactant(PROCESS, REACTANT),
                hasReferredElement(REACTANT, UPSTREAM_ENTITY_POOL),
                hasProduct(PROCESS, PRODUCT),
                hasReferredElement(PRODUCT, DOWNSTREAM_ENTITY_POOL)."""),
        docs="The single reactant->product hop in the production graph: the upstream entity pool is an element of a reactant and the downstream pool an element of a product of the same process. Passive voice (the process does the transforming, not the pool) keeps it language-neutral. Feeds the shared `isTransformedTo`/`isCyclicallyTransformedTo` cycle relations that gate transitive path extension; carries no influence kind itself.",
    ),
    Rule(
        identifier="paths_base:sbgn_pd:transitive_through_process",
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
    Rule(
        identifier="paths_base:composes_to",
        text=dedent("""\
            composesTo(positivelyInfluences, positivelyInfluences).
            composesTo(negativelyInfluences, negativelyInfluences).
            composesTo(modulates, modulates).
            composesTo(triggers, positivelyInfluences).
            composesTo(unknownPositivelyInfluences, unknownPositivelyInfluences).
            composesTo(unknownNegativelyInfluences, unknownNegativelyInfluences).
            composesTo(unknownModulates, unknownModulates).
            composesTo(unknownTriggers, unknownPositivelyInfluences)."""),
        docs="How an influence kind transforms across a process reactant->product hop (shared logic with CellDesigner; identity except triggering->positive and unknown_triggering->unknown_positive).",
    ),
)


_PATHS_BASE = RuleGroup(
    identifier="paths_base",
    profiles=_NON_CASQ_PROFILES,
    docs="Builds the kinded `propagatesInfluence(SOURCE_SPECIES, TARGET_SPECIES, INFLUENCE_KIND)` relation from PD reactions and modulation arcs. INFLUENCE_KIND is one of `positive`, `negative`, `triggering`, `modulation` and their `unknown_*` twins. Reaction modifiers and the matching species→species modulation arcs map to the *same* kind (e.g. a trigger modifier and a triggering arc both give `triggering`; catalysis and physical stimulation both give `positive`; inhibition gives `negative`). Reactant-chained transitivity extends paths through reactions, degrading `triggering`→`positive` (and `unknown_triggering`→`unknown_positive`) at each reaction hop via `composesTo`. Transitive extension is gated by the production-cycle relations: each variant defines the single reactant->product hop `isDirectlyTransformedTo/2`, the shared `isTransformedTo/2` is its transitive closure, and `isCyclicallyTransformedTo/2` marks the hops that lie inside a cycle so transitivity never propagates influence back around a production loop.",
    rules=(
        Rule(
            identifier="paths_base:is_transformed_to",
            text=dedent("""\
                isTransformedTo(UPSTREAM, DOWNSTREAM) :- isDirectlyTransformedTo(UPSTREAM, DOWNSTREAM).
                isTransformedTo(UPSTREAM, DOWNSTREAM) :-
                    isTransformedTo(UPSTREAM, INTERMEDIATE),
                    isDirectlyTransformedTo(INTERMEDIATE, DOWNSTREAM)."""),
            docs="Transitive closure of the single-hop `isDirectlyTransformedTo/2`: `isTransformedTo(UPSTREAM, DOWNSTREAM)` holds when DOWNSTREAM is reachable from UPSTREAM through one or more reactant->product hops (the direct hop included, as the base case). Language-agnostic; `isDirectlyTransformedTo/2` is supplied per language variant. Depends only on `isDirectlyTransformedTo`, never on `propagatesInfluence`, so the negation that gates transitivity stays stratified.",
        ),
        Rule(
            identifier="paths_base:is_cyclically_transformed_to",
            text=dedent("""\
                isCyclicallyTransformedTo(UPSTREAM, DOWNSTREAM) :-
                    isTransformedTo(UPSTREAM, DOWNSTREAM),
                    isTransformedTo(DOWNSTREAM, UPSTREAM)."""),
            docs="UPSTREAM and DOWNSTREAM are mutually reachable through the production graph -- both lie on a common cycle. At every transitivity guard site the hop in question is already a direct reactant->product edge, so requiring mutual reachability there is exactly equivalent to 'this hop is on a cycle'. The transitivity rules forbid extending a path across such a hop, blocking influence from leaking around production loops.",
        ),
    ),
    variants={
        "sbgn_pd": _PATHS_BASE_SBGN_PD,
        "celldesigner": (
        Rule(
            identifier="paths_base:catalyzer_to_product",
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
            identifier="paths_base:physical_stimulator_to_product",
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
            identifier="paths_base:trigger_to_product",
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
            identifier="paths_base:catalysis_modulation",
            text=dedent("""\
                propagatesInfluence(SOURCE_SPECIES, TARGET_SPECIES, positivelyInfluences) :-
                    catalysis(MODULATION),
                    hasSource(MODULATION, SOURCE_SPECIES),
                    hasTarget(MODULATION, TARGET_SPECIES)."""),
            docs="If a modulation arc is a catalysis, then there is a positive path from its source to its target. (The predicate is `catalysis`; an earlier `catalyzis` typo silently disabled this rule.)",
        ),
        Rule(
            identifier="paths_base:positive_influence_modulation",
            text=dedent("""\
                propagatesInfluence(SOURCE_SPECIES, TARGET_SPECIES, positivelyInfluences) :-
                    positiveInfluence(MODULATION),
                    hasSource(MODULATION, SOURCE_SPECIES),
                    hasTarget(MODULATION, TARGET_SPECIES)."""),
            docs="If a modulation is a positiveInfluence, then there is a positive path from its source to its target.",
        ),
        Rule(
            identifier="paths_base:physical_stimulation_modulation",
            text=dedent("""\
                propagatesInfluence(SOURCE_SPECIES, TARGET_SPECIES, positivelyInfluences) :-
                    physicalStimulation(MODULATION),
                    hasSource(MODULATION, SOURCE_SPECIES),
                    hasTarget(MODULATION, TARGET_SPECIES)."""),
            docs="If a modulation is a physicalStimulation, then there is a positive path from its source to its target.",
        ),
        Rule(
            identifier="paths_base:triggering_modulation",
            text=dedent("""\
                propagatesInfluence(SOURCE_SPECIES, TARGET_SPECIES, triggers) :-
                    triggering(MODULATION),
                    hasSource(MODULATION, SOURCE_SPECIES),
                    hasTarget(MODULATION, TARGET_SPECIES)."""),
            docs="If a modulation arc is a triggers, then there is a triggering path from its source to its target (a direct arc keeps the `triggering` kind).",
        ),
        Rule(
            identifier="paths_base:inhibitor_to_product",
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
            identifier="paths_base:modulator_to_product",
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
            identifier="paths_base:unknown_catalyzer_to_product",
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
            identifier="paths_base:unknown_inhibitor_to_product",
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
            identifier="paths_base:inhibition_modulation",
            text=dedent("""\
                propagatesInfluence(SOURCE_SPECIES, TARGET_SPECIES, negativelyInfluences) :-
                    inhibition(MODULATION),
                    hasSource(MODULATION, SOURCE_SPECIES),
                    hasTarget(MODULATION, TARGET_SPECIES)."""),
            docs="If a modulation is an inhibition, then there is a negative path from its source to its target.",
        ),
        Rule(
            identifier="paths_base:negative_influence_modulation",
            text=dedent("""\
                propagatesInfluence(SOURCE_SPECIES, TARGET_SPECIES, negativelyInfluences) :-
                    negativeInfluence(MODULATION),
                    hasSource(MODULATION, SOURCE_SPECIES),
                    hasTarget(MODULATION, TARGET_SPECIES)."""),
            docs="If a modulation is a negativeInfluence, then there is a negative path from its source to its target.",
        ),
        Rule(
            identifier="paths_base:modulation_modulation",
            text=dedent("""\
                propagatesInfluence(SOURCE_SPECIES, TARGET_SPECIES, modulates) :-
                    modulation(MODULATION),
                    not catalysis(MODULATION),
                    not physicalStimulation(MODULATION),
                    not inhibition(MODULATION),
                    not triggering(MODULATION),
                    not positiveInfluence(MODULATION),
                    not negativeInfluence(MODULATION),
                    hasSource(MODULATION, SOURCE_SPECIES),
                    hasTarget(MODULATION, TARGET_SPECIES)."""),
            docs="If a modulation arc is a *bare* modulation (the generic MODULATION arc — not one of the signed/typed subtypes), then there is a modulation path from its source to its target. The negations exclude the subtypes, which `modulation` is the umbrella over.",
        ),
        Rule(
            identifier="paths_base:unknown_catalysis_modulation",
            text=dedent("""\
                propagatesInfluence(SOURCE_SPECIES, TARGET_SPECIES, unknownPositivelyInfluences) :-
                    unknownCatalysis(MODULATION),
                    hasSource(MODULATION, SOURCE_SPECIES),
                    hasTarget(MODULATION, TARGET_SPECIES)."""),
            docs="If a modulation arc is an unknown catalysis, then there is an unknown-positive path from its source to its target.",
        ),
        Rule(
            identifier="paths_base:unknown_physical_stimulation_modulation",
            text=dedent("""\
                propagatesInfluence(SOURCE_SPECIES, TARGET_SPECIES, unknownPositivelyInfluences) :-
                    unknownPhysicalStimulation(MODULATION),
                    hasSource(MODULATION, SOURCE_SPECIES),
                    hasTarget(MODULATION, TARGET_SPECIES)."""),
            docs="If a modulation arc is an unknown physical stimulation, then there is an unknown-positive path from its source to its target.",
        ),
        Rule(
            identifier="paths_base:unknown_positive_influence_modulation",
            text=dedent("""\
                propagatesInfluence(SOURCE_SPECIES, TARGET_SPECIES, unknownPositivelyInfluences) :-
                    unknownPositiveInfluence(MODULATION),
                    hasSource(MODULATION, SOURCE_SPECIES),
                    hasTarget(MODULATION, TARGET_SPECIES)."""),
            docs="If a modulation arc is an unknown positive influence, then there is an unknown-positive path from its source to its target.",
        ),
        Rule(
            identifier="paths_base:unknown_inhibition_modulation",
            text=dedent("""\
                propagatesInfluence(SOURCE_SPECIES, TARGET_SPECIES, unknownNegativelyInfluences) :-
                    unknownInhibition(MODULATION),
                    hasSource(MODULATION, SOURCE_SPECIES),
                    hasTarget(MODULATION, TARGET_SPECIES)."""),
            docs="If a modulation arc is an unknown inhibition, then there is an unknown-negative path from its source to its target.",
        ),
        Rule(
            identifier="paths_base:unknown_negative_influence_modulation",
            text=dedent("""\
                propagatesInfluence(SOURCE_SPECIES, TARGET_SPECIES, unknownNegativelyInfluences) :-
                    unknownNegativeInfluence(MODULATION),
                    hasSource(MODULATION, SOURCE_SPECIES),
                    hasTarget(MODULATION, TARGET_SPECIES)."""),
            docs="If a modulation arc is an unknown negative influence, then there is an unknown-negative path from its source to its target.",
        ),
        Rule(
            identifier="paths_base:unknown_triggering_modulation",
            text=dedent("""\
                propagatesInfluence(SOURCE_SPECIES, TARGET_SPECIES, unknownTriggers) :-
                    unknownTriggering(MODULATION),
                    hasSource(MODULATION, SOURCE_SPECIES),
                    hasTarget(MODULATION, TARGET_SPECIES)."""),
            docs="If a modulation arc is an unknown triggers, then there is an unknown-triggering path from its source to its target (kept on the direct edge; it degrades to `unknown_positive` when composed through a reaction).",
        ),
        Rule(
            identifier="paths_base:unknown_modulation_modulation",
            text=dedent("""\
                propagatesInfluence(SOURCE_SPECIES, TARGET_SPECIES, unknownModulates) :-
                    unknownModulation(MODULATION),
                    not unknownCatalysis(MODULATION),
                    not unknownPhysicalStimulation(MODULATION),
                    not unknownInhibition(MODULATION),
                    not unknownTriggering(MODULATION),
                    not unknownPositiveInfluence(MODULATION),
                    not unknownNegativeInfluence(MODULATION),
                    hasSource(MODULATION, SOURCE_SPECIES),
                    hasTarget(MODULATION, TARGET_SPECIES)."""),
            docs="If a modulation arc is a *bare* unknown modulation (not one of the unknown subtypes), then there is an unknown-modulation path from its source to its target. The negations exclude the subtypes, which `unknownModulation` is the umbrella over.",
        ),
        Rule(
            identifier="paths_base:is_directly_transformed_to",
            text=dedent("""\
                isDirectlyTransformedTo(UPSTREAM_SPECIES, DOWNSTREAM_SPECIES) :-
                    reaction(REACTION),
                    hasReactant(REACTION, REACTANT),
                    hasReferredElement(REACTANT, UPSTREAM_SPECIES),
                    hasProduct(REACTION, PRODUCT),
                    hasReferredElement(PRODUCT, DOWNSTREAM_SPECIES)."""),
            docs="The single reactant->product hop in the production graph: the upstream species is referred to by a reactant and the downstream species by a product of the same reaction. Passive voice (the reaction does the transforming, not the species) keeps it language-neutral. Feeds the shared `isTransformedTo`/`isCyclicallyTransformedTo` cycle relations that gate transitive path extension; carries no influence kind itself.",
        ),
        Rule(
            identifier="paths_base:transitive_through_reaction",
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
        Rule(
            identifier="paths_base:composes_to",
            text=dedent("""\
                composesTo(positivelyInfluences, positivelyInfluences).
                composesTo(negativelyInfluences, negativelyInfluences).
                composesTo(modulates, modulates).
                composesTo(triggers, positivelyInfluences).
                composesTo(unknownPositivelyInfluences, unknownPositivelyInfluences).
                composesTo(unknownNegativelyInfluences, unknownNegativelyInfluences).
                composesTo(unknownModulates, unknownModulates).
                composesTo(unknownTriggers, unknownPositivelyInfluences)."""),
            docs="`composesTo(INCOMING_INFLUENCE_KIND, OUTGOING_INFLUENCE_KIND)`: how an influence kind transforms when a path is extended by one reaction reactant→product hop. Identity for every kind except `triggering` (a necessary-stimulation relationship is a property of the direct edge; composed through a reaction it weakens to a plain `positive` influence) and its unknown twin `unknown_triggering` → `unknown_positive`. Modulation composes like the signed kinds (stays `modulation`).",
        ),
        ),
    },
)

_PATHS_COMPLEX_TRAVERSAL = RuleGroup(
    identifier="paths_complex_traversal",
    profiles=frozenset({"normal_no_complex", "keep_species_no_complex"}),
    depends_on=frozenset({"paths_base"}),
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
    profiles=_NON_CASQ_PROFILES,
    docs="Non-casq derivation: emits `new(activity(KEY))` for every activity key, and the internal `influences(SOURCE_KEY, TARGET_KEY, INFLUENCE_KIND)` relation from kinded paths and from consumption-based reasoning (catalyzer/physicalStimulator/trigger negatively influence each consumed reactant, inhibitor positively influences each spared reactant — in both cases only reactants that are themselves activities; the unknown modifiers contribute the unknown twins). The internal `influences/3` relation is fanned out to the typed `new(...)` heads by the shared `influence_output` group.",
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
    # The consumption/sparing rules below are CellDesigner-only (they reason
    # over reaction modifiers); SBGN-PD has no analog yet, so its variant is
    # empty. This is a known gap to implement, not a deliberate design choice:
    # the same biology in CellDesigner vs SBGN-PD currently yields different AF
    # influences. The shared `activity` and `propagatesInfluence` rules above
    # carry both languages.
    variants={
        "sbgn_pd": (),
        "celldesigner": (
        Rule(
            identifier="influences_derivation:catalyzer_consumes_reactant",
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
            identifier="influences_derivation:physical_stimulator_consumes_reactant",
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
            identifier="influences_derivation:trigger_consumes_reactant",
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
            identifier="influences_derivation:inhibitor_spares_reactant",
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
            identifier="influences_derivation:unknown_catalyzer_consumes_reactant",
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
            identifier="influences_derivation:unknown_inhibitor_spares_reactant",
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
    profiles=_NON_CASQ_PROFILES | _CASQ_PROFILES,
    docs="Shared fan-out from the internal `influences(SOURCE, TARGET, INFLUENCE_KIND)` relation to the typed `new(...)` influence heads — one rule per kind. All pipelines (non-casq and casq) converge on `influences/3`; this group is the single place that turns a kind into its output predicate.",
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
#       `paths_base` rules already emit `propagatesInfluence(OPERATOR, TARGET, INFLUENCE_KIND)` (the source
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
# path unchanged. Not registered for casq (which keeps its own pipeline and
# never reads `hasActivity`/`propagatesInfluence`).
_GATES = RuleGroup(
    identifier="gates",
    profiles=_NON_CASQ_PROFILES,
    depends_on=frozenset({"activity_base", "paths_base"}),
    docs="Authored logical operators (CellDesigner `BooleanLogicGate`, SBGN-PD `LogicalOperator`): emits the operator node (`logicalOperator/2`, token-typed), its input edges (`logicalOperatorInput/2`, each input resolved through carrier/key), and the operator-sourced influence written straight into the internal `influences/3` relation by reusing the existing `propagatesInfluence(OPERATOR, TARGET, INFLUENCE_KIND)` closure. The widened influence `source` union (`predicates._INFLUENCE_SOURCE`) lets `influence_output` fan these out with no change. Provenance-agnostic and registered for the non-casq profiles only.",
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
                docs="One rule covering both Shape A (gate is a reaction modifier) and Shape B (gate is a modulation source): the `paths_base` rules already bind a `propagatesInfluence/3` whose source is the gate, so the gate only resolves its TARGET through carrier/key and writes the influence keyed by the operator. Inherits the transitive closure (Decision D1). The `booleanLogicGate` guard is essential -- without it every species `propagatesInfluence/3` source would be read as an operator key.",
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
                docs="SBGN-PD parallel of the CellDesigner operator-sourced influence rule (Shape B: the operator is a modulation source). Reuses the `propagatesInfluence(OPERATOR, TARGET, INFLUENCE_KIND)` closure emitted by `paths_base`, resolving only the target through carrier/key.",
            ),
        ),
    },
)


_CASQ_PARTICIPATION = RuleGroup(
    identifier="casq:participation",
    profiles=_CASQ_PROFILES,
    docs="CASQ helper relations describing how species participate in reactions (`activeParticipates` for reactant/modifier roles, `participates` adding products, `isProducedSpecies`, `isModifierSpecies` and `isModulationTarget`), used as conditions of the CASQ deletion rules.",
    rules=(
        Rule(
            identifier="casq:participation:active_from_reactant",
            text=dedent("""\
                activeParticipates(SPECIES, REACTION) :-
                    reaction(REACTION),
                    hasReactant(REACTION, REACTANT),
                    hasReferredElement(REACTANT, SPECIES)."""),
            docs="A species actively participates in a reaction if it is referred to by a reactant of that reaction.",
        ),
        Rule(
            identifier="casq:participation:active_from_modifier",
            text=dedent("""\
                activeParticipates(SPECIES, REACTION) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER),
                    hasReferredElement(MODIFIER, SPECIES)."""),
            docs="A species actively participates in a reaction if it is referred to by a modifier of that reaction.",
        ),
        Rule(
            identifier="casq:participation:inherits_active",
            text="participates(SPECIES, REACTION) :- activeParticipates(SPECIES, REACTION).",
            docs="Active participation implies participation.",
        ),
        Rule(
            identifier="casq:participation:from_product",
            text=dedent("""\
                participates(SPECIES, REACTION) :-
                    reaction(REACTION),
                    hasProduct(REACTION, PRODUCT),
                    hasReferredElement(PRODUCT, SPECIES)."""),
            docs="A species participates in a reaction if it is referred to by a product of that reaction.",
        ),
        Rule(
            identifier="casq:participation:is_produced",
            text=dedent("""\
                isProducedSpecies(SPECIES) :-
                    reaction(REACTION),
                    hasProduct(REACTION, PRODUCT),
                    hasReferredElement(PRODUCT, SPECIES)."""),
            docs="A species is produced if some reaction has a product referring to it.",
        ),
        Rule(
            identifier="casq:participation:is_modifier",
            text=dedent("""\
                isModifierSpecies(SPECIES) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER),
                    hasReferredElement(MODIFIER, SPECIES)."""),
            docs="A species is a modifier if some reaction has a modifier referring to it.",
        ),
        Rule(
            identifier="casq:participation:is_modulation_target",
            text=dedent("""\
                isModulationTarget(SPECIES) :-
                    modulation(MODULATION),
                    hasTarget(MODULATION, SPECIES)."""),
            docs="A species is a modulation target if some direct influence/modulation arc (`modulation/1` -- the umbrella over every typed arc: catalysis, inhibition, positiveInfluence, ...) points at it. Unlike a reaction-product influence into a deleted species (rewired via `bridgesToProduct` or blocked by the participation/`isProducedSpecies` guards), a modulation arc has no rewiring, so deleting its target would silently drop the influence. The deletion rules use `not isModulationTarget(...)` to refuse such deletions -- mirroring casq, whose deletions guard on the species having no incoming transitions (a CellDesigner influence arc is a reaction whose product is its target, so it counts as a transition).",
        ),
    ),
)

_CASQ_DELETE = RuleGroup(
    identifier="casq:delete",
    profiles=_CASQ_PROFILES,
    depends_on=frozenset({"casq:participation"}),
    docs="CASQ-style species pruning: rules 1-4 mark species for deletion based on heterodimer associations, name-preserving single-product reactions, and transports, mirroring the CASQ tool's removal heuristics.",
    rules=(
        Rule(
            identifier="casq:delete:rule_1",
            text=dedent("""\
                delete(RECEPTOR, rule_1) :-
                    heterodimerAssociation(REACTION),
                    hasReactant(REACTION, RECEPTOR_REACTANT),
                    hasReferredElement(RECEPTOR_REACTANT, RECEPTOR),
                    receptor(RECEPTOR),
                    not isModulationTarget(RECEPTOR),
                    hasReactant(REACTION, PARTNER_REACTANT),
                    hasReferredElement(PARTNER_REACTANT, PARTNER),
                    RECEPTOR != PARTNER,
                    #count{ REACTANT_SPECIES : hasReactant(REACTION, REACTANT), hasReferredElement(REACTANT, REACTANT_SPECIES) } = 2,
                    #count{ RECEPTOR_REACTION : participates(RECEPTOR, RECEPTOR_REACTION) } = 1,
                    #count{ PARTNER_REACTION : participates(PARTNER, PARTNER_REACTION) } = 1."""),
            docs="A receptor in a 2-reactant heterodimer association where receptor and partner each participate in only this reaction is deleted (rule_1), unless the receptor is the target of a modulation arc -- deleting it would silently drop that influence (which, unlike a reaction-product influence, is never rewired), so `not isModulationTarget(RECEPTOR)` blocks the deletion, mirroring casq's no-incoming-transitions guard.",
        ),
        Rule(
            identifier="casq:delete:rule_2",
            text=dedent("""\
                delete(SPECIES_1, rule_2) :-
                    heterodimerAssociation(REACTION),
                    hasReactant(REACTION, FIRST_REACTANT), hasReferredElement(FIRST_REACTANT, SPECIES_1),
                    hasReactant(REACTION, SECOND_REACTANT), hasReferredElement(SECOND_REACTANT, SPECIES_2),
                    SPECIES_1 != SPECIES_2,
                    not receptor(SPECIES_1),
                    not receptor(SPECIES_2),
                    not isModulationTarget(SPECIES_1),
                    not isModulationTarget(SPECIES_2),
                    #count{ REACTANT_SPECIES : hasReactant(REACTION, REACTANT), hasReferredElement(REACTANT, REACTANT_SPECIES) } = 2,
                    #count{ FIRST_SPECIES_REACTION : activeParticipates(SPECIES_1, FIRST_SPECIES_REACTION) } = 1,
                    #count{ SECOND_SPECIES_REACTION : activeParticipates(SPECIES_2, SECOND_SPECIES_REACTION) } = 1."""),
            docs="In a 2-reactant heterodimer association where neither reactant is a receptor and each actively participates only in this reaction, both reactant species are deleted (rule_2). The rule fires symmetrically for each side. The deletion is blocked if *either* reactant is the target of a modulation arc (`not isModulationTarget(SPECIES_1)`, `not isModulationTarget(SPECIES_2)`): such an influence is never rewired, so dropping either reactant would silently lose it -- and casq likewise refuses to delete the pair when either has incoming transitions.",
        ),
        Rule(
            identifier="casq:delete:rule_3",
            text=dedent("""\
                delete(REACTANT_SPECIES, rule_3) :-
                    reaction(REACTION),
                    hasReactant(REACTION, REACTANT), hasReferredElement(REACTANT, REACTANT_SPECIES),
                    hasProduct(REACTION, PRODUCT), hasReferredElement(PRODUCT, PRODUCT_SPECIES),
                    REACTANT_SPECIES != PRODUCT_SPECIES,
                    hasName(REACTANT_SPECIES, NAME), hasName(PRODUCT_SPECIES, NAME),
                    not isProducedSpecies(REACTANT_SPECIES),
                    not isModifierSpecies(REACTANT_SPECIES),
                    not isModulationTarget(REACTANT_SPECIES),
                    #count{ COUNTED_PRODUCT_SPECIES : hasProduct(REACTION, COUNTED_PRODUCT), hasReferredElement(COUNTED_PRODUCT, COUNTED_PRODUCT_SPECIES) } = 1,
                    #count{ CONSUMING_REACTION : hasReactant(CONSUMING_REACTION, CONSUMING_REACTANT), hasReferredElement(CONSUMING_REACTANT, REACTANT_SPECIES) } = 1."""),
            docs="In a single-product reaction where reactant and product share a name, the reactant is deleted (rule_3) if it is not produced anywhere else, never appears as a modifier, is not the target of a modulation arc (`not isModulationTarget` -- such an influence is never rewired and would be silently lost), and is consumed only by this reaction.",
        ),
        Rule(
            identifier="casq:delete:rule_4",
            text=dedent("""\
                delete(REACTANT_SPECIES, rule_4) :-
                    transport(REACTION),
                    hasReactant(REACTION, REACTANT), hasReferredElement(REACTANT, REACTANT_SPECIES),
                    hasProduct(REACTION, PRODUCT), hasReferredElement(PRODUCT, PRODUCT_SPECIES),
                    REACTANT_SPECIES != PRODUCT_SPECIES,
                    hasName(REACTANT_SPECIES, NAME), hasName(PRODUCT_SPECIES, NAME),
                    not isModulationTarget(REACTANT_SPECIES),
                    #count{ PARTICIPATED_REACTION : activeParticipates(REACTANT_SPECIES, PARTICIPATED_REACTION) } = 1,
                    #count{ COUNTED_PRODUCT_SPECIES : hasProduct(REACTION, COUNTED_PRODUCT), hasReferredElement(COUNTED_PRODUCT, COUNTED_PRODUCT_SPECIES) } = 1."""),
            docs="In a single-product transport where reactant and product share a name and the reactant actively participates only in this reaction, the reactant is deleted (rule_4), unless it is the target of a modulation arc (`not isModulationTarget` -- that influence has no rewiring and would be silently dropped).",
        ),
    ),
)

_CASQ_BRIDGED_PRODUCT = RuleGroup(
    identifier="casq:bridged_product",
    profiles=_CASQ_PROFILES,
    depends_on=frozenset({"casq:delete"}),
    docs="One-hop rewiring across species deleted by `rule_2` or `rule_4`: if reaction R1 produces a deleted species that reaction R2 consumes to make P, then R1 is treated as also producing P so influences can skip the deleted intermediate.",
    rules=(
        Rule(
            identifier="casq:bridged_product:from_rule_2",
            text=dedent("""\
                bridgesToProduct(REACTION_1, PRODUCT_SPECIES) :-
                    delete(DELETED_SPECIES, rule_2),
                    reaction(REACTION_1),
                    hasProduct(REACTION_1, REACTION_1_PRODUCT), hasReferredElement(REACTION_1_PRODUCT, DELETED_SPECIES),
                    reaction(REACTION_2),
                    hasReactant(REACTION_2, REACTION_2_REACTANT), hasReferredElement(REACTION_2_REACTANT, DELETED_SPECIES),
                    hasProduct(REACTION_2, REACTION_2_PRODUCT), hasReferredElement(REACTION_2_PRODUCT, PRODUCT_SPECIES)."""),
            docs="One-hop rewiring across a species deleted by rule_2: REACTION_1 produces DELETED_SPECIES and REACTION_2 consumes DELETED_SPECIES and produces PRODUCT_SPECIES, so REACTION_1 is treated as also producing PRODUCT_SPECIES.",
        ),
        Rule(
            identifier="casq:bridged_product:from_rule_4",
            text=dedent("""\
                bridgesToProduct(REACTION_1, PRODUCT_SPECIES) :-
                    delete(DELETED_SPECIES, rule_4),
                    reaction(REACTION_1),
                    hasProduct(REACTION_1, REACTION_1_PRODUCT), hasReferredElement(REACTION_1_PRODUCT, DELETED_SPECIES),
                    reaction(REACTION_2),
                    hasReactant(REACTION_2, REACTION_2_REACTANT), hasReferredElement(REACTION_2_REACTANT, DELETED_SPECIES),
                    hasProduct(REACTION_2, REACTION_2_PRODUCT), hasReferredElement(REACTION_2_PRODUCT, PRODUCT_SPECIES)."""),
            docs="One-hop rewiring across a species deleted by rule_4 (analog of casq:bridged_product:from_rule_2 for transport-driven deletes).",
        ),
    ),
)

_CASQ_ACTIVITY = RuleGroup(
    identifier="casq:activity",
    profiles=_CASQ_PROFILES,
    depends_on=frozenset({"casq:delete"}),
    docs="`casq` activity emission: every species resolves to its outermost top-level complex (recursively), and only surviving top-level entities become activities. A standalone top-level species resolves to itself and is keyed `keptSpeciesKey(SELF)`. A subunit -- at any nesting depth -- resolves to the outermost complex that contains it and is keyed by *that complex's* `keptSpeciesKey`; it is never emitted as its own activity. This mirrors casq, which collapses a subunit's participation onto its complex (a subunit is never a node in its own right). When the resolved top-level complex is deleted, the subunit has no surviving carrier and contributes no key, so its participation is dropped -- matching casq's deletion behaviour. Influence endpoints are resolved through `hasActivityKey` (in casq:influences), so an influence touching a subunit references its top-level complex's `keptSpeciesKey` key. The recursive resolution also fixes a casq bug: casq collapses only one nesting level and silently drops influences from more deeply nested subunits.",
    rules=(
        Rule(
            identifier="casq:activity:top_level_self",
            text=dedent("""\
                resolvesToTopLevel(SPECIES, SPECIES) :-
                    species(SPECIES),
                    not hasSubunit(_, SPECIES)."""),
            docs="A species that is not a subunit of any complex is its own top-level entity.",
        ),
        Rule(
            identifier="casq:activity:top_level_recursive",
            text=dedent("""\
                resolvesToTopLevel(SUBUNIT, TOPLEVEL) :-
                    hasSubunit(PARENT_COMPLEX, SUBUNIT),
                    resolvesToTopLevel(PARENT_COMPLEX, TOPLEVEL)."""),
            docs="A subunit resolves to the same top-level entity as its parent complex, recursively through nested complexes -- so a subunit at any depth resolves to its outermost complex.",
        ),
        Rule(
            identifier="casq:activity:key",
            text=dedent("""\
                hasActivityKey(SPECIES, keptSpeciesKey(TOPLEVEL)) :-
                    resolvesToTopLevel(SPECIES, TOPLEVEL),
                    not delete(TOPLEVEL, _)."""),
            docs="A species is keyed by the `keptSpeciesKey` activity of its surviving top-level complex (or of itself, when it is top-level). Subunits never get their own activity; when the top-level complex is deleted there is no key, so nothing it contains contributes.",
        ),
        Rule(
            identifier="casq:activity:emit",
            text="new(activity(KEY)) :- hasActivityKey(_, KEY).",
            docs="Every distinct activity key emits an activity node with that key.",
        ),
    ),
)

_CASQ_INFLUENCES = RuleGroup(
    identifier="casq:influences",
    profiles=_CASQ_PROFILES,
    depends_on=frozenset({"casq:bridged_product", "casq:activity"}),
    docs="Casq-specific influence emission, written into the internal `influences(SOURCE, TARGET, INFLUENCE_KIND)` relation (the shared `influence_output` group fans it out to the typed heads). Directly wires reactant/catalyzer/stimulator → product (positivelyInfluences), trigger → product (triggers), inhibitor → product (negativelyInfluences), the generic and unknown modifiers (modulation / unknown_*), and the species→species modulation arcs, each with a bridged-product variant to route across deleted intermediates. Because a bridged product crosses a reaction boundary, the bridged trigger variant degrades to positive (mirroring `composesTo`).",
    rules=(
        Rule(
            identifier="casq:influences:reactant_to_product",
            text=dedent("""\
                influences(SOURCE_KEY, TARGET_KEY, positivelyInfluences) :-
                    reaction(REACTION),
                    hasReactant(REACTION, REACTANT), hasReferredElement(REACTANT, SOURCE),
                    hasProduct(REACTION, PRODUCT), hasReferredElement(PRODUCT, TARGET),
                    hasActivityKey(SOURCE, SOURCE_KEY), hasActivityKey(TARGET, TARGET_KEY)."""),
            docs="A reactant of a reaction positively influences a product of that reaction.",
        ),
        Rule(
            identifier="casq:influences:reactant_to_bridged_product",
            text=dedent("""\
                influences(SOURCE_KEY, TARGET_KEY, positivelyInfluences) :-
                    reaction(REACTION),
                    hasReactant(REACTION, REACTANT), hasReferredElement(REACTANT, SOURCE),
                    bridgesToProduct(REACTION, TARGET),
                    hasActivityKey(SOURCE, SOURCE_KEY), hasActivityKey(TARGET, TARGET_KEY)."""),
            docs="Same as casq:influences:reactant_to_product but routed through a deleted intermediate via bridgesToProduct.",
        ),
        Rule(
            identifier="casq:influences:catalyzer_to_product",
            text=dedent("""\
                influences(SOURCE_KEY, TARGET_KEY, positivelyInfluences) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), catalyzer(MODIFIER),
                    hasReferredElement(MODIFIER, SOURCE),
                    hasProduct(REACTION, PRODUCT), hasReferredElement(PRODUCT, TARGET),
                    hasActivityKey(SOURCE, SOURCE_KEY), hasActivityKey(TARGET, TARGET_KEY)."""),
            docs="A catalyzer of a reaction positively influences a product of that reaction.",
        ),
        Rule(
            identifier="casq:influences:catalyzer_to_bridged_product",
            text=dedent("""\
                influences(SOURCE_KEY, TARGET_KEY, positivelyInfluences) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), catalyzer(MODIFIER),
                    hasReferredElement(MODIFIER, SOURCE),
                    bridgesToProduct(REACTION, TARGET),
                    hasActivityKey(SOURCE, SOURCE_KEY), hasActivityKey(TARGET, TARGET_KEY)."""),
            docs="Bridged variant of casq:influences:catalyzer_to_product.",
        ),
        Rule(
            identifier="casq:influences:physical_stimulator_to_product",
            text=dedent("""\
                influences(SOURCE_KEY, TARGET_KEY, positivelyInfluences) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), physicalStimulator(MODIFIER),
                    hasReferredElement(MODIFIER, SOURCE),
                    hasProduct(REACTION, PRODUCT), hasReferredElement(PRODUCT, TARGET),
                    hasActivityKey(SOURCE, SOURCE_KEY), hasActivityKey(TARGET, TARGET_KEY)."""),
            docs="A physical stimulator of a reaction positively influences a product of that reaction.",
        ),
        Rule(
            identifier="casq:influences:physical_stimulator_to_bridged_product",
            text=dedent("""\
                influences(SOURCE_KEY, TARGET_KEY, positivelyInfluences) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), physicalStimulator(MODIFIER),
                    hasReferredElement(MODIFIER, SOURCE),
                    bridgesToProduct(REACTION, TARGET),
                    hasActivityKey(SOURCE, SOURCE_KEY), hasActivityKey(TARGET, TARGET_KEY)."""),
            docs="Bridged variant of casq:influences:physical_stimulator_to_product.",
        ),
        Rule(
            identifier="casq:influences:trigger_to_product",
            text=dedent("""\
                influences(SOURCE_KEY, TARGET_KEY, triggers) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), trigger(MODIFIER),
                    hasReferredElement(MODIFIER, SOURCE),
                    hasProduct(REACTION, PRODUCT), hasReferredElement(PRODUCT, TARGET),
                    hasActivityKey(SOURCE, SOURCE_KEY), hasActivityKey(TARGET, TARGET_KEY)."""),
            docs="A trigger of a reaction triggers a product of that reaction (direct trigger→product, so it keeps the triggering kind).",
        ),
        Rule(
            identifier="casq:influences:trigger_to_bridged_product",
            text=dedent("""\
                influences(SOURCE_KEY, TARGET_KEY, positivelyInfluences) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), trigger(MODIFIER),
                    hasReferredElement(MODIFIER, SOURCE),
                    bridgesToProduct(REACTION, TARGET),
                    hasActivityKey(SOURCE, SOURCE_KEY), hasActivityKey(TARGET, TARGET_KEY)."""),
            docs="Bridged variant of casq:influences:trigger_to_product. A bridged product crosses a reaction boundary, so the triggering degrades to a plain positive influence (the casq analog of composesTo).",
        ),
        Rule(
            identifier="casq:influences:inhibitor_to_product",
            text=dedent("""\
                influences(SOURCE_KEY, TARGET_KEY, negativelyInfluences) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), inhibitor(MODIFIER),
                    hasReferredElement(MODIFIER, SOURCE),
                    hasProduct(REACTION, PRODUCT), hasReferredElement(PRODUCT, TARGET),
                    hasActivityKey(SOURCE, SOURCE_KEY), hasActivityKey(TARGET, TARGET_KEY)."""),
            docs="An inhibitor of a reaction negatively influences a product of that reaction.",
        ),
        Rule(
            identifier="casq:influences:inhibitor_to_bridged_product",
            text=dedent("""\
                influences(SOURCE_KEY, TARGET_KEY, negativelyInfluences) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), inhibitor(MODIFIER),
                    hasReferredElement(MODIFIER, SOURCE),
                    bridgesToProduct(REACTION, TARGET),
                    hasActivityKey(SOURCE, SOURCE_KEY), hasActivityKey(TARGET, TARGET_KEY)."""),
            docs="Bridged variant of casq:influences:inhibitor_to_product.",
        ),
        Rule(
            identifier="casq:influences:modulator_to_product",
            text=dedent("""\
                influences(SOURCE_KEY, TARGET_KEY, modulates) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), modulator(MODIFIER),
                    not physicalStimulator(MODIFIER),
                    not inhibitor(MODIFIER),
                    not trigger(MODIFIER),
                    hasReferredElement(MODIFIER, SOURCE),
                    hasProduct(REACTION, PRODUCT), hasReferredElement(PRODUCT, TARGET),
                    hasActivityKey(SOURCE, SOURCE_KEY), hasActivityKey(TARGET, TARGET_KEY)."""),
            docs="A bare modulator (generic MODULATION modifier) of a reaction modulates a product of that reaction.",
        ),
        Rule(
            identifier="casq:influences:modulator_to_bridged_product",
            text=dedent("""\
                influences(SOURCE_KEY, TARGET_KEY, modulates) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), modulator(MODIFIER),
                    not physicalStimulator(MODIFIER),
                    not inhibitor(MODIFIER),
                    not trigger(MODIFIER),
                    hasReferredElement(MODIFIER, SOURCE),
                    bridgesToProduct(REACTION, TARGET),
                    hasActivityKey(SOURCE, SOURCE_KEY), hasActivityKey(TARGET, TARGET_KEY)."""),
            docs="Bridged variant of casq:influences:modulator_to_product.",
        ),
        Rule(
            identifier="casq:influences:unknown_catalyzer_to_product",
            text=dedent("""\
                influences(SOURCE_KEY, TARGET_KEY, unknownPositivelyInfluences) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), unknownCatalyzer(MODIFIER),
                    hasReferredElement(MODIFIER, SOURCE),
                    hasProduct(REACTION, PRODUCT), hasReferredElement(PRODUCT, TARGET),
                    hasActivityKey(SOURCE, SOURCE_KEY), hasActivityKey(TARGET, TARGET_KEY)."""),
            docs="An unknown catalyzer of a reaction unknown-positively influences a product of that reaction.",
        ),
        Rule(
            identifier="casq:influences:unknown_catalyzer_to_bridged_product",
            text=dedent("""\
                influences(SOURCE_KEY, TARGET_KEY, unknownPositivelyInfluences) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), unknownCatalyzer(MODIFIER),
                    hasReferredElement(MODIFIER, SOURCE),
                    bridgesToProduct(REACTION, TARGET),
                    hasActivityKey(SOURCE, SOURCE_KEY), hasActivityKey(TARGET, TARGET_KEY)."""),
            docs="Bridged variant of casq:influences:unknown_catalyzer_to_product.",
        ),
        Rule(
            identifier="casq:influences:unknown_inhibitor_to_product",
            text=dedent("""\
                influences(SOURCE_KEY, TARGET_KEY, unknownNegativelyInfluences) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), unknownInhibitor(MODIFIER),
                    hasReferredElement(MODIFIER, SOURCE),
                    hasProduct(REACTION, PRODUCT), hasReferredElement(PRODUCT, TARGET),
                    hasActivityKey(SOURCE, SOURCE_KEY), hasActivityKey(TARGET, TARGET_KEY)."""),
            docs="An unknown inhibitor of a reaction unknown-negatively influences a product of that reaction.",
        ),
        Rule(
            identifier="casq:influences:unknown_inhibitor_to_bridged_product",
            text=dedent("""\
                influences(SOURCE_KEY, TARGET_KEY, unknownNegativelyInfluences) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), unknownInhibitor(MODIFIER),
                    hasReferredElement(MODIFIER, SOURCE),
                    bridgesToProduct(REACTION, TARGET),
                    hasActivityKey(SOURCE, SOURCE_KEY), hasActivityKey(TARGET, TARGET_KEY)."""),
            docs="Bridged variant of casq:influences:unknown_inhibitor_to_product.",
        ),
        Rule(
            identifier="casq:influences:catalysis_modulation",
            text=dedent("""\
                influences(SOURCE_KEY, TARGET_KEY, positivelyInfluences) :-
                    catalysis(MODULATION),
                    hasSource(MODULATION, SOURCE), hasTarget(MODULATION, TARGET),
                    hasActivityKey(SOURCE, SOURCE_KEY), hasActivityKey(TARGET, TARGET_KEY)."""),
            docs="A catalysis modulation arc emits a positive influence between its source and target activities. (The predicate is `catalysis`; an earlier `catalyzis` typo silently disabled this rule.)",
        ),
        Rule(
            identifier="casq:influences:positive_influence_modulation",
            text=dedent("""\
                influences(SOURCE_KEY, TARGET_KEY, positivelyInfluences) :-
                    positiveInfluence(MODULATION),
                    hasSource(MODULATION, SOURCE), hasTarget(MODULATION, TARGET),
                    hasActivityKey(SOURCE, SOURCE_KEY), hasActivityKey(TARGET, TARGET_KEY)."""),
            docs="A positiveInfluence modulation arc emits a positive influence between its source and target activities.",
        ),
        Rule(
            identifier="casq:influences:physical_stimulation_modulation",
            text=dedent("""\
                influences(SOURCE_KEY, TARGET_KEY, positivelyInfluences) :-
                    physicalStimulation(MODULATION),
                    hasSource(MODULATION, SOURCE), hasTarget(MODULATION, TARGET),
                    hasActivityKey(SOURCE, SOURCE_KEY), hasActivityKey(TARGET, TARGET_KEY)."""),
            docs="A physicalStimulation modulation arc emits a positive influence between its source and target activities.",
        ),
        Rule(
            identifier="casq:influences:triggering_modulation",
            text=dedent("""\
                influences(SOURCE_KEY, TARGET_KEY, triggers) :-
                    triggering(MODULATION),
                    hasSource(MODULATION, SOURCE), hasTarget(MODULATION, TARGET),
                    hasActivityKey(SOURCE, SOURCE_KEY), hasActivityKey(TARGET, TARGET_KEY)."""),
            docs="A triggering modulation arc emits a triggering influence between its source and target activities (direct arc keeps the triggering kind).",
        ),
        Rule(
            identifier="casq:influences:inhibition_modulation",
            text=dedent("""\
                influences(SOURCE_KEY, TARGET_KEY, negativelyInfluences) :-
                    inhibition(MODULATION),
                    hasSource(MODULATION, SOURCE), hasTarget(MODULATION, TARGET),
                    hasActivityKey(SOURCE, SOURCE_KEY), hasActivityKey(TARGET, TARGET_KEY)."""),
            docs="An inhibition modulation arc emits a negative influence between its source and target activities.",
        ),
        Rule(
            identifier="casq:influences:negative_influence_modulation",
            text=dedent("""\
                influences(SOURCE_KEY, TARGET_KEY, negativelyInfluences) :-
                    negativeInfluence(MODULATION),
                    hasSource(MODULATION, SOURCE), hasTarget(MODULATION, TARGET),
                    hasActivityKey(SOURCE, SOURCE_KEY), hasActivityKey(TARGET, TARGET_KEY)."""),
            docs="A negativeInfluence modulation arc emits a negative influence between its source and target activities.",
        ),
        Rule(
            identifier="casq:influences:modulation_modulation",
            text=dedent("""\
                influences(SOURCE_KEY, TARGET_KEY, modulates) :-
                    modulation(MODULATION),
                    not catalysis(MODULATION),
                    not physicalStimulation(MODULATION),
                    not inhibition(MODULATION),
                    not triggering(MODULATION),
                    not positiveInfluence(MODULATION),
                    not negativeInfluence(MODULATION),
                    hasSource(MODULATION, SOURCE), hasTarget(MODULATION, TARGET),
                    hasActivityKey(SOURCE, SOURCE_KEY), hasActivityKey(TARGET, TARGET_KEY)."""),
            docs="A bare modulation arc (the generic MODULATION arc, not one of the typed subtypes) emits a modulation influence between its source and target activities.",
        ),
        Rule(
            identifier="casq:influences:unknown_catalysis_modulation",
            text=dedent("""\
                influences(SOURCE_KEY, TARGET_KEY, unknownPositivelyInfluences) :-
                    unknownCatalysis(MODULATION),
                    hasSource(MODULATION, SOURCE), hasTarget(MODULATION, TARGET),
                    hasActivityKey(SOURCE, SOURCE_KEY), hasActivityKey(TARGET, TARGET_KEY)."""),
            docs="An unknown catalysis modulation arc emits an unknown-positive influence between its source and target activities.",
        ),
        Rule(
            identifier="casq:influences:unknown_physical_stimulation_modulation",
            text=dedent("""\
                influences(SOURCE_KEY, TARGET_KEY, unknownPositivelyInfluences) :-
                    unknownPhysicalStimulation(MODULATION),
                    hasSource(MODULATION, SOURCE), hasTarget(MODULATION, TARGET),
                    hasActivityKey(SOURCE, SOURCE_KEY), hasActivityKey(TARGET, TARGET_KEY)."""),
            docs="An unknown physical stimulation modulation arc emits an unknown-positive influence between its source and target activities.",
        ),
        Rule(
            identifier="casq:influences:unknown_positive_influence_modulation",
            text=dedent("""\
                influences(SOURCE_KEY, TARGET_KEY, unknownPositivelyInfluences) :-
                    unknownPositiveInfluence(MODULATION),
                    hasSource(MODULATION, SOURCE), hasTarget(MODULATION, TARGET),
                    hasActivityKey(SOURCE, SOURCE_KEY), hasActivityKey(TARGET, TARGET_KEY)."""),
            docs="An unknown positive influence modulation arc emits an unknown-positive influence between its source and target activities.",
        ),
        Rule(
            identifier="casq:influences:unknown_inhibition_modulation",
            text=dedent("""\
                influences(SOURCE_KEY, TARGET_KEY, unknownNegativelyInfluences) :-
                    unknownInhibition(MODULATION),
                    hasSource(MODULATION, SOURCE), hasTarget(MODULATION, TARGET),
                    hasActivityKey(SOURCE, SOURCE_KEY), hasActivityKey(TARGET, TARGET_KEY)."""),
            docs="An unknown inhibition modulation arc emits an unknown-negative influence between its source and target activities.",
        ),
        Rule(
            identifier="casq:influences:unknown_negative_influence_modulation",
            text=dedent("""\
                influences(SOURCE_KEY, TARGET_KEY, unknownNegativelyInfluences) :-
                    unknownNegativeInfluence(MODULATION),
                    hasSource(MODULATION, SOURCE), hasTarget(MODULATION, TARGET),
                    hasActivityKey(SOURCE, SOURCE_KEY), hasActivityKey(TARGET, TARGET_KEY)."""),
            docs="An unknown negative influence modulation arc emits an unknown-negative influence between its source and target activities.",
        ),
        Rule(
            identifier="casq:influences:unknown_triggering_modulation",
            text=dedent("""\
                influences(SOURCE_KEY, TARGET_KEY, unknownTriggers) :-
                    unknownTriggering(MODULATION),
                    hasSource(MODULATION, SOURCE), hasTarget(MODULATION, TARGET),
                    hasActivityKey(SOURCE, SOURCE_KEY), hasActivityKey(TARGET, TARGET_KEY)."""),
            docs="An unknown triggering modulation arc emits an unknown-triggering influence between its source and target activities.",
        ),
        Rule(
            identifier="casq:influences:unknown_modulation_modulation",
            text=dedent("""\
                influences(SOURCE_KEY, TARGET_KEY, unknownModulates) :-
                    unknownModulation(MODULATION),
                    not unknownCatalysis(MODULATION),
                    not unknownPhysicalStimulation(MODULATION),
                    not unknownInhibition(MODULATION),
                    not unknownTriggering(MODULATION),
                    not unknownPositiveInfluence(MODULATION),
                    not unknownNegativeInfluence(MODULATION),
                    hasSource(MODULATION, SOURCE), hasTarget(MODULATION, TARGET),
                    hasActivityKey(SOURCE, SOURCE_KEY), hasActivityKey(TARGET, TARGET_KEY)."""),
            docs="A bare unknown modulation arc (not one of the unknown subtypes) emits an unknown-modulation influence between its source and target activities.",
        ),
    ),
)


def _build_registry() -> RuleRegistry:
    registry = RuleRegistry()
    registry.register(
        [
            _ACTIVITY_BASE,
            _TOPOLOGY,
            _TOP_LEVEL,
            _PREPARATION_KEEP_SPECIES,
            _PREPARATION_KEEP_SPECIES_NO_COMPLEX,
            _PREPARATION_NO_COMPLEX,
            _PREPARATION_NORMAL,
            _PATHS_BASE,
            _PATHS_COMPLEX_TRAVERSAL,
            _INFLUENCES_DERIVATION,
            _INFLUENCE_OUTPUT,
            _GATES,
            _CASQ_PARTICIPATION,
            _CASQ_DELETE,
            _CASQ_BRIDGED_PRODUCT,
            _CASQ_ACTIVITY,
            _CASQ_INFLUENCES,
        ]
    )
    issues = registry.validate()
    if issues:
        raise RuntimeError(
            "pd2af rule registry is invalid:\n"
            + "\n".join(f"  - {issue}" for issue in issues)
        )
    return registry


def build_program(profile: str, language: str = "celldesigner") -> str:
    """Return the composed ASP program text for the given profile and
    input ``language``.

    Modes are aspcompose *profiles*; the input language is an aspcompose
    *variant*. Language-agnostic rules live in each group's ``rules`` and
    are emitted for every language; language-specific rules live in
    ``variants={"celldesigner": ..., "sbgn_pd": ...}`` and are selected
    here by ``resolve(variant=language)``.
    """
    registry = _build_registry()
    plan = CollectionPlan(registry)
    plan.add_profile(profile)
    return "\n".join(rule.text for rule in plan.resolve(variant=language))
