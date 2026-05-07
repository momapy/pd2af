"""ASP rule composition for the pd2af transformation modes.

The non-CASQ rules are organised in three layers:

* **topology** (shared across all non-CASQ profiles) — structural
  helpers: ``isSubunit``, ``topLevelAncestor``, ``hasActiveDescendant``,
  ``hasSomeCompartment``, ``hasSomeTemplate``, ``effectiveCompartment``.
* **preparation** (one rule group per non-CASQ mode) — emits
  ``activityCarrier(RAW_SPECIES, ACTIVITY_BEARER)`` and
  ``activityKey(ACTIVITY_BEARER, KEY)`` where ``KEY`` is one of three
  per-species wrappers: ``kept_species/1``, ``new_species_from_subunit/1``,
  or ``new_species_from_template/1``.
* **derivation** (shared across all non-CASQ profiles) — emits
  ``new(activity(KEY))`` and ``new(positivelyInfluences(...))`` /
  ``new(negativelyInfluences(...))`` from ``activityCarrier`` /
  ``activityKey``.

The ``casq`` profile keeps its own pipeline (deletion rules,
bridged-product rewiring, direct reactant→product / modifier→product
influence emission). It uses ``kept_species/1`` for every surviving
species's activity key, sharing only the predicate type with the
non-CASQ layers.
"""

from textwrap import dedent

from aspcompose import CollectionPlan, Rule, RuleGroup, RuleRegistry

_NON_CASQ_PROFILES = frozenset(
    {"normal", "no_complex", "keep_species", "keep_species_no_complex"}
)

_CASQ_PROFILES = frozenset({"casq"})

_ACTIVITY_BASE = RuleGroup(
    identifier="activity_base",
    profiles=_NON_CASQ_PROFILES | _CASQ_PROFILES,
    docs="Base rules deriving `hasActivity(SPECIES, REASON)` from PD signals: the explicit `hasActive` flag, an active structural state, the phenotype glyph, and being the source/modulator of a modulation arc.",
    rules=(
        Rule(
            identifier="activity_base:from_active_flag",
            text="hasActivity(SPECIES, active) :- species(SPECIES), hasActive(SPECIES, 1).",
            docs="If a species has its `hasActive` flag set to 1, then it has activity, with reason `active`.",
        ),
        Rule(
            identifier="activity_base:from_active_structural_state",
            text=dedent("""\
                hasActivity(SPECIES, structural_state_active) :-
                    species(SPECIES),
                    hasStructuralState(SPECIES, STRUCTURAL_STATE),
                    hasValue(STRUCTURAL_STATE, "active")."""),
            docs='If a species carries a structural state whose value is "active", then it has activity, with reason `structural_state_active`.',
        ),
        Rule(
            identifier="activity_base:from_phenotype",
            text="hasActivity(PHENOTYPE, phenotype) :- phenotype(PHENOTYPE).",
            docs="If a species is a phenotype, then it has activity, with reason `phenotype`.",
        ),
        Rule(
            identifier="activity_base:from_modulation_source",
            text=dedent("""\
                hasActivity(SOURCE, modulates(SOURCE, TARGET)) :-
                    species(SOURCE),
                    modulation(MODULATION),
                    hasSource(MODULATION, SOURCE),
                    hasTarget(MODULATION, TARGET)."""),
            docs="If a species is the source of a modulation with some target, then it has activity, with reason `modulates(source, target)`.",
        ),
        Rule(
            identifier="activity_base:from_reaction_modulator",
            text=dedent("""\
                hasActivity(SOURCE, modulates(SOURCE, TARGET)) :-
                    species(SOURCE),
                    modulator(MODULATOR),
                    hasReferredSpecies(MODULATOR, SOURCE),
                    hasModifier(TARGET, MODULATOR)."""),
            docs="If a species is referred to by a modulator that is a modifier of some target, then it has activity, with reason `modulates(source, target)`.",
        ),
    ),
)

_TOPOLOGY = RuleGroup(
    identifier="topology",
    profiles=_NON_CASQ_PROFILES,
    docs="Mode-agnostic structural helpers shared by all non-CASQ profiles: `isSubunit`, `hasActiveDescendant`, `hasSomeTemplate`.",
    rules=(
        Rule(
            identifier="topology:is_subunit",
            text="isSubunit(SUBUNIT) :- hasSubunit(_, SUBUNIT).",
            docs="A species is a subunit if it appears on the right side of any `hasSubunit` relation.",
        ),
        Rule(
            identifier="topology:has_active_descendant_direct",
            text=dedent("""\
                hasActiveDescendant(COMPLEX) :-
                    complex(COMPLEX),
                    hasSubunit(COMPLEX, SUBUNIT),
                    hasActivity(SUBUNIT, _)."""),
            docs="A complex has an active descendant if any direct subunit has activity.",
        ),
        Rule(
            identifier="topology:has_active_descendant_transitive",
            text=dedent("""\
                hasActiveDescendant(COMPLEX) :-
                    complex(COMPLEX),
                    hasSubunit(COMPLEX, NESTED_COMPLEX),
                    hasActiveDescendant(NESTED_COMPLEX)."""),
            docs="The `hasActiveDescendant` relation is transitive through complex containment.",
        ),
        Rule(
            identifier="topology:has_some_template",
            text="hasSomeTemplate(SPECIES) :- hasTemplate(SPECIES, _).",
            docs="A species has some template if it is linked to any template.",
        ),
    ),
)

_PREPARATION_KEEP_SPECIES = RuleGroup(
    identifier="preparation:keep_species",
    profiles=frozenset({"keep_species"}),
    depends_on=frozenset({"activity_base", "topology"}),
    docs="`keep-species` preparation: every species (subunit or top-level) with activity contributes its own `kept_species(SELF)` activity; influences stay on the original species (no rerouting). A complex with any active descendant additionally inherits an activity, so the assembly is also represented in the AF.",
    rules=(
        Rule(
            identifier="preparation:keep_species:complex_inherits_activity",
            text=dedent("""\
                hasActivity(COMPLEX, has_active_descendant) :-
                    complex(COMPLEX),
                    hasActiveDescendant(COMPLEX)."""),
            docs="A complex with any (transitive) active subunit is itself active in the AF view, with reason `has_active_descendant`.",
        ),
        Rule(
            identifier="preparation:keep_species:carrier",
            text="activityCarrier(SPECIES, SPECIES) :- species(SPECIES).",
            docs="Each species is its own activity carrier (no rerouting).",
        ),
        Rule(
            identifier="preparation:keep_species:key",
            text=dedent("""\
                activityKey(SPECIES, kept_species(SPECIES)) :-
                    hasActivity(SPECIES, _)."""),
            docs="Every species with activity is keyed by `kept_species(SELF)`.",
        ),
    ),
)

_PREPARATION_KEEP_SPECIES_NO_COMPLEX = RuleGroup(
    identifier="preparation:keep_species_no_complex",
    profiles=frozenset({"keep_species_no_complex"}),
    depends_on=frozenset({"activity_base", "topology"}),
    docs="`keep-species-no-complex` preparation: a complex with any active descendant is suppressed; every other species with activity contributes its own `kept_species(SELF)` activity. Subunits of suppressed complexes are promoted to top-level activities by reusing the subunit species itself; the layout strategy attaches them at top level.",
    rules=(
        Rule(
            identifier="preparation:keep_species_no_complex:suppressed",
            text=dedent("""\
                suppressed(COMPLEX) :-
                    complex(COMPLEX),
                    hasActiveDescendant(COMPLEX)."""),
            docs="A complex with any (transitive) active descendant is suppressed.",
        ),
        Rule(
            identifier="preparation:keep_species_no_complex:carrier",
            text="activityCarrier(SPECIES, SPECIES) :- species(SPECIES).",
            docs="Each species is its own activity carrier.",
        ),
        Rule(
            identifier="preparation:keep_species_no_complex:key",
            text=dedent("""\
                activityKey(SPECIES, kept_species(SPECIES)) :-
                    hasActivity(SPECIES, _),
                    not suppressed(SPECIES)."""),
            docs="A non-suppressed species with activity is keyed by `kept_species(SELF)` (subunits of suppressed complexes are reused as-is at top level).",
        ),
    ),
)

_PREPARATION_NO_COMPLEX = RuleGroup(
    identifier="preparation:no_complex",
    profiles=frozenset({"no_complex"}),
    depends_on=frozenset({"activity_base", "topology"}),
    docs="`no-complex` preparation: a complex with any active descendant is suppressed; templated active species are keyed by `new_species_from_template(SELF)` (proteoforms collapse via Python interning); templateless active species keep their identity via `kept_species(SELF)`.",
    rules=(
        Rule(
            identifier="preparation:no_complex:suppressed",
            text=dedent("""\
                suppressed(COMPLEX) :-
                    complex(COMPLEX),
                    hasActiveDescendant(COMPLEX)."""),
            docs="A complex with any (transitive) active descendant is suppressed.",
        ),
        Rule(
            identifier="preparation:no_complex:carrier",
            text="activityCarrier(SPECIES, SPECIES) :- species(SPECIES).",
            docs="Each species is its own activity carrier.",
        ),
        Rule(
            identifier="preparation:no_complex:key_templated",
            text=dedent("""\
                activityKey(SPECIES, new_species_from_template(SPECIES)) :-
                    hasActivity(SPECIES, _),
                    not suppressed(SPECIES),
                    hasTemplate(SPECIES, _)."""),
            docs="A non-suppressed templated active species is keyed by `new_species_from_template(SELF)`.",
        ),
        Rule(
            identifier="preparation:no_complex:key_templateless",
            text=dedent("""\
                activityKey(SPECIES, kept_species(SPECIES)) :-
                    hasActivity(SPECIES, _),
                    not suppressed(SPECIES),
                    not hasSomeTemplate(SPECIES)."""),
            docs="A non-suppressed templateless active species (phenotype, ion, simple molecule, surviving complex, or templateless subunit) is keyed by `kept_species(SELF)`.",
        ),
    ),
)

_PREPARATION_NORMAL = RuleGroup(
    identifier="preparation:normal",
    profiles=frozenset({"normal"}),
    depends_on=frozenset({"activity_base", "topology"}),
    docs="`normal` preparation: every species (subunit or top-level) with activity contributes its own activity (templated → `new_species_from_template(SELF)` so proteoforms collapse via Python interning; templateless → `kept_species(SELF)`); influences stay on the original species (no rerouting). A complex with any active descendant additionally inherits an activity (templateless, so `kept_species(SELF)`).",
    rules=(
        Rule(
            identifier="preparation:normal:complex_inherits_activity",
            text=dedent("""\
                hasActivity(COMPLEX, has_active_descendant) :-
                    complex(COMPLEX),
                    hasActiveDescendant(COMPLEX)."""),
            docs="A complex with any (transitive) active subunit is itself active in the AF view, with reason `has_active_descendant`.",
        ),
        Rule(
            identifier="preparation:normal:carrier",
            text="activityCarrier(SPECIES, SPECIES) :- species(SPECIES).",
            docs="Each species is its own activity carrier (no rerouting).",
        ),
        Rule(
            identifier="preparation:normal:key_templated",
            text=dedent("""\
                activityKey(SPECIES, new_species_from_template(SPECIES)) :-
                    hasActivity(SPECIES, _),
                    hasTemplate(SPECIES, _)."""),
            docs="A templated species with activity is keyed by `new_species_from_template(SELF)`.",
        ),
        Rule(
            identifier="preparation:normal:key_templateless",
            text=dedent("""\
                activityKey(SPECIES, kept_species(SPECIES)) :-
                    hasActivity(SPECIES, _),
                    not hasSomeTemplate(SPECIES)."""),
            docs="A templateless species with activity is keyed by `kept_species(SELF)`.",
        ),
    ),
)

_PATHS_BASE = RuleGroup(
    identifier="paths_base",
    profiles=_NON_CASQ_PROFILES,
    docs="Builds the signed `path(X, Y, SIGN)` relation from PD reactions and modulation arcs (catalyzer/stimulator/trigger/inhibitor → product, signed modulation arcs, and reactant-chained transitivity through reactions).",
    rules=(
        Rule(
            identifier="paths_base:catalyzer_to_product",
            text=dedent("""\
                path(START_SPECIES, END_SPECIES, positive) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER),
                    catalyzer(MODIFIER),
                    hasReferredSpecies(MODIFIER, START_SPECIES),
                    hasProduct(REACTION, PRODUCT),
                    hasReferredSpecies(PRODUCT, END_SPECIES)."""),
            docs="If a species is referred to by a catalyzer of a reaction and another species is referred to by a product of that reaction, then there is a positive path from the first to the second.",
        ),
        Rule(
            identifier="paths_base:physical_stimulator_to_product",
            text=dedent("""\
                path(START_SPECIES, END_SPECIES, positive) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER),
                    physicalStimulator(MODIFIER),
                    hasReferredSpecies(MODIFIER, START_SPECIES),
                    hasProduct(REACTION, PRODUCT),
                    hasReferredSpecies(PRODUCT, END_SPECIES)."""),
            docs="If a species is referred to by a physical stimulator of a reaction and another species is referred to by a product of that reaction, then there is a positive path from the first to the second.",
        ),
        Rule(
            identifier="paths_base:trigger_to_product",
            text=dedent("""\
                path(START_SPECIES, END_SPECIES, positive) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER),
                    trigger(MODIFIER),
                    hasReferredSpecies(MODIFIER, START_SPECIES),
                    hasProduct(REACTION, PRODUCT),
                    hasReferredSpecies(PRODUCT, END_SPECIES)."""),
            docs="If a species is referred to by a trigger of a reaction and another species is referred to by a product of that reaction, then there is a positive path from the first to the second.",
        ),
        Rule(
            identifier="paths_base:catalyzis_modulation",
            text=dedent("""\
                path(START_SPECIES, END_SPECIES, positive) :-
                    catalyzis(MODULATION),
                    hasSource(MODULATION, START_SPECIES),
                    hasTarget(MODULATION, END_SPECIES)."""),
            docs="If a modulation is a catalyzis, then there is a positive path from its source to its target.",
        ),
        Rule(
            identifier="paths_base:positive_influence_modulation",
            text=dedent("""\
                path(START_SPECIES, END_SPECIES, positive) :-
                    positiveInfluence(MODULATION),
                    hasSource(MODULATION, START_SPECIES),
                    hasTarget(MODULATION, END_SPECIES)."""),
            docs="If a modulation is a positiveInfluence, then there is a positive path from its source to its target.",
        ),
        Rule(
            identifier="paths_base:physical_stimulation_modulation",
            text=dedent("""\
                path(START_SPECIES, END_SPECIES, positive) :-
                    physicalStimulation(MODULATION),
                    hasSource(MODULATION, START_SPECIES),
                    hasTarget(MODULATION, END_SPECIES)."""),
            docs="If a modulation is a physicalStimulation, then there is a positive path from its source to its target.",
        ),
        Rule(
            identifier="paths_base:triggering_modulation",
            text=dedent("""\
                path(START_SPECIES, END_SPECIES, positive) :-
                    triggering(MODULATION),
                    hasSource(MODULATION, START_SPECIES),
                    hasTarget(MODULATION, END_SPECIES)."""),
            docs="If a modulation is a triggering, then there is a positive path from its source to its target.",
        ),
        Rule(
            identifier="paths_base:inhibitor_to_product",
            text=dedent("""\
                path(START_SPECIES, END_SPECIES, negative) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER),
                    inhibitor(MODIFIER),
                    hasReferredSpecies(MODIFIER, START_SPECIES),
                    hasProduct(REACTION, PRODUCT),
                    hasReferredSpecies(PRODUCT, END_SPECIES)."""),
            docs="If a species is referred to by an inhibitor of a reaction and another species is referred to by a product of that reaction, then there is a negative path from the first to the second.",
        ),
        Rule(
            identifier="paths_base:inhibition_modulation",
            text=dedent("""\
                path(START_SPECIES, END_SPECIES, negative) :-
                    inhibition(MODULATION),
                    hasSource(MODULATION, START_SPECIES),
                    hasTarget(MODULATION, END_SPECIES)."""),
            docs="If a modulation is an inhibition, then there is a negative path from its source to its target.",
        ),
        Rule(
            identifier="paths_base:negative_influence_modulation",
            text=dedent("""\
                path(START_SPECIES, END_SPECIES, negative) :-
                    negativeInfluence(MODULATION),
                    hasSource(MODULATION, START_SPECIES),
                    hasTarget(MODULATION, END_SPECIES)."""),
            docs="If a modulation is a negativeInfluence, then there is a negative path from its source to its target.",
        ),
        Rule(
            identifier="paths_base:transitive_through_reaction",
            text=dedent("""\
                path(START_SPECIES, END_SPECIES, SIGN) :-
                    path(START_SPECIES, INTERMEDIATE_SPECIES, SIGN),
                    reaction(REACTION),
                    hasReactant(REACTION, REACTANT),
                    hasReferredSpecies(REACTANT, INTERMEDIATE_SPECIES),
                    hasProduct(REACTION, PRODUCT),
                    hasReferredSpecies(PRODUCT, END_SPECIES)."""),
            docs="If there is a path of a given sign from a species to an intermediate species, and the intermediate species is referred to by a reactant of a reaction whose product refers to another species, then there is a path of that same sign from the first species to the second.",
        ),
    ),
)

_PATHS_COMPLEX_TRAVERSAL = RuleGroup(
    identifier="paths_complex_traversal",
    profiles=frozenset({"no_complex", "keep_species_no_complex"}),
    depends_on=frozenset({"paths_base"}),
    docs="Extends paths through complex containment for the no-complex profiles: a path touching a complex is propagated to/from each of its subunits so influences reach the surviving subunit activities.",
    rules=(
        Rule(
            identifier="paths_complex_traversal:into_subunits",
            text=dedent("""\
                path(START_SPECIES, SUBUNIT, SIGN) :-
                    path(START_SPECIES, END_SPECIES, SIGN),
                    complex(END_SPECIES),
                    hasSubunit(END_SPECIES, SUBUNIT)."""),
            docs="If there is a path of a given sign from a species to a complex, then there is a path of that same sign from the species to each subunit of the complex.",
        ),
        Rule(
            identifier="paths_complex_traversal:from_subunits",
            text=dedent("""\
                path(SUBUNIT, END_SPECIES, SIGN) :-
                    path(START_SPECIES, END_SPECIES, SIGN),
                    complex(START_SPECIES),
                    hasSubunit(START_SPECIES, SUBUNIT)."""),
            docs="If there is a path of a given sign from a complex to a species, then there is a path of that same sign from each subunit of the complex to the species.",
        ),
    ),
)

_INFLUENCES_DERIVATION = RuleGroup(
    identifier="influences_derivation",
    profiles=_NON_CASQ_PROFILES,
    docs="Profile-agnostic derivation: emits `new(activity(KEY))` for every activity key, and `new(positivelyInfluences/2)` / `new(negativelyInfluences/2)` from signed paths and from consumption-based reasoning (catalyzer/physicalStimulator/trigger negatively influence consumed reactant; inhibitor positively influences spared reactant).",
    rules=(
        Rule(
            identifier="influences_derivation:activity",
            text="new(activity(KEY)) :- activityKey(_, KEY).",
            docs="Every species with an activity key emits an activity node with that key.",
        ),
        Rule(
            identifier="influences_derivation:positive_path",
            text=dedent("""\
                new(positivelyInfluences(SOURCE_KEY, TARGET_KEY)) :-
                    path(RAW_SOURCE, RAW_TARGET, positive),
                    activityCarrier(RAW_SOURCE, ACTIVITY_SOURCE),
                    activityCarrier(RAW_TARGET, ACTIVITY_TARGET),
                    activityKey(ACTIVITY_SOURCE, SOURCE_KEY),
                    activityKey(ACTIVITY_TARGET, TARGET_KEY)."""),
            docs="A positive path between two raw species emits a positive influence between their activity-key images via their carriers.",
        ),
        Rule(
            identifier="influences_derivation:negative_path",
            text=dedent("""\
                new(negativelyInfluences(SOURCE_KEY, TARGET_KEY)) :-
                    path(RAW_SOURCE, RAW_TARGET, negative),
                    activityCarrier(RAW_SOURCE, ACTIVITY_SOURCE),
                    activityCarrier(RAW_TARGET, ACTIVITY_TARGET),
                    activityKey(ACTIVITY_SOURCE, SOURCE_KEY),
                    activityKey(ACTIVITY_TARGET, TARGET_KEY)."""),
            docs="A negative path between two raw species emits a negative influence between their activity-key images via their carriers.",
        ),
        Rule(
            identifier="influences_derivation:catalyzer_consumes_reactant",
            text=dedent("""\
                new(negativelyInfluences(SOURCE_KEY, TARGET_KEY)) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), catalyzer(MODIFIER),
                    hasReferredSpecies(MODIFIER, RAW_SOURCE),
                    hasReactant(REACTION, REACTANT),
                    hasReferredSpecies(REACTANT, RAW_TARGET),
                    activityCarrier(RAW_SOURCE, ACTIVITY_SOURCE),
                    activityCarrier(RAW_TARGET, ACTIVITY_TARGET),
                    activityKey(ACTIVITY_SOURCE, SOURCE_KEY),
                    activityKey(ACTIVITY_TARGET, TARGET_KEY)."""),
            docs="A catalyzer of a reaction negatively influences each of its reactants (consumption).",
        ),
        Rule(
            identifier="influences_derivation:physical_stimulator_consumes_reactant",
            text=dedent("""\
                new(negativelyInfluences(SOURCE_KEY, TARGET_KEY)) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), physicalStimulator(MODIFIER),
                    hasReferredSpecies(MODIFIER, RAW_SOURCE),
                    hasReactant(REACTION, REACTANT),
                    hasReferredSpecies(REACTANT, RAW_TARGET),
                    activityCarrier(RAW_SOURCE, ACTIVITY_SOURCE),
                    activityCarrier(RAW_TARGET, ACTIVITY_TARGET),
                    activityKey(ACTIVITY_SOURCE, SOURCE_KEY),
                    activityKey(ACTIVITY_TARGET, TARGET_KEY)."""),
            docs="A physical stimulator of a reaction negatively influences each of its reactants (consumption).",
        ),
        Rule(
            identifier="influences_derivation:trigger_consumes_reactant",
            text=dedent("""\
                new(negativelyInfluences(SOURCE_KEY, TARGET_KEY)) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), trigger(MODIFIER),
                    hasReferredSpecies(MODIFIER, RAW_SOURCE),
                    hasReactant(REACTION, REACTANT),
                    hasReferredSpecies(REACTANT, RAW_TARGET),
                    activityCarrier(RAW_SOURCE, ACTIVITY_SOURCE),
                    activityCarrier(RAW_TARGET, ACTIVITY_TARGET),
                    activityKey(ACTIVITY_SOURCE, SOURCE_KEY),
                    activityKey(ACTIVITY_TARGET, TARGET_KEY)."""),
            docs="A trigger of a reaction negatively influences each of its reactants (consumption).",
        ),
        Rule(
            identifier="influences_derivation:inhibitor_spares_reactant",
            text=dedent("""\
                new(positivelyInfluences(SOURCE_KEY, TARGET_KEY)) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), inhibitor(MODIFIER),
                    hasReferredSpecies(MODIFIER, RAW_SOURCE),
                    hasReactant(REACTION, REACTANT),
                    hasReferredSpecies(REACTANT, RAW_TARGET),
                    activityCarrier(RAW_SOURCE, ACTIVITY_SOURCE),
                    activityCarrier(RAW_TARGET, ACTIVITY_TARGET),
                    activityKey(ACTIVITY_SOURCE, SOURCE_KEY),
                    activityKey(ACTIVITY_TARGET, TARGET_KEY)."""),
            docs="An inhibitor of a reaction positively influences each of its reactants (sparing).",
        ),
    ),
)


_CASQ_PARTICIPATION = RuleGroup(
    identifier="casq:participation",
    profiles=_CASQ_PROFILES,
    docs="CASQ helper relations describing how species participate in reactions (`activeParticipates` for reactant/modifier roles, `participates` adding products, `isProducedSpecies` and `isModifierSpecies`), used as conditions of the CASQ deletion rules.",
    rules=(
        Rule(
            identifier="casq:participation:active_from_reactant",
            text=dedent("""\
                activeParticipates(SPECIES, REACTION) :-
                    reaction(REACTION),
                    hasReactant(REACTION, REACTANT),
                    hasReferredSpecies(REACTANT, SPECIES)."""),
            docs="A species actively participates in a reaction if it is referred to by a reactant of that reaction.",
        ),
        Rule(
            identifier="casq:participation:active_from_modifier",
            text=dedent("""\
                activeParticipates(SPECIES, REACTION) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER),
                    hasReferredSpecies(MODIFIER, SPECIES)."""),
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
                    hasReferredSpecies(PRODUCT, SPECIES)."""),
            docs="A species participates in a reaction if it is referred to by a product of that reaction.",
        ),
        Rule(
            identifier="casq:participation:is_produced",
            text=dedent("""\
                isProducedSpecies(SPECIES) :-
                    reaction(REACTION),
                    hasProduct(REACTION, PRODUCT),
                    hasReferredSpecies(PRODUCT, SPECIES)."""),
            docs="A species is produced if some reaction has a product referring to it.",
        ),
        Rule(
            identifier="casq:participation:is_modifier",
            text=dedent("""\
                isModifierSpecies(SPECIES) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER),
                    hasReferredSpecies(MODIFIER, SPECIES)."""),
            docs="A species is a modifier if some reaction has a modifier referring to it.",
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
                    hasReactant(REACTION, RT_RECEPTOR),
                    hasReferredSpecies(RT_RECEPTOR, RECEPTOR),
                    receptor(RECEPTOR),
                    hasReactant(REACTION, RT_PARTNER),
                    hasReferredSpecies(RT_PARTNER, PARTNER),
                    RECEPTOR != PARTNER,
                    #count{ X : hasReactant(REACTION, RT), hasReferredSpecies(RT, X) } = 2,
                    #count{ R1 : participates(RECEPTOR, R1) } = 1,
                    #count{ R2 : participates(PARTNER, R2) } = 1."""),
            docs="A receptor in a 2-reactant heterodimer association where receptor and partner each participate in only this reaction is deleted (rule_1).",
        ),
        Rule(
            identifier="casq:delete:rule_2",
            text=dedent("""\
                delete(SPECIES_1, rule_2) :-
                    heterodimerAssociation(REACTION),
                    hasReactant(REACTION, RT1), hasReferredSpecies(RT1, SPECIES_1),
                    hasReactant(REACTION, RT2), hasReferredSpecies(RT2, SPECIES_2),
                    SPECIES_1 != SPECIES_2,
                    not receptor(SPECIES_1),
                    not receptor(SPECIES_2),
                    #count{ X : hasReactant(REACTION, RT), hasReferredSpecies(RT, X) } = 2,
                    #count{ R1 : activeParticipates(SPECIES_1, R1) } = 1,
                    #count{ R2 : activeParticipates(SPECIES_2, R2) } = 1."""),
            docs="In a 2-reactant heterodimer association where neither reactant is a receptor and each actively participates only in this reaction, both reactant species are deleted (rule_2). The rule fires symmetrically for each side.",
        ),
        Rule(
            identifier="casq:delete:rule_3",
            text=dedent("""\
                delete(REACTANT, rule_3) :-
                    reaction(REACTION),
                    hasReactant(REACTION, RT), hasReferredSpecies(RT, REACTANT),
                    hasProduct(REACTION, P), hasReferredSpecies(P, PRODUCT),
                    REACTANT != PRODUCT,
                    hasName(REACTANT, NAME), hasName(PRODUCT, NAME),
                    not isProducedSpecies(REACTANT),
                    not isModifierSpecies(REACTANT),
                    #count{ X : hasProduct(REACTION, P2), hasReferredSpecies(P2, X) } = 1,
                    #count{ R : hasReactant(R, RT2), hasReferredSpecies(RT2, REACTANT) } = 1."""),
            docs="In a single-product reaction where reactant and product share a name, the reactant is deleted (rule_3) if it is not produced anywhere else, never appears as a modifier, and is consumed only by this reaction.",
        ),
        Rule(
            identifier="casq:delete:rule_4",
            text=dedent("""\
                delete(REACTANT, rule_4) :-
                    transport(REACTION),
                    hasReactant(REACTION, RT), hasReferredSpecies(RT, REACTANT),
                    hasProduct(REACTION, P), hasReferredSpecies(P, PRODUCT),
                    REACTANT != PRODUCT,
                    hasName(REACTANT, NAME), hasName(PRODUCT, NAME),
                    #count{ R : activeParticipates(REACTANT, R) } = 1,
                    #count{ X : hasProduct(REACTION, P2), hasReferredSpecies(P2, X) } = 1."""),
            docs="In a single-product transport where reactant and product share a name and the reactant actively participates only in this reaction, the reactant is deleted (rule_4).",
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
                bridgedProduct(REACTION_1, PRODUCT) :-
                    delete(SPECIES, rule_2),
                    reaction(REACTION_1),
                    hasProduct(REACTION_1, P1), hasReferredSpecies(P1, SPECIES),
                    reaction(REACTION_2),
                    hasReactant(REACTION_2, RT2), hasReferredSpecies(RT2, SPECIES),
                    hasProduct(REACTION_2, P2), hasReferredSpecies(P2, PRODUCT)."""),
            docs="One-hop rewiring across a species deleted by rule_2: REACTION_1 produces SPECIES and REACTION_2 consumes SPECIES and produces PRODUCT, so REACTION_1 is treated as also producing PRODUCT.",
        ),
        Rule(
            identifier="casq:bridged_product:from_rule_4",
            text=dedent("""\
                bridgedProduct(REACTION_1, PRODUCT) :-
                    delete(SPECIES, rule_4),
                    reaction(REACTION_1),
                    hasProduct(REACTION_1, P1), hasReferredSpecies(P1, SPECIES),
                    reaction(REACTION_2),
                    hasReactant(REACTION_2, RT2), hasReferredSpecies(RT2, SPECIES),
                    hasProduct(REACTION_2, P2), hasReferredSpecies(P2, PRODUCT)."""),
            docs="One-hop rewiring across a species deleted by rule_4 (analog of casq:bridged_product:from_rule_2 for transport-driven deletes).",
        ),
    ),
)

_CASQ_ACTIVITY = RuleGroup(
    identifier="casq:activity",
    profiles=_CASQ_PROFILES,
    depends_on=frozenset({"casq:delete"}),
    docs="`casq` activity emission: every PD species not marked for deletion contributes a `kept_species(SELF)` activity in the new AF map.",
    rules=(
        Rule(
            identifier="casq:activity:emit",
            text=dedent("""\
                new(activity(kept_species(SPECIES))) :-
                    species(SPECIES),
                    not delete(SPECIES, _)."""),
            docs="Every PD species that is not deleted emits a `kept_species(SELF)` activity in the AF.",
        ),
        Rule(
            identifier="casq:activity:contributes",
            text=dedent("""\
                contributesActivity(SPECIES) :-
                    species(SPECIES),
                    not delete(SPECIES, _)."""),
            docs="Every non-deleted species contributes an activity, used as a guard in the CASQ influence rules.",
        ),
    ),
)

_CASQ_INFLUENCES = RuleGroup(
    identifier="casq:influences",
    profiles=_CASQ_PROFILES,
    depends_on=frozenset({"casq:bridged_product", "casq:activity"}),
    docs="Casq-specific influence emission: directly wires reactant/catalyzer/stimulator/trigger → product (positive) and inhibitor → product (negative), plus modulation arcs, with a bridged-product variant for each modifier role to route across deleted intermediates.",
    rules=(
        Rule(
            identifier="casq:influences:reactant_to_product",
            text=dedent("""\
                new(positivelyInfluences(kept_species(SOURCE), kept_species(TARGET))) :-
                    reaction(REACTION),
                    hasReactant(REACTION, RT), hasReferredSpecies(RT, SOURCE),
                    hasProduct(REACTION, P), hasReferredSpecies(P, TARGET),
                    contributesActivity(SOURCE), contributesActivity(TARGET)."""),
            docs="A reactant of a reaction positively influences a product of that reaction.",
        ),
        Rule(
            identifier="casq:influences:reactant_to_bridged_product",
            text=dedent("""\
                new(positivelyInfluences(kept_species(SOURCE), kept_species(TARGET))) :-
                    reaction(REACTION),
                    hasReactant(REACTION, RT), hasReferredSpecies(RT, SOURCE),
                    bridgedProduct(REACTION, TARGET),
                    contributesActivity(SOURCE), contributesActivity(TARGET)."""),
            docs="Same as casq:influences:reactant_to_product but routed through a deleted intermediate via bridgedProduct.",
        ),
        Rule(
            identifier="casq:influences:catalyzer_to_product",
            text=dedent("""\
                new(positivelyInfluences(kept_species(SOURCE), kept_species(TARGET))) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), catalyzer(MODIFIER),
                    hasReferredSpecies(MODIFIER, SOURCE),
                    hasProduct(REACTION, P), hasReferredSpecies(P, TARGET),
                    contributesActivity(SOURCE), contributesActivity(TARGET)."""),
            docs="A catalyzer of a reaction positively influences a product of that reaction.",
        ),
        Rule(
            identifier="casq:influences:catalyzer_to_bridged_product",
            text=dedent("""\
                new(positivelyInfluences(kept_species(SOURCE), kept_species(TARGET))) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), catalyzer(MODIFIER),
                    hasReferredSpecies(MODIFIER, SOURCE),
                    bridgedProduct(REACTION, TARGET),
                    contributesActivity(SOURCE), contributesActivity(TARGET)."""),
            docs="Bridged variant of casq:influences:catalyzer_to_product.",
        ),
        Rule(
            identifier="casq:influences:physical_stimulator_to_product",
            text=dedent("""\
                new(positivelyInfluences(kept_species(SOURCE), kept_species(TARGET))) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), physicalStimulator(MODIFIER),
                    hasReferredSpecies(MODIFIER, SOURCE),
                    hasProduct(REACTION, P), hasReferredSpecies(P, TARGET),
                    contributesActivity(SOURCE), contributesActivity(TARGET)."""),
            docs="A physical stimulator of a reaction positively influences a product of that reaction.",
        ),
        Rule(
            identifier="casq:influences:physical_stimulator_to_bridged_product",
            text=dedent("""\
                new(positivelyInfluences(kept_species(SOURCE), kept_species(TARGET))) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), physicalStimulator(MODIFIER),
                    hasReferredSpecies(MODIFIER, SOURCE),
                    bridgedProduct(REACTION, TARGET),
                    contributesActivity(SOURCE), contributesActivity(TARGET)."""),
            docs="Bridged variant of casq:influences:physical_stimulator_to_product.",
        ),
        Rule(
            identifier="casq:influences:trigger_to_product",
            text=dedent("""\
                new(positivelyInfluences(kept_species(SOURCE), kept_species(TARGET))) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), trigger(MODIFIER),
                    hasReferredSpecies(MODIFIER, SOURCE),
                    hasProduct(REACTION, P), hasReferredSpecies(P, TARGET),
                    contributesActivity(SOURCE), contributesActivity(TARGET)."""),
            docs="A trigger of a reaction positively influences a product of that reaction.",
        ),
        Rule(
            identifier="casq:influences:trigger_to_bridged_product",
            text=dedent("""\
                new(positivelyInfluences(kept_species(SOURCE), kept_species(TARGET))) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), trigger(MODIFIER),
                    hasReferredSpecies(MODIFIER, SOURCE),
                    bridgedProduct(REACTION, TARGET),
                    contributesActivity(SOURCE), contributesActivity(TARGET)."""),
            docs="Bridged variant of casq:influences:trigger_to_product.",
        ),
        Rule(
            identifier="casq:influences:inhibitor_to_product",
            text=dedent("""\
                new(negativelyInfluences(kept_species(SOURCE), kept_species(TARGET))) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), inhibitor(MODIFIER),
                    hasReferredSpecies(MODIFIER, SOURCE),
                    hasProduct(REACTION, P), hasReferredSpecies(P, TARGET),
                    contributesActivity(SOURCE), contributesActivity(TARGET)."""),
            docs="An inhibitor of a reaction negatively influences a product of that reaction.",
        ),
        Rule(
            identifier="casq:influences:inhibitor_to_bridged_product",
            text=dedent("""\
                new(negativelyInfluences(kept_species(SOURCE), kept_species(TARGET))) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), inhibitor(MODIFIER),
                    hasReferredSpecies(MODIFIER, SOURCE),
                    bridgedProduct(REACTION, TARGET),
                    contributesActivity(SOURCE), contributesActivity(TARGET)."""),
            docs="Bridged variant of casq:influences:inhibitor_to_product.",
        ),
        Rule(
            identifier="casq:influences:catalyzis_modulation",
            text=dedent("""\
                new(positivelyInfluences(kept_species(SOURCE), kept_species(TARGET))) :-
                    catalyzis(MODULATION),
                    hasSource(MODULATION, SOURCE), hasTarget(MODULATION, TARGET),
                    contributesActivity(SOURCE), contributesActivity(TARGET)."""),
            docs="A catalyzis modulation arc emits a positive influence between its source and target activities.",
        ),
        Rule(
            identifier="casq:influences:positive_influence_modulation",
            text=dedent("""\
                new(positivelyInfluences(kept_species(SOURCE), kept_species(TARGET))) :-
                    positiveInfluence(MODULATION),
                    hasSource(MODULATION, SOURCE), hasTarget(MODULATION, TARGET),
                    contributesActivity(SOURCE), contributesActivity(TARGET)."""),
            docs="A positiveInfluence modulation arc emits a positive influence between its source and target activities.",
        ),
        Rule(
            identifier="casq:influences:physical_stimulation_modulation",
            text=dedent("""\
                new(positivelyInfluences(kept_species(SOURCE), kept_species(TARGET))) :-
                    physicalStimulation(MODULATION),
                    hasSource(MODULATION, SOURCE), hasTarget(MODULATION, TARGET),
                    contributesActivity(SOURCE), contributesActivity(TARGET)."""),
            docs="A physicalStimulation modulation arc emits a positive influence between its source and target activities.",
        ),
        Rule(
            identifier="casq:influences:triggering_modulation",
            text=dedent("""\
                new(positivelyInfluences(kept_species(SOURCE), kept_species(TARGET))) :-
                    triggering(MODULATION),
                    hasSource(MODULATION, SOURCE), hasTarget(MODULATION, TARGET),
                    contributesActivity(SOURCE), contributesActivity(TARGET)."""),
            docs="A triggering modulation arc emits a positive influence between its source and target activities.",
        ),
        Rule(
            identifier="casq:influences:inhibition_modulation",
            text=dedent("""\
                new(negativelyInfluences(kept_species(SOURCE), kept_species(TARGET))) :-
                    inhibition(MODULATION),
                    hasSource(MODULATION, SOURCE), hasTarget(MODULATION, TARGET),
                    contributesActivity(SOURCE), contributesActivity(TARGET)."""),
            docs="An inhibition modulation arc emits a negative influence between its source and target activities.",
        ),
        Rule(
            identifier="casq:influences:negative_influence_modulation",
            text=dedent("""\
                new(negativelyInfluences(kept_species(SOURCE), kept_species(TARGET))) :-
                    negativeInfluence(MODULATION),
                    hasSource(MODULATION, SOURCE), hasTarget(MODULATION, TARGET),
                    contributesActivity(SOURCE), contributesActivity(TARGET)."""),
            docs="A negativeInfluence modulation arc emits a negative influence between its source and target activities.",
        ),
    ),
)


def _build_registry() -> RuleRegistry:
    registry = RuleRegistry()
    registry.register(
        [
            _ACTIVITY_BASE,
            _TOPOLOGY,
            _PREPARATION_KEEP_SPECIES,
            _PREPARATION_KEEP_SPECIES_NO_COMPLEX,
            _PREPARATION_NO_COMPLEX,
            _PREPARATION_NORMAL,
            _PATHS_BASE,
            _PATHS_COMPLEX_TRAVERSAL,
            _INFLUENCES_DERIVATION,
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


def build_program(profile: str) -> str:
    """Return the composed ASP program text for the given profile."""
    registry = _build_registry()
    plan = CollectionPlan(registry)
    plan.add_profile(profile)
    return "\n".join(rule.text for rule in plan.resolve())
