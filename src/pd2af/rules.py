"""ASP rule composition for the pd2af transformation modes.

The rules are organized as aspcompose groups. Five profiles
(`normal`, `no_complex`, `keep_species`, `keep_species_no_complex`,
`casq`) compose the groups into the supported mode programs. The
first four profiles cover a 2x2 matrix on two orthogonal axes:

* species treatment — ``normal``/``no_complex`` merge proteoforms of
  the same template (and compartment) into a single activity (true
  PD->AF transform, the only style expressible in SBGN PD);
  ``keep_species``/``keep_species_no_complex`` keep each PD species
  as its own activity (CellDesigner-only).
* complex treatment — ``no_complex``/``keep_species_no_complex`` drop
  any complex that has an active subunit and route influences through
  its subunits; ``normal``/``keep_species`` keep complexes as their
  own activities and route influences through them.

The ``casq`` profile keeps each PD species as its own activity (like
``keep_species``), prunes species via CASQ-style deletion rules
(rule_1..rule_4) and emits influences directly from reaction
modifier/reactant -> product and from modulation arcs, with
single-hop rewiring across deleted intermediates.

All path rules use a single ``path(X, Y, SIGN)`` predicate so
transitivity and complex-subunit traversal can be written once and
parameterized over ``SIGN``.

Activity identity is split into two profile-specific relations so the
shared activity-derivation and influence rules remain
profile-agnostic:

* ``contributesActivity(SPECIES)`` — survives the profile filter and
  contributes an activity to the AF.
* ``activityKey(SPECIES, KEY)`` — this species's identity in the AF
  (a ``kept_species/1`` term in keep_species/keep_species_no_complex
  and casq; a ``kept_species/1`` or ``derived_proteoform_class/2``
  term in normal/no_complex).
"""

from textwrap import dedent

from aspcompose import CollectionPlan, Rule, RuleGroup, RuleRegistry

_ALL_PROFILES = frozenset(
    {"normal", "no_complex", "keep_species", "keep_species_no_complex"}
)

_CASQ_PROFILES = frozenset({"casq"})

_ACTIVITY_BASE = RuleGroup(
    identifier="activity_base",
    profiles=_ALL_PROFILES | _CASQ_PROFILES,
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

_HAS_ACTIVE_SUBUNIT = RuleGroup(
    identifier="has_active_subunit",
    profiles=frozenset({"no_complex", "keep_species_no_complex"}),
    depends_on=frozenset({"activity_base"}),
    docs="Flags a complex as having an active subunit when both the complex and one of its direct subunits have activity. Used by the no-complex profiles to drop such complexes in favor of their subunits.",
    rules=(
        Rule(
            identifier="has_active_subunit:detect",
            text=dedent("""\
                hasActiveSubunit(SPECIES) :-
                    complex(SPECIES),
                    hasActivity(SPECIES, _),
                    hasSubunit(SPECIES, SUBUNIT),
                    hasActivity(SUBUNIT, _)."""),
            docs="If a complex has activity and contains a subunit that also has activity, then it has an active subunit.",
        ),
    ),
)

_HAS_CONTRIBUTING_COMPLEX_ANCESTOR = RuleGroup(
    identifier="has_contributing_complex_ancestor",
    profiles=frozenset({"normal", "keep_species"}),
    depends_on=frozenset({"activity_base"}),
    docs="Computes the transitive `hasContributingComplexAncestor` relation, marking every (nested) subunit of an active complex so the keep-complex profiles can suppress those subunits in favor of the complex.",
    rules=(
        Rule(
            identifier="has_contributing_complex_ancestor:direct",
            text=dedent("""\
                hasContributingComplexAncestor(SUBUNIT) :-
                    complex(COMPLEX),
                    hasActivity(COMPLEX, _),
                    hasSubunit(COMPLEX, SUBUNIT)."""),
            docs="If a complex has activity, then each of its direct subunits has a contributing complex ancestor.",
        ),
        Rule(
            identifier="has_contributing_complex_ancestor:transitive",
            text=dedent("""\
                hasContributingComplexAncestor(SUBUNIT) :-
                    hasContributingComplexAncestor(COMPLEX),
                    hasSubunit(COMPLEX, SUBUNIT)."""),
            docs="The contributing-complex-ancestor relation is transitive through complex containment, so deeply nested subunits are also covered.",
        ),
    ),
)

_CONTRIBUTES_ACTIVITY_FLAT = RuleGroup(
    identifier="contributes_activity:flat",
    slot="contributes_activity",
    profiles=frozenset({"normal", "keep_species"}),
    depends_on=frozenset({"activity_base", "has_contributing_complex_ancestor"}),
    docs="`contributes_activity` slot for the keep-complex profiles: a species contributes an activity unless it sits inside a complex that itself contributes one (the complex stands in for its proteoform).",
    rules=(
        Rule(
            identifier="contributes_activity:flat:promote_unless_inside_contributing_complex",
            text=dedent("""\
                contributesActivity(SPECIES) :-
                    hasActivity(SPECIES, _),
                    not hasContributingComplexAncestor(SPECIES)."""),
            docs="A species with activity contributes an activity to the new AF unless it is a subunit (direct or nested) of a complex that itself contributes activity — the complex represents the active proteoform.",
        ),
    ),
)

_CONTRIBUTES_ACTIVITY_NO_ACTIVE_SUBUNITS = RuleGroup(
    identifier="contributes_activity:no_active_subunits",
    slot="contributes_activity",
    profiles=frozenset({"no_complex", "keep_species_no_complex"}),
    depends_on=frozenset({"activity_base", "has_active_subunit"}),
    docs="`contributes_activity` slot for the no-complex profiles: a species contributes an activity unless it is a complex that has an active subunit, in which case the subunit takes over.",
    rules=(
        Rule(
            identifier="contributes_activity:no_active_subunits:promote_unless_superseded",
            text=dedent("""\
                contributesActivity(SPECIES) :-
                    hasActivity(SPECIES, _),
                    not hasActiveSubunit(SPECIES)."""),
            docs="A species with activity contributes an activity to the new AF unless it is a complex that has an active subunit.",
        ),
    ),
)

_ACTIVITY_KEY_KEPT = RuleGroup(
    identifier="activity_key:kept",
    slot="activity_key",
    profiles=frozenset({"keep_species", "keep_species_no_complex", "casq"}),
    depends_on=frozenset({"activity_base"}),
    docs="`activity_key` slot for the keep-species and casq profiles: each PD species is its own AF identity (`kept_species(SPECIES)`), so no proteoform merging happens.",
    rules=(
        Rule(
            identifier="activity_key:kept:identity",
            text="activityKey(SPECIES, kept_species(SPECIES)) :- species(SPECIES).",
            docs="In keep_species, keep_species_no_complex, and casq modes, every species's activity key is itself.",
        ),
    ),
)

_ACTIVITY_KEY_MERGED = RuleGroup(
    identifier="activity_key:merged",
    slot="activity_key",
    profiles=frozenset({"normal", "no_complex"}),
    depends_on=frozenset({"activity_base"}),
    docs="`activity_key` slot for the merge profiles: templated species are keyed by `(template, compartment)` so proteoforms collapse into a single AF activity, while templateless species keep their own identity.",
    rules=(
        Rule(
            identifier="activity_key:merged:has_some_template",
            text="hasSomeTemplate(SPECIES) :- hasTemplate(SPECIES, _).",
            docs="A species has some template if it is linked to any template.",
        ),
        Rule(
            identifier="activity_key:merged:has_some_compartment",
            text="hasSomeCompartment(SPECIES) :- hasCompartment(SPECIES, _).",
            docs="A species has some compartment if it is linked to any compartment.",
        ),
        Rule(
            identifier="activity_key:merged:templated_with_compartment",
            text=dedent("""\
                activityKey(SPECIES, derived_proteoform_class(TEMPLATE, COMPARTMENT)) :-
                    hasTemplate(SPECIES, TEMPLATE),
                    hasCompartment(SPECIES, COMPARTMENT)."""),
            docs="A templated species in a compartment has a derived-proteoform-class key keyed by template and compartment.",
        ),
        Rule(
            identifier="activity_key:merged:templated_without_compartment",
            text=dedent("""\
                activityKey(SPECIES, derived_proteoform_class(TEMPLATE, no_compartment)) :-
                    hasTemplate(SPECIES, TEMPLATE),
                    not hasSomeCompartment(SPECIES)."""),
            docs="A templated species without a compartment has a derived-proteoform-class key with the sentinel `no_compartment`.",
        ),
        Rule(
            identifier="activity_key:merged:templateless",
            text=dedent("""\
                activityKey(SPECIES, kept_species(SPECIES)) :-
                    species(SPECIES),
                    not hasSomeTemplate(SPECIES)."""),
            docs="A templateless species (phenotype, ion, simple molecule, complex, etc.) keeps its own key when proteoforms are merged.",
        ),
    ),
)

_ACTIVITY_DERIVATION = RuleGroup(
    identifier="activity_derivation",
    profiles=_ALL_PROFILES | _CASQ_PROFILES,
    depends_on=frozenset({"contributes_activity", "activity_key"}),
    docs="Profile-agnostic step that emits a `new(activity(KEY))` fact for every species that contributes an activity, using whichever `activityKey` the active profile supplies.",
    rules=(
        Rule(
            identifier="activity_derivation:emit",
            text=dedent("""\
                new(activity(KEY)) :-
                    contributesActivity(SPECIES),
                    activityKey(SPECIES, KEY)."""),
            docs="If a species contributes an activity and has an activity key, then a new activity node with that key is emitted.",
        ),
    ),
)

_PATHS_BASE = RuleGroup(
    identifier="paths_base",
    profiles=_ALL_PROFILES,
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

_INFLUENCES_FROM_PATHS = RuleGroup(
    identifier="influences_from_paths",
    profiles=_ALL_PROFILES,
    depends_on=frozenset({"paths_base", "contributes_activity", "activity_key"}),
    docs="Emits `new(positivelyInfluences/2)` and `new(negativelyInfluences/2)` between activity keys whenever a signed path connects two species that both contribute activities.",
    rules=(
        Rule(
            identifier="influences_from_paths:positive",
            text=dedent("""\
                new(positivelyInfluences(SRC_KEY, TGT_KEY)) :-
                    path(SOURCE, TARGET, positive),
                    contributesActivity(SOURCE),
                    contributesActivity(TARGET),
                    activityKey(SOURCE, SRC_KEY),
                    activityKey(TARGET, TGT_KEY)."""),
            docs="If there is a positive path between two species and both contribute activities, then a positive influence between their activity keys is emitted.",
        ),
        Rule(
            identifier="influences_from_paths:negative",
            text=dedent("""\
                new(negativelyInfluences(SRC_KEY, TGT_KEY)) :-
                    path(SOURCE, TARGET, negative),
                    contributesActivity(SOURCE),
                    contributesActivity(TARGET),
                    activityKey(SOURCE, SRC_KEY),
                    activityKey(TARGET, TGT_KEY)."""),
            docs="If there is a negative path between two species and both contribute activities, then a negative influence between their activity keys is emitted.",
        ),
    ),
)

_INFLUENCES_CONSUMPTION = RuleGroup(
    identifier="influences_consumption",
    profiles=_ALL_PROFILES,
    depends_on=frozenset({"contributes_activity", "activity_key"}),
    docs="Emits influences between modifier and reactant of the same reaction: catalyzer/stimulator/trigger negatively influence the reactant they consume, while inhibitor positively influences the reactant it spares.",
    rules=(
        Rule(
            identifier="influences_consumption:catalyzer_consumes_reactant",
            text=dedent("""\
                new(negativelyInfluences(SRC_KEY, TGT_KEY)) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER),
                    catalyzer(MODIFIER),
                    hasReferredSpecies(MODIFIER, SOURCE),
                    hasReactant(REACTION, REACTANT),
                    hasReferredSpecies(REACTANT, TARGET),
                    contributesActivity(SOURCE),
                    contributesActivity(TARGET),
                    activityKey(SOURCE, SRC_KEY),
                    activityKey(TARGET, TGT_KEY)."""),
            docs="If a species is referred to by a catalyzer of a reaction and another species is referred to by a reactant of that reaction, and both contribute activities, then a negative influence between their activity keys is emitted.",
        ),
        Rule(
            identifier="influences_consumption:physical_stimulator_consumes_reactant",
            text=dedent("""\
                new(negativelyInfluences(SRC_KEY, TGT_KEY)) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER),
                    physicalStimulator(MODIFIER),
                    hasReferredSpecies(MODIFIER, SOURCE),
                    hasReactant(REACTION, REACTANT),
                    hasReferredSpecies(REACTANT, TARGET),
                    contributesActivity(SOURCE),
                    contributesActivity(TARGET),
                    activityKey(SOURCE, SRC_KEY),
                    activityKey(TARGET, TGT_KEY)."""),
            docs="If a species is referred to by a physical stimulator of a reaction and another species is referred to by a reactant of that reaction, and both contribute activities, then a negative influence between their activity keys is emitted.",
        ),
        Rule(
            identifier="influences_consumption:trigger_consumes_reactant",
            text=dedent("""\
                new(negativelyInfluences(SRC_KEY, TGT_KEY)) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER),
                    trigger(MODIFIER),
                    hasReferredSpecies(MODIFIER, SOURCE),
                    hasReactant(REACTION, REACTANT),
                    hasReferredSpecies(REACTANT, TARGET),
                    contributesActivity(SOURCE),
                    contributesActivity(TARGET),
                    activityKey(SOURCE, SRC_KEY),
                    activityKey(TARGET, TGT_KEY)."""),
            docs="If a species is referred to by a trigger of a reaction and another species is referred to by a reactant of that reaction, and both contribute activities, then a negative influence between their activity keys is emitted.",
        ),
        Rule(
            identifier="influences_consumption:inhibitor_spares_reactant",
            text=dedent("""\
                new(positivelyInfluences(SRC_KEY, TGT_KEY)) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER),
                    inhibitor(MODIFIER),
                    hasReferredSpecies(MODIFIER, SOURCE),
                    hasReactant(REACTION, REACTANT),
                    hasReferredSpecies(REACTANT, TARGET),
                    contributesActivity(SOURCE),
                    contributesActivity(TARGET),
                    activityKey(SOURCE, SRC_KEY),
                    activityKey(TARGET, TGT_KEY)."""),
            docs="If a species is referred to by an inhibitor of a reaction and another species is referred to by a reactant of that reaction, and both contribute activities, then a positive influence between their activity keys is emitted.",
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

_CASQ_CONTRIBUTES_ACTIVITY = RuleGroup(
    identifier="contributes_activity:casq",
    slot="contributes_activity",
    profiles=_CASQ_PROFILES,
    depends_on=frozenset({"casq:delete"}),
    docs="`contributes_activity` slot for the casq profile: every PD species not marked for deletion by any CASQ rule contributes an activity to the AF.",
    rules=(
        Rule(
            identifier="contributes_activity:casq:not_deleted",
            text=dedent("""\
                contributesActivity(SPECIES) :-
                    species(SPECIES),
                    not delete(SPECIES, _)."""),
            docs="Every PD species that is not deleted by any CASQ rule contributes an activity to the new AF.",
        ),
    ),
)

_CASQ_INFLUENCES = RuleGroup(
    identifier="casq:influences",
    profiles=_CASQ_PROFILES,
    depends_on=frozenset(
        {"casq:bridged_product", "contributes_activity", "activity_key"}
    ),
    docs="Casq-specific influence emission: directly wires reactant/catalyzer/stimulator/trigger → product (positive) and inhibitor → product (negative), plus modulation arcs, with a bridged-product variant for each modifier role to route across deleted intermediates.",
    rules=(
        Rule(
            identifier="casq:influences:reactant_to_product",
            text=dedent("""\
                new(positivelyInfluences(SRC_KEY, TGT_KEY)) :-
                    reaction(REACTION),
                    hasReactant(REACTION, RT), hasReferredSpecies(RT, SOURCE),
                    hasProduct(REACTION, P), hasReferredSpecies(P, TARGET),
                    contributesActivity(SOURCE), contributesActivity(TARGET),
                    activityKey(SOURCE, SRC_KEY), activityKey(TARGET, TGT_KEY)."""),
            docs="A reactant of a reaction positively influences a product of that reaction.",
        ),
        Rule(
            identifier="casq:influences:reactant_to_bridged_product",
            text=dedent("""\
                new(positivelyInfluences(SRC_KEY, TGT_KEY)) :-
                    reaction(REACTION),
                    hasReactant(REACTION, RT), hasReferredSpecies(RT, SOURCE),
                    bridgedProduct(REACTION, TARGET),
                    contributesActivity(SOURCE), contributesActivity(TARGET),
                    activityKey(SOURCE, SRC_KEY), activityKey(TARGET, TGT_KEY)."""),
            docs="Same as casq:influences:reactant_to_product but routed through a deleted intermediate via bridgedProduct.",
        ),
        Rule(
            identifier="casq:influences:catalyzer_to_product",
            text=dedent("""\
                new(positivelyInfluences(SRC_KEY, TGT_KEY)) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), catalyzer(MODIFIER),
                    hasReferredSpecies(MODIFIER, SOURCE),
                    hasProduct(REACTION, P), hasReferredSpecies(P, TARGET),
                    contributesActivity(SOURCE), contributesActivity(TARGET),
                    activityKey(SOURCE, SRC_KEY), activityKey(TARGET, TGT_KEY)."""),
            docs="A catalyzer of a reaction positively influences a product of that reaction.",
        ),
        Rule(
            identifier="casq:influences:catalyzer_to_bridged_product",
            text=dedent("""\
                new(positivelyInfluences(SRC_KEY, TGT_KEY)) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), catalyzer(MODIFIER),
                    hasReferredSpecies(MODIFIER, SOURCE),
                    bridgedProduct(REACTION, TARGET),
                    contributesActivity(SOURCE), contributesActivity(TARGET),
                    activityKey(SOURCE, SRC_KEY), activityKey(TARGET, TGT_KEY)."""),
            docs="Bridged variant of casq:influences:catalyzer_to_product.",
        ),
        Rule(
            identifier="casq:influences:physical_stimulator_to_product",
            text=dedent("""\
                new(positivelyInfluences(SRC_KEY, TGT_KEY)) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), physicalStimulator(MODIFIER),
                    hasReferredSpecies(MODIFIER, SOURCE),
                    hasProduct(REACTION, P), hasReferredSpecies(P, TARGET),
                    contributesActivity(SOURCE), contributesActivity(TARGET),
                    activityKey(SOURCE, SRC_KEY), activityKey(TARGET, TGT_KEY)."""),
            docs="A physical stimulator of a reaction positively influences a product of that reaction.",
        ),
        Rule(
            identifier="casq:influences:physical_stimulator_to_bridged_product",
            text=dedent("""\
                new(positivelyInfluences(SRC_KEY, TGT_KEY)) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), physicalStimulator(MODIFIER),
                    hasReferredSpecies(MODIFIER, SOURCE),
                    bridgedProduct(REACTION, TARGET),
                    contributesActivity(SOURCE), contributesActivity(TARGET),
                    activityKey(SOURCE, SRC_KEY), activityKey(TARGET, TGT_KEY)."""),
            docs="Bridged variant of casq:influences:physical_stimulator_to_product.",
        ),
        Rule(
            identifier="casq:influences:trigger_to_product",
            text=dedent("""\
                new(positivelyInfluences(SRC_KEY, TGT_KEY)) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), trigger(MODIFIER),
                    hasReferredSpecies(MODIFIER, SOURCE),
                    hasProduct(REACTION, P), hasReferredSpecies(P, TARGET),
                    contributesActivity(SOURCE), contributesActivity(TARGET),
                    activityKey(SOURCE, SRC_KEY), activityKey(TARGET, TGT_KEY)."""),
            docs="A trigger of a reaction positively influences a product of that reaction.",
        ),
        Rule(
            identifier="casq:influences:trigger_to_bridged_product",
            text=dedent("""\
                new(positivelyInfluences(SRC_KEY, TGT_KEY)) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), trigger(MODIFIER),
                    hasReferredSpecies(MODIFIER, SOURCE),
                    bridgedProduct(REACTION, TARGET),
                    contributesActivity(SOURCE), contributesActivity(TARGET),
                    activityKey(SOURCE, SRC_KEY), activityKey(TARGET, TGT_KEY)."""),
            docs="Bridged variant of casq:influences:trigger_to_product.",
        ),
        Rule(
            identifier="casq:influences:inhibitor_to_product",
            text=dedent("""\
                new(negativelyInfluences(SRC_KEY, TGT_KEY)) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), inhibitor(MODIFIER),
                    hasReferredSpecies(MODIFIER, SOURCE),
                    hasProduct(REACTION, P), hasReferredSpecies(P, TARGET),
                    contributesActivity(SOURCE), contributesActivity(TARGET),
                    activityKey(SOURCE, SRC_KEY), activityKey(TARGET, TGT_KEY)."""),
            docs="An inhibitor of a reaction negatively influences a product of that reaction.",
        ),
        Rule(
            identifier="casq:influences:inhibitor_to_bridged_product",
            text=dedent("""\
                new(negativelyInfluences(SRC_KEY, TGT_KEY)) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER), inhibitor(MODIFIER),
                    hasReferredSpecies(MODIFIER, SOURCE),
                    bridgedProduct(REACTION, TARGET),
                    contributesActivity(SOURCE), contributesActivity(TARGET),
                    activityKey(SOURCE, SRC_KEY), activityKey(TARGET, TGT_KEY)."""),
            docs="Bridged variant of casq:influences:inhibitor_to_product.",
        ),
        Rule(
            identifier="casq:influences:catalyzis_modulation",
            text=dedent("""\
                new(positivelyInfluences(SRC_KEY, TGT_KEY)) :-
                    catalyzis(MODULATION),
                    hasSource(MODULATION, SOURCE), hasTarget(MODULATION, TARGET),
                    contributesActivity(SOURCE), contributesActivity(TARGET),
                    activityKey(SOURCE, SRC_KEY), activityKey(TARGET, TGT_KEY)."""),
            docs="A catalyzis modulation arc emits a positive influence between its source and target activities.",
        ),
        Rule(
            identifier="casq:influences:positive_influence_modulation",
            text=dedent("""\
                new(positivelyInfluences(SRC_KEY, TGT_KEY)) :-
                    positiveInfluence(MODULATION),
                    hasSource(MODULATION, SOURCE), hasTarget(MODULATION, TARGET),
                    contributesActivity(SOURCE), contributesActivity(TARGET),
                    activityKey(SOURCE, SRC_KEY), activityKey(TARGET, TGT_KEY)."""),
            docs="A positiveInfluence modulation arc emits a positive influence between its source and target activities.",
        ),
        Rule(
            identifier="casq:influences:physical_stimulation_modulation",
            text=dedent("""\
                new(positivelyInfluences(SRC_KEY, TGT_KEY)) :-
                    physicalStimulation(MODULATION),
                    hasSource(MODULATION, SOURCE), hasTarget(MODULATION, TARGET),
                    contributesActivity(SOURCE), contributesActivity(TARGET),
                    activityKey(SOURCE, SRC_KEY), activityKey(TARGET, TGT_KEY)."""),
            docs="A physicalStimulation modulation arc emits a positive influence between its source and target activities.",
        ),
        Rule(
            identifier="casq:influences:triggering_modulation",
            text=dedent("""\
                new(positivelyInfluences(SRC_KEY, TGT_KEY)) :-
                    triggering(MODULATION),
                    hasSource(MODULATION, SOURCE), hasTarget(MODULATION, TARGET),
                    contributesActivity(SOURCE), contributesActivity(TARGET),
                    activityKey(SOURCE, SRC_KEY), activityKey(TARGET, TGT_KEY)."""),
            docs="A triggering modulation arc emits a positive influence between its source and target activities.",
        ),
        Rule(
            identifier="casq:influences:inhibition_modulation",
            text=dedent("""\
                new(negativelyInfluences(SRC_KEY, TGT_KEY)) :-
                    inhibition(MODULATION),
                    hasSource(MODULATION, SOURCE), hasTarget(MODULATION, TARGET),
                    contributesActivity(SOURCE), contributesActivity(TARGET),
                    activityKey(SOURCE, SRC_KEY), activityKey(TARGET, TGT_KEY)."""),
            docs="An inhibition modulation arc emits a negative influence between its source and target activities.",
        ),
        Rule(
            identifier="casq:influences:negative_influence_modulation",
            text=dedent("""\
                new(negativelyInfluences(SRC_KEY, TGT_KEY)) :-
                    negativeInfluence(MODULATION),
                    hasSource(MODULATION, SOURCE), hasTarget(MODULATION, TARGET),
                    contributesActivity(SOURCE), contributesActivity(TARGET),
                    activityKey(SOURCE, SRC_KEY), activityKey(TARGET, TGT_KEY)."""),
            docs="A negativeInfluence modulation arc emits a negative influence between its source and target activities.",
        ),
    ),
)


def _build_registry() -> RuleRegistry:
    registry = RuleRegistry()
    registry.register(
        [
            _ACTIVITY_BASE,
            _HAS_ACTIVE_SUBUNIT,
            _HAS_CONTRIBUTING_COMPLEX_ANCESTOR,
            _CONTRIBUTES_ACTIVITY_FLAT,
            _CONTRIBUTES_ACTIVITY_NO_ACTIVE_SUBUNITS,
            _ACTIVITY_KEY_KEPT,
            _ACTIVITY_KEY_MERGED,
            _ACTIVITY_DERIVATION,
            _PATHS_BASE,
            _PATHS_COMPLEX_TRAVERSAL,
            _INFLUENCES_FROM_PATHS,
            _INFLUENCES_CONSUMPTION,
            _CASQ_PARTICIPATION,
            _CASQ_DELETE,
            _CASQ_BRIDGED_PRODUCT,
            _CASQ_CONTRIBUTES_ACTIVITY,
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
