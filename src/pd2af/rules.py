"""ASP rule composition for the pd2af and pd2af-no-complex modes.

The rules are organized as aspcompose groups. Two profiles (`default` and
`no_complex`) compose the groups into the two supported mode programs. All
path rules use a single `path(X, Y, SIGN)` predicate so transitivity and
complex-subunit traversal can be written once and parameterized over `SIGN`.
"""

from textwrap import dedent

from aspcompose import CollectionPlan, Rule, RuleGroup, RuleRegistry

_ACTIVITY_BASE = RuleGroup(
    identifier="activity_base",
    profiles=frozenset({"default", "no_complex"}),
    rules=(
        Rule(
            identifier="activity_base:from_active_flag",
            text="hasActivity(SPECIES, active) :- species(SPECIES), hasActive(SPECIES, 1).",
            documentation="If a species has its `hasActive` flag set to 1, then it has activity, with reason `active`.",
        ),
        Rule(
            identifier="activity_base:from_active_structural_state",
            text=dedent("""\
                hasActivity(SPECIES, structural_state_active) :-
                    species(SPECIES),
                    hasStructuralState(SPECIES, STRUCTURAL_STATE),
                    hasValue(STRUCTURAL_STATE, "active")."""),
            documentation='If a species carries a structural state whose value is "active", then it has activity, with reason `structural_state_active`.',
        ),
        Rule(
            identifier="activity_base:from_phenotype",
            text="hasActivity(PHENOTYPE, phenotype) :- phenotype(PHENOTYPE).",
            documentation="If a species is a phenotype, then it has activity, with reason `phenotype`.",
        ),
        Rule(
            identifier="activity_base:from_modulation_source",
            text=dedent("""\
                hasActivity(SOURCE, modulates(SOURCE, TARGET)) :-
                    species(SOURCE),
                    modulation(MODULATION),
                    hasSource(MODULATION, SOURCE),
                    hasTarget(MODULATION, TARGET)."""),
            documentation="If a species is the source of a modulation with some target, then it has activity, with reason `modulates(source, target)`.",
        ),
        Rule(
            identifier="activity_base:from_reaction_modulator",
            text=dedent("""\
                hasActivity(SOURCE, modulates(SOURCE, TARGET)) :-
                    species(SOURCE),
                    modulator(MODULATOR),
                    hasReferredSpecies(MODULATOR, SOURCE),
                    hasModifier(TARGET, MODULATOR)."""),
            documentation="If a species is referred to by a modulator that is a modifier of some target, then it has activity, with reason `modulates(source, target)`.",
        ),
    ),
)

_ACTIVITY_DERIVATION_FLAT = RuleGroup(
    identifier="activity_derivation:flat",
    slot="activity_derivation",
    profiles=frozenset({"default"}),
    depends_on=frozenset({"activity_base"}),
    rules=(
        Rule(
            identifier="activity_derivation:flat:promote_all",
            text="new(activity(SPECIES)) :- hasActivity(SPECIES, _).",
            documentation="If a species has activity, then it is a new activity node.",
        ),
    ),
)

_ACTIVITY_DERIVATION_NO_ACTIVE_SUBUNITS = RuleGroup(
    identifier="activity_derivation:no_active_subunits",
    slot="activity_derivation",
    profiles=frozenset({"no_complex"}),
    depends_on=frozenset({"activity_base"}),
    rules=(
        Rule(
            identifier="activity_derivation:no_active_subunits:detect_active_subunit",
            text=dedent("""\
                hasActiveSubunit(SPECIES) :-
                    complex(SPECIES),
                    hasActivity(SPECIES, _),
                    hasSubunit(SPECIES, SUBUNIT),
                    hasActivity(SUBUNIT, _)."""),
            documentation="If a complex has activity and contains a subunit that also has activity, then it has an active subunit.",
        ),
        Rule(
            identifier="activity_derivation:no_active_subunits:promote_unless_superseded",
            text=dedent("""\
                new(activity(SPECIES)) :-
                    hasActivity(SPECIES, _),
                    not hasActiveSubunit(SPECIES)."""),
            documentation="If a species has activity and does not have an active subunit, then it is a new activity node.",
        ),
    ),
)

_PATHS_BASE = RuleGroup(
    identifier="paths_base",
    profiles=frozenset({"default", "no_complex"}),
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
            documentation="If a species is referred to by a catalyzer of a reaction and another species is referred to by a product of that reaction, then there is a positive path from the first to the second.",
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
            documentation="If a species is referred to by a physical stimulator of a reaction and another species is referred to by a product of that reaction, then there is a positive path from the first to the second.",
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
            documentation="If a species is referred to by a trigger of a reaction and another species is referred to by a product of that reaction, then there is a positive path from the first to the second.",
        ),
        Rule(
            identifier="paths_base:catalyzis_modulation",
            text=dedent("""\
                path(START_SPECIES, END_SPECIES, positive) :-
                    catalyzis(MODULATION),
                    hasSource(MODULATION, START_SPECIES),
                    hasTarget(MODULATION, END_SPECIES)."""),
            documentation="If a modulation is a catalyzis, then there is a positive path from its source to its target.",
        ),
        Rule(
            identifier="paths_base:positive_influence_modulation",
            text=dedent("""\
                path(START_SPECIES, END_SPECIES, positive) :-
                    positiveInfluence(MODULATION),
                    hasSource(MODULATION, START_SPECIES),
                    hasTarget(MODULATION, END_SPECIES)."""),
            documentation="If a modulation is a positiveInfluence, then there is a positive path from its source to its target.",
        ),
        Rule(
            identifier="paths_base:physical_stimulation_modulation",
            text=dedent("""\
                path(START_SPECIES, END_SPECIES, positive) :-
                    physicalStimulation(MODULATION),
                    hasSource(MODULATION, START_SPECIES),
                    hasTarget(MODULATION, END_SPECIES)."""),
            documentation="If a modulation is a physicalStimulation, then there is a positive path from its source to its target.",
        ),
        Rule(
            identifier="paths_base:triggering_modulation",
            text=dedent("""\
                path(START_SPECIES, END_SPECIES, positive) :-
                    triggering(MODULATION),
                    hasSource(MODULATION, START_SPECIES),
                    hasTarget(MODULATION, END_SPECIES)."""),
            documentation="If a modulation is a triggering, then there is a positive path from its source to its target.",
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
            documentation="If a species is referred to by an inhibitor of a reaction and another species is referred to by a product of that reaction, then there is a negative path from the first to the second.",
        ),
        Rule(
            identifier="paths_base:inhibition_modulation",
            text=dedent("""\
                path(START_SPECIES, END_SPECIES, negative) :-
                    inhibition(MODULATION),
                    hasSource(MODULATION, START_SPECIES),
                    hasTarget(MODULATION, END_SPECIES)."""),
            documentation="If a modulation is an inhibition, then there is a negative path from its source to its target.",
        ),
        Rule(
            identifier="paths_base:negative_influence_modulation",
            text=dedent("""\
                path(START_SPECIES, END_SPECIES, negative) :-
                    negativeInfluence(MODULATION),
                    hasSource(MODULATION, START_SPECIES),
                    hasTarget(MODULATION, END_SPECIES)."""),
            documentation="If a modulation is a negativeInfluence, then there is a negative path from its source to its target.",
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
            documentation="If there is a path of a given sign from a species to an intermediate species, and the intermediate species is referred to by a reactant of a reaction whose product refers to another species, then there is a path of that same sign from the first species to the second.",
        ),
    ),
)

_PATHS_COMPLEX_TRAVERSAL = RuleGroup(
    identifier="paths_complex_traversal",
    profiles=frozenset({"no_complex"}),
    depends_on=frozenset({"paths_base"}),
    rules=(
        Rule(
            identifier="paths_complex_traversal:into_subunits",
            text=dedent("""\
                path(START_SPECIES, SUBUNIT, SIGN) :-
                    path(START_SPECIES, END_SPECIES, SIGN),
                    complex(END_SPECIES),
                    hasSubunit(END_SPECIES, SUBUNIT)."""),
            documentation="If there is a path of a given sign from a species to a complex, then there is a path of that same sign from the species to each subunit of the complex.",
        ),
        Rule(
            identifier="paths_complex_traversal:from_subunits",
            text=dedent("""\
                path(SUBUNIT, END_SPECIES, SIGN) :-
                    path(START_SPECIES, END_SPECIES, SIGN),
                    complex(START_SPECIES),
                    hasSubunit(START_SPECIES, SUBUNIT)."""),
            documentation="If there is a path of a given sign from a complex to a species, then there is a path of that same sign from each subunit of the complex to the species.",
        ),
    ),
)

_INFLUENCES_FROM_PATHS = RuleGroup(
    identifier="influences_from_paths",
    profiles=frozenset({"default", "no_complex"}),
    depends_on=frozenset({"paths_base", "activity_derivation"}),
    rules=(
        Rule(
            identifier="influences_from_paths:positive",
            text=dedent("""\
                new(positivelyInfluences(START_SPECIES, END_SPECIES)) :-
                    path(START_SPECIES, END_SPECIES, positive),
                    new(activity(START_SPECIES)),
                    new(activity(END_SPECIES))."""),
            documentation="If there is a positive path between two species and both are activity nodes, then the first positively influences the second.",
        ),
        Rule(
            identifier="influences_from_paths:negative",
            text=dedent("""\
                new(negativelyInfluences(START_SPECIES, END_SPECIES)) :-
                    path(START_SPECIES, END_SPECIES, negative),
                    new(activity(START_SPECIES)),
                    new(activity(END_SPECIES))."""),
            documentation="If there is a negative path between two species and both are activity nodes, then the first negatively influences the second.",
        ),
    ),
)

_INFLUENCES_CONSUMPTION = RuleGroup(
    identifier="influences_consumption",
    profiles=frozenset({"default", "no_complex"}),
    depends_on=frozenset({"activity_derivation"}),
    rules=(
        Rule(
            identifier="influences_consumption:catalyzer_consumes_reactant",
            text=dedent("""\
                new(negativelyInfluences(START_SPECIES, END_SPECIES)) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER),
                    catalyzer(MODIFIER),
                    hasReferredSpecies(MODIFIER, START_SPECIES),
                    hasReactant(REACTION, REACTANT),
                    hasReferredSpecies(REACTANT, END_SPECIES),
                    new(activity(START_SPECIES)),
                    new(activity(END_SPECIES))."""),
            documentation="If a species is referred to by a catalyzer of a reaction, another species is referred to by a reactant of that reaction, and both are activity nodes, then the first negatively influences the second.",
        ),
        Rule(
            identifier="influences_consumption:physical_stimulator_consumes_reactant",
            text=dedent("""\
                new(negativelyInfluences(START_SPECIES, END_SPECIES)) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER),
                    physicalStimulator(MODIFIER),
                    hasReferredSpecies(MODIFIER, START_SPECIES),
                    hasReactant(REACTION, REACTANT),
                    hasReferredSpecies(REACTANT, END_SPECIES),
                    new(activity(START_SPECIES)),
                    new(activity(END_SPECIES))."""),
            documentation="If a species is referred to by a physical stimulator of a reaction, another species is referred to by a reactant of that reaction, and both are activity nodes, then the first negatively influences the second.",
        ),
        Rule(
            identifier="influences_consumption:trigger_consumes_reactant",
            text=dedent("""\
                new(negativelyInfluences(START_SPECIES, END_SPECIES)) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER),
                    trigger(MODIFIER),
                    hasReferredSpecies(MODIFIER, START_SPECIES),
                    hasReactant(REACTION, REACTANT),
                    hasReferredSpecies(REACTANT, END_SPECIES),
                    new(activity(START_SPECIES)),
                    new(activity(END_SPECIES))."""),
            documentation="If a species is referred to by a trigger of a reaction, another species is referred to by a reactant of that reaction, and both are activity nodes, then the first negatively influences the second.",
        ),
        Rule(
            identifier="influences_consumption:inhibitor_spares_reactant",
            text=dedent("""\
                new(positivelyInfluences(START_SPECIES, END_SPECIES)) :-
                    reaction(REACTION),
                    hasModifier(REACTION, MODIFIER),
                    inhibitor(MODIFIER),
                    hasReferredSpecies(MODIFIER, START_SPECIES),
                    hasReactant(REACTION, REACTANT),
                    hasReferredSpecies(REACTANT, END_SPECIES),
                    new(activity(START_SPECIES)),
                    new(activity(END_SPECIES))."""),
            documentation="If a species is referred to by an inhibitor of a reaction, another species is referred to by a reactant of that reaction, and both are activity nodes, then the first positively influences the second.",
        ),
    ),
)


def _build_registry() -> RuleRegistry:
    registry = RuleRegistry()
    registry.register(
        [
            _ACTIVITY_BASE,
            _ACTIVITY_DERIVATION_FLAT,
            _ACTIVITY_DERIVATION_NO_ACTIVE_SUBUNITS,
            _PATHS_BASE,
            _PATHS_COMPLEX_TRAVERSAL,
            _INFLUENCES_FROM_PATHS,
            _INFLUENCES_CONSUMPTION,
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
