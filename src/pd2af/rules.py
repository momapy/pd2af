"""ASP rule composition for the pd2af transformation modes.

The rules are organized as aspcompose groups. Four profiles
(`normal`, `no_complex`, `keep_species`, `keep_species_no_complex`)
compose the groups into the four supported mode programs. The
profiles cover a 2x2 matrix on two orthogonal axes:

* species treatment — ``normal``/``no_complex`` merge proteoforms of
  the same template (and compartment) into a single activity (true
  PD->AF transform, the only style expressible in SBGN PD);
  ``keep_species``/``keep_species_no_complex`` keep each PD species
  as its own activity (CellDesigner-only).
* complex treatment — ``no_complex``/``keep_species_no_complex`` drop
  any complex that has an active subunit and route influences through
  its subunits; ``normal``/``keep_species`` keep complexes as their
  own activities and route influences through them.

All path rules use a single ``path(X, Y, SIGN)`` predicate so
transitivity and complex-subunit traversal can be written once and
parameterized over ``SIGN``.

Activity identity is split into two profile-specific relations so the
shared activity-derivation and influence rules remain
profile-agnostic:

* ``contributesActivity(SPECIES)`` — survives the profile filter and
  contributes an activity to the AF.
* ``activityKey(SPECIES, KEY)`` — this species's identity in the AF
  (a ``kept_species/1`` term in keep_species/keep_species_no_complex;
  a ``kept_species/1`` or ``derived_proteoform_class/2`` term in
  normal/no_complex).
"""

from textwrap import dedent

from aspcompose import CollectionPlan, Rule, RuleGroup, RuleRegistry

_ALL_PROFILES = frozenset(
    {"normal", "no_complex", "keep_species", "keep_species_no_complex"}
)

_ACTIVITY_BASE = RuleGroup(
    identifier="activity_base",
    profiles=_ALL_PROFILES,
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

_HAS_ACTIVE_SUBUNIT = RuleGroup(
    identifier="has_active_subunit",
    profiles=frozenset({"no_complex", "keep_species_no_complex"}),
    depends_on=frozenset({"activity_base"}),
    rules=(
        Rule(
            identifier="has_active_subunit:detect",
            text=dedent("""\
                hasActiveSubunit(SPECIES) :-
                    complex(SPECIES),
                    hasActivity(SPECIES, _),
                    hasSubunit(SPECIES, SUBUNIT),
                    hasActivity(SUBUNIT, _)."""),
            documentation="If a complex has activity and contains a subunit that also has activity, then it has an active subunit.",
        ),
    ),
)

_HAS_CONTRIBUTING_COMPLEX_ANCESTOR = RuleGroup(
    identifier="has_contributing_complex_ancestor",
    profiles=frozenset({"normal", "keep_species"}),
    depends_on=frozenset({"activity_base"}),
    rules=(
        Rule(
            identifier="has_contributing_complex_ancestor:direct",
            text=dedent("""\
                hasContributingComplexAncestor(SUBUNIT) :-
                    complex(COMPLEX),
                    hasActivity(COMPLEX, _),
                    hasSubunit(COMPLEX, SUBUNIT)."""),
            documentation="If a complex has activity, then each of its direct subunits has a contributing complex ancestor.",
        ),
        Rule(
            identifier="has_contributing_complex_ancestor:transitive",
            text=dedent("""\
                hasContributingComplexAncestor(SUBUNIT) :-
                    hasContributingComplexAncestor(COMPLEX),
                    hasSubunit(COMPLEX, SUBUNIT)."""),
            documentation="The contributing-complex-ancestor relation is transitive through complex containment, so deeply nested subunits are also covered.",
        ),
    ),
)

_CONTRIBUTES_ACTIVITY_FLAT = RuleGroup(
    identifier="contributes_activity:flat",
    slot="contributes_activity",
    profiles=frozenset({"normal", "keep_species"}),
    depends_on=frozenset({"activity_base", "has_contributing_complex_ancestor"}),
    rules=(
        Rule(
            identifier="contributes_activity:flat:promote_unless_inside_contributing_complex",
            text=dedent("""\
                contributesActivity(SPECIES) :-
                    hasActivity(SPECIES, _),
                    not hasContributingComplexAncestor(SPECIES)."""),
            documentation="A species with activity contributes an activity to the new AF unless it is a subunit (direct or nested) of a complex that itself contributes activity — the complex represents the active proteoform.",
        ),
    ),
)

_CONTRIBUTES_ACTIVITY_NO_ACTIVE_SUBUNITS = RuleGroup(
    identifier="contributes_activity:no_active_subunits",
    slot="contributes_activity",
    profiles=frozenset({"no_complex", "keep_species_no_complex"}),
    depends_on=frozenset({"activity_base", "has_active_subunit"}),
    rules=(
        Rule(
            identifier="contributes_activity:no_active_subunits:promote_unless_superseded",
            text=dedent("""\
                contributesActivity(SPECIES) :-
                    hasActivity(SPECIES, _),
                    not hasActiveSubunit(SPECIES)."""),
            documentation="A species with activity contributes an activity to the new AF unless it is a complex that has an active subunit.",
        ),
    ),
)

_ACTIVITY_KEY_KEPT = RuleGroup(
    identifier="activity_key:kept",
    slot="activity_key",
    profiles=frozenset({"keep_species", "keep_species_no_complex"}),
    depends_on=frozenset({"activity_base"}),
    rules=(
        Rule(
            identifier="activity_key:kept:identity",
            text="activityKey(SPECIES, kept_species(SPECIES)) :- species(SPECIES).",
            documentation="In keep_species and keep_species_no_complex modes, every species's activity key is itself.",
        ),
    ),
)

_ACTIVITY_KEY_MERGED = RuleGroup(
    identifier="activity_key:merged",
    slot="activity_key",
    profiles=frozenset({"normal", "no_complex"}),
    depends_on=frozenset({"activity_base"}),
    rules=(
        Rule(
            identifier="activity_key:merged:has_some_template",
            text="hasSomeTemplate(SPECIES) :- hasTemplate(SPECIES, _).",
            documentation="A species has some template if it is linked to any template.",
        ),
        Rule(
            identifier="activity_key:merged:has_some_compartment",
            text="hasSomeCompartment(SPECIES) :- hasCompartment(SPECIES, _).",
            documentation="A species has some compartment if it is linked to any compartment.",
        ),
        Rule(
            identifier="activity_key:merged:templated_with_compartment",
            text=dedent("""\
                activityKey(SPECIES, derived_proteoform_class(TEMPLATE, COMPARTMENT)) :-
                    hasTemplate(SPECIES, TEMPLATE),
                    hasCompartment(SPECIES, COMPARTMENT)."""),
            documentation="A templated species in a compartment has a derived-proteoform-class key keyed by template and compartment.",
        ),
        Rule(
            identifier="activity_key:merged:templated_without_compartment",
            text=dedent("""\
                activityKey(SPECIES, derived_proteoform_class(TEMPLATE, no_compartment)) :-
                    hasTemplate(SPECIES, TEMPLATE),
                    not hasSomeCompartment(SPECIES)."""),
            documentation="A templated species without a compartment has a derived-proteoform-class key with the sentinel `no_compartment`.",
        ),
        Rule(
            identifier="activity_key:merged:templateless",
            text=dedent("""\
                activityKey(SPECIES, kept_species(SPECIES)) :-
                    species(SPECIES),
                    not hasSomeTemplate(SPECIES)."""),
            documentation="A templateless species (phenotype, ion, simple molecule, complex, etc.) keeps its own key when proteoforms are merged.",
        ),
    ),
)

_ACTIVITY_DERIVATION = RuleGroup(
    identifier="activity_derivation",
    profiles=_ALL_PROFILES,
    depends_on=frozenset({"contributes_activity", "activity_key"}),
    rules=(
        Rule(
            identifier="activity_derivation:emit",
            text=dedent("""\
                new(activity(KEY)) :-
                    contributesActivity(SPECIES),
                    activityKey(SPECIES, KEY)."""),
            documentation="If a species contributes an activity and has an activity key, then a new activity node with that key is emitted.",
        ),
    ),
)

_PATHS_BASE = RuleGroup(
    identifier="paths_base",
    profiles=_ALL_PROFILES,
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
    profiles=frozenset({"no_complex", "keep_species_no_complex"}),
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
    profiles=_ALL_PROFILES,
    depends_on=frozenset({"paths_base", "contributes_activity", "activity_key"}),
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
            documentation="If there is a positive path between two species and both contribute activities, then a positive influence between their activity keys is emitted.",
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
            documentation="If there is a negative path between two species and both contribute activities, then a negative influence between their activity keys is emitted.",
        ),
    ),
)

_INFLUENCES_CONSUMPTION = RuleGroup(
    identifier="influences_consumption",
    profiles=_ALL_PROFILES,
    depends_on=frozenset({"contributes_activity", "activity_key"}),
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
            documentation="If a species is referred to by a catalyzer of a reaction and another species is referred to by a reactant of that reaction, and both contribute activities, then a negative influence between their activity keys is emitted.",
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
            documentation="If a species is referred to by a physical stimulator of a reaction and another species is referred to by a reactant of that reaction, and both contribute activities, then a negative influence between their activity keys is emitted.",
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
            documentation="If a species is referred to by a trigger of a reaction and another species is referred to by a reactant of that reaction, and both contribute activities, then a negative influence between their activity keys is emitted.",
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
            documentation="If a species is referred to by an inhibitor of a reaction and another species is referred to by a reactant of that reaction, and both contribute activities, then a positive influence between their activity keys is emitted.",
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
