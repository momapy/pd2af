import pytest

import pd2af.rules


_PROFILES = (
    "normal",
    "normal_no_complex",
    "keep_species",
    "keep_species_no_complex",
    "keep_reactions",
    "casq",
)


class TestBuildProgram:
    @pytest.mark.parametrize("profile", _PROFILES)
    def test_profile_returns_program_text(self, profile):
        program = pd2af.rules.build_program(profile)
        assert isinstance(program, str)
        assert len(program) > 0

    @pytest.mark.parametrize("profile", _PROFILES)
    def test_profile_has_activity_and_influence_rules(self, profile):
        program = pd2af.rules.build_program(profile)
        assert "hasActivity" in program
        assert "hasActivityKey" in program
        assert "new(activity(KEY))" in program
        assert "new(positivelyInfluences" in program
        assert "new(negativelyInfluences" in program

    @pytest.mark.parametrize("profile", _PROFILES)
    def test_profile_emits_all_typed_influence_heads(self, profile):
        program = pd2af.rules.build_program(profile)
        for head in (
            "new(positivelyInfluences",
            "new(negativelyInfluences",
            "new(modulates",
            "new(triggers",
            "new(unknownPositivelyInfluences",
            "new(unknownNegativelyInfluences",
            "new(unknownModulates",
            "new(unknownTriggers",
        ):
            assert head in program, (profile, head)

    @pytest.mark.parametrize("profile", _PROFILES)
    def test_profile_fans_out_internal_influences_relation(self, profile):
        program = pd2af.rules.build_program(profile)
        # The typed heads are derived from the internal influences/3 relation.
        assert "influences(SOURCE, TARGET, positivelyInfluences)" in program
        assert "influences(SOURCE_KEY, TARGET_KEY," in program

    def test_path_inference_profiles_carry_kind_through_composes_to(self):
        for profile in ("normal", "normal_no_complex", "keep_species", "keep_species_no_complex"):
            program = pd2af.rules.build_program(profile)
            assert "composesTo(triggers, positivelyInfluences)" in program
            assert "composesTo(INCOMING_INFLUENCE_KIND, OUTGOING_INFLUENCE_KIND)" in program

    def test_normal_no_complex_variants_promote_active_subunits(self):
        for profile in ("normal_no_complex", "keep_species_no_complex"):
            program = pd2af.rules.build_program(profile)
            assert "promotedSubunitKey" in program

    def test_keep_complex_variants_omit_subunit_promotion(self):
        for profile in ("normal", "keep_species"):
            program = pd2af.rules.build_program(profile)
            assert "promotedSubunitKey" not in program

    def test_normal_no_complex_variants_include_complex_traversal(self):
        for profile in ("normal_no_complex", "keep_species_no_complex"):
            program = pd2af.rules.build_program(profile)
            assert "propagatesInfluence(SOURCE, SUBUNIT, INFLUENCE_KIND)" in program
            assert "propagatesInfluence(SUBUNIT, TARGET, INFLUENCE_KIND)" in program

    def test_keep_complex_variants_omit_complex_traversal(self):
        for profile in ("normal", "keep_species"):
            program = pd2af.rules.build_program(profile)
            assert "propagatesInfluence(SOURCE, SUBUNIT, INFLUENCE_KIND)" not in program
            assert "propagatesInfluence(SUBUNIT, TARGET, INFLUENCE_KIND)" not in program

    def test_complex_keeping_profiles_key_by_top_level(self):
        # `normal`/`keep-species` route every species (including subunits) to
        # its top-level complex via the shared `resolvesToTopLevel` relation, so a
        # subunit never gets its own key.
        for profile in ("normal", "keep_species"):
            program = pd2af.rules.build_program(profile)
            assert "resolvesToTopLevel(SPECIES, TOPLEVEL)" in program
            assert "hasActivityKey(SPECIES, keptSpeciesKey(TOPLEVEL))" in program

    def test_no_profile_uses_new_species_from_template(self):
        for profile in _PROFILES:
            program = pd2af.rules.build_program(profile)
            assert "new_species_from_template" not in program

    def test_keep_species_profiles_use_kept_species_only(self):
        for profile in ("keep_species", "keep_species_no_complex"):
            program = pd2af.rules.build_program(profile)
            assert "new_species_from_template" not in program
            assert "keptSpeciesKey" in program

    def test_unknown_profile_raises(self):
        with pytest.raises(ValueError):
            pd2af.rules.build_program("does-not-exist")


class TestCycleAwareTransitivity:
    """Transitive influence extension is gated by the production-cycle
    relations so a source feeding a production cycle does not leak influence
    back around the loop. Only the path-inference profiles extend paths, so only
    they carry the guard; casq has neither the relations nor the guard, and
    `keep_reactions` defines the relations (they come with `paths:core`) but
    never extends a path, so it has no guard either."""

    _PATH_INFERENCE = (
        "normal",
        "normal_no_complex",
        "keep_species",
        "keep_species_no_complex",
    )

    @pytest.mark.parametrize("profile", _PATH_INFERENCE)
    @pytest.mark.parametrize("language", ("celldesigner", "sbgn_pd"))
    def test_path_inference_profiles_define_cycle_relations(self, profile, language):
        program = pd2af.rules.build_program(profile, language=language)
        assert "isDirectlyTransformedTo" in program
        assert "isTransformedTo" in program
        assert "isCyclicallyTransformedTo" in program

    @pytest.mark.parametrize("profile", _PATH_INFERENCE)
    def test_celldesigner_transitivity_guards_against_cycles(self, profile):
        program = pd2af.rules.build_program(profile, language="celldesigner")
        assert (
            "not isCyclicallyTransformedTo(INTERMEDIATE_SPECIES, TARGET_SPECIES)"
            in program
        )

    @pytest.mark.parametrize("profile", _PATH_INFERENCE)
    def test_sbgn_pd_transitivity_guards_against_cycles(self, profile):
        program = pd2af.rules.build_program(profile, language="sbgn_pd")
        assert (
            "not isCyclicallyTransformedTo(INTERMEDIATE_ENTITY_POOL, TARGET_ENTITY_POOL)"
            in program
        )

    def test_casq_omits_cycle_relations(self):
        program = pd2af.rules.build_program("casq")
        assert "isDirectlyTransformedTo" not in program
        assert "isTransformedTo" not in program
        assert "isCyclicallyTransformedTo" not in program


class TestKeepReactionsProfile:
    """`keep_reactions` keeps the PD topology itself: every species is an
    activity and every reaction is rendered as reactant->product positive
    influences, so it carries the single-hop modulation influences but none of
    the inference layers."""

    def test_every_species_is_an_activity_candidate(self):
        program = pd2af.rules.build_program("keep_reactions")
        assert (
            "hasActivityCandidate(SPECIES, isSpecies) :- species(SPECIES)." in program
        )

    def test_reactant_positively_influences_product(self):
        program = pd2af.rules.build_program("keep_reactions")
        assert (
            "propagatesInfluence(REACTANT_SPECIES, PRODUCT_SPECIES, positivelyInfluences)"
            in program
        )

    def test_keeps_direct_modulation_influences(self):
        program = pd2af.rules.build_program("keep_reactions")
        # the modulation-arc rule and its kind table, plus the modifier->product
        # rules, all come from `paths:core`/`modulation_kind`.
        assert "hasModulationKind(MODULATION, triggers) :- triggering(MODULATION)." in program
        assert "hasSource(MODULATION, SOURCE_SPECIES)" in program
        assert "catalyzer(MODIFIER)" in program

    def test_omits_multi_hop_chaining(self):
        program = pd2af.rules.build_program("keep_reactions")
        assert "not isCyclicallyTransformedTo" not in program
        assert "composesTo(INCOMING_INFLUENCE_KIND, OUTGOING_INFLUENCE_KIND)" not in program

    def test_omits_consumption_and_sparing(self):
        program = pd2af.rules.build_program("keep_reactions")
        # the consumption/sparing rules are the only ones that make a reaction's
        # reactant the *target* of a modifier's influence.
        assert "hasReferredElement(REACTANT, TARGET_SPECIES)" not in program
        registry = pd2af.rules.build_registry()
        assert "influences_consumption" not in registry.profiles["keep_reactions"]

    def test_keys_by_top_level_and_keeps_complexes(self):
        program = pd2af.rules.build_program("keep_reactions")
        assert "hasActivityKey(SPECIES, keptSpeciesKey(TOPLEVEL))" in program
        assert "promotedSubunitKey" not in program

    def test_omits_inert_activity_feature_groups(self):
        registry = pd2af.rules.build_registry()
        included = registry.profiles["keep_reactions"]
        assert "activity:core" in included
        for group_id in (
            "activity:phenotype",
            "activity:active_marker",
            "activity:modulation_source",
            "activity:gate_input",
        ):
            assert group_id not in included


class TestInfluencesConsumptionGroup:
    """The consumption/sparing reasoning is a group of its own, carried by the
    path-inference profiles only and excludable on its own."""

    @pytest.mark.parametrize(
        "profile",
        ("normal", "normal_no_complex", "keep_species", "keep_species_no_complex"),
    )
    def test_path_inference_profiles_include_it(self, profile):
        program = pd2af.rules.build_program(profile, "celldesigner")
        assert "hasReferredElement(REACTANT, TARGET_SPECIES)" in program
        registry = pd2af.rules.build_registry()
        assert "influences_consumption" in registry.profiles[profile]

    @pytest.mark.parametrize("profile", ("keep_reactions", "casq"))
    def test_other_profiles_omit_it(self, profile):
        registry = pd2af.rules.build_registry()
        assert "influences_consumption" not in registry.profiles[profile]

    def test_excluding_it_keeps_the_rest_of_the_derivation(self):
        full = pd2af.rules.build_program("keep_species", "celldesigner")
        pruned = pd2af.rules.build_program(
            "keep_species", "celldesigner", exclude_groups=("influences_consumption",)
        )
        assert "new(activity(KEY)) :- hasActivityKey(_, KEY)." in pruned
        assert len(pruned.splitlines()) < len(full.splitlines())


class TestMergedProfilesSbgnPdVariant:
    """`normal`/`normal_no_complex` must emit working rules for SBGN-PD input:
    entity-pool carriers (not the CellDesigner `species` carrier). The
    templated/mergeable gating is gone -- keys are structural roles only and
    PTM stripping happens at the build stage."""

    @pytest.mark.parametrize("profile", ("normal", "normal_no_complex"))
    def test_sbgn_pd_variant_uses_entity_pool_carrier(self, profile):
        program = pd2af.rules.build_program(profile, language="sbgn_pd")
        assert "new_species_from_template" not in program
        assert "isMergeableEntity" not in program
        assert (
            "hasActivityCarrier(ENTITY_POOL, ENTITY_POOL) :- entityPool(ENTITY_POOL)." in program
        )
        # The CellDesigner species carrier must NOT leak into the SBGN-PD program
        # (its absence is exactly the carrier bug this variant fixes).
        assert (
            "hasActivityCarrier(SPECIES, SPECIES) :- species(SPECIES)."
            not in program
        )

    @pytest.mark.parametrize("profile", ("normal", "normal_no_complex"))
    def test_celldesigner_variant_uses_species_carrier(self, profile):
        program = pd2af.rules.build_program(profile, language="celldesigner")
        assert (
            "hasActivityCarrier(SPECIES, SPECIES) :- species(SPECIES)." in program
        )
        assert "isMergeableEntity" not in program
        assert "new_species_from_template" not in program


def test_registry_validates():
    registry = pd2af.rules.build_registry()
    assert registry is not None


def test_registry_registers_all_profiles():
    registry = pd2af.rules.build_registry()
    assert set(registry.profiles) == set(_PROFILES)


class TestExcludeGroups:
    def test_exclude_phenotype_drops_only_that_rule(self):
        full = pd2af.rules.build_program("keep_species", "celldesigner")
        pruned = pd2af.rules.build_program(
            "keep_species", "celldesigner", exclude_groups=("activity:phenotype",)
        )
        phenotype_rule = (
            "hasActivityCandidate(PHENOTYPE, isPhenotype) :- phenotype(PHENOTYPE)."
        )
        assert phenotype_rule in full
        assert phenotype_rule not in pruned
        assert set(pruned.splitlines()) == set(full.splitlines()) - {phenotype_rule}

    def test_exclude_paths_chaining_keeps_single_hop_only(self):
        pruned = pd2af.rules.build_program(
            "normal", "celldesigner", exclude_groups=("paths:chaining",)
        )
        # the multi-hop rule is gone, but the cycle relation it consumed stays.
        assert "not isCyclicallyTransformedTo" not in pruned
        assert "isCyclicallyTransformedTo(UPSTREAM, DOWNSTREAM)" in pruned

    def test_exclude_mandatory_group_raises_friendly_error(self):
        with pytest.raises(ValueError) as excinfo:
            pd2af.rules.build_program(
                "keep_species", "celldesigner", exclude_groups=("activity:core",)
            )
        message = str(excinfo.value)
        assert "cannot exclude 'activity:core'" in message
        assert "required by" in message

    def test_exclude_unregistered_group_raises(self):
        with pytest.raises(ValueError):
            pd2af.rules.build_program(
                "keep_species", "celldesigner", exclude_groups=("does:not:exist",)
            )

    def test_disable_rule_drops_one_table_entry(self):
        catalysis = (
            "hasModulationKind(MODULATION, positivelyInfluences) :- "
            "catalysis(MODULATION)."
        )
        full = pd2af.rules.build_program("keep_species", "celldesigner")
        pruned = pd2af.rules.build_program(
            "keep_species",
            "celldesigner",
            exclude_rules=("modulation_kind:celldesigner:catalysis",),
        )
        assert catalysis in full
        assert catalysis not in pruned

    def test_disable_unknown_rule_raises(self):
        with pytest.raises(ValueError):
            pd2af.rules.build_program(
                "keep_species", "celldesigner", exclude_rules=("no:such:rule",)
            )

    @pytest.mark.parametrize("profile", _PROFILES)
    def test_excludable_groups_are_dependency_leaves(self, profile):
        registry = pd2af.rules.build_registry()
        excludable, mandatory = pd2af.rules.get_excludable_groups(profile)
        included = registry.profiles[profile]
        assert excludable | mandatory == set(included)
        assert not (excludable & mandatory)
        for group_id in included:
            dependents = {
                other
                for other in included
                if group_id in registry.groups[other].depends_on
            }
            if dependents:
                assert group_id in mandatory
            else:
                assert group_id in excludable


class TestCasqProfile:
    def test_casq_includes_delete_and_bridged_product(self):
        program = pd2af.rules.build_program("casq")
        assert "delete(" in program
        assert "rule_1" in program
        assert "rule_2" in program
        assert "rule_3" in program
        assert "rule_4" in program
        assert "bridgesToProduct" in program

    def test_casq_excludes_path_and_complex_traversal(self):
        program = pd2af.rules.build_program("casq")
        assert "propagatesInfluence(" not in program
        assert "propagatesInfluence(SUBUNIT, TARGET, INFLUENCE_KIND)" not in program
        assert "hasContributingComplexAncestor" not in program
        assert "promotedSubunitKey" not in program

    def test_casq_uses_kept_species_keys_only(self):
        program = pd2af.rules.build_program("casq")
        assert "keptSpeciesKey" in program
        assert "new_species_from_template" not in program

    def test_casq_does_not_emit_inhibitor_spares_reactant(self):
        # casq.lp has no rule analogous to inhibitor-spares-reactant; that
        # belongs to the path-based modes.
        program = pd2af.rules.build_program("casq")
        assert "inhibitor_spares_reactant" not in program
