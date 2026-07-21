import pytest

import pd2af.rules


_PROFILES = (
    "normal",
    "normal_no_complex",
    "keep_species",
    "keep_species_no_complex",
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

    def test_non_casq_profiles_carry_kind_through_composes_to(self):
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
    back around the loop. The relations exist in every non-casq profile (which
    has `path`) and are absent from casq (which has none)."""

    _NON_CASQ = ("normal", "normal_no_complex", "keep_species", "keep_species_no_complex")

    @pytest.mark.parametrize("profile", _NON_CASQ)
    @pytest.mark.parametrize("language", ("celldesigner", "sbgn_pd"))
    def test_non_casq_profiles_define_cycle_relations(self, profile, language):
        program = pd2af.rules.build_program(profile, language=language)
        assert "isDirectlyTransformedTo" in program
        assert "isTransformedTo" in program
        assert "isCyclicallyTransformedTo" in program

    @pytest.mark.parametrize("profile", _NON_CASQ)
    def test_celldesigner_transitivity_guards_against_cycles(self, profile):
        program = pd2af.rules.build_program(profile, language="celldesigner")
        assert (
            "not isCyclicallyTransformedTo(INTERMEDIATE_SPECIES, TARGET_SPECIES)"
            in program
        )

    @pytest.mark.parametrize("profile", _NON_CASQ)
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
