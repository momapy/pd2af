import aspcompose
import pytest

import pd2af.modes
import pd2af.asp.rules

from tests._helpers import build_program_for_mode_name


_MODES = (
    "normal",
    "no-complex",
    "keep-reactions",
)


def _group_ids(mode_name):
    """The rule groups a mode is made of, as the mode itself declares them."""
    return frozenset(pd2af.modes.get_transformation_mode(mode_name).rule_group_ids)


class TestBuildProgram:
    @pytest.mark.parametrize("mode", _MODES)
    def test_mode_returns_program_text(self, mode):
        program = build_program_for_mode_name(mode)
        assert isinstance(program, str)
        assert len(program) > 0

    @pytest.mark.parametrize("mode", _MODES)
    def test_mode_has_activity_and_influence_rules(self, mode):
        program = build_program_for_mode_name(mode)
        assert "hasActivity" in program
        assert "hasActivityKey" in program
        assert "new(activity(KEY))" in program
        assert "new(positivelyInfluences" in program
        assert "new(negativelyInfluences" in program

    @pytest.mark.parametrize("mode", _MODES)
    def test_mode_emits_all_typed_influence_heads(self, mode):
        program = build_program_for_mode_name(mode)
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
            assert head in program, (mode, head)

    @pytest.mark.parametrize("mode", _MODES)
    def test_mode_fans_out_internal_influences_relation(self, mode):
        program = build_program_for_mode_name(mode)
        # The typed heads are derived from the internal influences/3 relation.
        assert "influences(SOURCE, TARGET, positivelyInfluences)" in program
        assert "influences(SOURCE_KEY, TARGET_KEY," in program

    def test_path_inference_modes_carry_kind_through_composes_to(self):
        for mode in ("normal", "no-complex"):
            program = build_program_for_mode_name(mode)
            assert "composesTo(triggers, positivelyInfluences)" in program
            assert (
                "composesTo(INCOMING_INFLUENCE_KIND, OUTGOING_INFLUENCE_KIND)"
                in program
            )

    def test_no_complex_promotes_active_subunits(self):
        program = build_program_for_mode_name("no-complex")
        assert "promotedSubunitKey" in program

    def test_normal_omits_subunit_promotion(self):
        program = build_program_for_mode_name("normal")
        assert "promotedSubunitKey" not in program

    def test_no_complex_includes_complex_traversal(self):
        for mode in ("no-complex",):
            program = build_program_for_mode_name(mode)
            assert "propagatesInfluence(SOURCE, SUBUNIT, INFLUENCE_KIND)" in program
            assert "propagatesInfluence(SUBUNIT, TARGET, INFLUENCE_KIND)" in program

    def test_keep_complex_modes_omit_complex_traversal(self):
        for mode in ("normal", "keep-reactions"):
            program = build_program_for_mode_name(mode)
            assert "propagatesInfluence(SOURCE, SUBUNIT, INFLUENCE_KIND)" not in program
            assert "propagatesInfluence(SUBUNIT, TARGET, INFLUENCE_KIND)" not in program

    def test_complex_keeping_modes_key_by_top_level(self):
        # `normal`/`keep-reactions` route every species (including subunits) to
        # its top-level complex via the shared `resolvesToTopLevel` relation, so a
        # subunit never gets its own key.
        for mode in ("normal", "keep-reactions"):
            program = build_program_for_mode_name(mode)
            assert "resolvesToTopLevel(SPECIES, TOPLEVEL)" in program
            assert "hasActivityKey(SPECIES, keptSpeciesKey(TOPLEVEL))" in program

    def test_no_mode_uses_new_species_from_template(self):
        for mode in _MODES:
            program = build_program_for_mode_name(mode)
            assert "new_species_from_template" not in program

    def test_every_mode_uses_kept_species_only(self):
        for mode in _MODES:
            program = build_program_for_mode_name(mode)
            assert "new_species_from_template" not in program
            assert "keptSpeciesKey" in program

    def test_unknown_mode_raises(self):
        with pytest.raises(ValueError):
            build_program_for_mode_name("does-not-exist")


class TestCycleAwareTransitivity:
    """Transitive influence extension is gated by the production-cycle
    relations so a source feeding a production cycle does not leak influence
    back around the loop. Only the path-inference modes extend paths, so only
    they carry the guard; `keep-reactions` defines the relations (they come
    with `paths:core`) but never extends a path, so it has no guard."""

    _PATH_INFERENCE = (
        "normal",
        "no-complex",
    )

    @pytest.mark.parametrize("mode", _PATH_INFERENCE)
    @pytest.mark.parametrize("language", tuple(pd2af.modes.LANGUAGES))
    def test_path_inference_modes_define_cycle_relations(self, mode, language):
        program = build_program_for_mode_name(mode, language=language)
        assert "isDirectlyTransformedTo" in program
        assert "isTransformedTo" in program
        assert "isCyclicallyTransformedTo" in program

    @pytest.mark.parametrize("mode", _PATH_INFERENCE)
    def test_celldesigner_transitivity_guards_against_cycles(self, mode):
        program = build_program_for_mode_name(mode, language="celldesigner")
        assert (
            "not isCyclicallyTransformedTo(INTERMEDIATE_SPECIES, TARGET_SPECIES)"
            in program
        )

    @pytest.mark.parametrize("mode", _PATH_INFERENCE)
    def test_sbgn_pd_transitivity_guards_against_cycles(self, mode):
        program = build_program_for_mode_name(mode, language="sbgn_pd")
        assert (
            "not isCyclicallyTransformedTo(INTERMEDIATE_ENTITY_POOL, TARGET_ENTITY_POOL)"
            in program
        )

    def test_sbgn_pd_operator_input_resolves_through_operator_key(self):
        # An SBGN-PD operator input whose referred element is another logical
        # operator resolves through that operator's key, guarded on the input
        # operator emitting a node.
        program = build_program_for_mode_name("normal", language="sbgn_pd")
        assert (
            "new(logicalOperatorInput(logicalOperatorKey(OPERATOR), "
            "logicalOperatorKey(INPUT_OPERATOR)))" in program
        )
        assert "new(logicalOperator(logicalOperatorKey(INPUT_OPERATOR), _))" in program

    def test_celldesigner_has_no_operator_input_rule(self):
        # A BooleanLogicGateInput always refers to a species, so the CD
        # variant has no operator-referred input edge.
        program = build_program_for_mode_name("normal", language="celldesigner")
        assert "logicalOperatorKey(INPUT_OPERATOR)" not in program


class TestKeepReactionsMode:
    """`keep-reactions` keeps the PD topology itself: every species is an
    activity and every reaction is rendered as reactant->product positive
    influences, so it carries the single-hop modulation influences but none of
    the inference layers."""

    def test_every_species_is_an_activity_candidate(self):
        program = build_program_for_mode_name("keep-reactions")
        assert (
            "hasActivityCandidate(SPECIES, isSpecies) :- species(SPECIES)." in program
        )

    def test_reactant_positively_influences_product(self):
        program = build_program_for_mode_name("keep-reactions")
        assert (
            "propagatesInfluence(REACTANT_SPECIES, PRODUCT_SPECIES, positivelyInfluences)"
            in program
        )

    def test_keeps_direct_modulation_influences(self):
        program = build_program_for_mode_name("keep-reactions")
        # the modulation-arc rule and its kind table, plus the modifier->product
        # rules, all come from `paths:core`/`influences:kind`.
        assert (
            "hasInfluenceKind(MODULATION, triggers) :- triggering(MODULATION)."
            in program
        )
        assert "hasSource(MODULATION, SOURCE_SPECIES)" in program
        assert "catalyzer(MODIFIER)" in program

    def test_omits_multi_hop_chaining(self):
        program = build_program_for_mode_name("keep-reactions")
        assert "not isCyclicallyTransformedTo" not in program
        assert (
            "composesTo(INCOMING_INFLUENCE_KIND, OUTGOING_INFLUENCE_KIND)"
            not in program
        )

    def test_omits_consumption_and_sparing(self):
        program = build_program_for_mode_name("keep-reactions")
        # the consumption/sparing rules are the only ones that make a reaction's
        # reactant the *target* of a modifier's influence.
        assert "hasReferredElement(REACTANT, TARGET_SPECIES)" not in program
        assert "influences:consumption" not in _group_ids("keep-reactions")

    def test_keys_by_top_level_and_keeps_complexes(self):
        program = build_program_for_mode_name("keep-reactions")
        assert "hasActivityKey(SPECIES, keptSpeciesKey(TOPLEVEL))" in program
        assert "promotedSubunitKey" not in program

    def test_omits_inert_activity_feature_groups(self):
        included = _group_ids("keep-reactions")
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
    path-inference modes only and excludable on its own."""

    @pytest.mark.parametrize(
        "mode",
        ("normal", "no-complex"),
    )
    def test_path_inference_modes_include_it(self, mode):
        program = build_program_for_mode_name(mode, "celldesigner")
        assert "hasReferredElement(REACTANT, TARGET_SPECIES)" in program
        assert "influences:consumption" in _group_ids(mode)

    def test_keep_reactions_omits_it(self):
        assert "influences:consumption" not in _group_ids("keep-reactions")

    def test_excluding_it_keeps_the_rest_of_the_derivation(self):
        full = build_program_for_mode_name("normal", "celldesigner")
        pruned = build_program_for_mode_name(
            "normal", "celldesigner", exclude_groups=("influences:consumption",)
        )
        assert "new(activity(KEY)) :- hasActivityKey(_, KEY)." in pruned
        assert len(pruned.splitlines()) < len(full.splitlines())

    def test_sbgn_pd_variant_walks_modulation_to_reactant(self):
        # The SBGN-PD variant follows a modulation arc onto a process down to
        # its reactants: stimulation consumes (negative), inhibition spares.
        program = build_program_for_mode_name("normal", "sbgn_pd")
        assert "stimulation(MODULATION)" in program
        assert "influences(SOURCE_KEY, TARGET_KEY, negativelyInfluences)" in program
        assert "inhibition(MODULATION)" in program
        assert "influences(SOURCE_KEY, TARGET_KEY, positivelyInfluences)" in program
        assert "hasReactant(PROCESS, REACTANT)" in program
        assert "hasSource(MODULATION, SOURCE_ENTITY_POOL)" in program


class TestSbgnPdVariant:
    """`normal`/`no-complex` must emit working rules for SBGN-PD input:
    entity-pool carriers (not the CellDesigner `species` carrier). Keys are
    structural roles only and PTM stripping happens at the build stage."""

    @pytest.mark.parametrize("mode", ("normal", "no-complex"))
    def test_sbgn_pd_variant_uses_entity_pool_carrier(self, mode):
        program = build_program_for_mode_name(mode, language="sbgn_pd")
        assert "new_species_from_template" not in program
        assert "isMergeableEntity" not in program
        assert (
            "hasActivityCarrier(ENTITY_POOL, ENTITY_POOL) :- entityPool(ENTITY_POOL)."
            in program
        )
        # The CellDesigner species carrier must NOT leak into the SBGN-PD program
        # (its absence is exactly the carrier bug this variant fixes).
        assert (
            "hasActivityCarrier(SPECIES, SPECIES) :- species(SPECIES)." not in program
        )

    @pytest.mark.parametrize("mode", ("normal", "no-complex"))
    def test_celldesigner_variant_uses_species_carrier(self, mode):
        program = build_program_for_mode_name(mode, language="celldesigner")
        assert "hasActivityCarrier(SPECIES, SPECIES) :- species(SPECIES)." in program
        assert "isMergeableEntity" not in program
        assert "new_species_from_template" not in program


def test_registry_validates():
    registry = pd2af.asp.rules.build_registry()
    assert registry is not None


def test_registry_holds_every_group_the_modes_name():
    registry = pd2af.asp.rules.build_registry()
    for mode_name in _MODES:
        for group_id in _group_ids(mode_name):
            assert group_id in registry.groups


class TestExcludeGroups:
    def test_exclude_phenotype_drops_only_that_rule(self):
        full = build_program_for_mode_name("normal", "celldesigner")
        pruned = build_program_for_mode_name(
            "normal", "celldesigner", exclude_groups=("activity:phenotype",)
        )
        phenotype_rule = (
            "hasActivityCandidate(PHENOTYPE, isPhenotype) :- phenotype(PHENOTYPE)."
        )
        assert phenotype_rule in full
        assert phenotype_rule not in pruned
        assert set(pruned.splitlines()) == set(full.splitlines()) - {phenotype_rule}

    def test_exclude_paths_chaining_keeps_single_hop_only(self):
        pruned = build_program_for_mode_name(
            "normal", "celldesigner", exclude_groups=("paths:chaining",)
        )
        # the multi-hop rule is gone, but the cycle relation it consumed stays.
        assert "not isCyclicallyTransformedTo" not in pruned
        assert "isCyclicallyTransformedTo(UPSTREAM, DOWNSTREAM)" in pruned

    def test_exclude_mandatory_group_raises_friendly_error(self):
        with pytest.raises(ValueError) as excinfo:
            build_program_for_mode_name(
                "normal", "celldesigner", exclude_groups=("activity:core",)
            )
        message = str(excinfo.value)
        assert "cannot exclude 'activity:core'" in message
        assert "required by" in message

    def test_exclude_unregistered_group_raises(self):
        with pytest.raises(ValueError):
            build_program_for_mode_name(
                "normal", "celldesigner", exclude_groups=("does:not:exist",)
            )

    def test_exclude_rule_drops_one_table_entry(self):
        catalysis = (
            "hasInfluenceKind(MODULATION, positivelyInfluences) :- "
            "catalysis(MODULATION)."
        )
        full = build_program_for_mode_name("normal", "celldesigner")
        pruned = build_program_for_mode_name(
            "normal",
            "celldesigner",
            exclude_rules=("influences:kind:celldesigner:catalysis",),
        )
        assert catalysis in full
        assert catalysis not in pruned

    def test_exclude_unknown_rule_raises(self):
        with pytest.raises(ValueError):
            build_program_for_mode_name(
                "normal", "celldesigner", exclude_rules=("no:such:rule",)
            )

    @pytest.mark.parametrize("mode", _MODES)
    def test_excludable_groups_are_dependency_leaves(self, mode):
        registry = pd2af.asp.rules.build_registry()
        excludable, mandatory = pd2af.asp.rules.get_excludable_groups(mode)
        included = _group_ids(mode)
        assert excludable | mandatory == set(included)
        assert not (excludable & mandatory)
        for group_id in included:
            # A group is depended on either by id or through the slot it fills.
            depended_on_names = {group_id}
            if registry.groups[group_id].slot:
                depended_on_names.add(registry.groups[group_id].slot)
            dependents = {
                other
                for other in included
                if depended_on_names & registry.groups[other].depends_on
            }
            if dependents:
                assert group_id in mandatory
            else:
                assert group_id in excludable


class TestPreparationSlot:
    """`preparation` is an aspcompose slot: every mode must fill it exactly
    once, and the consumers of the activity-key predicates depend on the slot
    rather than on a particular filler."""

    @pytest.mark.parametrize("mode", _MODES)
    def test_every_mode_fills_the_slot(self, mode):
        registry = pd2af.asp.rules.build_registry()
        fillers = registry.slots["preparation"] & _group_ids(mode)
        assert len(fillers) == 1

    def test_excluding_the_filler_raises_rather_than_emitting_no_activities(self):
        with pytest.raises(ValueError) as excinfo:
            build_program_for_mode_name(
                "normal", "celldesigner", exclude_groups=("preparation:complex",)
            )
        assert "slot 'preparation'" in str(excinfo.value)

    @pytest.mark.parametrize("mode", _MODES)
    def test_the_filler_precedes_its_consumers(self, mode):
        program = build_program_for_mode_name(mode)
        assert program.index("hasActivityKey(SPECIES") < program.index(
            "new(activity(KEY)) :- hasActivityKey(_, KEY)."
        )


class TestModeExtensionPoint:
    """A mode contributed through the entry point composes exactly like a
    built-in one: its `rule_group_references` name registered groups and its
    `rule_group_definitions` are registered alongside them."""

    def test_contributed_mode_composes_references_and_definitions(self, monkeypatch):
        contributed_group = aspcompose.RuleGroup(
            identifier="contributed:emit",
            depends_on=frozenset({"influences:output"}),
            docs="A throwaway group contributed by a test-only mode.",
            rules=(
                aspcompose.Rule(
                    identifier="contributed:emit:everything",
                    text="influences(a, b, positivelyInfluences).",
                    docs="A single ground influence, so the program is checkable.",
                ),
            ),
        )
        contributed_mode = pd2af.modes.TransformationMode(
            name="contributed",
            docs="a test-only mode contributed through the entry point",
            rule_group_references=("influences:output",),
            rule_group_definitions=(contributed_group,),
        )

        class _FakeEntryPoint:
            name = "contributed"
            value = "tests.test_rules:contributed_mode"

            @staticmethod
            def load():
                return contributed_mode

        monkeypatch.setattr(
            pd2af.modes.importlib.metadata,
            "entry_points",
            lambda group: [_FakeEntryPoint],
        )
        pd2af.modes.get_transformation_modes.cache_clear()
        try:
            assert "contributed" in pd2af.modes.get_transformation_modes()
            program = build_program_for_mode_name("contributed")
            assert "influences(a, b, positivelyInfluences)." in program
            assert "new(positivelyInfluences(SOURCE, TARGET))" in program
        finally:
            pd2af.modes.get_transformation_modes.cache_clear()

    def test_dangling_group_reference_names_the_mode_and_the_group(self, monkeypatch):
        contributed_mode = pd2af.modes.TransformationMode(
            name="dangling",
            docs="a test-only mode naming a group nothing registers",
            rule_group_references=("influences:output", "no_such_group"),
        )

        class _FakeEntryPoint:
            name = "dangling"
            value = "tests.test_rules:contributed_mode"

            @staticmethod
            def load():
                return contributed_mode

        monkeypatch.setattr(
            pd2af.modes.importlib.metadata,
            "entry_points",
            lambda group: [_FakeEntryPoint],
        )
        pd2af.modes.get_transformation_modes.cache_clear()
        try:
            with pytest.raises(RuntimeError) as raised:
                pd2af.asp.rules.build_registry()
            message = str(raised.value)
            assert "dangling" in message
            assert "no_such_group" in message
            assert "rule_group_references" in message
        finally:
            pd2af.modes.get_transformation_modes.cache_clear()


class TestRuleVariantsKeyOnLanguages:
    """Rule variants are keyed on the input-language tokens, so the registry
    can never grow a variant for a language pd2af does not define."""

    def test_every_variant_key_is_a_registered_language(self):
        registry = pd2af.asp.rules.build_registry()
        variant_keys = {
            variant_key
            for group in registry.groups.values()
            for variant_key in group.variants
        }
        assert variant_keys
        assert variant_keys <= set(pd2af.modes.LANGUAGES)
