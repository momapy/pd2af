import pytest

import pd2af.rules


class TestBuildProgram:
    def test_default_profile_returns_program_text(self):
        program = pd2af.rules.build_program("default")
        assert isinstance(program, str)
        assert len(program) > 0

    def test_no_complex_profile_returns_program_text(self):
        program = pd2af.rules.build_program("no_complex")
        assert isinstance(program, str)
        assert len(program) > 0

    def test_pure_af_profile_returns_program_text(self):
        program = pd2af.rules.build_program("pure_af")
        assert isinstance(program, str)
        assert len(program) > 0

    def test_default_profile_has_activity_and_influence_rules(self):
        program = pd2af.rules.build_program("default")
        assert "hasActivity" in program
        assert "contributesActivity" in program
        assert "activityKey" in program
        assert "new(activity(KEY))" in program
        assert "new(positivelyInfluences" in program
        assert "new(negativelyInfluences" in program

    def test_no_complex_profile_includes_complex_traversal_rules(self):
        program = pd2af.rules.build_program("no_complex")
        assert "hasActiveSubunit" in program
        # default profile should not include the no-complex-only rules
        default_program = pd2af.rules.build_program("default")
        assert "hasActiveSubunit" not in default_program

    def test_default_profile_does_not_include_complex_traversal(self):
        program = pd2af.rules.build_program("default")
        # complex traversal subunit rules are unique to no_complex / pure_af
        assert "hasSubunit(END_SPECIES, SUBUNIT)" not in program
        assert "hasSubunit(START_SPECIES, SUBUNIT)" not in program

    def test_pure_af_profile_includes_derived_proteoform_class_rules(self):
        program = pd2af.rules.build_program("pure_af")
        assert "derived_proteoform_class" in program
        assert "no_compartment" in program
        assert "hasActiveSubunit" in program
        # complex traversal carried over
        assert "hasSubunit(END_SPECIES, SUBUNIT)" in program

    def test_default_and_no_complex_use_kept_species_only(self):
        for profile in ("default", "no_complex"):
            program = pd2af.rules.build_program(profile)
            assert "derived_proteoform_class" not in program
            assert "kept_species" in program

    def test_unknown_profile_raises(self):
        with pytest.raises(ValueError):
            pd2af.rules.build_program("does-not-exist")


def test_registry_validates():
    # _build_registry raises if validation fails; calling it is the assertion.
    registry = pd2af.rules._build_registry()
    assert registry is not None


def test_contributes_activity_slot_has_two_fillers():
    registry = pd2af.rules._build_registry()
    fillers = registry.slots["contributes_activity"]
    assert len(fillers) == 2
    assert "contributes_activity:flat" in fillers
    assert "contributes_activity:no_active_subunits" in fillers


def test_activity_key_slot_has_two_fillers():
    registry = pd2af.rules._build_registry()
    fillers = registry.slots["activity_key"]
    assert len(fillers) == 2
    assert "activity_key:kept" in fillers
    assert "activity_key:pure_af" in fillers
