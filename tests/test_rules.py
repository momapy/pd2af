import pytest

import pd2af.rules


_PROFILES = ("normal", "no_complex", "keep_species", "keep_species_no_complex")


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
        assert "contributesActivity" in program
        assert "activityKey" in program
        assert "new(activity(KEY))" in program
        assert "new(positivelyInfluences" in program
        assert "new(negativelyInfluences" in program

    def test_no_complex_variants_include_active_subunit_rules(self):
        for profile in ("no_complex", "keep_species_no_complex"):
            program = pd2af.rules.build_program(profile)
            assert "hasActiveSubunit" in program

    def test_keep_complex_variants_omit_active_subunit_rules(self):
        for profile in ("normal", "keep_species"):
            program = pd2af.rules.build_program(profile)
            assert "hasActiveSubunit" not in program

    def test_no_complex_variants_include_complex_traversal(self):
        for profile in ("no_complex", "keep_species_no_complex"):
            program = pd2af.rules.build_program(profile)
            assert "hasSubunit(END_SPECIES, SUBUNIT)" in program
            assert "hasSubunit(START_SPECIES, SUBUNIT)" in program

    def test_keep_complex_variants_omit_complex_traversal(self):
        for profile in ("normal", "keep_species"):
            program = pd2af.rules.build_program(profile)
            assert "hasSubunit(END_SPECIES, SUBUNIT)" not in program
            assert "hasSubunit(START_SPECIES, SUBUNIT)" not in program

    def test_merged_profiles_use_derived_proteoform_class(self):
        for profile in ("normal", "no_complex"):
            program = pd2af.rules.build_program(profile)
            assert "derived_proteoform_class" in program
            assert "no_compartment" in program

    def test_keep_species_profiles_use_kept_species_only(self):
        for profile in ("keep_species", "keep_species_no_complex"):
            program = pd2af.rules.build_program(profile)
            assert "derived_proteoform_class" not in program
            assert "kept_species" in program

    def test_unknown_profile_raises(self):
        with pytest.raises(ValueError):
            pd2af.rules.build_program("does-not-exist")


def test_registry_validates():
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
    assert "activity_key:merged" in fillers
