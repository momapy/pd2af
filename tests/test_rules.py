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

    def test_default_profile_has_activity_and_influence_rules(self):
        program = pd2af.rules.build_program("default")
        assert "hasActivity" in program
        assert "new(activity(SPECIES))" in program
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
        # complex traversal subunit rules are unique to no_complex
        assert "hasSubunit(END_SPECIES, SUBUNIT)" not in program
        assert "hasSubunit(START_SPECIES, SUBUNIT)" not in program

    def test_unknown_profile_raises(self):
        with pytest.raises(ValueError):
            pd2af.rules.build_program("does-not-exist")


def test_registry_validates():
    # _build_registry raises if validation fails; calling it is the assertion.
    registry = pd2af.rules._build_registry()
    assert registry is not None
