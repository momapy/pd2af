import os
import types

import pytest

import momapy.celldesigner
import momapy.io.core
from momapy.sbml.model import BQBiol, RDFAnnotation

import pd2af
import pd2af.modes

from tests._helpers import (
    BINDING_ACTIVATION_MAPS_DIR,
    MAPS_DIR,
    has_dot_binary,
    modulation_set,
    read_cd_map,
    species_names,
)


def _assert_recursively_stripped(species_iterable):
    """Every species (and subunit, recursively) carries no PTM decoration and
    keeps a real, reader-resolvable id (no synthesized prefix)."""
    for species in species_iterable:
        assert not species.id_.startswith("new_species_from_template")
        assert getattr(species, "homomultimer", 1) == 1
        assert not getattr(species, "structural_states", frozenset())
        assert not getattr(species, "modifications", frozenset())
        template = getattr(species, "template", None)
        if template is not None:
            assert not getattr(template, "modification_residues", frozenset())
            assert not getattr(template, "regions", frozenset())
        _assert_recursively_stripped(getattr(species, "subunits", ()) or ())


@pytest.fixture(scope="module")
def out_keep_species(example_cd_map):
    return pd2af.transform(
        example_cd_map, mode="normal", keep_species=True, layout_mode="plain"
    ).obj


@pytest.fixture(scope="module")
def out_keep_species_no_complex(example_cd_map):
    return pd2af.transform(
        example_cd_map, mode="no-complex", keep_species=True, layout_mode="plain"
    ).obj


@pytest.fixture(scope="module")
def out_keep_reactions(example_cd_map):
    return pd2af.transform(
        example_cd_map, mode="keep-reactions", layout_mode="plain"
    ).obj


@pytest.fixture(scope="module")
def out_normal_no_complex(example_cd_map):
    if not has_dot_binary():
        pytest.skip("graphviz `dot` binary not on PATH")
    return pd2af.transform(example_cd_map, mode="no-complex", layout_mode="auto").obj


@pytest.fixture(scope="module")
def out_normal(example_cd_map):
    if not has_dot_binary():
        pytest.skip("graphviz `dot` binary not on PATH")
    return pd2af.transform(example_cd_map, mode="normal", layout_mode="auto").obj


class TestTransformExampleKeepSpeciesMode:
    """Golden test: example.xml under `--keep-species` and plain layout."""

    def test_returns_celldesigner_map(self, out_keep_species):
        assert isinstance(out_keep_species, momapy.celldesigner.CellDesignerMap)

    def test_preserves_celldesigner_model_type(self, out_keep_species):
        assert isinstance(out_keep_species.model, momapy.celldesigner.CellDesignerModel)

    def test_expected_active_species(self, out_keep_species):
        # Active subunit C is subsumed into its containing complex D under
        # `normal` (which keeps complexes), so it does not appear as a
        # top-level activity. A is an activity because binding it activates C.
        assert species_names(out_keep_species.model) == [
            "A",
            "B",
            "D",
            "E",
            "F",
            "G",
        ]

    def test_expected_modulations(self, out_keep_species):
        assert modulation_set(out_keep_species.model) == {
            # A activates C by binding, and C is carried by D
            ("PositiveInfluence", "A", "D"),
            # B catalyzes the reaction producing A, and consumes A
            ("PositiveInfluence", "B", "A"),
            ("NegativeInfluence", "B", "A"),
            ("PositiveInfluence", "B", "D"),
            ("PositiveInfluence", "D", "F"),
            ("NegativeInfluence", "B", "E"),
            ("NegativeInfluence", "G", "B"),
        }

    def test_layout_present_with_plain_mode(self, out_keep_species):
        assert out_keep_species.layout is not None
        assert len(out_keep_species.layout.layout_elements) > 0


class TestTransformExampleKeepSpeciesNoComplexMode:
    """Golden test: example.xml under `no-complex` with `--keep-species`."""

    def test_complex_subunits_replace_complex(self, out_keep_species_no_complex):
        # Complex D contains subunits A and C; under `no-complex`,
        # D drops out and influences route through its active subunits
        # instead. Only C has activity, so D->F becomes C->F and B->D becomes
        # B->C.
        names = species_names(out_keep_species_no_complex.model)
        assert "D" not in names
        assert "C" in names

    def test_expected_modulations(self, out_keep_species_no_complex):
        assert modulation_set(out_keep_species_no_complex.model) == {
            ("PositiveInfluence", "A", "C"),
            ("PositiveInfluence", "B", "A"),
            ("NegativeInfluence", "B", "A"),
            ("PositiveInfluence", "B", "C"),
            ("PositiveInfluence", "C", "F"),
            ("NegativeInfluence", "B", "E"),
            ("NegativeInfluence", "G", "B"),
        }


class TestTransformExampleKeepReactionsMode:
    """Golden test: example.xml under keep-reactions mode and plain layout."""

    def test_every_species_becomes_an_activity(self, out_keep_reactions):
        # All eight top-level species, unlike `normal` which keeps only the
        # five that carry an activity signal. Subunits A and C stay folded into
        # complex D. Two species are named E: the plain form s7 and its active
        # form s7_active.
        assert species_names(out_keep_reactions.model) == [
            "A",
            "B",
            "C",
            "D",
            "E",
            "E",
            "F",
            "G",
        ]

    def test_expected_modulations(self, out_keep_reactions):
        assert modulation_set(out_keep_reactions.model) == {
            # every (reactant, product) pair of every reaction, as a positive
            # influence -- re1 has the same species on both sides (a self
            # influence), re2 is the heterodimer association, re3 is s7 ->
            # s7_active
            ("PositiveInfluence", "A", "A"),
            ("PositiveInfluence", "A", "D"),
            ("PositiveInfluence", "C", "D"),
            ("PositiveInfluence", "E", "E"),
            # the direct modulation influences: two reaction modifiers and the
            # one modulation arc
            ("PositiveInfluence", "B", "A"),
            ("NegativeInfluence", "B", "E"),
            ("PositiveInfluence", "D", "F"),
        }

    def test_reaction_without_products_contributes_nothing(self, out_keep_reactions):
        # re6 (B -> external sink, stimulated by G) has no product, so it yields
        # neither a triggering from B nor a modifier influence from G. G is its
        # only participant that would otherwise become an influence source.
        source_names = {
            modulation.source.name
            for modulation in out_keep_reactions.model.modulations
        }
        assert "G" not in source_names

    def test_omits_inferred_influences(self, out_keep_reactions):
        modulations = modulation_set(out_keep_reactions.model)
        # B -> D needs multi-hop chaining, G -> B needs consumption reasoning;
        # `normal` derives both, keep-reactions derives neither.
        assert ("PositiveInfluence", "B", "D") not in modulations
        assert ("NegativeInfluence", "G", "B") not in modulations

    def test_layout_present_with_plain_mode(self, out_keep_reactions):
        assert out_keep_reactions.layout is not None
        assert len(out_keep_reactions.layout.layout_elements) > 0


class TestTransformExampleNoComplexMode:
    """Golden test: example.xml under `no-complex` (merged forms, dropped
    complexes)."""

    def test_returns_celldesigner_map(self, out_normal_no_complex):
        assert isinstance(out_normal_no_complex, momapy.celldesigner.CellDesignerMap)

    def test_complex_drops_out(self, out_normal_no_complex):
        # Complex D has an active subunit; `no-complex` drops D in favour of
        # the active subunit C.
        names = species_names(out_normal_no_complex.model)
        assert "D" not in names

    def test_synthetic_species_have_clean_template(self, out_normal_no_complex):
        for species in out_normal_no_complex.model.species:
            template = getattr(species, "template", None)
            if template is None:
                continue
            if not template.id_.startswith("merged_template__"):
                continue
            if hasattr(template, "modification_residues"):
                assert template.modification_residues == frozenset()
            if hasattr(template, "regions"):
                assert template.regions == frozenset()

    def test_synthetic_species_have_no_states_or_modifications(
        self, out_normal_no_complex
    ):
        for species in out_normal_no_complex.model.species:
            if not species.id_.startswith("merged__"):
                continue
            assert species.homomultimer == 1
            structural_states = getattr(species, "structural_states", None)
            if structural_states is not None:
                assert structural_states == frozenset()
            modifications = getattr(species, "modifications", None)
            if modifications is not None:
                assert modifications == frozenset()


class TestTransformExampleNormalMode:
    """Golden test: example.xml under normal mode (merge proteoforms, keep
    complexes)."""

    def test_returns_celldesigner_map(self, out_normal):
        assert isinstance(out_normal, momapy.celldesigner.CellDesignerMap)

    def test_complex_is_kept(self, out_normal):
        # Complex D is kept as its own activity. Its active subunit C is no
        # longer a separate top-level activity -- it is a structural component
        # of D, so any influence it carries routes to D (same as `normal`).
        names = species_names(out_normal.model)
        assert "D" in names
        assert "C" not in names

    def test_complex_routes_through_itself(self, out_normal):
        # Influences involving the complex go through the complex (keptSpeciesKey)
        # rather than through subunits.
        mods = modulation_set(out_normal.model)
        assert ("PositiveInfluence", "B", "D") in mods
        assert ("PositiveInfluence", "D", "F") in mods

    def test_complex_kept_with_kept_species_id(self, out_normal):
        # The complex survives with its original id (not a synthesized one).
        complex_species = [
            s
            for s in out_normal.model.species
            if isinstance(s, momapy.celldesigner.Complex)
        ]
        assert complex_species
        for s in complex_species:
            assert not s.id_.startswith("new_species_from_template")

    def test_merged_species_carry_no_decorations(self, out_normal):
        # Merged modes strip every PTM decoration (recursively, incl. subunits)
        # and keep each species' real id -- there is no synthesized id prefix.
        _assert_recursively_stripped(out_normal.model.species)


class TestMergedStripping:
    """The merged reading strips all PTM decorations (recursively);
    `--keep-species` keeps them. Exercised on a committed map that carries real
    modifications and structural states -- example.xml is decoration-free, so it cannot prove the
    strip on its own."""

    @pytest.fixture(scope="class")
    def rich_map(self):
        return read_cd_map(os.path.join(MAPS_DIR, "Apoptosis_pathway.xml"))

    def test_normal_strips_all_decorations(self, rich_map):
        out = pd2af.transform(rich_map, mode="normal", layout_mode=None).obj
        _assert_recursively_stripped(out.model.species)

    def test_no_complex_strips_all_decorations(self, rich_map):
        out = pd2af.transform(rich_map, mode="no-complex", layout_mode=None).obj
        _assert_recursively_stripped(out.model.species)

    @pytest.mark.parametrize("mode", ["normal", "no-complex"])
    def test_keep_species_retains_decorations(self, rich_map, mode):
        out = pd2af.transform(
            rich_map, mode=mode, keep_species=True, layout_mode=None
        ).obj
        decorated = [
            s
            for s in out.model.species
            if getattr(s, "modifications", None)
            or getattr(s, "structural_states", None)
        ]
        assert decorated  # `keep_species` preserves PTM decorations

    # A map with curated residue-modification, structural-state and active-border
    # glyphs, so the auto-layout strip below has every kind of stripped PD layout
    # decoration to remove.
    @pytest.fixture(scope="class")
    def decorated_layout_map(self):
        return read_cd_map(os.path.join(MAPS_DIR, "FOXO3_activity.xml"))

    _ACTIVE_LAYOUT_CLASSES = tuple(
        getattr(momapy.celldesigner, name)
        for name in dir(momapy.celldesigner)
        if name.endswith("ActiveLayout")
    )

    @classmethod
    def _walk_layout(cls, layout):
        yield layout
        for child in getattr(layout, "layout_elements", ()) or ():
            yield from cls._walk_layout(child)

    @classmethod
    def _count_decoration_glyphs(cls, layout):
        # PTM decorations + active-border siblings.
        decoration_classes = (
            momapy.celldesigner.ModificationLayout,
            momapy.celldesigner.StructuralStateLayout,
        ) + cls._ACTIVE_LAYOUT_CLASSES
        return sum(
            isinstance(node, decoration_classes) for node in cls._walk_layout(layout)
        )

    @pytest.mark.skipif(
        not has_dot_binary(), reason="graphviz `dot` binary not on PATH"
    )
    def test_normal_auto_layout_strips_decoration_glyphs(self, decorated_layout_map):
        # The model strip is not enough: the auto layout reuses the input species
        # glyph, whose decoration / active-border sub-glyphs would otherwise
        # render on a stripped species.
        out = pd2af.transform(
            decorated_layout_map, mode="normal", layout_mode="auto"
        ).obj
        assert self._count_decoration_glyphs(out.layout) == 0
        # the real species glyphs (and subunit structure) must still be there
        protein_glyphs = [
            node
            for node in out.layout.descendants()
            if isinstance(node, momapy.celldesigner.GenericProteinLayout)
        ]
        assert protein_glyphs

    def test_keep_species_layout_retains_decoration_glyphs(self, decorated_layout_map):
        out = pd2af.transform(
            decorated_layout_map, mode="normal", keep_species=True, layout_mode="plain"
        ).obj
        assert self._count_decoration_glyphs(out.layout) > 0


class TestTransformLayoutModes:
    def test_no_layout_mode_yields_map_without_layout(self, example_cd_map):
        out = pd2af.transform(
            example_cd_map, mode="normal", keep_species=True, layout_mode=None
        ).obj
        assert isinstance(out, momapy.celldesigner.CellDesignerMap)
        assert out.layout is None
        assert len(out.model.species) > 0

    def test_plain_mode_attaches_layout(self, example_cd_map):
        out = pd2af.transform(
            example_cd_map, mode="normal", keep_species=True, layout_mode="plain"
        ).obj
        assert out.layout is not None
        assert len(out.layout.layout_elements) > 0

    @pytest.mark.skipif(
        not has_dot_binary(), reason="graphviz `dot` binary not on PATH"
    )
    def test_auto_mode_runs_with_dot(self, example_cd_map):
        out = pd2af.transform(
            example_cd_map, mode="normal", keep_species=True, layout_mode="auto"
        ).obj
        assert out.layout is not None
        assert len(out.layout.layout_elements) > 0


class TestTransformSetActive:
    """`set_active` marks elements active by id_, surfacing them as
    activities even when the map gives them no structural activity signal."""

    def test_set_active_surfaces_non_active_species(self, example_cd_map):
        # Species A (id `s1`) is a bare, non-active species: absent from the
        # baseline `keep_species` activities once binding activation, which
        # makes it an activity, is excluded.
        exclude_groups = ("influences:binding_activation",)
        baseline = pd2af.transform(
            example_cd_map,
            mode="normal",
            keep_species=True,
            layout_mode="plain",
            exclude_groups=exclude_groups,
        ).obj
        assert "A" not in species_names(baseline.model)
        with_active = pd2af.transform(
            example_cd_map,
            mode="normal",
            keep_species=True,
            layout_mode="plain",
            set_active=["s1"],
            exclude_groups=exclude_groups,
        ).obj
        assert "A" in species_names(with_active.model)

    def test_unknown_set_active_id_raises(self, example_cd_map):
        with pytest.raises(ValueError):
            pd2af.transform(
                example_cd_map,
                mode="normal",
                keep_species=True,
                layout_mode="plain",
                set_active=["not-an-id"],
            )


class TestTransformSetInactive:
    """`set_inactive` suppresses activities by id_, dropping species that the
    map would otherwise surface as active."""

    def test_set_inactive_suppresses_default_active_species(self, example_cd_map):
        # Species B (id `s2`) is active by default in the `keep_species` output.
        baseline = pd2af.transform(
            example_cd_map, mode="normal", keep_species=True, layout_mode="plain"
        ).obj
        assert "B" in species_names(baseline.model)
        with_inactive = pd2af.transform(
            example_cd_map,
            mode="normal",
            keep_species=True,
            layout_mode="plain",
            set_inactive=["s2"],
        ).obj
        assert "B" not in species_names(with_inactive.model)

    def test_conflicting_set_active_and_set_inactive_raises(self, example_cd_map):
        with pytest.raises(ValueError):
            pd2af.transform(
                example_cd_map,
                mode="normal",
                keep_species=True,
                layout_mode="plain",
                set_active=["s1"],
                set_inactive=["s1"],
            )

    def test_unknown_set_inactive_id_raises(self, example_cd_map):
        with pytest.raises(ValueError):
            pd2af.transform(
                example_cd_map,
                mode="normal",
                keep_species=True,
                layout_mode="plain",
                set_inactive=["not-an-id"],
            )


class TestTransformSetAllActive:
    """`set_all_active` turns every top-level species into an activity;
    per-id `set_inactive` still overrides it and the two global toggles are
    mutually exclusive."""

    def test_set_all_active_surfaces_bare_species(self, example_cd_map):
        with_all_active = pd2af.transform(
            example_cd_map,
            mode="normal",
            keep_species=True,
            layout_mode="plain",
            set_all_active=True,
            exclude_groups=("influences:binding_activation",),
        ).obj
        # A is bare (never active by default, binding activation excluded);
        # the toggle surfaces it.
        assert "A" in species_names(with_all_active.model)

    def test_set_inactive_overrides_set_all_active(self, example_cd_map):
        with_override = pd2af.transform(
            example_cd_map,
            mode="normal",
            keep_species=True,
            layout_mode="plain",
            set_all_active=True,
            set_inactive=["s2"],
        ).obj
        assert "B" not in species_names(with_override.model)

    def test_both_global_toggles_raises(self, example_cd_map):
        with pytest.raises(ValueError):
            pd2af.transform(
                example_cd_map,
                mode="normal",
                keep_species=True,
                layout_mode="plain",
                set_all_active=True,
                set_all_inactive=True,
            )


class TestTransformSetAllInactive:
    """`set_all_inactive` suppresses every activity; per-id `set_active` still
    forces the named ids back on."""

    def test_set_all_inactive_drops_every_species(self, example_cd_map):
        with_all_inactive = pd2af.transform(
            example_cd_map,
            mode="normal",
            keep_species=True,
            layout_mode="plain",
            set_all_inactive=True,
        ).obj
        assert not with_all_inactive.model.species

    def test_set_active_overrides_set_all_inactive(self, example_cd_map):
        with_override = pd2af.transform(
            example_cd_map,
            mode="normal",
            keep_species=True,
            layout_mode="plain",
            set_all_inactive=True,
            set_active=["s1"],
        ).obj
        assert species_names(with_override.model) == ["A"]


class TestTransformErrors:
    def test_unknown_mode_raises_value_error(self, example_cd_map):
        with pytest.raises(ValueError):
            pd2af.transform(example_cd_map, mode="not-a-mode", layout_mode="plain")

    def test_unknown_layout_mode_raises_value_error(self, example_cd_map):
        with pytest.raises(ValueError):
            pd2af.transform(
                example_cd_map,
                mode="normal",
                keep_species=True,
                layout_mode="not-a-layout-mode",
            )

    def test_unknown_influence_pairing_raises_value_error(self, example_cd_map):
        with pytest.raises(ValueError):
            pd2af.transform(
                example_cd_map,
                mode="normal",
                keep_species=True,
                layout_mode="plain",
                influence_pairing="not-a-pairing",
            )

    def test_layout_incompatibility_message_uses_tokens(self, example_cd_map):
        with pytest.raises(ValueError) as error:
            pd2af.transform(example_cd_map, mode="normal", layout_mode="plain")
        message = str(error.value)
        assert "'celldesigner'" in message
        assert "'plain'" in message
        assert "'dot'" in message
        assert "LayoutMode." not in message
        assert "Language." not in message

    def test_language_incompatibility_message_uses_tokens(self, sbgn_example_map):
        with pytest.raises(ValueError) as error:
            pd2af.transform(sbgn_example_map, mode="keep-reactions")
        message = str(error.value)
        assert "'sbgn_pd'" in message
        assert "'celldesigner'" in message
        assert "Language." not in message

    def test_model_input_layout_message_uses_tokens(self, example_cd_map):
        with pytest.raises(ValueError) as error:
            pd2af.transform(
                example_cd_map.model,
                mode="normal",
                keep_species=True,
                layout_mode="plain",
            )
        assert "'plain'" in str(error.value)
        assert "LayoutMode." not in str(error.value)

    @pytest.mark.parametrize("mode", ["normal", "no-complex"])
    @pytest.mark.parametrize("layout_mode", ["plain", "overlay"])
    def test_merged_forms_reject_input_derived_layout(
        self, example_cd_map, mode, layout_mode
    ):
        with pytest.raises(ValueError) as error:
            pd2af.transform(example_cd_map, mode=mode, layout_mode=layout_mode)
        assert "merged forms" in str(error.value)

    @pytest.mark.parametrize("mode", ["normal", "no-complex"])
    @pytest.mark.parametrize("layout_mode", ["plain", "overlay"])
    def test_keep_species_accepts_non_auto_layout(
        self, example_cd_map, mode, layout_mode
    ):
        out = pd2af.transform(
            example_cd_map, mode=mode, keep_species=True, layout_mode=layout_mode
        ).obj
        assert out.layout is not None


class TestTransformDropCompartments:
    """`drop_compartments` merges every compartment into the default one.

    Exercised on a committed map whose neuron and astrocyte compartments hold
    many of the same species -- example.xml declares a single compartment, so it
    cannot show species merging across compartments.
    """

    ALL_MODES = [
        "normal",
        "no-complex",
        "keep-reactions",
    ]

    @pytest.fixture(scope="class")
    def multi_compartment_map(self):
        return read_cd_map(os.path.join(MAPS_DIR, "Glycolysis.xml"))

    @pytest.mark.parametrize("mode", ALL_MODES)
    def test_only_the_default_compartment_survives(self, multi_compartment_map, mode):
        out = pd2af.transform(
            multi_compartment_map, mode=mode, layout_mode=None, drop_compartments=True
        ).obj
        assert len(multi_compartment_map.model.compartments) > 1
        (compartment,) = out.model.compartments
        assert compartment.id_ == "default"

    @pytest.mark.parametrize("mode", ALL_MODES)
    def test_every_species_is_in_the_default_compartment(
        self, multi_compartment_map, mode
    ):
        out = pd2af.transform(
            multi_compartment_map, mode=mode, layout_mode=None, drop_compartments=True
        ).obj
        (compartment,) = out.model.compartments
        assert out.model.species
        assert all(species.compartment is compartment for species in out.model.species)

    @pytest.mark.parametrize("mode", ALL_MODES)
    def test_species_differing_only_by_compartment_merge(
        self, multi_compartment_map, mode
    ):
        kept = pd2af.transform(multi_compartment_map, mode=mode, layout_mode=None).obj
        merged = pd2af.transform(
            multi_compartment_map, mode=mode, layout_mode=None, drop_compartments=True
        ).obj
        assert len(merged.model.species) < len(kept.model.species)

    @pytest.mark.parametrize("mode", ALL_MODES)
    def test_influences_that_become_equal_merge(self, multi_compartment_map, mode):
        kept = pd2af.transform(multi_compartment_map, mode=mode, layout_mode=None).obj
        merged = pd2af.transform(
            multi_compartment_map, mode=mode, layout_mode=None, drop_compartments=True
        ).obj
        assert len(merged.model.modulations) < len(kept.model.modulations)

    def test_compartment_free_map_is_unaffected(self, example_cd_map):
        kept = pd2af.transform(
            example_cd_map, mode="normal", keep_species=True, layout_mode=None
        ).obj
        merged = pd2af.transform(
            example_cd_map,
            mode="normal",
            keep_species=True,
            layout_mode=None,
            drop_compartments=True,
        ).obj
        assert species_names(merged.model) == species_names(kept.model)
        assert modulation_set(merged.model) == modulation_set(kept.model)

    @pytest.mark.parametrize("mode", ALL_MODES)
    def test_output_round_trips(self, tmp_path, multi_compartment_map, mode):
        if not has_dot_binary():
            pytest.skip("graphviz `dot` not available")
        out = pd2af.transform(
            multi_compartment_map, mode=mode, layout_mode="dot", drop_compartments=True
        ).obj
        path = str(tmp_path / "drop_compartments.xml")
        momapy.io.core.write(out, path, writer="celldesigner")
        back = momapy.io.core.read(path).obj
        assert len(back.model.species) == len(out.model.species)
        assert [c.id_ for c in back.model.compartments] == ["default"]

    @pytest.mark.parametrize("mode", ["normal", "keep-reactions"])
    @pytest.mark.parametrize("layout_mode", ["plain", "overlay"])
    def test_input_derived_layout_is_rejected(self, example_cd_map, mode, layout_mode):
        with pytest.raises(ValueError) as error:
            pd2af.transform(
                example_cd_map,
                mode=mode,
                keep_species=True,
                layout_mode=layout_mode,
                drop_compartments=True,
            )
        assert "merged compartments" in str(error.value)


class TestTransformIsPure:
    """transform should not mutate the input map's model."""

    def test_original_model_unchanged(self, example_map_path):
        cd_map = read_cd_map(example_map_path)
        before_species = sorted(s.id_ for s in cd_map.model.species)
        before_reactions = sorted(r.id_ for r in cd_map.model.reactions)
        pd2af.transform(cd_map, mode="normal", keep_species=True, layout_mode=None)
        after_species = sorted(s.id_ for s in cd_map.model.species)
        after_reactions = sorted(r.id_ for r in cd_map.model.reactions)
        assert before_species == after_species
        assert before_reactions == after_reactions

    def test_normal_no_complex_does_not_mutate_template(self, example_map_path):
        if not has_dot_binary():
            pytest.skip("graphviz `dot` binary not on PATH")
        cd_map = read_cd_map(example_map_path)
        before_template_ids = sorted(t.id_ for t in cd_map.model.species_templates)
        before_template_residues = {
            t.id_: getattr(t, "modification_residues", None)
            for t in cd_map.model.species_templates
        }
        pd2af.transform(cd_map, mode="no-complex", layout_mode="auto")
        after_template_ids = sorted(t.id_ for t in cd_map.model.species_templates)
        after_template_residues = {
            t.id_: getattr(t, "modification_residues", None)
            for t in cd_map.model.species_templates
        }
        assert before_template_ids == after_template_ids
        assert before_template_residues == after_template_residues


def _species_and_subunits(species_iterable):
    """Every species of the iterable plus its subunits, recursively."""
    elements = set()
    for species in species_iterable:
        elements.add(species)
        elements |= _species_and_subunits(getattr(species, "subunits", ()) or ())
    return elements


class TestProvenance:
    """TransformerResult.output_element_to_input_elements maps each output AF element to the input
    elements it derives from; every key is a real element of the output model
    (a species, one of its subunits, a gate or a compartment) and the inverse
    round-trips."""

    def test_provenance_keys_are_output_model_elements(self, example_cd_map):
        result = pd2af.transform(
            example_cd_map, mode="normal", keep_species=True, layout_mode=None
        )
        model = result.obj.model
        output_elements = (
            _species_and_subunits(model.species)
            | set(model.boolean_logic_gates)
            | set(model.compartments)
        )
        assert (
            result.output_element_to_input_elements
        )  # at least one output has a traced source
        for output_element in result.output_element_to_input_elements:
            assert output_element in output_elements

    def test_inverse_round_trips(self, example_cd_map):
        result = pd2af.transform(
            example_cd_map, mode="normal", keep_species=True, layout_mode=None
        )
        for (
            output_element,
            input_elements,
        ) in result.output_element_to_input_elements.items():
            for input_element in input_elements:
                assert (
                    output_element
                    in result.output_element_to_input_elements.inverse[
                        id(input_element)
                    ]
                )


class TestTransformModelInput:
    """`transform` accepts a bare model (not a full map): it returns the
    transformed model, forces `layout_mode` to None, and rejects any explicit
    non-None layout mode."""

    def test_returns_celldesigner_model(self, example_cd_map):
        result = pd2af.transform(example_cd_map.model, mode="normal", keep_species=True)
        assert isinstance(result.obj, momapy.celldesigner.CellDesignerModel)

    def test_model_output_matches_map_output(self, example_cd_map):
        from_model = pd2af.transform(
            example_cd_map.model, mode="normal", keep_species=True
        ).obj
        from_map = pd2af.transform(
            example_cd_map, mode="normal", keep_species=True, layout_mode=None
        ).obj
        assert species_names(from_model) == species_names(from_map.model)

    def test_provenance_available_for_model_input(self, example_cd_map):
        result = pd2af.transform(example_cd_map.model, mode="normal", keep_species=True)
        output_elements = (
            _species_and_subunits(result.obj.species)
            | set(result.obj.boolean_logic_gates)
            | set(result.obj.compartments)
        )
        assert result.output_element_to_input_elements
        for output_element in result.output_element_to_input_elements:
            assert output_element in output_elements

    @pytest.mark.parametrize("layout_mode", ["dot", "plain", "overlay"])
    def test_rejects_explicit_layout_mode(self, example_cd_map, layout_mode):
        with pytest.raises(ValueError):
            pd2af.transform(
                example_cd_map.model,
                mode="normal",
                keep_species=True,
                layout_mode=layout_mode,
            )

    @pytest.mark.parametrize("layout_mode", ["auto", None])
    def test_accepts_auto_or_none_layout_mode(self, example_cd_map, layout_mode):
        result = pd2af.transform(
            example_cd_map.model,
            mode="normal",
            keep_species=True,
            layout_mode=layout_mode,
        )
        assert isinstance(result.obj, momapy.celldesigner.CellDesignerModel)


class TestAnnotationCarry:
    """Input RDF annotations and notes are carried onto the output elements
    they derive from, through the transform's provenance, and survive a
    round-trip write/read."""

    @pytest.fixture(scope="class")
    def annotated_reader_result(self):
        # JNK_pathway carries RDF annotations and notes on top-level species.
        return momapy.io.core.read(os.path.join(MAPS_DIR, "JNK_pathway.xml"))

    def test_keep_species_carries_species_annotations(self, annotated_reader_result):
        reader_result = annotated_reader_result
        result = pd2af.transform(
            reader_result.obj,
            mode="normal",
            keep_species=True,
            layout_mode=None,
            element_to_annotations=reader_result.element_to_annotations,
            element_to_notes=reader_result.element_to_notes,
        )
        assert result.element_to_annotations  # something was carried
        # every carried bucket is exactly the union of its provenance sources'
        for output_element, annotations in result.element_to_annotations.items():
            expected = frozenset()
            for source in result.output_element_to_input_elements.get(
                output_element, ()
            ):
                expected |= reader_result.element_to_annotations.get(
                    source, frozenset()
                )
            assert annotations == expected
        # at least one carried element is a top-level species (not just the map)
        species = set(result.obj.model.species)
        assert any(
            output_element in species
            for output_element in result.element_to_annotations
        )

    def test_carries_compartment_annotation(self, annotated_reader_result):
        reader_result = annotated_reader_result
        input_compartments = list(reader_result.obj.model.compartments)
        assert input_compartments
        target_compartment = input_compartments[0]
        annotation = RDFAnnotation(
            qualifier=BQBiol.IS,
            resources=frozenset(["urn:miriam:go:GO:0005737"]),
        )
        # Inject a compartment annotation the source map does not carry.
        element_to_annotations = dict(reader_result.element_to_annotations or {})
        element_to_annotations[target_compartment] = frozenset([annotation])
        result = pd2af.transform(
            reader_result.obj,
            mode="normal",
            keep_species=True,
            layout_mode=None,
            element_to_annotations=element_to_annotations,
            element_to_notes=reader_result.element_to_notes,
        )
        # the output compartment (content-equal to the input) carries it
        assert annotation in result.element_to_annotations.get(
            target_compartment, frozenset()
        )

    def test_round_trip_preserves_a_species_annotation(
        self, annotated_reader_result, tmp_path
    ):
        reader_result = annotated_reader_result
        result = pd2af.transform(
            reader_result.obj,
            mode="normal",
            keep_species=True,
            layout_mode="plain",
            element_to_annotations=reader_result.element_to_annotations,
            element_to_notes=reader_result.element_to_notes,
        )
        species = set(result.obj.model.species)
        carried_species = [
            output_element
            for output_element in result.element_to_annotations
            if output_element in species
            and result.element_to_annotations[output_element]
        ]
        assert carried_species
        target_species = carried_species[0]
        expected_annotations = result.element_to_annotations[target_species]
        path = str(tmp_path / "out.xml")
        momapy.io.core.write(
            result.obj,
            path,
            writer="celldesigner",
            element_to_annotations=result.element_to_annotations,
            element_to_notes=result.element_to_notes,
        )
        reread = momapy.io.core.read(path)  # must not raise
        assert (reread.element_to_annotations or {}).get(
            target_species
        ) == expected_annotations


class TestTransformationOptions:
    """An option is resolved from the user's value, then the mode's, then the
    vocabulary's default -- so a mode suggests and never overrules."""

    def test_keep_reactions_keeps_the_species_apart_by_default(self, example_cd_map):
        kept = pd2af.transform(
            example_cd_map, mode="keep-reactions", layout_mode=None
        ).obj
        merged = pd2af.transform(
            example_cd_map,
            mode="keep-reactions",
            keep_species=False,
            layout_mode=None,
        ).obj
        assert len(merged.model.species) < len(kept.model.species)

    def test_explicit_value_overrules_the_mode(self, example_cd_map):
        explicit = pd2af.transform(
            example_cd_map,
            mode="keep-reactions",
            keep_species=False,
            layout_mode=None,
        ).obj
        _assert_recursively_stripped(explicit.model.species)

    def test_every_default_option_is_a_known_option(self):
        for mode in pd2af.modes.get_transformation_modes().values():
            for name in mode.default_options:
                assert name in pd2af.modes.TRANSFORMATION_OPTIONS

    def test_unknown_default_option_fails_to_load(self):
        mode = pd2af.modes.TransformationMode(
            name="bogus-option-mode",
            docs="names an option that does not exist",
            default_options=types.MappingProxyType({"not_an_option": True}),
        )
        with pytest.raises(ValueError) as error:
            pd2af.modes._check_default_options_are_known(mode)
        assert "not_an_option" in str(error.value)


class TestOptionsDecideTheLayoutModes:
    @pytest.mark.parametrize("layout_mode", ["plain", "overlay"])
    def test_keep_species_puts_the_input_derived_modes_back_on_offer(
        self, example_cd_map, layout_mode
    ):
        out = pd2af.transform(
            example_cd_map, mode="normal", keep_species=True, layout_mode=layout_mode
        ).obj
        assert out.layout is not None

    @pytest.mark.parametrize("layout_mode", ["plain", "overlay"])
    def test_drop_compartments_takes_them_away_again(self, example_cd_map, layout_mode):
        with pytest.raises(ValueError):
            pd2af.transform(
                example_cd_map,
                mode="normal",
                keep_species=True,
                drop_compartments=True,
                layout_mode=layout_mode,
            )

    def test_message_names_merged_forms(self, example_cd_map):
        with pytest.raises(ValueError) as error:
            pd2af.transform(example_cd_map, mode="normal", layout_mode="plain")
        assert "merged forms on 'celldesigner' input" in str(error.value)

    def test_message_names_merged_compartments(self, example_cd_map):
        with pytest.raises(ValueError) as error:
            pd2af.transform(
                example_cd_map,
                mode="normal",
                keep_species=True,
                drop_compartments=True,
                layout_mode="plain",
            )
        assert "merged compartments on 'celldesigner' input" in str(error.value)

    def test_message_names_both(self, example_cd_map):
        with pytest.raises(ValueError) as error:
            pd2af.transform(
                example_cd_map,
                mode="normal",
                drop_compartments=True,
                layout_mode="plain",
            )
        assert "merged forms and merged compartments" in str(error.value)

    def test_message_names_no_merging_when_none_applies(self, sbgn_example_map):
        with pytest.raises(ValueError) as error:
            pd2af.transform(
                sbgn_example_map,
                mode="normal",
                keep_species=True,
                layout_mode="overlay",
            )
        message = str(error.value)
        assert message.startswith("'sbgn_pd' input supports layout_mode")
        assert "merged" not in message


def _binding_activation_map(name):
    return read_cd_map(os.path.join(BINDING_ACTIVATION_MAPS_DIR, f"{name}.xml"))


class TestBindingActivation:
    """`L + R -> L:R` with `R` drawn active only inside the complex: `L`
    activated `R` by binding, so `L` is an activity that positively influences
    the activity carrying `R` (the complex in `normal`, the promoted `R` in
    `no-complex`). Each case is a small map in `tests/maps/binding_activation/`,
    drawn in both languages."""

    def _transform(self, name, mode, **kwargs):
        out = pd2af.transform(
            _binding_activation_map(name), mode=mode, layout_mode=None, **kwargs
        ).obj
        return species_names(out.model), modulation_set(out.model)

    # M -> L, L + R -> L:R (R active only in the complex), L:R + X -> L:R:X.
    def test_normal_targets_the_complex(self):
        names, modulations = self._transform("ligand_receptor", "normal")
        assert names == ["L", "L:R", "L:R:X", "M"]
        assert ("PositiveInfluence", "L", "L:R") in modulations

    def test_no_complex_targets_the_promoted_subunit(self):
        names, modulations = self._transform("ligand_receptor", "no-complex")
        assert names == ["L", "L:R:X", "M", "R"]
        assert ("PositiveInfluence", "L", "R") in modulations

    @pytest.mark.parametrize("mode", ("normal", "no-complex"))
    def test_excluding_the_group_restores_the_plain_reading(self, mode):
        names, modulations = self._transform(
            "ligand_receptor", mode, exclude_groups=("influences:binding_activation",)
        )
        assert "L" not in names
        assert not {one for one in modulations if one[1] == "L"}

    @pytest.mark.parametrize("mode", ("normal", "no-complex"))
    def test_upstream_path_reaches_the_activator(self, mode):
        _, modulations = self._transform("ligand_receptor", mode)
        assert ("PositiveInfluence", "M", "L") in modulations
        assert ("PositiveInfluence", "M", "L:R:X") in modulations

    def test_upstream_path_reaches_the_promoted_subunit(self):
        _, modulations = self._transform("ligand_receptor", "no-complex")
        assert ("PositiveInfluence", "M", "R") in modulations

    @pytest.mark.parametrize("mode", ("normal", "no-complex"))
    def test_binding_edge_does_not_chain_downstream(self, mode):
        _, modulations = self._transform("ligand_receptor", mode)
        assert ("PositiveInfluence", "L", "L:R:X") not in modulations

    @pytest.mark.parametrize("mode", ("normal", "no-complex"))
    def test_suppressed_subunit_draws_no_edge(self, mode):
        names, modulations = self._transform(
            "ligand_receptor", mode, set_inactive=["lr_r_active"]
        )
        assert "L" not in names
        assert not {one for one in modulations if one[1] == "L"}

    # Ras:GTP (Ras active) + Raf -> Ras:GTP:Raf (Ras and Raf active).
    def test_active_recruiter_is_a_source(self):
        _, modulations = self._transform("active_recruiter", "normal")
        assert modulations == {("PositiveInfluence", "Ras:GTP", "Ras:GTP:Raf")}
        names, modulations = self._transform("active_recruiter", "no-complex")
        assert names == ["Raf", "Ras"]
        assert modulations == {("PositiveInfluence", "Ras", "Raf")}

    # A + B -> A:B, both active in the complex.
    @pytest.mark.parametrize("mode", ("normal", "no-complex"))
    def test_two_newly_active_reactants_draw_no_edge(self, mode):
        names, modulations = self._transform("mutual", mode)
        assert names == (["A:B"] if mode == "normal" else ["A", "B"])
        assert modulations == set()

    def test_explicitly_active_reactant_is_not_newly_activated(self):
        names, modulations = self._transform("mutual", "no-complex", set_active=["a"])
        assert names == ["A", "B"]
        assert modulations == {("PositiveInfluence", "A", "B")}

    # L + R:S (R inactive) -> L:(R:S) (R active, nested one level down).
    def test_nested_subunits(self):
        _, modulations = self._transform("nested", "normal")
        assert modulations == {("PositiveInfluence", "L", "L:(R:S)")}
        names, modulations = self._transform("nested", "no-complex")
        assert names == ["L", "R"]
        assert modulations == {("PositiveInfluence", "L", "R")}

    # cAMP (a simple molecule, so no template) + PKA -> cAMP:PKA.
    def test_template_free_species_match_by_class_and_name(self):
        _, modulations = self._transform("template_free", "normal")
        assert modulations == {("PositiveInfluence", "PKA", "cAMP:PKA")}
        _, modulations = self._transform("template_free", "no-complex")
        assert modulations == {("PositiveInfluence", "PKA", "cAMP")}

    # X the simple molecule and X the protein share a name, not an entity.
    @pytest.mark.parametrize("mode", ("normal", "no-complex"))
    def test_equal_names_of_different_classes_do_not_match(self, mode):
        names, modulations = self._transform("different_kinds", mode)
        assert "P" not in names
        assert modulations == set()

    # R -> R:R, both subunits active.
    @pytest.mark.parametrize("mode", ("normal", "no-complex"))
    def test_same_entity_is_never_its_own_activator(self, mode):
        _, modulations = self._transform("homodimer", mode)
        assert modulations == set()

    # L + R <-> L:R (R active).
    @pytest.mark.parametrize("mode", ("normal", "no-complex"))
    def test_reversible_reaction_is_read_in_its_drawn_direction(self, mode):
        _, modulations = self._transform("reversible", mode)
        assert modulations == {
            ("PositiveInfluence", "L", "L:R" if mode == "normal" else "R")
        }

    # L (membrane) + R (cytosol) -> L:R (R active); R also catalyzes Y -> Z.
    @pytest.mark.parametrize("mode", ("normal", "no-complex"))
    @pytest.mark.parametrize("keep_species", (False, True))
    @pytest.mark.parametrize("drop_compartments", (False, True))
    def test_no_self_influence(self, mode, keep_species, drop_compartments):
        names, modulations = self._transform(
            "merge",
            mode,
            keep_species=keep_species,
            drop_compartments=drop_compartments,
        )
        assert (
            "PositiveInfluence",
            "L",
            "L:R" if mode == "normal" else "R",
        ) in modulations
        assert ("PositiveInfluence", "R", "Z") in modulations
        assert not {one for one in modulations if one[1] == one[2]}
        if mode == "no-complex" and not keep_species and drop_compartments:
            # The promoted R merges with the R that catalyzes the other reaction.
            assert names == ["L", "R", "Z"]

    @pytest.mark.parametrize("mode", ("normal", "no-complex"))
    @pytest.mark.parametrize("name", ("ligand_receptor", "nested", "merge"))
    def test_output_round_trips(self, tmp_path, name, mode):
        out = pd2af.transform(
            _binding_activation_map(name),
            mode=mode,
            keep_species=True,
            layout_mode="plain",
        ).obj
        path = os.path.join(tmp_path, f"{name}_{mode}.xml")
        momapy.io.core.write(out, path, writer="celldesigner")
        momapy.io.core.read(path)
