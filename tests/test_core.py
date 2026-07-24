import os

import pytest

import momapy.celldesigner
import momapy.io.core
from momapy.sbml.model import BQBiol, RDFAnnotation

import pd2af

from tests._helpers import (
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
        example_cd_map, mode="keep-species", layout_mode="plain"
    ).obj


@pytest.fixture(scope="module")
def out_keep_species_no_complex(example_cd_map):
    return pd2af.transform(
        example_cd_map, mode="keep-species-no-complex", layout_mode="plain"
    ).obj


@pytest.fixture(scope="module")
def out_normal_no_complex(example_cd_map):
    if not has_dot_binary():
        pytest.skip("graphviz `dot` binary not on PATH")
    return pd2af.transform(
        example_cd_map, mode="normal-no-complex", layout_mode="auto"
    ).obj


@pytest.fixture(scope="module")
def out_normal(example_cd_map):
    if not has_dot_binary():
        pytest.skip("graphviz `dot` binary not on PATH")
    return pd2af.transform(example_cd_map, mode="normal", layout_mode="auto").obj


class TestTransformExampleKeepSpeciesMode:
    """Golden test: example.xml under keep-species mode and plain layout."""

    def test_returns_celldesigner_map(self, out_keep_species):
        assert isinstance(out_keep_species, momapy.celldesigner.CellDesignerMap)

    def test_preserves_celldesigner_model_type(self, out_keep_species):
        assert isinstance(
            out_keep_species.model, momapy.celldesigner.CellDesignerModel
        )

    def test_expected_active_species(self, out_keep_species):
        # Active subunit C is subsumed into its containing complex D under
        # keep-species (which keeps complexes), so it does not appear as a
        # top-level activity.
        assert species_names(out_keep_species.model) == [
            "B", "D", "E", "F", "G",
        ]

    def test_expected_modulations(self, out_keep_species):
        assert modulation_set(out_keep_species.model) == {
            ("PositiveInfluence", "B", "D"),
            ("PositiveInfluence", "D", "F"),
            ("NegativeInfluence", "B", "E"),
            ("NegativeInfluence", "G", "B"),
        }

    def test_layout_present_with_plain_mode(self, out_keep_species):
        assert out_keep_species.layout is not None
        assert len(out_keep_species.layout.layout_elements) > 0


class TestTransformExampleKeepSpeciesNoComplexMode:
    """Golden test: example.xml under keep-species-no-complex mode."""

    def test_complex_subunits_replace_complex(self, out_keep_species_no_complex):
        # Complex D contains subunits A and C; under keep-species-no-complex,
        # D drops out and influences route through its active subunits
        # instead. Only C has activity, so D->F becomes C->F and B->D becomes
        # B->C.
        names = species_names(out_keep_species_no_complex.model)
        assert "D" not in names
        assert "C" in names

    def test_expected_modulations(self, out_keep_species_no_complex):
        assert modulation_set(out_keep_species_no_complex.model) == {
            ("PositiveInfluence", "B", "C"),
            ("PositiveInfluence", "C", "F"),
            ("NegativeInfluence", "B", "E"),
            ("NegativeInfluence", "G", "B"),
        }


class TestTransformExampleNoComplexMode:
    """Golden test: example.xml under normal-no-complex mode (merge proteoforms,
    drop complexes)."""

    def test_returns_celldesigner_map(self, out_normal_no_complex):
        assert isinstance(out_normal_no_complex, momapy.celldesigner.CellDesignerMap)

    def test_complex_drops_out(self, out_normal_no_complex):
        # Complex D has an active subunit; normal-no-complex drops D in favour of
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
        # of D, so any influence it carries routes to D (same as keep-species).
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


class TestMergedModeStripping:
    """The merged modes strip all PTM decorations (recursively); keep-species
    keeps them. Exercised on a committed map that carries real modifications and
    structural states -- example.xml is decoration-free, so it cannot prove the
    strip on its own."""

    @pytest.fixture(scope="class")
    def rich_map(self):
        return read_cd_map(os.path.join(MAPS_DIR, "Apoptosis_pathway.xml"))

    def test_normal_strips_all_decorations(self, rich_map):
        out = pd2af.transform(rich_map, mode="normal", layout_mode=None).obj
        _assert_recursively_stripped(out.model.species)

    def test_normal_no_complex_strips_all_decorations(self, rich_map):
        out = pd2af.transform(
            rich_map, mode="normal-no-complex", layout_mode=None
        ).obj
        _assert_recursively_stripped(out.model.species)

    def test_keep_species_retains_decorations(self, rich_map):
        out = pd2af.transform(rich_map, mode="keep-species", layout_mode=None).obj
        decorated = [
            s
            for s in out.model.species
            if getattr(s, "modifications", None)
            or getattr(s, "structural_states", None)
        ]
        assert decorated  # keep-species preserves PTM decorations

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

    def test_keep_species_layout_retains_decoration_glyphs(
        self, decorated_layout_map
    ):
        out = pd2af.transform(
            decorated_layout_map, mode="keep-species", layout_mode="plain"
        ).obj
        assert self._count_decoration_glyphs(out.layout) > 0


class TestTransformLayoutModes:
    def test_no_layout_mode_yields_map_without_layout(self, example_cd_map):
        out = pd2af.transform(
            example_cd_map, mode="keep-species", layout_mode=None
        ).obj
        assert isinstance(out, momapy.celldesigner.CellDesignerMap)
        assert out.layout is None
        assert len(out.model.species) > 0

    def test_plain_mode_attaches_layout(self, example_cd_map):
        out = pd2af.transform(
            example_cd_map, mode="keep-species", layout_mode="plain"
        ).obj
        assert out.layout is not None
        assert len(out.layout.layout_elements) > 0

    @pytest.mark.skipif(
        not has_dot_binary(), reason="graphviz `dot` binary not on PATH"
    )
    def test_auto_mode_runs_with_dot(self, example_cd_map):
        out = pd2af.transform(
            example_cd_map, mode="keep-species", layout_mode="auto"
        ).obj
        assert out.layout is not None
        assert len(out.layout.layout_elements) > 0


class TestTransformSetActive:
    """`set_active` marks elements active by id_, surfacing them as
    activities even when the map gives them no structural activity signal."""

    def test_set_active_surfaces_non_active_species(self, example_cd_map):
        # Species A (id `s1`) is a bare, non-active species: absent from the
        # baseline keep-species activities.
        baseline = pd2af.transform(
            example_cd_map, mode="keep-species", layout_mode="plain"
        ).obj
        assert "A" not in species_names(baseline.model)
        with_active = pd2af.transform(
            example_cd_map,
            mode="keep-species",
            layout_mode="plain",
            set_active=["s1"],
        ).obj
        assert "A" in species_names(with_active.model)

    def test_unknown_set_active_id_raises(self, example_cd_map):
        with pytest.raises(ValueError):
            pd2af.transform(
                example_cd_map,
                mode="keep-species",
                layout_mode="plain",
                set_active=["not-an-id"],
            )


class TestTransformSetInactive:
    """`set_inactive` suppresses activities by id_, dropping species that the
    map would otherwise surface as active."""

    def test_set_inactive_suppresses_default_active_species(
        self, example_cd_map
    ):
        # Species B (id `s2`) is active by default in the keep-species output.
        baseline = pd2af.transform(
            example_cd_map, mode="keep-species", layout_mode="plain"
        ).obj
        assert "B" in species_names(baseline.model)
        with_inactive = pd2af.transform(
            example_cd_map,
            mode="keep-species",
            layout_mode="plain",
            set_inactive=["s2"],
        ).obj
        assert "B" not in species_names(with_inactive.model)

    def test_conflicting_set_active_and_set_inactive_raises(
        self, example_cd_map
    ):
        with pytest.raises(ValueError):
            pd2af.transform(
                example_cd_map,
                mode="keep-species",
                layout_mode="plain",
                set_active=["s1"],
                set_inactive=["s1"],
            )

    def test_unknown_set_inactive_id_raises(self, example_cd_map):
        with pytest.raises(ValueError):
            pd2af.transform(
                example_cd_map,
                mode="keep-species",
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
            mode="keep-species",
            layout_mode="plain",
            set_all_active=True,
        ).obj
        # A is bare (never active by default); the toggle surfaces it.
        assert "A" in species_names(with_all_active.model)

    def test_set_inactive_overrides_set_all_active(self, example_cd_map):
        with_override = pd2af.transform(
            example_cd_map,
            mode="keep-species",
            layout_mode="plain",
            set_all_active=True,
            set_inactive=["s2"],
        ).obj
        assert "B" not in species_names(with_override.model)

    def test_both_global_toggles_raises(self, example_cd_map):
        with pytest.raises(ValueError):
            pd2af.transform(
                example_cd_map,
                mode="keep-species",
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
            mode="keep-species",
            layout_mode="plain",
            set_all_inactive=True,
        ).obj
        assert not with_all_inactive.model.species

    def test_set_active_overrides_set_all_inactive(self, example_cd_map):
        with_override = pd2af.transform(
            example_cd_map,
            mode="keep-species",
            layout_mode="plain",
            set_all_inactive=True,
            set_active=["s1"],
        ).obj
        assert species_names(with_override.model) == ["A"]


class TestTransformErrors:
    def test_unknown_mode_raises_value_error(self, example_cd_map):
        with pytest.raises(ValueError):
            pd2af.transform(
                example_cd_map, mode="not-a-mode", layout_mode="plain"
            )

    @pytest.mark.parametrize("mode", ["normal", "normal-no-complex"])
    @pytest.mark.parametrize("layout_mode", ["plain", "overlay"])
    def test_merged_modes_reject_input_derived_layout(
        self, example_cd_map, mode, layout_mode
    ):
        with pytest.raises(ValueError):
            pd2af.transform(example_cd_map, mode=mode, layout_mode=layout_mode)

    @pytest.mark.parametrize(
        "mode", ["keep-species", "keep-species-no-complex"]
    )
    @pytest.mark.parametrize("layout_mode", ["plain", "overlay"])
    def test_keep_species_modes_accept_non_auto_layout(
        self, example_cd_map, mode, layout_mode
    ):
        out = pd2af.transform(
            example_cd_map, mode=mode, layout_mode=layout_mode
        ).obj
        assert out.layout is not None


class TestTransformIsPure:
    """transform should not mutate the input map's model."""

    def test_original_model_unchanged(self, example_map_path):
        cd_map = read_cd_map(example_map_path)
        before_species = sorted(s.id_ for s in cd_map.model.species)
        before_reactions = sorted(r.id_ for r in cd_map.model.reactions)
        pd2af.transform(cd_map, mode="keep-species", layout_mode=None)
        after_species = sorted(s.id_ for s in cd_map.model.species)
        after_reactions = sorted(r.id_ for r in cd_map.model.reactions)
        assert before_species == after_species
        assert before_reactions == after_reactions

    def test_normal_no_complex_does_not_mutate_template(self, example_map_path):
        if not has_dot_binary():
            pytest.skip("graphviz `dot` binary not on PATH")
        cd_map = read_cd_map(example_map_path)
        before_template_ids = sorted(
            t.id_ for t in cd_map.model.species_templates
        )
        before_template_residues = {
            t.id_: getattr(t, "modification_residues", None)
            for t in cd_map.model.species_templates
        }
        pd2af.transform(cd_map, mode="normal-no-complex", layout_mode="auto")
        after_template_ids = sorted(
            t.id_ for t in cd_map.model.species_templates
        )
        after_template_residues = {
            t.id_: getattr(t, "modification_residues", None)
            for t in cd_map.model.species_templates
        }
        assert before_template_ids == after_template_ids
        assert before_template_residues == after_template_residues


class TestProvenance:
    """TransformerResult.provenance maps each output AF element to the input
    elements it derives from; every key is a real element of the output model
    and the inverse round-trips."""

    def test_provenance_keys_are_output_model_elements(self, example_cd_map):
        result = pd2af.transform(
            example_cd_map, mode="keep-species", layout_mode=None
        )
        model = result.obj.model
        output_elements = (
            set(model.species)
            | set(model.boolean_logic_gates)
            | set(model.compartments)
        )
        assert result.provenance  # at least one output has a traced source
        for output_element in result.provenance:
            assert output_element in output_elements

    def test_inverse_round_trips(self, example_cd_map):
        result = pd2af.transform(
            example_cd_map, mode="keep-species", layout_mode=None
        )
        for output_element, input_elements in result.provenance.items():
            for input_element in input_elements:
                assert (
                    output_element
                    in result.provenance.inverse[id(input_element)]
                )


class TestTransformModelInput:
    """`transform` accepts a bare model (not a full map): it returns the
    transformed model, forces `layout_mode` to None, and rejects any explicit
    non-None layout mode."""

    def test_returns_celldesigner_model(self, example_cd_map):
        result = pd2af.transform(example_cd_map.model, mode="keep-species")
        assert isinstance(result.obj, momapy.celldesigner.CellDesignerModel)

    def test_model_output_matches_map_output(self, example_cd_map):
        from_model = pd2af.transform(
            example_cd_map.model, mode="keep-species"
        ).obj
        from_map = pd2af.transform(
            example_cd_map, mode="keep-species", layout_mode=None
        ).obj
        assert species_names(from_model) == species_names(from_map.model)

    def test_provenance_available_for_model_input(self, example_cd_map):
        result = pd2af.transform(example_cd_map.model, mode="keep-species")
        output_elements = (
            set(result.obj.species)
            | set(result.obj.boolean_logic_gates)
            | set(result.obj.compartments)
        )
        assert result.provenance
        for output_element in result.provenance:
            assert output_element in output_elements

    @pytest.mark.parametrize("layout_mode", ["dot", "plain", "overlay"])
    def test_rejects_explicit_layout_mode(self, example_cd_map, layout_mode):
        with pytest.raises(ValueError):
            pd2af.transform(
                example_cd_map.model, mode="keep-species", layout_mode=layout_mode
            )

    @pytest.mark.parametrize("layout_mode", ["auto", None])
    def test_accepts_auto_or_none_layout_mode(self, example_cd_map, layout_mode):
        result = pd2af.transform(
            example_cd_map.model, mode="keep-species", layout_mode=layout_mode
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

    def test_keep_species_carries_species_annotations(
        self, annotated_reader_result
    ):
        reader_result = annotated_reader_result
        result = pd2af.transform(
            reader_result.obj,
            mode="keep-species",
            layout_mode=None,
            element_to_annotations=reader_result.element_to_annotations,
            element_to_notes=reader_result.element_to_notes,
        )
        assert result.element_to_annotations  # something was carried
        # every carried bucket is exactly the union of its provenance sources'
        for output_element, annotations in result.element_to_annotations.items():
            expected = frozenset()
            for source in result.provenance.get(output_element, ()):
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
            mode="keep-species",
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
            mode="keep-species",
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
