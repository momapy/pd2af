"""Unit tests for the annotation/note remap utility.

``carry_annotations_through_provenance`` re-keys the input-keyed side-tables
onto the output elements through the transform's ``output -> inputs`` provenance,
unioning every source's metadata onto their shared output.
"""

import os

import pytest

import momapy.io.core
import momapy.utils
from momapy.sbml.model import BQBiol, RDFAnnotation

import pd2af
from pd2af.annotations import carry_annotations_through_provenance

from tests._helpers import MAPS_DIR


def _annotation(resource):
    return RDFAnnotation(qualifier=BQBiol.IS, resources=frozenset([resource]))


def _collect_subunits(species_iterable):
    """Every subunit of the species, recursively."""
    subunits = set()
    for species in species_iterable:
        for subunit in getattr(species, "subunits", ()) or ():
            subunits.add(subunit)
            subunits |= _collect_subunits(getattr(subunit, "subunits", ()) or ())
    return subunits


class TestCarryAnnotationsThroughProvenance:
    def test_unions_sources_onto_shared_output(self):
        input_one, input_two, output = object(), object(), object()
        provenance = momapy.utils.FrozenIdentityMultiDict(
            {output: [input_one, input_two]}
        )
        annotation_one, annotation_two = _annotation("urn:one"), _annotation(
            "urn:two"
        )
        note_one, note_two = "<body>one</body>", "<body>two</body>"
        output_annotations, output_notes = carry_annotations_through_provenance(
            provenance,
            {
                input_one: frozenset([annotation_one]),
                input_two: frozenset([annotation_two]),
            },
            {
                input_one: frozenset([note_one]),
                input_two: frozenset([note_two]),
            },
        )
        assert output_annotations[output] == frozenset(
            [annotation_one, annotation_two]
        )
        assert output_notes[output] == frozenset([note_one, note_two])

    def test_empty_and_none_inputs_return_empty(self):
        provenance = momapy.utils.FrozenIdentityMultiDict(
            {object(): [object()]}
        )
        assert carry_annotations_through_provenance(
            provenance, None, None
        ) == ({}, {})
        assert carry_annotations_through_provenance(None, {}, {}) == ({}, {})

    def test_map_level_metadata_rekeyed_onto_output_map(self):
        input_map, output_map = object(), object()
        annotation, note = _annotation("urn:map"), "<body>map</body>"
        output_annotations, output_notes = carry_annotations_through_provenance(
            momapy.utils.FrozenIdentityMultiDict({}),
            {input_map: frozenset([annotation])},
            {input_map: frozenset([note])},
            input_map=input_map,
            output_map=output_map,
        )
        assert output_annotations[output_map] == frozenset([annotation])
        assert output_notes[output_map] == frozenset([note])
        # keyed by the OUTPUT map, never the input map
        assert input_map not in output_annotations
        assert input_map not in output_notes


class TestSubunitAnnotationCarry:
    """Annotations on the subunits of a complex reach the output subunits that
    stand for them: subunits are never activity keys, so they reach provenance
    only through the subunit-tree pairing, and CellDesigner writes their RDF
    inside ``<celldesigner:listOfIncludedSpecies>``."""

    @pytest.fixture(scope="class")
    def annotated_reader_result(self):
        # Electron_Transport_Chain_disruption carries uniprot annotations on the
        # subunits of its complexes, not only on top-level species.
        return momapy.io.core.read(
            os.path.join(MAPS_DIR, "Electron_Transport_Chain_disruption.xml")
        )

    @pytest.mark.parametrize(
        "mode",
        ["keep-reactions", "keep-species", "normal", "normal-no-complex"],
    )
    def test_subunit_annotations_carry(self, annotated_reader_result, mode):
        reader_result = annotated_reader_result
        input_subunits = _collect_subunits(reader_result.obj.model.species)
        annotated_input_subunits = {
            subunit
            for subunit in input_subunits
            if reader_result.element_to_annotations.get(subunit)
        }
        assert annotated_input_subunits
        result = pd2af.transform(
            reader_result.obj,
            mode=mode,
            layout_mode=None,
            element_to_annotations=reader_result.element_to_annotations,
            element_to_notes=reader_result.element_to_notes,
        )
        output_subunits = _collect_subunits(result.obj.model.species)
        # a subunit that survives into the output carries its input annotations
        surviving_subunits = annotated_input_subunits & output_subunits
        assert surviving_subunits
        for subunit in surviving_subunits:
            assert reader_result.element_to_annotations[
                subunit
            ] <= result.element_to_annotations.get(subunit, frozenset())
        # and the output subunits are genuine provenance keys
        assert any(
            subunit in result.provenance for subunit in surviving_subunits
        )

    def test_round_trip_preserves_included_species_annotations(
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
        path = str(tmp_path / "out.xml")
        momapy.io.core.write(
            result.obj,
            path,
            writer="celldesigner",
            element_to_annotations=result.element_to_annotations,
            element_to_notes=result.element_to_notes,
        )
        with open(path) as written_file:
            written = written_file.read()
        opening_tag = "<celldesigner:listOfIncludedSpecies"
        closing_tag = "</celldesigner:listOfIncludedSpecies>"
        start = written.index(opening_tag)
        included_species_block = written[start : written.index(closing_tag)]
        assert "urn:miriam:uniprot" in included_species_block
        reread = momapy.io.core.read(path)  # must not raise
        reread_subunits = _collect_subunits(reread.obj.model.species)
        assert any(
            (reread.element_to_annotations or {}).get(subunit)
            for subunit in reread_subunits
        )
