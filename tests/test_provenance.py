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
from pd2af.building.provenance import carry_annotations_through_provenance

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
        annotation_one, annotation_two = _annotation("urn:one"), _annotation("urn:two")
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
        assert output_annotations[output] == frozenset([annotation_one, annotation_two])
        assert output_notes[output] == frozenset([note_one, note_two])

    def test_empty_and_none_inputs_return_empty(self):
        provenance = momapy.utils.FrozenIdentityMultiDict({object(): [object()]})
        assert carry_annotations_through_provenance(provenance, None, None) == ({}, {})
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
            subunit in result.output_element_to_input_elements
            for subunit in surviving_subunits
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


class TestNoCompartmentAnnotationCarry:
    """An activity merged from several compartments gathers the annotations of
    every species that collapsed into it, the way a merged proteoform does."""

    @pytest.fixture(scope="class")
    def annotated_reader_result(self):
        # Glycolysis holds the same annotated species in its neuron and
        # astrocyte compartments.
        return momapy.io.core.read(os.path.join(MAPS_DIR, "Glycolysis.xml"))

    @pytest.mark.parametrize("mode", ["keep-species", "normal"])
    def test_merged_activity_unions_its_sources(self, annotated_reader_result, mode):
        reader_result = annotated_reader_result
        result = pd2af.transform(
            reader_result.obj,
            mode=mode,
            layout_mode=None,
            no_compartment=True,
            element_to_annotations=reader_result.element_to_annotations,
            element_to_notes=reader_result.element_to_notes,
        )
        merged = {
            output_element: input_elements
            for output_element, input_elements in (
                result.output_element_to_input_elements.items()
            )
            if len(input_elements) > 1 and output_element in result.obj.model.species
        }
        assert merged
        annotated_merged = {
            output_element: input_elements
            for output_element, input_elements in merged.items()
            if any(
                reader_result.element_to_annotations.get(input_element)
                for input_element in input_elements
            )
        }
        assert annotated_merged
        for output_element, input_elements in annotated_merged.items():
            expected = frozenset().union(
                *(
                    reader_result.element_to_annotations.get(input_element, frozenset())
                    for input_element in input_elements
                )
            )
            assert expected <= result.element_to_annotations.get(
                output_element, frozenset()
            )

    def test_removed_compartments_leave_provenance(self, annotated_reader_result):
        reader_result = annotated_reader_result
        result = pd2af.transform(
            reader_result.obj,
            mode="keep-species",
            layout_mode=None,
            no_compartment=True,
            element_to_annotations=reader_result.element_to_annotations,
            element_to_notes=reader_result.element_to_notes,
        )
        (compartment,) = result.obj.model.compartments
        compartment_keys = [
            output_element
            for output_element in result.output_element_to_input_elements
            if output_element in result.obj.model.compartments
        ]
        assert compartment_keys == [compartment]


class TestDroppedOperatorProvenance:
    """Provenance covers only the elements actually included in the output: an
    operator that sources no influence never reaches the model, so it must not
    appear as a provenance key; a surviving operator does, carrying its gate's
    origin."""

    def test_srr_dropped_gates_are_absent_from_provenance(self):
        srr_map = momapy.io.core.read(os.path.join(MAPS_DIR, "SRR_signaling.xml")).obj
        result = pd2af.transform(srr_map, mode="keep-species", layout_mode=None)
        assert len(result.obj.model.boolean_logic_gates) == 0
        gate_class_names = {"AndGate", "OrGate", "NotGate", "UnknownGate"}
        assert not any(
            type(output_element).__name__ in gate_class_names
            for output_element in result.output_element_to_input_elements.keys()
        )

    def test_creb_surviving_gate_keeps_provenance(self):
        creb_map = momapy.io.core.read(os.path.join(MAPS_DIR, "CREB_activity.xml")).obj
        result = pd2af.transform(creb_map, mode="keep-species", layout_mode=None)
        (gate,) = result.obj.model.boolean_logic_gates
        gate_inputs = result.output_element_to_input_elements.get(gate)
        assert gate_inputs is not None
        gate_ids = {input_gate.id_ for input_gate in creb_map.model.boolean_logic_gates}
        assert {input_element.id_ for input_element in gate_inputs} & gate_ids
