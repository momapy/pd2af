"""Unit tests for the annotation/note remap utility.

``carry_annotations_through_provenance`` re-keys the input-keyed side-tables
onto the output elements through the transform's ``output -> inputs`` provenance,
unioning every source's metadata onto their shared output.
"""

import momapy.utils
from momapy.sbml.model import BQBiol, RDFAnnotation

from pd2af.annotations import carry_annotations_through_provenance


def _annotation(resource):
    return RDFAnnotation(qualifier=BQBiol.IS, resources=frozenset([resource]))


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
