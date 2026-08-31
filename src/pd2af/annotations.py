"""Carry input annotations and notes onto the transformed output elements.

momapy stores annotations and notes not on model elements but in side-tables on
the ``ReaderResult`` (``element_to_annotations`` / ``element_to_notes``, each a
``Mapping[model_element -> frozenset]``). The transform builds new output
elements, so those input-keyed side-tables no longer address anything in the
output map. This module re-keys them onto the output elements through the
transform's provenance -- which maps each output element to the input elements
it derives from -- unioning the annotations/notes of every merged source onto
their shared output element.
"""

import typing

import momapy.utils


def carry_annotations_through_provenance(
    provenance: momapy.utils.FrozenIdentityMultiDict,
    input_element_to_annotations: dict | None,
    input_element_to_notes: dict | None,
    input_map: typing.Any = None,
    output_map: typing.Any = None,
) -> tuple[dict, dict]:
    """Re-key input annotation/note side-tables onto the output elements.

    Args:
        provenance: the transform's ``FrozenIdentityMultiDict`` mapping each
            output element to the ``frozenset`` of input elements it derives
            from (``TransformerResult.provenance``).
        input_element_to_annotations: ``Mapping[input_element -> frozenset]`` of
            RDF annotations from the input ``ReaderResult`` (or ``None``).
        input_element_to_notes: ``Mapping[input_element -> frozenset[str]]`` of
            notes from the input ``ReaderResult`` (or ``None``).
        input_map: the input map, whose own annotations/notes are keyed by the
            map object itself (the map is not a provenance element); pass it
            with ``output_map`` to carry the map-level metadata.
        output_map: the output map -- the key the writer looks the map-level
            annotations up under.

    Returns:
        An ``(output_element_to_annotations, output_element_to_notes)`` tuple of
        plain dicts keyed by output model element with ``frozenset`` values,
        ready to hand to the momapy writer. For every output element the
        annotations/notes of all its provenance input sources are unioned, so a
        merged activity gathers the metadata of every proteoform that collapsed
        into it.
    """
    input_element_to_annotations = input_element_to_annotations or {}
    input_element_to_notes = input_element_to_notes or {}
    output_element_to_annotations = {}
    output_element_to_notes = {}

    def union_metadata_onto_output(
        output_element: typing.Any, input_element: typing.Any
    ) -> None:
        annotations = input_element_to_annotations.get(input_element)
        if annotations:
            output_element_to_annotations[output_element] = (
                output_element_to_annotations.get(output_element, frozenset())
                | annotations
            )
        notes = input_element_to_notes.get(input_element)
        if notes:
            output_element_to_notes[output_element] = (
                output_element_to_notes.get(output_element, frozenset()) | notes
            )

    for output_element, input_elements in (provenance or {}).items():
        for input_element in input_elements:
            union_metadata_onto_output(output_element, input_element)

    # The input map's own annotations/notes are keyed by the map object, which
    # is not a provenance element; re-key them onto the output map object (what
    # the writer looks the map-level metadata up under).
    if input_map is not None and output_map is not None:
        union_metadata_onto_output(output_map, input_map)

    return output_element_to_annotations, output_element_to_notes
