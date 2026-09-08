"""Where each output element comes from, and the metadata it inherits.

Provenance maps each output element to the input elements it derives from;
:func:`make_provenance_from_context` builds it from the slots the model pass
filled, and :func:`carry_annotations_through_provenance` reads it back to
re-key the input annotations and notes onto the output.

momapy stores annotations and notes not on model elements but in side-tables on
the ``ReaderResult`` (``element_to_annotations`` / ``element_to_notes``, each a
``Mapping[model_element -> frozenset]``). The transform builds new output
elements, so those input-keyed side-tables no longer address anything in the
output map. This module re-keys them onto the output elements through the
transform's provenance -- which maps each output element to the input elements
it derives from -- unioning the annotations/notes of every merged source onto
their shared output element.
"""

import collections.abc
import typing

import momapy.utils

import pd2af.building.context


def record_provenance_for_subunit_trees(
    output_species: typing.Any,
    input_species: typing.Any,
    input_model_element_to_canonical_model_element: dict,
    record_pair: collections.abc.Callable[[typing.Any, typing.Any], None],
) -> None:
    """Pair the subunits of an ``(output_species, input_species)`` pair.

    Each pairing is recorded, recursing to arbitrary depth for nested complexes.
    Subunits are never activity keys -- the ASP ``topLevel`` relation resolves a
    subunit at any depth to its outermost complex -- so they reach provenance
    only through this walk. Pairing cannot be positional (``subunits`` is a
    ``frozenset``) nor by identity (the merged modes rebuild every subunit while
    stripping, and the kept modes may keep a content-equal twin from another
    complex). Each input subunit resolves to an output subunit through
    ``input_model_element_to_canonical_model_element``, the
    ``id(input) -> canonical`` map the model pass records while stripping and
    promoting, and falls back to content-equality (momapy model elements are
    frozen dataclasses excluding ``id_`` from ``__eq__``/``__hash__``, so a dict
    keyed by subunit is a content index). The canonical is accepted only when it
    is a subunit of the paired output complex, keeping every recorded key a
    genuine part of that complex.

    A no-op for elements without subunits: SBGN-AF activities, logical
    operators, gates and compartments.
    """
    output_subunits = getattr(output_species, "subunits", None)
    input_subunits = getattr(input_species, "subunits", None)
    if not output_subunits or not input_subunits:
        return
    output_subunit_by_content = {
        output_subunit: output_subunit for output_subunit in output_subunits
    }
    for input_subunit in input_subunits:
        canonical_subunit = input_model_element_to_canonical_model_element.get(
            id(input_subunit)
        )
        output_subunit = None
        if canonical_subunit is not None:
            output_subunit = output_subunit_by_content.get(canonical_subunit)
        if output_subunit is None:
            output_subunit = output_subunit_by_content.get(input_subunit)
        if output_subunit is None:
            continue
        record_pair(output_subunit, input_subunit)
        record_provenance_for_subunit_trees(
            output_subunit,
            input_subunit,
            input_model_element_to_canonical_model_element,
            record_pair,
        )


def make_provenance_from_context(
    context: pd2af.building.context.BuilderContext,
) -> momapy.utils.FrozenIdentityMultiDict:
    """Build the output-element -> input-elements provenance mapping.

    Provenance answers "where did this output come from": it maps each output
    AF element to the ``frozenset`` of input PD/CD elements it derives from. The
    relation is many-to-one -- the merged modes content-dedup several input
    proteoforms (or compartments) into one output -- so the origin direction is
    the genuinely multi-valued one, keyed by output.

    Walks the ``key_to_*`` dicts the model pass populates for *every* activity
    and operator atom (the ``*_emissions`` lists only record the first element
    per identity and would miss the many-to-one dedup pairings). Each activity
    key carries the input element's clingo id in ``key.species`` and each
    operator key in ``key.gate``; both resolve back to the input model element
    through ``context.clingo_id_to_model_element``. Compartments are folded in
    from ``context.compartment_emissions`` (already input/output element pairs).

    The subunits of a complex are provenance keys too, paired by
    :func:`record_provenance_for_subunit_trees` from the species pairs above so
    that annotations and notes carried on a subunit reach the output subunit
    that stands for it. Provenance keys are therefore the output species (or
    activities), their subunits at any depth, the gates (or operators) and the
    compartments.

    Returns a :class:`momapy.utils.FrozenIdentityMultiDict` whose forward maps
    each output element to the ``frozenset`` of input elements it derives from,
    and whose ``.inverse`` (keyed by ``id(input_element)``) recovers the output.
    """
    output_element_to_input_elements = {}

    def record_pair(output_element: typing.Any, input_element: typing.Any) -> None:
        output_element_to_input_elements.setdefault(output_element, set()).add(
            input_element
        )

    def record_provenance(clingo_id: str, output_element: typing.Any) -> None:
        if output_element is None:
            return
        input_element = context.clingo_id_to_model_element[clingo_id]
        record_pair(output_element, input_element)
        record_provenance_for_subunit_trees(
            output_element,
            input_element,
            context.input_model_element_to_canonical_model_element,
            record_pair,
        )

    for activity_key, output_element in context.key_to_activity.items():
        record_provenance(activity_key.species, output_element)
    # Only the operators that reached the model (the influence-sourcing ones)
    # have provenance: the others are absent from the output map, so a key
    # pointing at them would be a dangling provenance entry. The keys of
    # content-equal operators that merged into one canonical instance still
    # walk, so every merged origin is recorded.
    emitted_operators = {id(operator) for operator, _ in context.operator_emissions}
    for operator_key, output_operator in context.key_to_operator.items():
        if id(output_operator) not in emitted_operators:
            continue
        record_provenance(operator_key.gate, output_operator)
    for input_compartment, output_compartment in context.compartment_emissions:
        record_pair(output_compartment, input_compartment)

    return momapy.utils.FrozenIdentityMultiDict(
        {
            output_element: frozenset(input_elements)
            for output_element, input_elements in (
                output_element_to_input_elements.items()
            )
        }
    )


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
