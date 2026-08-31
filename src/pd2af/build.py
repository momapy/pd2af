"""Coordinator: drive the two-phase :class:`pd2af.context.BuilderContext` pipeline.

Pass 1 (:func:`pd2af.celldesigner.building_model.make_and_add_model`) walks clingo
activity / influence atoms and populates ``context.model`` with canonical,
content-deduped compartments, templates, species and modulations.
References between elements are wired to canonical instances at
construction time.

Pass 2 (:func:`pd2af.celldesigner.building_layout.make_and_add_layout`) -- skipped
entirely when ``layout_mode is None`` -- populates ``context.layout`` and
``context.layout_model_mapping``, branching on ``layout_mode``.

``build_map`` creates the :class:`pd2af.context.BuilderContext`, runs the two
passes, then assembles the final map from the three context slots and returns
it with the provenance mapping. :func:`pd2af.core.transform` -- its only caller
-- builds the public :class:`pd2af.core.TransformerResult` from that pair.
"""

import momapy.builder
import momapy.celldesigner
import momapy.sbgn.af
import momapy.utils

import pd2af.celldesigner.building_layout
import pd2af.celldesigner.building_model
import pd2af.context
import pd2af.languages
import pd2af.modes
import pd2af.sbgn.building_layout
import pd2af.sbgn.building_model
import pd2af.utils


def build_map(
    map_,
    layout_mode,
    clingo_model,
    clingo_id_to_model_element,
    influence_pairing="cross",
    *,
    mode: pd2af.modes.TransformationMode,
):
    context = pd2af.context.BuilderContext(
        input_map=map_,
        layout_mode=layout_mode,
        clingo_id_to_model_element=clingo_id_to_model_element,
        influence_pairing=influence_pairing,
        mode=mode,
    )
    language = pd2af.languages.get_language_from_map_or_model(map_)
    if language == pd2af.languages.SBGN_PD:
        pd2af.sbgn.building_model.make_and_add_model(context, clingo_model)
        if layout_mode is not None:
            pd2af.sbgn.building_layout.make_and_add_layout(context)
        map_builder_class = momapy.builder.get_or_make_builder_cls(
            momapy.sbgn.af.SBGNAFMap
        )
        # SBGN-AF logical operators have input/output connectors: rank their
        # logic-arc inputs upstream (reversed_arc_classes) and re-attach the
        # operator arcs to the connector tips after graphviz repositions
        # (operator_arc_resolver). Both hooks are no-ops on operator-free maps,
        # so non-operator SBGN output is unchanged.
        auto_layout_arguments = {
            "compartment_layout_classes": (momapy.sbgn.af.CompartmentLayout,),
            "reversed_arc_classes": (momapy.sbgn.af.LogicArcLayout,),
            "operator_arc_resolver": (
                pd2af.sbgn.building_layout.resolve_operator_arc_segments
            ),
        }
    else:
        pd2af.celldesigner.building_model.make_and_add_model(context, clingo_model)
        if layout_mode is not None:
            pd2af.celldesigner.building_layout.make_and_add_layout(context)
        map_builder_class = momapy.builder.get_or_make_builder_cls(
            momapy.celldesigner.CellDesignerMap
        )
        auto_layout_arguments = {}

    map_builder = map_builder_class(
        model=context.model,
        layout=context.layout,
        layout_model_mapping=context.layout_model_mapping,
    )
    new_map = momapy.builder.object_from_builder(map_builder)

    # the "dot" mode (graphviz) repositions an already-built layout. The
    # compartment-layout classes differ per language (see
    # pd2af.utils.make_auto_layout).
    if layout_mode == "dot":
        new_map = pd2af.utils.make_auto_layout(new_map, **auto_layout_arguments)

    return new_map, make_provenance_from_context(context)


def record_provenance_for_subunit_trees(
    output_species,
    input_species,
    input_model_element_to_canonical_model_element,
    record_pair,
):
    """Pair the subunits of an ``(output_species, input_species)`` pair and
    record each pairing, recursing to arbitrary depth for nested complexes.

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


def make_provenance_from_context(context):
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

    def record_pair(output_element, input_element):
        output_element_to_input_elements.setdefault(output_element, set()).add(
            input_element
        )

    def record_provenance(clingo_id, output_element):
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
    for operator_key, output_element in context.key_to_operator.items():
        record_provenance(operator_key.gate, output_element)
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
