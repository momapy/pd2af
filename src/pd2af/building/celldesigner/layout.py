"""Build the AF layout, reusing the input map's glyphs where possible.

``make_and_add_layout`` is the layout pass: it walks the model built by
:mod:`pd2af.building.celldesigner.model` and populates ``context.layout`` and
``context.layout_model_mapping``, branching on ``layout_mode``. The leaf
primitives below -- synthetic nodes, modulation arcs, mapping helpers, and
background cloning for overlay -- do the per-element construction.

:mod:`pd2af.core` invokes this pass through a
:class:`pd2af.building.context.BuilderContext` after
the model pass.
"""

import dataclasses
import math
import typing

import momapy.builder
import momapy.celldesigner
import momapy.core.layout
import momapy.core.mapping
import momapy.geometry

import pd2af.building.celldesigner.model
import pd2af.building.context
import pd2af.building.layout

_GATE_CLASS_TO_LAYOUT_CLASS = {
    momapy.celldesigner.AndGate: momapy.celldesigner.AndGateLayout,
    momapy.celldesigner.OrGate: momapy.celldesigner.OrGateLayout,
    momapy.celldesigner.NotGate: momapy.celldesigner.NotGateLayout,
    momapy.celldesigner.UnknownGate: momapy.celldesigner.UnknownGateLayout,
}


_SPECIES_CLASS_TO_LAYOUT_CLASS = {
    momapy.celldesigner.GenericProtein: momapy.celldesigner.GenericProteinLayout,
    momapy.celldesigner.TruncatedProtein: momapy.celldesigner.TruncatedProteinLayout,
    momapy.celldesigner.Receptor: momapy.celldesigner.ReceptorLayout,
    momapy.celldesigner.IonChannel: momapy.celldesigner.IonChannelLayout,
    momapy.celldesigner.Gene: momapy.celldesigner.GeneLayout,
    momapy.celldesigner.RNA: momapy.celldesigner.RNALayout,
    momapy.celldesigner.AntisenseRNA: momapy.celldesigner.AntisenseRNALayout,
    momapy.celldesigner.Phenotype: momapy.celldesigner.PhenotypeLayout,
    momapy.celldesigner.Ion: momapy.celldesigner.IonLayout,
    momapy.celldesigner.SimpleMolecule: momapy.celldesigner.SimpleMoleculeLayout,
    momapy.celldesigner.Drug: momapy.celldesigner.DrugLayout,
    momapy.celldesigner.Unknown: momapy.celldesigner.UnknownLayout,
    momapy.celldesigner.Complex: momapy.celldesigner.ComplexLayout,
}

# Species-glyph child sub-glyphs the merged reading strips: the PTM decorations
# (ModificationLayout / StructuralStateLayout) and the active-border
# sibling (a `*ActiveLayout` the reader appends to a species glyph's
# `layout_elements` for an active species). Dropping all of them makes a stripped
# species render plain, matching its decoration-free model (whose `active` flag
# is cleared in `make_stripped_species`, so the two stay in sync).
_STRIPPABLE_SPECIES_DECORATION_CLASSES = (
    momapy.celldesigner.ModificationLayout,
    momapy.celldesigner.StructuralStateLayout,
    momapy.celldesigner.GenericProteinActiveLayout,
    momapy.celldesigner.IonChannelActiveLayout,
    momapy.celldesigner.ComplexActiveLayout,
    momapy.celldesigner.SimpleMoleculeActiveLayout,
    momapy.celldesigner.IonActiveLayout,
    momapy.celldesigner.UnknownActiveLayout,
    momapy.celldesigner.DegradedActiveLayout,
    momapy.celldesigner.GeneActiveLayout,
    momapy.celldesigner.PhenotypeActiveLayout,
    momapy.celldesigner.RNAActiveLayout,
    momapy.celldesigner.AntisenseRNAActiveLayout,
    momapy.celldesigner.TruncatedProteinActiveLayout,
    momapy.celldesigner.ReceptorActiveLayout,
    momapy.celldesigner.DrugActiveLayout,
)


# Influence (modulation) model class -> its arc layout class. NegativeInfluence
# and UnknownNegativeInfluence have no own ``*Layout``; they reuse the inhibition
# arc layouts (mirrors the celldesigner reader, which maps NEGATIVE_INFLUENCE ->
# InhibitionLayout and UNKNOWN_NEGATIVE_INFLUENCE -> UnknownInhibitionLayout).
_MODULATION_CLASS_TO_LAYOUT_CLASS = {
    momapy.celldesigner.PositiveInfluence: momapy.celldesigner.PositiveInfluenceLayout,
    momapy.celldesigner.NegativeInfluence: momapy.celldesigner.InhibitionLayout,
    momapy.celldesigner.Modulation: momapy.celldesigner.ModulationLayout,
    momapy.celldesigner.Triggering: momapy.celldesigner.TriggeringLayout,
    momapy.celldesigner.UnknownPositiveInfluence: momapy.celldesigner.UnknownPositiveInfluenceLayout,
    momapy.celldesigner.UnknownNegativeInfluence: momapy.celldesigner.UnknownInhibitionLayout,
    momapy.celldesigner.UnknownModulation: momapy.celldesigner.UnknownModulationLayout,
    momapy.celldesigner.UnknownTriggering: momapy.celldesigner.UnknownTriggeringLayout,
}


def new_layout_and_mapping_builders() -> typing.Any:
    """A fresh, empty (layout builder, layout-model mapping builder) pair."""
    layout_builder_class = momapy.builder.get_or_make_builder_cls(
        momapy.core.layout.Layout
    )
    return (
        layout_builder_class(),
        momapy.core.mapping.LayoutModelMappingBuilder(),
    )


def clone_layout_pruning_foreground(
    input_layout_element: typing.Any, foreground_ids: set[str], object_to_builder: dict
) -> typing.Any:
    """Clone an input layout element as a builder, to serve as dimmed background.

    Any descendant subtree whose ``id_`` is in ``foreground_ids``
    (already drawn in the foreground) is omitted, and ``None`` is returned when
    the element itself belongs to the foreground.

    Clones share ``object_to_builder`` so structure shared across the input
    layout stays shared once cloned, and they are builders -- not the input
    frozen objects -- so dimming the background never mutates the input
    map's layout.
    """
    if input_layout_element.id_ in foreground_ids:
        return None
    clone = momapy.builder.builder_from_object(
        input_layout_element, object_to_builder=object_to_builder
    )
    _prune_foreground_from_clone(
        input_layout_element, clone, foreground_ids, object_to_builder
    )
    return clone


def _prune_foreground_from_clone(
    input_layout_element: typing.Any,
    clone: typing.Any,
    foreground_ids: set[str],
    object_to_builder: dict,
):
    # The only foreground glyphs that can sit *inside* a background clone are
    # subunits nested via `layout_elements` (e.g. promoted subunits of a
    # dissolved complex). `builder_from_object` already cloned the whole
    # subtree into `object_to_builder`; rebuild `layout_elements` keeping only
    # the clones whose input original is not in the foreground, recursing so
    # deeper nestings are pruned too.
    input_subunits = getattr(input_layout_element, "layout_elements", None)
    if not input_subunits:
        return
    surviving_clones = []
    for input_subunit in input_subunits:
        if input_subunit.id_ in foreground_ids:
            continue
        subunit_clone = object_to_builder[id(input_subunit)]
        _prune_foreground_from_clone(
            input_subunit, subunit_clone, foreground_ids, object_to_builder
        )
        surviving_clones.append(subunit_clone)
    clone.layout_elements = surviving_clones


def make_decoration_stripped_layout(
    input_layout: typing.Any, original_to_stripped: dict[int, typing.Any] | None
) -> typing.Any:
    """Return a frozen copy of ``input_layout`` stripped of its decoration glyphs.

    Every PTM-decoration and active-border sub-glyph is dropped, at any depth,
    and ``id(original) -> stripped`` recorded for every kept element.

    The merged reading (`keep_species` off) strips the model's species of their
    decorations; this drops the matching layout glyphs
    (``ModificationLayout`` / ``StructuralStateLayout`` and the ``*ActiveLayout``
    active-border sibling -- see ``_STRIPPABLE_SPECIES_DECORATION_CLASSES``) so a
    stripped species renders plain. Subunit glyphs of kept complexes are
    preserved (and recursively stripped), so the real glyph size and subunit
    structure -- which the dot layout deliberately keeps -- survive.

    Frozen objects throughout (via :func:`dataclasses.replace`), never builders,
    so the result pickles and round-trips like the reused input glyphs of the
    non-stripping path; a subtree that needs no change is returned unchanged
    (shared with the input). ``original_to_stripped`` lets the caller map the
    surviving glyphs back to their model elements
    (:func:`add_mappings_for_layout_and_descendants`).
    """
    children = getattr(input_layout, "layout_elements", ()) or ()
    stripped_children = []
    changed = False
    for child in children:
        if isinstance(child, _STRIPPABLE_SPECIES_DECORATION_CLASSES):
            changed = True
            continue
        stripped_child = make_decoration_stripped_layout(child, original_to_stripped)
        if stripped_child is not child:
            changed = True
        stripped_children.append(stripped_child)
    if changed:
        stripped = dataclasses.replace(
            input_layout, layout_elements=tuple(stripped_children)
        )
    else:
        stripped = input_layout
    original_to_stripped[id(input_layout)] = stripped
    return stripped


def make_synthetic_layout(species: typing.Any, index: int) -> typing.Any:
    """Build a placeholder Node for ``species``.

    ``index`` seeds the position so two synthetic layouts for
    content-equal species (same class, same name) are themselves
    content-distinct. Layout dataclasses have ``compare=False`` on
    ``id_``, so identical content collapses to one mapping key — and
    the layout-model mapping then drops every synthetic but the first.
    The position itself is throwaway: ``make_auto_layout`` repositions
    every node before render.
    """
    layout_class = _SPECIES_CLASS_TO_LAYOUT_CLASS.get(type(species))
    if layout_class is None:
        raise ValueError(
            f"no default layout class registered for species type "
            f"{type(species).__name__}"
        )
    position = momapy.geometry.Point(float(index), 0.0)
    label = momapy.core.layout.TextLayout(
        text=getattr(species, "name", "") or "",
        position=position,
    )
    return layout_class(position=position, label=label)


def make_modulation_arc(
    context: pd2af.building.context.BuilderContext,
    modulation: typing.Any,
    source_layout: typing.Any,
    target_layout: typing.Any,
) -> typing.Any:
    """Build the arc layout for ``modulation`` between two node layouts.

    In the ``dot`` mode the segments are placeholders that graphviz overwrites;
    otherwise they are computed from the two layouts' real geometry.
    """
    arc_class = _MODULATION_CLASS_TO_LAYOUT_CLASS[type(modulation)]
    if context.layout_mode == "dot":
        segments = pd2af.building.layout.PLACEHOLDER_ARC_SEGMENTS
    else:
        segments = pd2af.building.layout.make_arc_segments_from_source_and_target(
            source_layout, target_layout
        )
    return arc_class(
        source=source_layout,
        target=target_layout,
        segments=tuple(segments),
    )


def make_synthetic_gate_layout(gate: typing.Any, index: int) -> typing.Any:
    """Build a placeholder gate node for ``gate``.

    Used by the ``dot`` mode when the input gate has no curated layout.
    ``index`` seeds the position so two content-equal synthetic gates stay
    content-distinct; ``make_auto_layout`` repositions every node before render.
    """
    layout_class = _GATE_CLASS_TO_LAYOUT_CLASS.get(type(gate))
    if layout_class is None:
        raise ValueError(
            f"no default layout class registered for gate type {type(gate).__name__}"
        )
    position = momapy.geometry.Point(float(index), 0.0)
    return layout_class(position=position)


def make_logic_arc(
    context: pd2af.building.context.BuilderContext,
    gate_layout: typing.Any,
    input_layout: typing.Any,
) -> typing.Any:
    """Build a ``LogicArcLayout`` from a gate to one of its input species.

    The CellDesigner writer locates a gate's inputs by scanning for logic arcs
    whose ``source`` is the gate layout, so the arc runs gate -> input species
    (the CellDesigner convention), not input -> gate.
    """
    if context.layout_mode == "dot":
        segments = pd2af.building.layout.PLACEHOLDER_ARC_SEGMENTS
    else:
        segments = pd2af.building.layout.make_arc_segments_from_source_and_target(
            gate_layout, input_layout
        )
    return momapy.celldesigner.LogicArcLayout(
        source=gate_layout,
        target=input_layout,
        segments=tuple(segments),
    )


def add_mappings_for_layout_and_descendants(
    context: pd2af.building.context.BuilderContext,
    input_layout: typing.Any,
    original_to_stripped: dict[int, typing.Any] | None = None,
):
    """Map an input glyph and its descendants to their model elements.

    With ``original_to_stripped``, the glyphs actually placed are the stripped
    clones of :func:`make_decoration_stripped_layout`: a clone may be a fresh
    ``dataclasses.replace`` copy that is not a key of the input layout-model
    mapping (keyed by the input frozen objects), so we walk the input originals
    and follow each to its stripped counterpart. Dropped decoration glyphs are
    absent from that map, so they never get mapped -- and being absent from the
    layout, never render.
    """
    input_layout_model_mapping = context.input_map.layout_model_mapping
    for original in [input_layout] + list(input_layout.descendants()):
        placed = original
        if original_to_stripped is not None:
            placed = original_to_stripped.get(id(original))
            if placed is None:
                continue
        model_value = input_layout_model_mapping.get_mapping(original)
        if model_value is None:
            continue
        model_value = context.input_model_element_to_canonical_model_element.get(
            id(model_value), model_value
        )
        context.layout_model_mapping.add_mapping(placed, model_value)


def add_modulation_mapping(
    mapping_builder: typing.Any,
    arc: typing.Any,
    source_layout: typing.Any,
    target_layout: typing.Any,
    modulation: typing.Any,
):
    """Map the arc, together with its two endpoint clusters, to ``modulation``.

    An endpoint that is already part of a cluster (a species drawn as several
    glyphs) contributes that whole cluster, so the mapping stays consistent with
    the one the endpoint's own layout was registered under.
    """
    mapping_builder.add_mapping(
        frozenset([arc])
        | pd2af.building.layout.get_key_of_layout_element(
            mapping_builder, source_layout
        )
        | pd2af.building.layout.get_key_of_layout_element(
            mapping_builder, target_layout
        ),
        modulation,
        representative=arc,
    )


# ---------------------------------------------------------------------------
# Pass 2 -- layout construction
# ---------------------------------------------------------------------------


def make_and_add_layout(context: pd2af.building.context.BuilderContext):
    """Build ``context.layout`` and ``context.layout_model_mapping`` (pass 2)."""
    context.layout, context.layout_model_mapping = new_layout_and_mapping_builders()

    for compartment in pd2af.building.celldesigner.model.compartments_outermost_first(
        context.model.compartments
    ):
        _make_and_add_compartment_layout(context, compartment)
    # Only compartments have been added so far, so this counts exactly the
    # foreground compartment layouts -- the insertion point for the overlay
    # background below.
    compartment_count = len(context.layout.layout_elements)
    for _key_class, species, input_species in context.activity_emissions:
        _make_and_add_species_layout(context, species, input_species)
    for gate, input_gate in context.operator_emissions:
        _make_and_add_gate_modulation_layouts(context, gate, input_gate)
    for modulation in context.model.modulations:
        if not isinstance(modulation.source, momapy.celldesigner.BooleanLogicGate):
            _make_and_add_modulation_layout(context, modulation)

    # Overlay = the plain foreground built above + the input map's remaining
    # glyphs cloned in as dimmed, unmapped background. The background carries
    # PD context for rendering only; being unmapped, the model-driven writer
    # drops it, so overlay round-trips identically to plain.
    if context.layout_mode == "overlay":
        foreground = list(context.layout.layout_elements)
        _add_dimmed_background(context, foreground, compartment_count)
        context.layout = pd2af.building.layout.highlight_layout_elements(
            foreground, context.layout
        )
    pd2af.building.layout.harmonize_root_layout(context.layout)


def _make_and_add_compartment_layout(
    context: pd2af.building.context.BuilderContext, compartment: typing.Any
):
    input_layouts = pd2af.building.layout.get_layout_elements_sorted_by_id(
        context.input_map.layout_model_mapping, compartment
    )
    if not input_layouts:
        return
    context.layout.layout_elements.extend(input_layouts)
    for input_layout in input_layouts:
        add_mappings_for_layout_and_descendants(context, input_layout)


def _get_input_container(
    context: pd2af.building.context.BuilderContext, input_layout: typing.Any
) -> typing.Any:
    for top_level_layout in context.input_map.layout.layout_elements:
        if top_level_layout is input_layout:
            return input_layout.container
        for descendant in top_level_layout.descendants():
            if descendant is input_layout:
                return top_level_layout.container
    return None


def _get_mapped_container(
    context: pd2af.building.context.BuilderContext, container: typing.Any
) -> typing.Any:
    if container is None:
        return None
    if isinstance(
        context.layout_model_mapping.get_mapping(container),
        momapy.celldesigner.Compartment,
    ):
        return container
    return None


def _add_species_layouts(
    context: pd2af.building.context.BuilderContext,
    species: typing.Any,
    input_layouts: typing.Any,
):
    placed_layouts = []
    for input_layout in input_layouts:
        original_to_placed = {}
        if context.keep_species:
            placed_layout = input_layout
            for original in [input_layout] + list(input_layout.descendants()):
                original_to_placed[id(original)] = original
        else:
            placed_layout = make_decoration_stripped_layout(
                input_layout, original_to_placed
            )
        container = _get_mapped_container(
            context, _get_input_container(context, input_layout)
        )
        if container is not placed_layout.container:
            placed_layout = dataclasses.replace(placed_layout, container=container)
            original_to_placed[id(input_layout)] = placed_layout
        add_mappings_for_layout_and_descendants(
            context, input_layout, original_to_placed
        )
        placed_layouts.append(placed_layout)
    context.layout.layout_elements.extend(placed_layouts)
    context.model_element_to_layout_elements[id(species)] = tuple(placed_layouts)


def _make_and_add_species_layout(
    context: pd2af.building.context.BuilderContext,
    species: typing.Any,
    input_species: typing.Any,
):
    input_layouts = (
        pd2af.building.layout.get_layout_elements_sorted_by_id(
            context.input_map.layout_model_mapping, input_species
        )
        if input_species is not None
        else None
    )

    if input_layouts:
        if context.layout_mode == "dot" and len(input_layouts) > 1:
            # Auto repositions every node, so a cloned species' extra glyphs
            # carry no spatial meaning -- and the cross-product pairing would
            # multiply each influence arc N*M. Keep one glyph (its real size and
            # subunit structure) and let graphviz place it. plain/overlay keep
            # every clone, where the curated positions are meaningful.
            input_layouts = input_layouts[:1]
        _add_species_layouts(context, species, input_layouts)
    elif context.layout_mode == "dot":
        synthetic_layout = make_synthetic_layout(species, context.synthetic_index)
        context.synthetic_index += 1
        context.layout.layout_elements.append(synthetic_layout)
        context.layout_model_mapping.add_mapping(synthetic_layout, species)
        context.model_element_to_layout_elements[id(species)] = (synthetic_layout,)


def _make_and_add_gate_modulation_layouts(
    context: pd2af.building.context.BuilderContext,
    gate: typing.Any,
    input_gate: typing.Any,
):
    """Draw the modulations ``gate`` is the source of, each arc from its own gate copy.

    A CellDesigner gate has one output, so a gate glyph is the source of exactly
    one modulation arc. When the gate has several arcs, each copy is placed on
    the line from the centroid of its inputs toward its target, at the distance
    the input gate glyph is from that centroid but no further than half way.
    """
    input_gate_layouts = (
        pd2af.building.layout.get_layout_elements_sorted_by_id(
            context.input_map.layout_model_mapping, input_gate
        )
        if input_gate is not None
        else []
    )
    if not input_gate_layouts:
        if context.layout_mode != "dot":
            return
        input_gate_layouts = [make_synthetic_gate_layout(gate, context.synthetic_index)]
        context.synthetic_index += 1
    prefer_nearest = context.influence_pairing == "nearest" and context.layout_mode in (
        "plain",
        "overlay",
    )
    modulation_target_layout_pairs = []
    for modulation in context.model.modulations:
        if modulation.source is not gate:
            continue
        target_layouts = context.model_element_to_layout_elements.get(
            id(modulation.target)
        )
        if not target_layouts:
            continue
        if prefer_nearest:
            _, target_layout = pd2af.building.layout.nearest_layout_pair(
                input_gate_layouts, target_layouts
            )
            target_layouts = [target_layout]
        for target_layout in target_layouts:
            modulation_target_layout_pairs.append((modulation, target_layout))
    used_positions = []
    for copy_index, (modulation, target_layout) in enumerate(
        modulation_target_layout_pairs
    ):
        gate_input_layout_pairs = _get_gate_input_layout_pairs(
            context, gate, target_layout
        )
        reference_layout = pd2af.building.layout.get_nearest_layout_element(
            input_gate_layouts, target_layout.position
        )
        if len(modulation_target_layout_pairs) == 1:
            gate_layout = reference_layout
        else:
            if context.layout_mode == "dot":
                position = momapy.geometry.Point(float(context.synthetic_index), 0.0)
                context.synthetic_index += 1
            else:
                position = _make_gate_copy_position(
                    reference_layout,
                    gate_input_layout_pairs,
                    target_layout,
                    used_positions,
                )
            used_positions.append(position)
            gate_layout = dataclasses.replace(
                reference_layout,
                id_=f"{reference_layout.id_}_{copy_index}",
                position=position,
            )
        context.layout.layout_elements.append(gate_layout)
        _add_logic_arcs_and_gate_mapping(
            context, gate, gate_layout, gate_input_layout_pairs
        )
        arc = make_modulation_arc(context, modulation, gate_layout, target_layout)
        context.layout.layout_elements.append(arc)
        add_modulation_mapping(
            context.layout_model_mapping,
            arc,
            gate_layout,
            target_layout,
            modulation,
        )


def _get_gate_input_layout_pairs(
    context: pd2af.building.context.BuilderContext,
    gate: typing.Any,
    target_layout: typing.Any,
) -> list:
    """Return, for each input of ``gate``, its glyph nearest ``target_layout``."""
    gate_input_layout_pairs = []
    for gate_input in sorted(gate.inputs, key=lambda gate_input: gate_input.id_):
        species_layouts = context.model_element_to_layout_elements.get(
            id(gate_input.referred_element)
        )
        if not species_layouts:
            continue
        gate_input_layout_pairs.append(
            (
                gate_input,
                pd2af.building.layout.get_nearest_layout_element(
                    species_layouts, target_layout.position
                ),
            )
        )
    return gate_input_layout_pairs


def _make_gate_copy_position(
    reference_layout: typing.Any,
    gate_input_layout_pairs: list,
    target_layout: typing.Any,
    used_positions: list,
) -> typing.Any:
    """Return the position of one gate copy among several.

    It is on the line from the centroid of the inputs toward the target, at the
    distance ``reference_layout`` is from that centroid but no further than half
    way to the target, moved right by one gate width until no other copy of the
    same gate is there.
    """
    position = reference_layout.position
    if gate_input_layout_pairs:
        input_layouts = [layout for _, layout in gate_input_layout_pairs]
        centroid = momapy.geometry.Point(
            sum(layout.position.x for layout in input_layouts) / len(input_layouts),
            sum(layout.position.y for layout in input_layouts) / len(input_layouts),
        )
        target_distance = math.dist(
            centroid.to_tuple(), target_layout.position.to_tuple()
        )
        if target_distance > 0:
            reference_distance = math.dist(
                centroid.to_tuple(), reference_layout.position.to_tuple()
            )
            ratio = min(reference_distance / target_distance, 0.5)
            position = momapy.geometry.Point(
                centroid.x + (target_layout.position.x - centroid.x) * ratio,
                centroid.y + (target_layout.position.y - centroid.y) * ratio,
            )
    while position in used_positions:
        position = momapy.geometry.Point(
            position.x + reference_layout.width, position.y
        )
    return position


def _add_logic_arcs_and_gate_mapping(
    context: pd2af.building.context.BuilderContext,
    gate: typing.Any,
    gate_layout: typing.Any,
    gate_input_layout_pairs: list,
):
    """Draw a logic arc from ``gate_layout`` to each input, and map them.

    Each logic arc is mapped to its gate input together with the key of the
    input glyph, and the gate copy to the gate together with the keys of its
    logic arcs, the copy being the representative.
    """
    gate_key = frozenset([gate_layout])
    for gate_input, input_layout in gate_input_layout_pairs:
        logic_arc = make_logic_arc(context, gate_layout, input_layout)
        context.layout.layout_elements.append(logic_arc)
        logic_arc_key = frozenset(
            [logic_arc]
        ) | pd2af.building.layout.get_key_of_layout_element(
            context.layout_model_mapping, input_layout
        )
        context.layout_model_mapping.add_mapping(
            logic_arc_key, gate_input, representative=logic_arc
        )
        gate_key |= logic_arc_key
    context.layout_model_mapping.add_mapping(gate_key, gate, representative=gate_layout)


def _make_and_add_modulation_layout(
    context: pd2af.building.context.BuilderContext, modulation: typing.Any
):
    source_layouts = context.model_element_to_layout_elements.get(id(modulation.source))
    target_layouts = context.model_element_to_layout_elements.get(id(modulation.target))
    if not source_layouts or not target_layouts:
        return
    # "nearest" collapses the N*M fan-out to the single closest pair, but only
    # where positions are real (plain/overlay); in dot they are throwaway
    # placeholders that graphviz overwrites, so the cross product is kept.
    prefer_nearest = context.influence_pairing == "nearest" and context.layout_mode in (
        "plain",
        "overlay",
    )
    for source_layout, target_layout in pd2af.building.layout.influence_layout_pairs(
        source_layouts, target_layouts, prefer_nearest
    ):
        arc = make_modulation_arc(context, modulation, source_layout, target_layout)
        context.layout.layout_elements.append(arc)
        add_modulation_mapping(
            context.layout_model_mapping,
            arc,
            source_layout,
            target_layout,
            modulation,
        )


def _add_dimmed_background(
    context: pd2af.building.context.BuilderContext,
    foreground: typing.Any,
    insert_index: int,
):
    """Clone the input layout's remaining glyphs into ``context.layout``.

    They go in as unmapped background, for the dimmer to grey out.

    ``foreground`` is the set of layout elements built by the plain path
    above (input objects, or copies keeping their ``id_``). Any input subtree
    already represented there -- a top-level glyph reused verbatim, or a
    promoted subunit lifted out of a dissolved complex -- is pruned from the
    clones, so the background never duplicates a foreground glyph nor
    collides with its ``id_`` in the dimming selector.

    momapy draws ``layout_elements`` in list order (later elements paint on
    top), so the clones are inserted at ``insert_index`` -- the count of
    foreground compartments, which are the leading elements -- to land *after*
    them but *before* the rest of the foreground. This three-layer z-order is
    required: compartments render an opaque white interior that would hide the
    background if it sat behind them, while the foreground species and freshly
    built influence arcs must stay on top of the background.
    """
    foreground_ids = set()
    for layout_element in foreground:
        foreground_ids.add(layout_element.id_)
        for descendant in layout_element.descendants():
            foreground_ids.add(descendant.id_)
    background_clones = []
    for input_layout_element in context.input_map.layout.layout_elements:
        background_clone = clone_layout_pruning_foreground(
            input_layout_element, foreground_ids, context.object_to_builder
        )
        if background_clone is not None:
            if hasattr(background_clone, "container"):
                background_clone.container = None
            background_clones.append(background_clone)
    context.layout.layout_elements[insert_index:insert_index] = background_clones
