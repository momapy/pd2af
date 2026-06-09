"""Build the AF layout, reusing the input map's glyphs where possible.

``make_and_add_layout`` is the layout pass: it walks the model built by
:mod:`pd2af.celldesigner.building_model` and populates ``context.layout`` and
``context.layout_model_mapping``, branching on ``layout_mode``. The leaf
primitives below -- synthetic nodes, modulation arcs, mapping helpers, and
background cloning for overlay -- do the per-element construction.

:mod:`pd2af.build` owns the ``BuilderContext`` and invokes this pass after
the model pass.
"""

import momapy.builder
import momapy.celldesigner
import momapy.core.layout
import momapy.core.mapping
import momapy.geometry

import pd2af.celldesigner.building_model
import pd2af.utils


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


def new_layout_and_mapping_builders():
    layout_builder_class = momapy.builder.get_or_make_builder_cls(
        momapy.core.layout.Layout
    )
    return (
        layout_builder_class(),
        momapy.core.mapping.LayoutModelMappingBuilder(),
    )


def clone_layout_pruning_foreground(
    input_layout_element, foreground_ids, object_to_builder
):
    """Clone an input layout element as a builder for use as dimmed
    background, omitting any descendant subtree whose input original is in
    ``foreground_ids`` (already drawn in the foreground). Returns ``None``
    when the element itself belongs to the foreground.

    Clones share ``object_to_builder`` so structure shared across the input
    layout stays shared once cloned, and they are builders -- not the input
    frozen objects -- so dimming the background never mutates the input
    map's layout.
    """
    if id(input_layout_element) in foreground_ids:
        return None
    clone = momapy.builder.builder_from_object(
        input_layout_element, object_to_builder=object_to_builder
    )
    _prune_foreground_from_clone(
        input_layout_element, clone, foreground_ids, object_to_builder
    )
    return clone


def _prune_foreground_from_clone(
    input_layout_element, clone, foreground_ids, object_to_builder
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
        if id(input_subunit) in foreground_ids:
            continue
        subunit_clone = object_to_builder[id(input_subunit)]
        _prune_foreground_from_clone(
            input_subunit, subunit_clone, foreground_ids, object_to_builder
        )
        surviving_clones.append(subunit_clone)
    clone.layout_elements = surviving_clones


def make_synthetic_layout(species, index):
    """Build a placeholder Node for ``species``.

    ``index`` seeds the position so two synthetic layouts for
    content-equal species (same class, same name) are themselves
    content-distinct. Layout dataclasses have ``compare=False`` on
    ``id_``, so identical content collapses to one mapping key — and
    the layout-model mapping then drops every synthetic but the first.
    The position itself is throwaway: ``auto_layout`` repositions every
    node before render.
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


def make_modulation_arc(modulation, source_layout, target_layout):
    arc_class = _MODULATION_CLASS_TO_LAYOUT_CLASS[type(modulation)]
    if source_layout is target_layout:
        start_point = source_layout.anchor_point("north_north_west")
        end_point = source_layout.anchor_point("north_north_east")
    else:
        start_point = source_layout.own_border(target_layout.center())
        end_point = target_layout.own_border(source_layout.center())
        if start_point is None:
            start_point = source_layout.north_west()
        if end_point is None:
            end_point = target_layout.north_east()
    segment = momapy.geometry.Segment(start_point, end_point)
    return arc_class(
        source=source_layout,
        target=target_layout,
        segments=(segment,),
    )


def add_mappings_for_layout_and_descendants(
    input_layout_model_mapping,
    layout_element,
    mapping_builder,
    input_model_element_to_canonical_model_element=None,
):
    elements = [layout_element] + list(layout_element.descendants())
    for element in elements:
        if element in input_layout_model_mapping:
            model_value = input_layout_model_mapping[element]
            if input_model_element_to_canonical_model_element is not None:
                model_value = input_model_element_to_canonical_model_element.get(
                    id(model_value), model_value
                )
            mapping_builder.add_mapping(element, model_value)


def add_modulation_mapping(
    mapping_builder, arc, source_layout, target_layout, modulation
):
    source_key = mapping_builder._singleton_to_key.get(source_layout)
    source_cluster = (
        source_key if source_key is not None else frozenset([source_layout])
    )
    target_key = mapping_builder._singleton_to_key.get(target_layout)
    target_cluster = (
        target_key if target_key is not None else frozenset([target_layout])
    )
    mapping_builder.add_mapping(
        frozenset([arc]) | source_cluster | target_cluster,
        modulation,
        anchor=arc,
    )


# ---------------------------------------------------------------------------
# Pass 2 -- layout construction
# ---------------------------------------------------------------------------


def make_and_add_layout(context):
    context.layout, context.layout_model_mapping = new_layout_and_mapping_builders()

    for compartment in pd2af.celldesigner.building_model.compartments_outermost_first(
        context.model.compartments
    ):
        _make_and_add_compartment_layout(context, compartment)
    # Only compartments have been added so far, so this counts exactly the
    # foreground compartment layouts -- the insertion point for the overlay
    # background below.
    compartment_count = len(context.layout.layout_elements)
    for _key_class, species, input_species in context.species_emissions:
        _make_and_add_species_layout(context, species, input_species)
    for modulation in context.model.modulations:
        _make_and_add_modulation_layout(context, modulation)

    # Overlay = the plain foreground built above + the input map's remaining
    # glyphs cloned in as dimmed, unmapped background. The background carries
    # PD context for rendering only; being unmapped, the model-driven writer
    # drops it, so overlay round-trips identically to plain.
    if context.layout_mode == "overlay":
        foreground = list(context.layout.layout_elements)
        _add_dimmed_background(context, foreground, compartment_count)
        context.layout = pd2af.utils.highlight_layout_elements(
            foreground, context.layout
        )
    pd2af.utils.harmonize_root_layout(context.layout)


def _make_and_add_compartment_layout(context, compartment):
    input_layouts = context.input_map.layout_model_mapping.get_mapping(compartment)
    if not input_layouts:
        return
    context.layout.layout_elements.extend(input_layouts)
    for input_layout in input_layouts:
        add_mappings_for_layout_and_descendants(
            context.input_map.layout_model_mapping,
            input_layout,
            context.layout_model_mapping,
            input_model_element_to_canonical_model_element=(
                context.input_model_element_to_canonical_model_element
            ),
        )


def _make_and_add_species_layout(context, species, input_species):
    input_layouts = (
        context.input_map.layout_model_mapping.get_mapping(input_species)
        if input_species is not None
        else None
    )

    if input_layouts:
        context.layout.layout_elements.extend(input_layouts)
        for input_layout in input_layouts:
            add_mappings_for_layout_and_descendants(
                context.input_map.layout_model_mapping,
                input_layout,
                context.layout_model_mapping,
                input_model_element_to_canonical_model_element=(
                    context.input_model_element_to_canonical_model_element
                ),
            )
        context.model_element_to_layout_elements[id(species)] = tuple(input_layouts)
    elif context.layout_mode == "auto":
        synthetic_layout = make_synthetic_layout(species, context.synthetic_index)
        context.synthetic_index += 1
        context.layout.layout_elements.append(synthetic_layout)
        context.layout_model_mapping.add_mapping(synthetic_layout, species)
        context.model_element_to_layout_elements[id(species)] = (synthetic_layout,)


def _make_and_add_modulation_layout(context, modulation):
    source_layouts = context.model_element_to_layout_elements.get(
        id(modulation.source)
    )
    target_layouts = context.model_element_to_layout_elements.get(
        id(modulation.target)
    )
    if not source_layouts or not target_layouts:
        return
    # "nearest" collapses the N*M fan-out to the single closest pair, but only
    # where positions are real (plain/overlay); in auto they are throwaway
    # placeholders that graphviz overwrites, so the cross product is kept.
    prefer_nearest = (
        context.influence_pairing == "nearest"
        and context.layout_mode in ("plain", "overlay")
    )
    for source_layout, target_layout in pd2af.utils.influence_layout_pairs(
        source_layouts, target_layouts, prefer_nearest
    ):
        arc = make_modulation_arc(modulation, source_layout, target_layout)
        context.layout.layout_elements.append(arc)
        add_modulation_mapping(
            context.layout_model_mapping,
            arc,
            source_layout,
            target_layout,
            modulation,
        )


def _add_dimmed_background(context, foreground, insert_index):
    """Clone the input layout's remaining glyphs into ``context.layout`` as
    unmapped background, for the dimmer to grey out.

    ``foreground`` is the set of layout elements built by the plain path
    above (input objects shared with the input map). Any input subtree
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
        foreground_ids.add(id(layout_element))
        for descendant in layout_element.descendants():
            foreground_ids.add(id(descendant))
    background_clones = []
    for input_layout_element in context.input_map.layout.layout_elements:
        background_clone = clone_layout_pruning_foreground(
            input_layout_element, foreground_ids, context.object_to_builder
        )
        if background_clone is not None:
            background_clones.append(background_clone)
    context.layout.layout_elements[insert_index:insert_index] = background_clones
