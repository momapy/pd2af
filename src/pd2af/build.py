"""Coordinator: drive the canonicalisation + builder-append pipeline.

Two passes. Pass 1 (canonicalisation) walks clingo activity / influence
atoms and produces canonical, content-deduped lists of compartments,
templates, species, and modulations. References between elements are
wired to canonical instances at construction time, so Pass 2 never
encounters a non-canonical reference — no post-pass remap is needed.

Pass 2 (builder-append + layout) writes the canonical elements into a
fresh model builder and, in the same loop, constructs the layout tree
and the layout-model mapping inline; ``layout_mode`` drives an inline
branch (``None`` / ``"plain"`` / ``"auto"`` / ``"overlay"``) calling
leaf primitives from :mod:`pd2af.layouts`.

Pass 1 follows the layered pipeline order from
``plans/model-layout-split.md``:

  (1) ingredients pass
  (2) canonical compartments
  (3) canonical templates (kept first, then stripped)
  (4) canonical species (kept layers first, then synthesized)
  (5) canonical modulations

Each layer's references point at already-canonical instances by the
time the next layer reads them, so no cross-layer remap is needed.
"""

import dataclasses
import itertools

import momapy.builder
import momapy.celldesigner

import pd2af.layouts
import pd2af.model
import pd2af.predicates
import pd2af.utils


_SPECIES_LAYER_ORDER = (
    pd2af.predicates.kept_species,
    pd2af.predicates.kept_subunit,
    pd2af.predicates.promoted_subunit,
    pd2af.predicates.new_species_from_template,
)

_KEPT_KEY_CLASSES = (
    pd2af.predicates.kept_species,
    pd2af.predicates.kept_subunit,
    pd2af.predicates.promoted_subunit,
)


@dataclasses.dataclass(frozen=True)
class _Ingredients:
    """Output of the ingredients pass.

    ``activity_atoms_by_key_class`` groups activity atoms in layer order.
    ``influence_atoms`` is the raw list of influence atoms; the
    modulation layer dispatches on each one's predicate type.
    """

    activity_atoms_by_key_class: dict
    influence_atoms: list


def build_map(map_, layout_mode, clingo_model, clingo_id_to_model_element):
    ingredients = _collect_ingredients(clingo_model)
    subunit_to_top_level = pd2af.model.build_subunit_to_top_level(map_)
    cache = {}

    compartments = _canonicalise_compartments(
        ingredients, clingo_id_to_model_element, subunit_to_top_level
    )
    templates = _canonicalise_templates(
        ingredients, clingo_id_to_model_element, cache
    )
    species_emissions, key_to_species = _canonicalise_species(
        ingredients, clingo_id_to_model_element, subunit_to_top_level, cache
    )
    modulations = _canonicalise_modulations(
        ingredients, key_to_species, cache
    )

    return _build_and_drive_layout(
        map_,
        layout_mode,
        compartments,
        templates,
        species_emissions,
        modulations,
    )


# ---------------------------------------------------------------------------
# (1) ingredients pass


def _collect_ingredients(clingo_model):
    activity_atoms_by_key_class = {
        key_class: [] for key_class in _SPECIES_LAYER_ORDER
    }
    influence_atoms = []
    for atom in clingo_model.query(pd2af.predicates.new).all():
        payload = atom.object_
        if isinstance(payload, pd2af.predicates.activity):
            activity_atoms_by_key_class[type(payload.key)].append(payload)
        elif isinstance(
            payload,
            (
                pd2af.predicates.positivelyInfluences,
                pd2af.predicates.negativelyInfluences,
            ),
        ):
            influence_atoms.append(payload)
    return _Ingredients(
        activity_atoms_by_key_class=activity_atoms_by_key_class,
        influence_atoms=influence_atoms,
    )


# ---------------------------------------------------------------------------
# (2) canonical compartments


def _canonicalise_compartments(
    ingredients, clingo_id_to_model_element, subunit_to_top_level
):
    immediate_compartments = set()
    for key_class in _SPECIES_LAYER_ORDER:
        for atom in ingredients.activity_atoms_by_key_class[key_class]:
            input_species = clingo_id_to_model_element[atom.key.species]
            compartment = _compartment_for_input_species(
                input_species, subunit_to_top_level
            )
            if compartment is not None:
                immediate_compartments.add(compartment)
    return pd2af.model.compartments_outermost_first(
        pd2af.model.collect_ancestor_compartments(immediate_compartments)
    )


def _compartment_for_input_species(input_species, subunit_to_top_level):
    if getattr(input_species, "compartment", None) is not None:
        return input_species.compartment
    return pd2af.model.get_parent_complex_compartment(
        input_species, subunit_to_top_level
    )


# ---------------------------------------------------------------------------
# (3) canonical templates


def _canonicalise_templates(ingredients, clingo_id_to_model_element, cache):
    canonical_templates = []
    seen_template_identities = set()

    def record(template):
        if id(template) in seen_template_identities:
            return
        seen_template_identities.add(id(template))
        canonical_templates.append(template)

    # Kept templates first: walk each kept-key input species's subunit
    # tree and register every template it carries by identity. This
    # primes the cache so that a content-equal stripped candidate
    # collapses onto its kept twin in the next phase.
    for key_class in _KEPT_KEY_CLASSES:
        for atom in ingredients.activity_atoms_by_key_class[key_class]:
            input_species = clingo_id_to_model_element[atom.key.species]
            for input_template in _walk_templates(input_species):
                canonical = pd2af.model.register_or_reuse(
                    input_template, cache
                )
                record(canonical)

    # Stripped templates: one per content cell. Built from each
    # new_species_from_template atom's input template; the species
    # layer reads back the canonical via the cache.
    for atom in ingredients.activity_atoms_by_key_class[
        pd2af.predicates.new_species_from_template
    ]:
        input_species = clingo_id_to_model_element[atom.key.species]
        input_template = input_species.template
        if input_template is None:
            raise ValueError(
                f"new_species_from_template key references species "
                f"{input_species.id_!r} which has no template"
            )
        canonical = pd2af.model.get_or_make_stripped_template(
            input_template, cache
        )
        record(canonical)

    return canonical_templates


def _walk_templates(species):
    template = getattr(species, "template", None)
    if template is not None:
        yield template
    for subunit in getattr(species, "subunits", ()) or ():
        yield from _walk_templates(subunit)


# ---------------------------------------------------------------------------
# (4) canonical species


def _canonicalise_species(
    ingredients, clingo_id_to_model_element, subunit_to_top_level, cache
):
    species_emissions = []
    key_to_species = {}
    seen_species_identities = set()
    for key_class in _SPECIES_LAYER_ORDER:
        for atom in ingredients.activity_atoms_by_key_class[key_class]:
            species = _resolve_activity_key(
                atom.key,
                clingo_id_to_model_element,
                subunit_to_top_level,
                cache,
            )
            key_to_species[atom.key] = species
            if id(species) in seen_species_identities:
                continue
            seen_species_identities.add(id(species))
            species_emissions.append((key_class, species))
    return species_emissions, key_to_species


def _resolve_activity_key(
    key, clingo_id_to_model_element, subunit_to_top_level, cache
):
    if isinstance(key, _KEPT_KEY_CLASSES):
        input_species = clingo_id_to_model_element[key.species]
        return pd2af.model.get_or_make_kept_species(input_species, cache)
    if isinstance(key, pd2af.predicates.new_species_from_template):
        input_species = clingo_id_to_model_element[key.species]
        compartment = _compartment_for_input_species(
            input_species, subunit_to_top_level
        )
        stripped_template = pd2af.model.get_or_make_stripped_template(
            input_species.template, cache
        )
        return pd2af.model.get_or_make_synthesized_species(
            input_species, stripped_template, compartment, cache
        )
    raise ValueError(f"unknown activity key wrapper {type(key).__name__}")


# ---------------------------------------------------------------------------
# (5) canonical modulations


def _canonicalise_modulations(ingredients, key_to_species, cache):
    modulations = []
    seen_modulation_identities = set()
    for atom in ingredients.influence_atoms:
        modulation_class = pd2af.predicates.predicate_to_model_element_class[
            type(atom)
        ]
        source = key_to_species[atom.source]
        target = key_to_species[atom.target]
        canonical = pd2af.model.get_or_make_modulation(
            modulation_class, source, target, cache
        )
        if id(canonical) in seen_modulation_identities:
            continue
        seen_modulation_identities.add(id(canonical))
        modulations.append(canonical)
    return modulations


# ---------------------------------------------------------------------------
# Pass 2 — builder-append + layout construction


@dataclasses.dataclass
class _LayoutState:
    map_: object
    layout_builder: object = None
    mapping_builder: object = None
    model_element_to_layout_elements: dict = dataclasses.field(
        default_factory=dict
    )
    object_to_builder: dict = dataclasses.field(default_factory=dict)
    kept: set = dataclasses.field(default_factory=set)
    synthetic_index: int = 0


def _build_and_drive_layout(
    map_,
    layout_mode,
    compartments,
    templates,
    species_emissions,
    modulations,
):
    model_builder = momapy.builder.get_or_make_builder_cls(
        momapy.celldesigner.CellDesignerModel
    )()
    state = _LayoutState(map_=map_)

    # Phase 0 — set up layout + mapping builders.
    if layout_mode in ("plain", "auto"):
        state.layout_builder, state.mapping_builder = (
            pd2af.layouts.new_layout_and_mapping_builders()
        )
    elif layout_mode == "overlay":
        state.layout_builder, state.mapping_builder = (
            pd2af.layouts.clone_input_layout_and_mapping(
                map_, state.object_to_builder
            )
        )

    # Phase 1 — compartments.
    for compartment in compartments:
        model_builder.compartments.add(compartment)
        if layout_mode is None:
            continue
        input_layouts = map_.layout_model_mapping.get_mapping(compartment)
        if not input_layouts:
            continue
        if layout_mode in ("plain", "auto"):
            state.layout_builder.layout_elements.extend(input_layouts)
            for input_layout in input_layouts:
                pd2af.layouts.copy_subtree_mappings(
                    map_.layout_model_mapping,
                    input_layout,
                    state.mapping_builder,
                )
        elif layout_mode == "overlay":
            for input_layout in input_layouts:
                clone = state.object_to_builder.get(id(input_layout))
                if clone is None:
                    raise ValueError(
                        f"compartment layout {input_layout!r} has no clone in "
                        f"object_to_builder"
                    )
                state.kept.add(clone)

    # Phase 2 — templates.
    for template in templates:
        model_builder.species_templates.add(template)

    # Phase 3 — species.
    for key_class, species in species_emissions:
        is_kept_subunit = key_class is pd2af.predicates.kept_subunit
        if not is_kept_subunit:
            model_builder.species.add(species)
        if layout_mode is None:
            continue

        has_input = key_class in _KEPT_KEY_CLASSES
        input_layouts = (
            map_.layout_model_mapping.get_mapping(species)
            if has_input
            else None
        )

        if layout_mode in ("plain", "auto"):
            if input_layouts:
                if not is_kept_subunit:
                    state.layout_builder.layout_elements.extend(input_layouts)
                    for input_layout in input_layouts:
                        pd2af.layouts.copy_subtree_mappings(
                            map_.layout_model_mapping,
                            input_layout,
                            state.mapping_builder,
                        )
                state.model_element_to_layout_elements[id(species)] = tuple(
                    input_layouts
                )
            elif layout_mode == "auto" and not is_kept_subunit:
                synthetic_layout = pd2af.layouts.make_synthetic_layout(
                    species, state.synthetic_index
                )
                state.synthetic_index += 1
                state.layout_builder.layout_elements.append(synthetic_layout)
                state.mapping_builder.add_mapping(synthetic_layout, species)
                state.model_element_to_layout_elements[id(species)] = (
                    synthetic_layout,
                )
        elif layout_mode == "overlay":
            if input_layouts:
                for input_layout in input_layouts:
                    clone = state.object_to_builder.get(id(input_layout))
                    if clone is None:
                        raise ValueError(
                            f"species layout {input_layout!r} has no clone in "
                            f"object_to_builder"
                        )
                    state.kept.add(clone)
                state.model_element_to_layout_elements[id(species)] = tuple(
                    input_layouts
                )

    # Phase 4 — modulations.
    for modulation in modulations:
        model_builder.modulations.add(modulation)
        if layout_mode is None:
            continue
        source_layouts = state.model_element_to_layout_elements.get(
            id(modulation.source)
        )
        target_layouts = state.model_element_to_layout_elements.get(
            id(modulation.target)
        )
        if not source_layouts or not target_layouts:
            continue
        for source_layout, target_layout in itertools.product(
            source_layouts, target_layouts
        ):
            if layout_mode in ("plain", "auto"):
                arc = pd2af.layouts.make_modulation_arc(
                    modulation, source_layout, target_layout
                )
                state.layout_builder.layout_elements.append(arc)
                pd2af.layouts.add_modulation_mapping(
                    state.mapping_builder,
                    arc,
                    source_layout,
                    target_layout,
                    modulation,
                )
            elif layout_mode == "overlay":
                # source_builder / target_builder are clones already in
                # state.kept (added by Phase 1/3 against the shared
                # object_to_builder). Only the new arc_builder needs to be
                # marked kept here.
                arc_builder, source_builder, target_builder = (
                    pd2af.layouts.make_overlay_modulation_arc(
                        state, modulation, source_layout, target_layout
                    )
                )
                state.layout_builder.layout_elements.append(arc_builder)
                state.kept.add(arc_builder)
                pd2af.layouts.add_modulation_mapping(
                    state.mapping_builder,
                    arc_builder,
                    source_builder,
                    target_builder,
                    modulation,
                )

    # Phase 5 — finish.
    if layout_mode is None:
        return pd2af.layouts.finalize_map(model_builder, None, None)

    if layout_mode == "overlay":
        state.layout_builder = pd2af.utils.highlight_layout_elements(
            state.kept, state.layout_builder
        )

    pd2af.utils.harmonize_root_layout(state.layout_builder)
    new_map = pd2af.layouts.finalize_map(
        model_builder, state.layout_builder, state.mapping_builder
    )

    if layout_mode == "auto":
        new_map = pd2af.utils.auto_layout(new_map)

    return new_map
