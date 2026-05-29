"""Coordinator: drive the two-phase BuilderContext pipeline.

Pass 1 (``_make_and_add_model``) walks clingo activity / influence atoms
and populates ``context.model`` with canonical, content-deduped
compartments, templates, species and modulations. References between
elements are wired to canonical instances at construction time.

Pass 2 (``_make_and_add_layout``) -- skipped entirely when
``layout_mode is None`` -- populates ``context.layout`` and
``context.layout_model_mapping``, branching on ``layout_mode`` to call
leaf primitives from :mod:`pd2af.layouts`.

The orchestrator (``build_map``) assembles the final CellDesignerMap
from the three context slots and returns the frozen map.
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


@dataclasses.dataclass
class BuilderContext:
    # --- inputs ---
    input_map: object
    layout_mode: str | None
    clingo_id_to_model_element: dict

    # --- outputs being built ---
    model: object = None
    layout: object = None
    layout_model_mapping: object = None

    # --- Pass-1 -> Pass-2 handoff ---
    species_emissions: list = dataclasses.field(default_factory=list)

    # --- Pass-1 scratch ---
    cache: dict = dataclasses.field(default_factory=dict)
    subunit_to_top_level: dict = None
    activity_atoms_by_key_class: dict = dataclasses.field(default_factory=dict)
    influence_atoms: list = dataclasses.field(default_factory=list)
    key_to_species: dict = dataclasses.field(default_factory=dict)

    # --- Pass-2 scratch ---
    model_element_to_layout_elements: dict = dataclasses.field(default_factory=dict)
    object_to_builder: dict = dataclasses.field(default_factory=dict)
    kept: set = dataclasses.field(default_factory=set)
    synthetic_index: int = 0


def build_map(map_, layout_mode, clingo_model, clingo_id_to_model_element):
    context = BuilderContext(
        input_map=map_,
        layout_mode=layout_mode,
        clingo_id_to_model_element=clingo_id_to_model_element,
    )
    _make_and_add_model(context, clingo_model)
    if layout_mode is not None:
        _make_and_add_layout(context)

    map_builder_class = momapy.builder.get_or_make_builder_cls(
        momapy.celldesigner.CellDesignerMap
    )
    map_builder = map_builder_class(
        model=context.model,
        layout=context.layout,
        layout_model_mapping=context.layout_model_mapping,
    )
    new_map = momapy.builder.object_from_builder(map_builder)

    if layout_mode == "auto":
        new_map = pd2af.utils.auto_layout(new_map)
    return new_map


# ---------------------------------------------------------------------------
# Pass 1 -- model construction
# ---------------------------------------------------------------------------


def _make_and_add_model(context, clingo_model):
    context.model = momapy.builder.get_or_make_builder_cls(
        momapy.celldesigner.CellDesignerModel
    )()
    context.subunit_to_top_level = pd2af.model.build_subunit_to_top_level(
        context.input_map
    )
    _collect_ingredients(context, clingo_model)
    _make_and_add_compartments(context)
    _make_and_add_templates(context)
    _make_and_add_species(context)
    _make_and_add_modulations(context)


def _collect_ingredients(context, clingo_model):
    context.activity_atoms_by_key_class = {
        key_class: [] for key_class in _SPECIES_LAYER_ORDER
    }
    for atom in clingo_model.query(pd2af.predicates.new).all():
        payload = atom.object_
        if isinstance(payload, pd2af.predicates.activity):
            context.activity_atoms_by_key_class[type(payload.key)].append(payload)
        elif isinstance(
            payload,
            (
                pd2af.predicates.positivelyInfluences,
                pd2af.predicates.negativelyInfluences,
            ),
        ):
            context.influence_atoms.append(payload)


def _make_and_add_compartments(context):
    immediate_compartments = set()
    for key_class in _SPECIES_LAYER_ORDER:
        for atom in context.activity_atoms_by_key_class[key_class]:
            input_species = context.clingo_id_to_model_element[atom.key.species]
            compartment = _compartment_for_input_species(context, input_species)
            if compartment is not None:
                immediate_compartments.add(compartment)
    for compartment in pd2af.model.compartments_outermost_first(
        pd2af.model.collect_ancestor_compartments(immediate_compartments)
    ):
        context.model.compartments.add(compartment)


def _compartment_for_input_species(context, input_species):
    if getattr(input_species, "compartment", None) is not None:
        return input_species.compartment
    return pd2af.model.get_parent_complex_compartment(
        input_species, context.subunit_to_top_level
    )


def _make_and_add_templates(context):
    seen_template_identities = set()

    def record(template):
        if id(template) in seen_template_identities:
            return
        seen_template_identities.add(id(template))
        context.model.species_templates.add(template)

    # Kept templates first: walk each kept-key input species's subunit
    # tree and register every template it carries by identity. This
    # primes the cache so that a content-equal stripped candidate
    # collapses onto its kept twin in the next phase.
    for key_class in _KEPT_KEY_CLASSES:
        for atom in context.activity_atoms_by_key_class[key_class]:
            input_species = context.clingo_id_to_model_element[atom.key.species]
            for input_template in _walk_templates(input_species):
                canonical = pd2af.model.register_or_reuse(
                    input_template, context.cache
                )
                record(canonical)

    # Stripped templates: one per content cell. Built from each
    # new_species_from_template atom's input template; the species
    # layer reads back the canonical via the cache.
    for atom in context.activity_atoms_by_key_class[
        pd2af.predicates.new_species_from_template
    ]:
        input_species = context.clingo_id_to_model_element[atom.key.species]
        input_template = input_species.template
        if input_template is None:
            raise ValueError(
                f"new_species_from_template key references species "
                f"{input_species.id_!r} which has no template"
            )
        canonical = pd2af.model.get_or_make_stripped_template(
            input_template, context.cache
        )
        record(canonical)


def _walk_templates(species):
    template = getattr(species, "template", None)
    if template is not None:
        yield template
    for subunit in getattr(species, "subunits", ()) or ():
        yield from _walk_templates(subunit)


def _make_and_add_species(context):
    seen_species_identities = set()
    for key_class in _SPECIES_LAYER_ORDER:
        for atom in context.activity_atoms_by_key_class[key_class]:
            species = _resolve_activity_key(context, atom.key)
            context.key_to_species[atom.key] = species
            if id(species) in seen_species_identities:
                continue
            seen_species_identities.add(id(species))
            context.species_emissions.append((key_class, species))
            if key_class is not pd2af.predicates.kept_subunit:
                context.model.species.add(species)


def _resolve_activity_key(context, key):
    if isinstance(key, _KEPT_KEY_CLASSES):
        input_species = context.clingo_id_to_model_element[key.species]
        return pd2af.model.get_or_make_kept_species(input_species, context.cache)
    if isinstance(key, pd2af.predicates.new_species_from_template):
        input_species = context.clingo_id_to_model_element[key.species]
        compartment = _compartment_for_input_species(context, input_species)
        stripped_template = pd2af.model.get_or_make_stripped_template(
            input_species.template, context.cache
        )
        return pd2af.model.get_or_make_synthesized_species(
            input_species, stripped_template, compartment, context.cache
        )
    raise ValueError(f"unknown activity key wrapper {type(key).__name__}")


def _make_and_add_modulations(context):
    seen_modulation_identities = set()
    for atom in context.influence_atoms:
        modulation_class = pd2af.predicates.predicate_to_model_element_class[
            type(atom)
        ]
        source = context.key_to_species[atom.source]
        target = context.key_to_species[atom.target]
        canonical = pd2af.model.get_or_make_modulation(
            modulation_class, source, target, context.cache
        )
        if id(canonical) in seen_modulation_identities:
            continue
        seen_modulation_identities.add(id(canonical))
        context.model.modulations.add(canonical)


# ---------------------------------------------------------------------------
# Pass 2 -- layout construction
# ---------------------------------------------------------------------------


def _make_and_add_layout(context):
    if context.layout_mode in ("plain", "auto"):
        context.layout, context.layout_model_mapping = (
            pd2af.layouts.new_layout_and_mapping_builders()
        )
    elif context.layout_mode == "overlay":
        context.layout, context.layout_model_mapping = (
            pd2af.layouts.clone_input_layout_and_mapping(
                context.input_map, context.object_to_builder
            )
        )

    for compartment in pd2af.model.compartments_outermost_first(
        context.model.compartments
    ):
        _make_and_add_compartment_layout(context, compartment)
    for key_class, species in context.species_emissions:
        _make_and_add_species_layout(context, key_class, species)
    for modulation in context.model.modulations:
        _make_and_add_modulation_layout(context, modulation)

    pd2af.utils.harmonize_root_layout(context.layout)
    if context.layout_mode == "overlay":
        context.layout = pd2af.utils.highlight_layout_elements(
            context.kept, context.layout
        )


def _make_and_add_compartment_layout(context, compartment):
    input_layouts = context.input_map.layout_model_mapping.get_mapping(compartment)
    if not input_layouts:
        return
    if context.layout_mode in ("plain", "auto"):
        context.layout.layout_elements.extend(input_layouts)
        for input_layout in input_layouts:
            pd2af.layouts.copy_subtree_mappings(
                context.input_map.layout_model_mapping,
                input_layout,
                context.layout_model_mapping,
            )
    elif context.layout_mode == "overlay":
        for input_layout in input_layouts:
            clone = context.object_to_builder.get(id(input_layout))
            if clone is None:
                raise ValueError(
                    f"compartment layout {input_layout!r} has no clone in "
                    f"object_to_builder"
                )
            context.kept.add(clone)


def _make_and_add_species_layout(context, key_class, species):
    is_kept_subunit = key_class is pd2af.predicates.kept_subunit
    has_input = key_class in _KEPT_KEY_CLASSES
    input_layouts = (
        context.input_map.layout_model_mapping.get_mapping(species)
        if has_input
        else None
    )

    if context.layout_mode in ("plain", "auto"):
        if input_layouts:
            if not is_kept_subunit:
                context.layout.layout_elements.extend(input_layouts)
                for input_layout in input_layouts:
                    pd2af.layouts.copy_subtree_mappings(
                        context.input_map.layout_model_mapping,
                        input_layout,
                        context.layout_model_mapping,
                    )
            context.model_element_to_layout_elements[id(species)] = tuple(
                input_layouts
            )
        elif context.layout_mode == "auto" and not is_kept_subunit:
            synthetic_layout = pd2af.layouts.make_synthetic_layout(
                species, context.synthetic_index
            )
            context.synthetic_index += 1
            context.layout.layout_elements.append(synthetic_layout)
            context.layout_model_mapping.add_mapping(synthetic_layout, species)
            context.model_element_to_layout_elements[id(species)] = (
                synthetic_layout,
            )
    elif context.layout_mode == "overlay":
        if input_layouts:
            for input_layout in input_layouts:
                clone = context.object_to_builder.get(id(input_layout))
                if clone is None:
                    raise ValueError(
                        f"species layout {input_layout!r} has no clone in "
                        f"object_to_builder"
                    )
                context.kept.add(clone)
            context.model_element_to_layout_elements[id(species)] = tuple(
                input_layouts
            )


def _make_and_add_modulation_layout(context, modulation):
    source_layouts = context.model_element_to_layout_elements.get(
        id(modulation.source)
    )
    target_layouts = context.model_element_to_layout_elements.get(
        id(modulation.target)
    )
    if not source_layouts or not target_layouts:
        return
    for source_layout, target_layout in itertools.product(
        source_layouts, target_layouts
    ):
        if context.layout_mode in ("plain", "auto"):
            arc = pd2af.layouts.make_modulation_arc(
                modulation, source_layout, target_layout
            )
            context.layout.layout_elements.append(arc)
            pd2af.layouts.add_modulation_mapping(
                context.layout_model_mapping,
                arc,
                source_layout,
                target_layout,
                modulation,
            )
        elif context.layout_mode == "overlay":
            arc_builder, source_builder, target_builder = (
                pd2af.layouts.make_overlay_modulation_arc(
                    context, modulation, source_layout, target_layout
                )
            )
            context.layout.layout_elements.append(arc_builder)
            context.kept.add(arc_builder)
            pd2af.layouts.add_modulation_mapping(
                context.layout_model_mapping,
                arc_builder,
                source_builder,
                target_builder,
                modulation,
            )
