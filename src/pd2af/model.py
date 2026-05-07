"""Resolve clingo facts into AF model elements.

Pure function — does not import :mod:`pd2af.solver`. Takes the solver's
outputs (a clingo model and the id-to-model-element registry built when
facts were generated) plus the original input map (used to walk
ancestor-compartment lookups for synthesized subunit-promoted species)
and returns a :class:`Resolution` describing the new model:
compartments, templates, species, and modulations to place.

The lockstep guarantee with :mod:`pd2af.layouts` is that the same
``Species`` instance returned in :attr:`Resolution.species` is the one
the layout strategy uses as the layout-model mapping value. Identity
is shared end-to-end, which is what the CellDesigner writer's
identity-keyed lookups (e.g. ``<listOfSpeciesAliases>``) need.
"""

import dataclasses

import momapy.celldesigner

import pd2af.predicates


_FLAT_COMPLEX_MODES = frozenset({"normal", "no-complex"})

_NO_COMPARTMENT_SENTINEL = "no_compartment"

_TEMPLATE_MERGE_PREFIX = "new_species_from_template__"

_STRIPPED_TEMPLATE_PREFIX = "merged_template__"

_TEMPLATE_TO_SPECIES_CLASS = {
    momapy.celldesigner.GenericProteinTemplate: momapy.celldesigner.GenericProtein,
    momapy.celldesigner.TruncatedProteinTemplate: momapy.celldesigner.TruncatedProtein,
    momapy.celldesigner.ReceptorTemplate: momapy.celldesigner.Receptor,
    momapy.celldesigner.IonChannelTemplate: momapy.celldesigner.IonChannel,
    momapy.celldesigner.GeneTemplate: momapy.celldesigner.Gene,
    momapy.celldesigner.RNATemplate: momapy.celldesigner.RNA,
    momapy.celldesigner.AntisenseRNATemplate: momapy.celldesigner.AntisenseRNA,
}


@dataclasses.dataclass(frozen=True)
class Resolution:
    """Plan for building the new AF map.

    ``species`` entries are ``(species, existing_species_or_None,
    is_subunit)`` triples. ``existing_species`` is the input-map species
    whose layout is reusable for the new species, or ``None`` when the
    species was synthesised (merged proteoform / promoted subunit) or
    its layout subtree no longer matches (flattened complex).
    ``is_subunit`` is ``True`` for subunits of kept (non-flattened)
    complexes — those go into the layout but not into
    ``model_builder.species``.
    """

    compartments: list
    templates: list
    species: list
    modulations: list


def resolve(clingo_model, clingo_id_to_model_element, mode, input_map):
    flatten_complexes = mode in _FLAT_COMPLEX_MODES
    activity_atoms = _get_activity_atoms(clingo_model)
    influence_atoms = _get_influence_atoms(clingo_model)
    key_to_resolution = {}
    intern_table = {}
    stripped_template_cache = {}
    subunit_to_top_level = _build_subunit_to_top_level(input_map)
    resolved_species = []
    for atom in activity_atoms:
        species, existing_species = _resolve_activity_key(
            atom.key,
            clingo_id_to_model_element,
            key_to_resolution,
            intern_table,
            stripped_template_cache,
            subunit_to_top_level,
            flatten_complexes=flatten_complexes,
        )
        resolved_species.append((species, existing_species))
    immediate_compartments = {
        species.compartment
        for species, _ in resolved_species
        if getattr(species, "compartment", None) is not None
    }
    compartments = _compartments_outermost_first(
        _collect_ancestor_compartments(immediate_compartments)
    )
    species_iterable = (species for species, _ in resolved_species)
    templates = list(_collect_templates_from_species(species_iterable))
    species_entries, canonical_species = _walk_species(
        resolved_species, flatten_complexes=flatten_complexes
    )
    modulations = []
    for atom in influence_atoms:
        modulation_class = pd2af.predicates.predicate_to_model_element_class[
            type(atom)
        ]
        source_species, _ = _resolve_activity_key(
            atom.source,
            clingo_id_to_model_element,
            key_to_resolution,
            intern_table,
            stripped_template_cache,
            subunit_to_top_level,
            flatten_complexes=flatten_complexes,
        )
        target_species, _ = _resolve_activity_key(
            atom.target,
            clingo_id_to_model_element,
            key_to_resolution,
            intern_table,
            stripped_template_cache,
            subunit_to_top_level,
            flatten_complexes=flatten_complexes,
        )
        source_species = canonical_species.get(source_species, source_species)
        target_species = canonical_species.get(target_species, target_species)
        modulations.append(
            modulation_class(source=source_species, target=target_species)
        )
    return Resolution(
        compartments=compartments,
        templates=templates,
        species=species_entries,
        modulations=modulations,
    )


def _get_activity_atoms(clingo_model):
    return [
        atom.object_
        for atom in clingo_model.query(pd2af.predicates.new).all()
        if isinstance(atom.object_, pd2af.predicates.activity)
    ]


def _get_influence_atoms(clingo_model):
    return [
        atom.object_
        for atom in clingo_model.query(pd2af.predicates.new).all()
        if isinstance(
            atom.object_,
            (
                pd2af.predicates.positivelyInfluences,
                pd2af.predicates.negativelyInfluences,
            ),
        )
    ]


def _resolve_activity_key(
    key,
    clingo_id_to_model_element,
    key_to_resolution,
    intern_table,
    stripped_template_cache,
    subunit_to_top_level,
    *,
    flatten_complexes,
):
    cached = key_to_resolution.get(key)
    if cached is not None:
        return cached
    if isinstance(key, pd2af.predicates.kept_species):
        input_species = clingo_id_to_model_element[key.species]
        if flatten_complexes:
            species = _flatten_complex(input_species)
            # Flattening drops subunits, so the input layout subtree
            # (carrying subunit/state/modification glyphs) no longer
            # matches the model. Treat as having no usable input layout.
            existing_species = (
                input_species if species is input_species else None
            )
        else:
            species = input_species
            existing_species = input_species
    elif isinstance(key, pd2af.predicates.new_species_from_template):
        input_species = clingo_id_to_model_element[key.species]
        candidate = _build_from_template(
            input_species, subunit_to_top_level, stripped_template_cache
        )
        species = _intern_species(candidate, intern_table)
        existing_species = None
    else:
        raise ValueError(f"unknown activity key wrapper {type(key).__name__}")
    resolution = (species, existing_species)
    key_to_resolution[key] = resolution
    return resolution


def _intern_species(species, intern_table):
    """Content-dedup species across all activity-key branches; smaller
    ``id_`` wins on collision. Mirrors
    ``momapy.io.utils.register_model_element`` so two content-equal
    activity sources resolve to a single identity end-to-end (model
    set, modulation source/target, layout-model mapping value)."""
    existing = intern_table.get(species)
    if existing is None:
        intern_table[species] = species
        return species
    species_id = species.id_ or ""
    existing_id = existing.id_ or ""
    if species_id < existing_id:
        intern_table[species] = species
        return species
    return existing


def _walk_species(resolved_species, *, flatten_complexes):
    """Linearise the species list to the order the layout strategy
    consumes: top-level first, then subunits of kept complexes.

    Content-based dedup with a ``canonical_species`` map. In two
    passes:

    1. **Subunits-of-kept-complexes first.** Every subunit reachable
       through a resolved complex's subtree is registered into
       ``canonical_species`` and emitted as an ``is_subunit=True``
       entry. This pins the canonical identity of each proteoform to
       the subunit-inside-complex.
    2. **Top-level entries.** A resolved species that is content-equal
       to an already-registered subunit (e.g. a subunit promoted from
       a suppressed complex) is skipped — modulation source/target
       are remapped through ``canonical_species`` so arcs land on the
       subunit's alias inside the kept complex (where the layout is
       structurally correct).
    """
    species_entries = []
    canonical_species = {}
    if not flatten_complexes:
        for species, _existing in resolved_species:
            for subunit in getattr(species, "subunits", ()) or ():
                if subunit in canonical_species:
                    continue
                canonical_species[subunit] = subunit
                species_entries.append((subunit, subunit, True))
    for species, existing_species in resolved_species:
        if species in canonical_species:
            continue
        canonical_species[species] = species
        species_entries.append((species, existing_species, False))
    return species_entries, canonical_species


def _flatten_complex(species):
    if not isinstance(species, momapy.celldesigner.Complex):
        return species
    if not getattr(species, "subunits", None):
        return species
    return dataclasses.replace(species, subunits=frozenset())


def _stripped_template_for(template, cache):
    cached = cache.get(id(template))
    if cached is not None:
        return cached
    fields_to_clear = {}
    if hasattr(template, "modification_residues"):
        fields_to_clear["modification_residues"] = frozenset()
    if hasattr(template, "regions"):
        fields_to_clear["regions"] = frozenset()
    stripped = dataclasses.replace(
        template,
        id_=f"{_STRIPPED_TEMPLATE_PREFIX}{template.id_}",
        **fields_to_clear,
    )
    cache[id(template)] = stripped
    return stripped


def _build_from_template(species, subunit_to_top_level, stripped_template_cache):
    template = species.template
    if template is None:
        raise ValueError(
            f"new_species_from_template key references species "
            f"{species.id_!r} which has no template"
        )
    species_class = _TEMPLATE_TO_SPECIES_CLASS.get(type(template))
    if species_class is None:
        raise ValueError(
            f"cannot synthesize a templated species for template class "
            f"{type(template).__name__}"
        )
    stripped_template = _stripped_template_for(template, stripped_template_cache)
    compartment = _ancestor_compartment(species, subunit_to_top_level)
    compartment_id = (
        compartment.id_ if compartment is not None else _NO_COMPARTMENT_SENTINEL
    )
    return species_class(
        id_=f"{_TEMPLATE_MERGE_PREFIX}{template.id_}__{compartment_id}",
        name=template.name,
        template=stripped_template,
        compartment=compartment,
    )


def _build_subunit_to_top_level(input_map):
    """Map every subunit (transitively) to its top-level species in the
    input map. Top-level species map to themselves so the lookup is
    total over all model species.
    """
    subunit_to_top_level = {}

    def walk(species, top_level):
        subunit_to_top_level[id(species)] = top_level
        for subunit in getattr(species, "subunits", ()) or ():
            walk(subunit, top_level)

    for species in input_map.model.species:
        walk(species, species)
    return subunit_to_top_level


def _ancestor_compartment(species, subunit_to_top_level):
    if getattr(species, "compartment", None) is not None:
        return species.compartment
    top_level = subunit_to_top_level.get(id(species), species)
    return getattr(top_level, "compartment", None)


def _compartments_outermost_first(compartments):
    def depth(compartment):
        result = 0
        seen = set()
        current = compartment.outside
        while current is not None and id(current) not in seen:
            seen.add(id(current))
            result += 1
            current = current.outside
        return result
    return sorted(compartments, key=depth)


def _collect_ancestor_compartments(compartments):
    expanded = set(compartments)
    frontier = expanded
    while True:
        next_frontier = set(
            compartment.outside
            for compartment in frontier
            if compartment.outside is not None
            and compartment.outside not in expanded
        )
        if not next_frontier:
            break
        expanded |= next_frontier
        frontier = next_frontier
    return expanded


def _collect_templates_from_species(species_iterable):
    collected = set()

    def visit(species):
        template = getattr(species, "template", None)
        if template is not None:
            collected.add(template)
        for subunit in getattr(species, "subunits", ()) or ():
            visit(subunit)

    for species in species_iterable:
        visit(species)
    return collected
