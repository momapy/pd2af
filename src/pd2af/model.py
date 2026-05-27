"""Construct AF model elements from clingo activity keys.

Stateless ``get_or_make_*`` helpers. Each helper interns the constructed
element through a shared content-keyed cache (``register_or_reuse``),
so that two content-equal elements collapse to a single Python identity
end-to-end. The coordinator in :mod:`pd2af.build` owns the cache and
orchestrates the layer-ordered canonicalisation.

The canonicity policy is "first-registered wins". Combined with the
coordinator's layer order (kept_species → kept_subunit →
promoted_subunit → new_species_from_template), this guarantees that
a kept input-map element is always the canonical instance for its
content class, never displaced by a freshly synthesized one.
"""

import dataclasses

import momapy.celldesigner


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


def register_or_reuse(element, cache):
    """Intern ``element`` by content in ``cache``. First-registered wins:
    if a content-equal element is already cached, return it; otherwise
    record ``element`` as the canonical instance and return it.
    """
    existing = cache.get(element)
    if existing is not None:
        return existing
    cache[element] = element
    return element


def get_or_make_kept_species(input_species, cache):
    """Canonical species for a ``kept_species``/``kept_subunit``/
    ``promoted_subunit`` key. The input species is the canonical
    instance — register it so later content-equal candidates collapse
    onto it.
    """
    return register_or_reuse(input_species, cache)


def get_or_make_stripped_template(input_template, cache):
    """Strip proteoform decorations from ``input_template`` and intern
    by content. Two distinct input templates that strip to the same
    content yield a single canonical stripped template.

    Kept templates must already be registered in the cache (per pipeline
    step (3), kept templates are collected before stripped ones are
    built). When a stripped candidate is content-equal to a kept
    template, ``register_or_reuse`` returns the kept canonical and no
    duplicate enters the cache.
    """
    fields_to_clear = {}
    if hasattr(input_template, "modification_residues"):
        fields_to_clear["modification_residues"] = frozenset()
    if hasattr(input_template, "regions"):
        fields_to_clear["regions"] = frozenset()
    candidate = dataclasses.replace(
        input_template,
        id_=f"{_STRIPPED_TEMPLATE_PREFIX}{input_template.id_}",
        **fields_to_clear,
    )
    return register_or_reuse(candidate, cache)


def get_or_make_synthesized_species(
    input_species, stripped_template, compartment, cache
):
    """Build a synthesized species from a stripped template and a
    compartment. Two ``new_species_from_template`` keys whose (template,
    compartment) cells coincide yield a single canonical species.
    """
    template = input_species.template
    if template is None:
        raise ValueError(
            f"new_species_from_template key references species "
            f"{input_species.id_!r} which has no template"
        )
    species_class = _TEMPLATE_TO_SPECIES_CLASS.get(type(template))
    if species_class is None:
        raise ValueError(
            f"cannot synthesize a templated species for template class "
            f"{type(template).__name__}"
        )
    compartment_id = (
        compartment.id_ if compartment is not None else _NO_COMPARTMENT_SENTINEL
    )
    candidate = species_class(
        id_=f"{_TEMPLATE_MERGE_PREFIX}{template.id_}__{compartment_id}",
        name=template.name,
        template=stripped_template,
        compartment=compartment,
    )
    return register_or_reuse(candidate, cache)


def get_or_make_modulation(modulation_class, source, target, cache):
    """Build a modulation and intern by content."""
    candidate = modulation_class(source=source, target=target)
    return register_or_reuse(candidate, cache)


def get_parent_complex_compartment(subunit, subunit_to_top_level):
    """Resolve a subunit's effective compartment by walking to its
    containing top-level species. Subunits typically have no
    ``compartment`` attribute of their own.
    """
    if getattr(subunit, "compartment", None) is not None:
        return subunit.compartment
    top_level = subunit_to_top_level.get(id(subunit), subunit)
    return getattr(top_level, "compartment", None)


def build_subunit_to_top_level(input_map):
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


def compartments_outermost_first(compartments):
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


def collect_ancestor_compartments(compartments):
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


def collect_templates_from_species(species_iterable):
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
