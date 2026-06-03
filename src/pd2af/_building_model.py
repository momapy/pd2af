"""Build the AF model from clingo activity / influence atoms.

``make_and_add_model`` is the model pass: it walks the activity atoms in
layer order (kept_species → kept_subunit → promoted_subunit →
new_species_from_template) and populates ``context.model`` with canonical,
content-deduped compartments, templates, species and modulations.

The stateless ``get_or_make_*`` leaf helpers do the actual element
construction. Each interns the constructed element through a shared
content-keyed cache (``register_or_reuse``), so two content-equal elements
collapse to a single Python identity end-to-end. The canonicity policy is
"first-registered wins"; combined with the layer order above, a kept
input-map element is always the canonical instance for its content class,
never displaced by a freshly synthesized one.

:mod:`pd2af.build` owns the ``BuilderContext`` and drives this pass, then
the layout pass.
"""

import dataclasses

import momapy.builder
import momapy.celldesigner

import pd2af.predicates


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


def get_or_make_kept_species_or_subunit(input_species, cache):
    """Canonical species for a ``kept_species`` or ``kept_subunit`` key.
    The input species is the canonical instance — register it so later
    content-equal candidates collapse onto it.
    """
    return register_or_reuse(input_species, cache)


def get_or_make_promoted_subunit_species(
    input_subunit,
    subunit_to_top_level,
    cache,
    input_model_element_to_canonical_model_element,
):
    """Canonical species for a ``promoted_subunit`` key. CellDesigner
    subunits carry ``compartment=None`` (inherited from the parent
    complex); when promoted to top-level they need the parent's
    compartment, otherwise the writer substitutes a synthetic ``default``
    compartment that the reader picks up, breaking content-eq on
    round-trip.

    Two distinct input subunits promoted to the same compartment
    collapse via the cache; if a content-equal kept species is already
    cached (from the kept_species layer), the corrected version
    collapses onto that.

    Records ``id(input_subunit) -> canonical_species`` in
    ``input_model_element_to_canonical_model_element`` whenever the
    canonical differs from the input — Pass 2's layout mapping copies
    use that dict to substitute stale value references that still point
    at the input subunit.
    """
    parent_compartment = get_parent_complex_compartment(
        input_subunit, subunit_to_top_level
    )
    if parent_compartment is None:
        return register_or_reuse(input_subunit, cache)
    corrected_species = dataclasses.replace(
        input_subunit, compartment=parent_compartment
    )
    canonical_species = register_or_reuse(corrected_species, cache)
    if canonical_species is not input_subunit:
        input_model_element_to_canonical_model_element[id(input_subunit)] = (
            canonical_species
        )
    return canonical_species


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


# ---------------------------------------------------------------------------
# Pass 1 -- model construction
# ---------------------------------------------------------------------------


def make_and_add_model(context, clingo_model):
    context.model = momapy.builder.get_or_make_builder_cls(
        momapy.celldesigner.CellDesignerModel
    )()
    context.subunit_to_top_level = build_subunit_to_top_level(context.input_map)
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
    for compartment in compartments_outermost_first(
        collect_ancestor_compartments(immediate_compartments)
    ):
        context.model.compartments.add(compartment)


def _compartment_for_input_species(context, input_species):
    if getattr(input_species, "compartment", None) is not None:
        return input_species.compartment
    return get_parent_complex_compartment(
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
                canonical = register_or_reuse(input_template, context.cache)
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
        canonical = get_or_make_stripped_template(input_template, context.cache)
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
            input_species = (
                context.clingo_id_to_model_element[atom.key.species]
                if key_class in _KEPT_KEY_CLASSES
                else None
            )
            context.species_emissions.append((key_class, species, input_species))
            if key_class is not pd2af.predicates.kept_subunit:
                context.model.species.add(species)


def _resolve_activity_key(context, key):
    if isinstance(key, pd2af.predicates.promoted_subunit):
        input_species = context.clingo_id_to_model_element[key.species]
        return get_or_make_promoted_subunit_species(
            input_species,
            context.subunit_to_top_level,
            context.cache,
            context.input_model_element_to_canonical_model_element,
        )
    if isinstance(
        key, (pd2af.predicates.kept_species, pd2af.predicates.kept_subunit)
    ):
        input_species = context.clingo_id_to_model_element[key.species]
        return get_or_make_kept_species_or_subunit(input_species, context.cache)
    if isinstance(key, pd2af.predicates.new_species_from_template):
        input_species = context.clingo_id_to_model_element[key.species]
        compartment = _compartment_for_input_species(context, input_species)
        stripped_template = get_or_make_stripped_template(
            input_species.template, context.cache
        )
        return get_or_make_synthesized_species(
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
        canonical = get_or_make_modulation(
            modulation_class, source, target, context.cache
        )
        if id(canonical) in seen_modulation_identities:
            continue
        seen_modulation_identities.add(id(canonical))
        context.model.modulations.add(canonical)
