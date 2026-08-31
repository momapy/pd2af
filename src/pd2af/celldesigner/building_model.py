"""Build the AF model from clingo activity / influence atoms.

``make_and_add_model`` is the model pass: it walks the activity atoms in
layer order (keptSpeciesKey → promotedSubunitKey) and populates ``context.model``
with canonical, content-deduped compartments, templates, species and
modulations. In the merged modes (``normal``/``normal-no-complex``) each activity's
species is stripped of its PTM decorations (recursively, including subunits)
by ``get_or_make_stripped_species``; the other modes reuse input species by
identity.

The stateless ``get_or_make_*`` leaf helpers do the actual element
construction. Each interns the constructed element through a shared
content-keyed cache (``register_or_reuse``), so two content-equal elements
collapse to a single Python identity end-to-end. The canonicity policy is
"first-registered wins"; combined with the layer order above, a kept
input-map element is always the canonical instance for its content class,
never displaced by a freshly stripped one.

:mod:`pd2af.build` owns the ``BuilderContext`` and drives this pass, then
the layout pass.
"""

import dataclasses

import momapy.builder
import momapy.celldesigner

import pd2af.building_model
import pd2af.predicates
from pd2af.utils import add_model_element_if_new, register_or_reuse


_STRIPPED_TEMPLATE_PREFIX = "merged_template__"


_SPECIES_LAYER_ORDER = (
    pd2af.predicates.keptSpeciesKey,
    pd2af.predicates.promotedSubunitKey,
)

# Operator-type token (``logicalOperator.type_``) -> CellDesigner gate class.
# The NOT token is ``not_`` because bare ``not`` is a reserved clingo keyword.
_OPERATOR_TYPE_TO_GATE_CLASS = {
    "and": momapy.celldesigner.AndGate,
    "or": momapy.celldesigner.OrGate,
    "not_": momapy.celldesigner.NotGate,
    "unknown": momapy.celldesigner.UnknownGate,
}


def get_or_make_kept_species_key_or_subunit(input_species, cache):
    """Canonical species for a ``keptSpeciesKey`` key (non-merged modes).
    The input species is the canonical instance — register it so later
    content-equal candidates collapse onto it.
    """
    return register_or_reuse(input_species, cache)


def get_or_make_promoted_subunit_key_species(
    input_subunit,
    subunit_to_top_level,
    cache,
    input_model_element_to_canonical_model_element,
):
    """Canonical species for a ``promotedSubunitKey`` key. CellDesigner
    subunits carry ``compartment=None`` (inherited from the parent
    complex); when promoted to top-level they need the parent's
    compartment, otherwise the writer substitutes a synthetic ``default``
    compartment that the reader picks up, breaking content-eq on
    round-trip.

    Two distinct input subunits promoted to the same compartment
    collapse via the cache; if a content-equal kept species is already
    cached (from the keptSpeciesKey layer), the corrected version
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


def get_or_make_stripped_species(
    input_species,
    compartment,
    cache,
    input_model_element_to_canonical_model_element,
):
    """Return a decoration-free canonical species for ``input_species`` (the
    merged modes ``normal``/``normal-no-complex``). Clears every post-translational
    decoration -- ``active`` (-> False), ``homomultimer`` (-> 1),
    ``modifications``, ``structural_states``, and, via a stripped template,
    ``modification_residues``/``regions`` -- sets the effective ``compartment``,
    and strips ``subunits`` recursively (a ``frozenset`` collapses subunits that
    strip to equal content). Built with ``dataclasses.replace`` so the concrete
    class and the input's real, reader-resolvable ``id_`` are preserved;
    content-equal stripped twins still intern to one canonical because ``id_``
    is ``compare=False``. Records ``id(input) -> canonical`` when they differ so
    Pass 2 can substitute stale references (mirrors
    ``get_or_make_promoted_subunit_key_species``).
    """
    fields = {}
    if hasattr(input_species, "active"):
        fields["active"] = False
    if hasattr(input_species, "homomultimer"):
        fields["homomultimer"] = 1
    if hasattr(input_species, "modifications"):
        fields["modifications"] = frozenset()
    if hasattr(input_species, "structural_states"):
        fields["structural_states"] = frozenset()
    if hasattr(input_species, "compartment"):
        fields["compartment"] = compartment
    template = getattr(input_species, "template", None)
    if template is not None:
        fields["template"] = get_or_make_stripped_template(template, cache)
    subunits = getattr(input_species, "subunits", None)
    if subunits:
        fields["subunits"] = frozenset(
            get_or_make_stripped_species(
                subunit,
                getattr(subunit, "compartment", None),
                cache,
                input_model_element_to_canonical_model_element,
            )
            for subunit in subunits
        )
    candidate = dataclasses.replace(input_species, **fields)
    canonical = register_or_reuse(candidate, cache)
    if canonical is not input_species:
        input_model_element_to_canonical_model_element[id(input_species)] = (
            canonical
        )
    return canonical


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


# ---------------------------------------------------------------------------
# Pass 1 -- model construction
# ---------------------------------------------------------------------------


def make_and_add_model(context, clingo_model):
    context.model = momapy.builder.get_or_make_builder_cls(
        momapy.celldesigner.CellDesignerModel
    )()
    context.subunit_to_top_level = pd2af.building_model.build_subunit_to_top_level(
        context.input_map.model.species
    )
    pd2af.building_model.collect_atoms(context, clingo_model)
    _make_and_add_compartments(context)
    _make_and_add_templates(context)
    _make_and_add_species(context)
    pd2af.building_model.make_and_add_operators(
        context,
        _OPERATOR_TYPE_TO_GATE_CLASS,
        momapy.celldesigner.BooleanLogicGateInput,
        context.model.boolean_logic_gates,
    )
    _make_and_add_modulations(context)


def _activity_atoms_in_layer_order(context):
    """The activity atoms sorted by the layer their key class belongs to, so a
    kept species is always registered before a promoted subunit."""
    return sorted(
        context.activity_atoms,
        key=lambda atom: _SPECIES_LAYER_ORDER.index(type(atom.key)),
    )


def _make_and_add_compartments(context):
    immediate_compartments = set()
    for atom in _activity_atoms_in_layer_order(context):
        input_species = context.clingo_id_to_model_element[atom.key.species]
        compartment = _compartment_for_input_species(context, input_species)
        if compartment is not None:
            immediate_compartments.add(compartment)
    for compartment in compartments_outermost_first(
        collect_ancestor_compartments(immediate_compartments)
    ):
        context.model.compartments.add(compartment)
        # Compartments carry over by identity, so the output compartment is the
        # input one; record it for provenance so its annotations/notes carry.
        context.compartment_emissions.append((compartment, compartment))


def _compartment_for_input_species(context, input_species):
    if getattr(input_species, "compartment", None) is not None:
        return input_species.compartment
    return get_parent_complex_compartment(
        input_species, context.subunit_to_top_level
    )


def _make_and_add_templates(context):
    strip = context.mode.merges_proteoforms
    seen_template_identities = set()

    # Register the templates each activity carries (walking subunit trees). In
    # the merged modes the registered template is the *stripped* one, so the
    # stripped species built in the species pass finds its canonical template
    # already present (same `context.cache`); otherwise the input template is
    # registered verbatim. This runs before `_make_and_add_species`.
    for atom in _activity_atoms_in_layer_order(context):
        input_species = context.clingo_id_to_model_element[atom.key.species]
        for input_template in _walk_templates(input_species):
            if strip:
                canonical = get_or_make_stripped_template(
                    input_template, context.cache
                )
            else:
                canonical = register_or_reuse(input_template, context.cache)
            add_model_element_if_new(
                context.model.species_templates,
                canonical,
                seen_template_identities,
            )


def _walk_templates(species):
    template = getattr(species, "template", None)
    if template is not None:
        yield template
    for subunit in getattr(species, "subunits", ()) or ():
        yield from _walk_templates(subunit)


def _make_and_add_species(context):
    seen_species_identities = set()
    for atom in _activity_atoms_in_layer_order(context):
        species = _resolve_activity_key(context, atom.key)
        context.key_to_activity[atom.key] = species
        if add_model_element_if_new(
            context.model.species, species, seen_species_identities
        ):
            input_species = context.clingo_id_to_model_element[
                atom.key.species
            ]
            context.activity_emissions.append(
                (type(atom.key), species, input_species)
            )


def _resolve_activity_key(context, key):
    """Resolve an activity key to its output species. In the merged modes
    (``normal``/``normal-no-complex``) every species is stripped of its PTM decorations
    (recursively, including subunits) so content-equal proteoforms collapse;
    otherwise the input species is reused by identity (``keptSpeciesKey``) or
    promoted with a corrected compartment (``promotedSubunitKey``)."""
    input_species = context.clingo_id_to_model_element[key.species]
    if context.mode.merges_proteoforms:
        compartment = _compartment_for_input_species(context, input_species)
        return get_or_make_stripped_species(
            input_species,
            compartment,
            context.cache,
            context.input_model_element_to_canonical_model_element,
        )
    if isinstance(key, pd2af.predicates.promotedSubunitKey):
        return get_or_make_promoted_subunit_key_species(
            input_species,
            context.subunit_to_top_level,
            context.cache,
            context.input_model_element_to_canonical_model_element,
        )
    if isinstance(key, pd2af.predicates.keptSpeciesKey):
        return get_or_make_kept_species_key_or_subunit(input_species, context.cache)
    raise ValueError(f"unknown activity key wrapper {type(key).__name__}")


def _make_and_add_modulations(context):
    seen_modulation_identities = set()
    for atom in context.influence_atoms:
        modulation_class = pd2af.predicates.predicate_to_model_element_class[
            type(atom)
        ]
        source = pd2af.building_model.resolve_influence_source(
            context, atom.source
        )
        if source is None:
            continue
        target = context.key_to_activity[atom.target]
        canonical = get_or_make_modulation(
            modulation_class, source, target, context.cache
        )
        add_model_element_if_new(
            context.model.modulations, canonical, seen_modulation_identities
        )
