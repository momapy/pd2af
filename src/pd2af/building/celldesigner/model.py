"""Build the AF model from clingo activity / influence atoms.

``make_and_add_model`` is the model pass: it walks the activity atoms in
layer order (keptSpeciesKey → promotedSubunitKey) and populates ``context.model``
with canonical, content-deduped compartments, templates, species and
modulations. When the forms of the same base entity merge (``keep_species``
off) each activity's species is stripped of its PTM decorations (recursively,
including subunits) by ``get_or_make_stripped_species``; otherwise the input
species is kept, placed in its effective compartment by
``get_or_make_species_in_compartment``.

Which compartment a species goes in is ``_compartment_for_input_species``: its
own, its complex's when it is a promoted subunit, or -- with the
``drop_compartments`` option -- the input map's default compartment, whatever
the species' own. Every compartment then merges into that one, and species and
influences that become content-equal collapse through the shared cache.

The stateless ``get_or_make_*`` leaf helpers do the actual element
construction. Each interns the constructed element through a shared
content-keyed cache (``register_or_reuse``), so two content-equal elements
collapse to a single Python identity end-to-end. The canonicity policy is
"first-registered wins"; combined with the layer order above, a kept
input-map element is always the canonical instance for its content class,
never displaced by a freshly stripped one.

:mod:`pd2af.core` drives this pass through a
:class:`pd2af.building.context.BuilderContext`, then
the layout pass.
"""

import dataclasses
import typing

import momapy.builder
import momapy.celldesigner

import pd2af.asp.predicates
import pd2af.building.context
import pd2af.building.model
from pd2af.building.model import add_model_element_if_new, register_or_reuse

_STRIPPED_TEMPLATE_PREFIX = "merged_template__"

# The id CellDesigner gives the compartment a species with no compartment of its
# own belongs to, and the one every species goes in under `drop_compartments`.
_DEFAULT_COMPARTMENT_ID = "default"


_SPECIES_LAYER_ORDER = (
    pd2af.asp.predicates.keptSpeciesKey,
    pd2af.asp.predicates.promotedSubunitKey,
)

# Operator-type token (``logicalOperator.type_``) -> CellDesigner gate class.
# The NOT token is ``not_`` because bare ``not`` is a reserved clingo keyword.
_OPERATOR_TYPE_TO_GATE_CLASS = {
    "and": momapy.celldesigner.AndGate,
    "or": momapy.celldesigner.OrGate,
    "not_": momapy.celldesigner.NotGate,
    "unknown": momapy.celldesigner.UnknownGate,
}


def get_or_make_species_in_compartment(
    input_species: typing.Any,
    compartment: typing.Any,
    cache: dict,
    input_model_element_to_canonical_model_element: dict,
) -> typing.Any:
    """Canonical species for ``input_species`` placed in ``compartment``.

    Used with ``keep_species`` on, which keeps each PD species as its own
    activity, for both activity-key kinds. A species already in
    ``compartment`` — the usual ``keptSpeciesKey`` case — is its own canonical
    instance, so it is registered as is and later
    content-equal candidates collapse onto it. Otherwise the species is rebuilt
    in ``compartment``: a promoted subunit carries ``compartment=None``
    (inherited from the parent complex) and needs its complex's compartment
    once it stands at top level, and the ``drop_compartments`` option puts
    every species in the default compartment so that species differing only by
    compartment collapse.

    Two input species placed in the same compartment collapse via the cache; a
    rebuilt species that is content-equal to an already-cached kept one
    collapses onto that.

    Records ``id(input_species) -> canonical`` in
    ``input_model_element_to_canonical_model_element`` whenever the canonical
    differs from the input — Pass 2's layout mapping copies use that dict to
    substitute stale value references that still point at the input species.
    """
    if getattr(input_species, "compartment", None) is compartment:
        return register_or_reuse(input_species, cache)
    candidate = dataclasses.replace(input_species, compartment=compartment)
    canonical = register_or_reuse(candidate, cache)
    if canonical is not input_species:
        input_model_element_to_canonical_model_element[id(input_species)] = canonical
    return canonical


def get_or_make_stripped_template(
    input_template: typing.Any, cache: dict
) -> typing.Any:
    """Strip the PTM decorations from ``input_template`` and intern by content.

    Two distinct input templates that strip to the same content yield a single
    canonical stripped template.

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
    input_species: typing.Any,
    compartment: typing.Any,
    cache: dict,
    input_model_element_to_canonical_model_element: dict,
) -> typing.Any:
    """Return a decoration-free canonical species for ``input_species``.

    Used with ``keep_species`` off. The stripped species is built by
    :func:`make_stripped_species` and interned by content, so content-equal
    stripped twins collapse to one canonical (``id_`` is ``compare=False``).
    Only the species itself is interned: its subunits are built fresh for it,
    so two complexes never share a subunit object.

    Records ``id(input) -> canonical`` when they differ so Pass 2 can
    substitute stale references (mirrors
    ``get_or_make_species_in_compartment``), and records each input subunit, at
    any depth, against the subunit of the canonical that stands for it
    (:func:`record_canonical_subunits`).
    """
    input_subunit_to_stripped_subunit = {}
    candidate = make_stripped_species(
        input_species, compartment, cache, input_subunit_to_stripped_subunit
    )
    canonical = register_or_reuse(candidate, cache)
    if canonical is not input_species:
        input_model_element_to_canonical_model_element[id(input_species)] = canonical
    record_canonical_subunits(
        input_species,
        canonical,
        input_subunit_to_stripped_subunit,
        input_model_element_to_canonical_model_element,
    )
    return canonical


def make_stripped_species(
    input_species: typing.Any,
    compartment: typing.Any,
    cache: dict,
    input_subunit_to_stripped_subunit: dict,
) -> typing.Any:
    """Build a decoration-free copy of ``input_species``, without interning it.

    Clears every post-translational decoration -- ``active`` (-> False),
    ``homomultimer`` (-> 1), ``modifications``, ``structural_state``, and, via
    a stripped template, ``modification_residues``/``regions`` -- sets the
    effective ``compartment``, and strips ``subunits`` recursively (a
    ``frozenset`` collapses subunits that strip to equal content). Built with
    ``dataclasses.replace`` so the concrete class and the input's real,
    reader-resolvable ``id_`` are preserved. The template is interned through
    ``cache``; the subunits are not, and ``id(input_subunit) -> stripped`` is
    recorded in ``input_subunit_to_stripped_subunit`` for every subunit at any
    depth.
    """
    fields = {}
    if hasattr(input_species, "active"):
        fields["active"] = False
    if hasattr(input_species, "homomultimer"):
        fields["homomultimer"] = 1
    if hasattr(input_species, "modifications"):
        fields["modifications"] = frozenset()
    if hasattr(input_species, "structural_state"):
        fields["structural_state"] = None
    if hasattr(input_species, "compartment"):
        fields["compartment"] = compartment
    template = getattr(input_species, "template", None)
    if template is not None:
        fields["template"] = get_or_make_stripped_template(template, cache)
    subunits = getattr(input_species, "subunits", None)
    if subunits:
        stripped_subunits = []
        for subunit in subunits:
            stripped_subunit = make_stripped_species(
                subunit,
                getattr(subunit, "compartment", None),
                cache,
                input_subunit_to_stripped_subunit,
            )
            input_subunit_to_stripped_subunit[id(subunit)] = stripped_subunit
            stripped_subunits.append(stripped_subunit)
        fields["subunits"] = frozenset(stripped_subunits)
    return dataclasses.replace(input_species, **fields)


def record_canonical_subunits(
    input_species: typing.Any,
    canonical_species: typing.Any,
    input_subunit_to_stripped_subunit: dict,
    input_model_element_to_canonical_model_element: dict,
) -> None:
    """Record each input subunit against the subunit of ``canonical_species`` it became.

    The stripped copy :func:`make_stripped_species` built for an input subunit
    is not always the one in the model: two subunits of one complex that strip
    to equal content collapse in the ``frozenset``, and a complex that collapses
    onto an earlier content-equal one keeps none of its own copies. Each copy is
    content-equal to a subunit of the canonical, so that subunit is looked up by
    content and recorded, recursing into nested complexes.
    """
    canonical_subunit_by_content = {
        canonical_subunit: canonical_subunit
        for canonical_subunit in getattr(canonical_species, "subunits", None) or ()
    }
    for input_subunit in getattr(input_species, "subunits", None) or ():
        stripped_subunit = input_subunit_to_stripped_subunit[id(input_subunit)]
        canonical_subunit = canonical_subunit_by_content[stripped_subunit]
        input_model_element_to_canonical_model_element[id(input_subunit)] = (
            canonical_subunit
        )
        record_canonical_subunits(
            input_subunit,
            canonical_subunit,
            input_subunit_to_stripped_subunit,
            input_model_element_to_canonical_model_element,
        )


def get_or_make_modulation(
    modulation_class: type, source: typing.Any, target: typing.Any, cache: dict
) -> typing.Any:
    """Build a modulation and intern by content."""
    candidate = modulation_class(source=source, target=target)
    return register_or_reuse(candidate, cache)


def get_parent_complex_compartment(
    subunit: typing.Any, subunit_to_top_level: dict[int, typing.Any]
) -> typing.Any:
    """Resolve a subunit's effective compartment.

    Walks to its containing top-level species: subunits typically have no
    ``compartment`` attribute of their own.
    """
    if getattr(subunit, "compartment", None) is not None:
        return subunit.compartment
    top_level = subunit_to_top_level.get(id(subunit), subunit)
    return getattr(top_level, "compartment", None)


def find_default_compartment(input_model: typing.Any) -> typing.Any:
    """The default compartment of ``input_model``, or ``None`` if it has none.

    CellDesigner declares a compartment with id ``default`` in every map, and
    it is the one the writer points a species with no compartment at. It is the
    compartment every activity goes in under the ``drop_compartments`` option,
    so the output still declares the compartment its species refer to.
    """
    for compartment in input_model.compartments:
        if compartment.id_ == _DEFAULT_COMPARTMENT_ID:
            return compartment
    return None


def compartments_outermost_first(compartments: typing.Any) -> typing.Any:
    """``compartments`` sorted by nesting depth, so a container precedes its content."""

    def depth(compartment: typing.Any) -> typing.Any:
        result = 0
        seen = set()
        current = compartment.outside
        while current is not None and id(current) not in seen:
            seen.add(id(current))
            result += 1
            current = current.outside
        return result

    return sorted(compartments, key=depth)


def collect_ancestor_compartments(compartments: typing.Any) -> typing.Any:
    """``compartments`` plus every compartment they are nested in, transitively."""
    expanded = set(compartments)
    frontier = expanded
    while True:
        next_frontier = {
            compartment.outside
            for compartment in frontier
            if compartment.outside is not None and compartment.outside not in expanded
        }
        if not next_frontier:
            break
        expanded |= next_frontier
        frontier = next_frontier
    return expanded


# ---------------------------------------------------------------------------
# Pass 1 -- model construction
# ---------------------------------------------------------------------------


def make_and_add_model(
    context: pd2af.building.context.BuilderContext, clingo_model: typing.Any
):
    """Build ``context.model`` from the clingo atoms (pass 1)."""
    context.model = momapy.builder.get_or_make_builder_cls(
        momapy.celldesigner.CellDesignerModel
    )()
    context.subunit_to_top_level = pd2af.building.model.build_subunit_to_top_level(
        context.input_map.model.species
    )
    if context.drop_compartments:
        context.default_compartment = find_default_compartment(context.input_map.model)
    pd2af.building.model.collect_atoms(context, clingo_model)
    _make_and_add_compartments(context)
    _make_and_add_templates(context)
    _make_and_add_species(context)
    pd2af.building.model.make_and_add_operators(
        context,
        _OPERATOR_TYPE_TO_GATE_CLASS,
        momapy.celldesigner.BooleanLogicGateInput,
        context.model.boolean_logic_gates,
    )
    _make_and_add_modulations(context)


def _activity_atoms_in_layer_order(
    context: pd2af.building.context.BuilderContext,
) -> typing.Any:
    """The activity atoms sorted by the layer their key class belongs to.

    A kept species is thus always registered before a promoted subunit.
    """
    return sorted(
        context.activity_atoms,
        key=lambda atom: _SPECIES_LAYER_ORDER.index(type(atom.key)),
    )


def _make_and_add_compartments(context: pd2af.building.context.BuilderContext):
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


def _compartment_for_input_species(
    context: pd2af.building.context.BuilderContext, input_species: typing.Any
) -> typing.Any:
    """The compartment ``input_species`` becomes an activity in.

    Its own, or its top-level complex's when it is a subunit. With the
    ``drop_compartments`` option it is the default compartment whatever the
    species' own, so the compartment pass collects that one alone and species
    differing only by compartment merge.
    """
    if context.drop_compartments:
        return context.default_compartment
    if getattr(input_species, "compartment", None) is not None:
        return input_species.compartment
    return get_parent_complex_compartment(input_species, context.subunit_to_top_level)


def _make_and_add_templates(context: pd2af.building.context.BuilderContext):
    strip = not context.keep_species
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
                canonical = get_or_make_stripped_template(input_template, context.cache)
            else:
                canonical = register_or_reuse(input_template, context.cache)
            add_model_element_if_new(
                context.model.species_templates,
                canonical,
                seen_template_identities,
            )


def _walk_templates(species: typing.Any):
    template = getattr(species, "template", None)
    if template is not None:
        yield template
    for subunit in getattr(species, "subunits", ()) or ():
        yield from _walk_templates(subunit)


def _make_and_add_species(context: pd2af.building.context.BuilderContext):
    seen_species_identities = set()
    for atom in _activity_atoms_in_layer_order(context):
        species = _resolve_activity_key(context, atom.key)
        context.key_to_activity[atom.key] = species
        if add_model_element_if_new(
            context.model.species, species, seen_species_identities
        ):
            input_species = context.clingo_id_to_model_element[atom.key.species]
            context.activity_emissions.append((type(atom.key), species, input_species))


def _resolve_activity_key(
    context: pd2af.building.context.BuilderContext, key: typing.Any
) -> typing.Any:
    """Resolve an activity key to its output species.

    With ``keep_species`` off every species is stripped of its PTM decorations
    (recursively, including subunits) so the content-equal forms of the same
    base entity collapse; otherwise the input species is kept. Either way it goes in the
    compartment ``_compartment_for_input_species`` gives it, which is what
    lifts a promoted subunit into its complex's compartment and what puts every
    species in the default one under the ``drop_compartments`` option. A key class
    outside ``_SPECIES_LAYER_ORDER`` never reaches here:
    ``_activity_atoms_in_layer_order`` raises on it first.
    """
    input_species = context.clingo_id_to_model_element[key.species]
    compartment = _compartment_for_input_species(context, input_species)
    if not context.keep_species:
        return get_or_make_stripped_species(
            input_species,
            compartment,
            context.cache,
            context.input_model_element_to_canonical_model_element,
        )
    return get_or_make_species_in_compartment(
        input_species,
        compartment,
        context.cache,
        context.input_model_element_to_canonical_model_element,
    )


def _make_and_add_modulations(context: pd2af.building.context.BuilderContext):
    seen_modulation_identities = set()
    for atom in context.influence_atoms:
        modulation_class = pd2af.asp.predicates.predicate_to_model_element_class[
            type(atom)
        ]
        source = pd2af.building.model.resolve_influence_source(context, atom.source)
        if source is None:
            continue
        target = context.key_to_activity[atom.target]
        canonical = get_or_make_modulation(
            modulation_class, source, target, context.cache
        )
        add_model_element_if_new(
            context.model.modulations, canonical, seen_modulation_identities
        )
