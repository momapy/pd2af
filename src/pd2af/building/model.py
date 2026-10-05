"""The parts of the model pass that are the same for either language.

:mod:`pd2af.building.celldesigner.model` and :mod:`pd2af.building.sbgn.model`
build their own compartments, species and influences, but read the clingo atoms,
resolve an influence source and walk the input subunit trees the same way. The logical-operator pass is
shared too: the CellDesigner builder makes ``BooleanLogicGate`` objects and the
SBGN-AF builder ``LogicalOperator`` objects, and otherwise the two do the
same thing -- read the operator atoms, resolve each operator's inputs through
``key_to_activity``, and add to the model those operators that actually source
an influence. Each builder passes in the classes it wants and the model
collection to add to.
"""

import collections.abc
import typing

import pd2af.asp.predicates
import pd2af.building.context


def register_or_reuse(element: typing.Any, cache: dict) -> typing.Any:
    """Intern ``element`` by content in ``cache``.

    First-registered wins: if a content-equal element is already cached, return
    it; otherwise record ``element`` as the canonical instance and return it.
    """
    existing = cache.get(element)
    if existing is not None:
        return existing
    cache[element] = element
    return element


def add_model_element_if_new(
    collection: typing.Any, model_element: typing.Any, seen_identities: set[int]
) -> bool:
    """Append ``model_element`` to ``collection`` unless it was already appended.

    Identities already added are tracked in ``seen_identities``.

    ``model_element`` is assumed to already be the canonical instance (e.g. the
    result of :func:`register_or_reuse`); this only guards against adding the
    same identity twice. Returns ``True`` if it was added this call, ``False``
    if it was a duplicate.
    """
    if id(model_element) in seen_identities:
        return False
    seen_identities.add(id(model_element))
    collection.add(model_element)
    return True


def collect_atoms(
    context: pd2af.building.context.BuilderContext, clingo_model: typing.Any
) -> None:
    """Sort the ``new(...)`` atoms into the context scratch lists by payload type.

    Activities, influences, logical operators and operator inputs each get their
    own list.
    """
    for atom in clingo_model.query(pd2af.asp.predicates.new).all():
        payload = atom.object_
        if isinstance(payload, pd2af.asp.predicates.activity):
            context.activity_atoms.append(payload)
        elif isinstance(payload, pd2af.asp.predicates.INFLUENCE_PREDICATES):
            context.influence_atoms.append(payload)
        elif isinstance(payload, pd2af.asp.predicates.logicalOperator):
            context.operator_atoms.append(payload)
        elif isinstance(payload, pd2af.asp.predicates.logicalOperatorInput):
            context.operator_input_atoms.append(payload)


def resolve_influence_source(
    context: pd2af.building.context.BuilderContext, source_key: typing.Any
) -> typing.Any:
    """Resolve an influence ``source`` key to its model element.

    A logical operator resolves through ``key_to_operator``, any activity key
    through ``key_to_activity``. ``None`` when the gate or the activity was not
    built, so the caller skips the edge.
    """
    if isinstance(source_key, pd2af.asp.predicates.logicalOperatorKey):
        return context.key_to_operator.get(source_key)
    return context.key_to_activity.get(source_key)


def build_subunit_to_top_level(
    top_level_elements: collections.abc.Iterable[typing.Any],
) -> dict[int, typing.Any]:
    """Map ``id(subunit)`` -> the top-level element it belongs to.

    Every subunit of ``top_level_elements`` is covered, at any depth.
    A subunit is a structural component, never an independent activity, so both
    builders attribute it to its outermost element: the CellDesigner builder
    takes the top-level species' compartment, the SBGN-AF builder the top-level
    entity pool's.
    """
    subunit_to_top_level = {}

    def walk(element: typing.Any, top_level_element: typing.Any) -> None:
        for subunit in getattr(element, "subunits", None) or ():
            subunit_to_top_level[id(subunit)] = top_level_element
            walk(subunit, top_level_element)

    for element in top_level_elements:
        walk(element, element)
    return subunit_to_top_level


def make_and_add_operators(
    context: pd2af.building.context.BuilderContext,
    operator_type_to_class: dict[str, type],
    operator_input_class: type,
    model_operators: typing.Any,
) -> None:
    """Build every authored logical operator, adding the used ones to the model.

    An operator reaches ``model_operators`` when it sources an influence or is
    reachable from such an operator through input edges (SBGN-PD nesting: the
    inner operators of an influence-sourcing AND are part of the output). It is
    added only when it has at least one resolved input: an operator whose
    inputs all failed to resolve (every referred element carries no activity
    key, e.g. every input suppressed) is dropped entirely -- it emits neither a
    node nor a sourced influence, and its consumer keeps only its remaining
    resolved inputs. An operator with some resolved inputs keeps exactly those
    inputs; the unresolved ones are skipped.

    Every surviving operator is built into ``context.key_to_operator``, so the
    influence pass can resolve an operator source. The added operators are
    appended to ``context.operator_emissions`` in dependency order, inner
    operators first, so the layout pass finds each input's glyph already
    registered. The CellDesigner writer emits a gate only through its
    modulation, so dropping unused gates keeps the output round-trip-safe.
    """
    inputs_by_operator = {}
    for input_atom in context.operator_input_atoms:
        inputs_by_operator.setdefault(input_atom.operator, []).append(input_atom.input)
    used_operator_keys = {
        atom.source
        for atom in context.influence_atoms
        if isinstance(atom.source, pd2af.asp.predicates.logicalOperatorKey)
    }
    atoms_by_key = {atom.key: atom for atom in context.operator_atoms}
    seen_operator_identities = set()
    building = set()
    emissions_by_key = {}
    for atom in context.operator_atoms:
        operator = _get_or_make_operator_recursive(
            context,
            atom.key,
            atoms_by_key,
            inputs_by_operator,
            operator_type_to_class,
            operator_input_class,
            building,
        )
        if operator is None:
            continue
        if atom.key in used_operator_keys:
            emissions_by_key[atom.key] = operator
    emitted = set()

    def _emit(operator_key: typing.Any) -> None:
        if operator_key in emitted:
            return
        emitted.add(operator_key)
        for input_key in inputs_by_operator.get(operator_key, ()):
            if isinstance(input_key, pd2af.asp.predicates.logicalOperatorKey):
                _emit(input_key)
        operator = context.key_to_operator.get(operator_key)
        if operator is None:
            return
        if add_model_element_if_new(
            model_operators, operator, seen_operator_identities
        ):
            input_operator = context.clingo_id_to_model_element[operator_key.gate]
            context.operator_emissions.append((operator, input_operator))

    for operator_key in emissions_by_key:
        _emit(operator_key)


def _get_or_make_operator_recursive(
    context: pd2af.building.context.BuilderContext,
    operator_key: typing.Any,
    atoms_by_key: dict,
    inputs_by_operator: dict,
    operator_type_to_class: dict[str, type],
    operator_input_class: type,
    building: set,
) -> typing.Any:
    """Resolve an operator key to its model element, building on demand.

    An operator's inputs may themselves be operators (SBGN-PD nesting), whose
    model elements may not be built yet; those resolve recursively through
    this same function. ``key_to_operator`` memoises, and ``building`` refuses
    to re-enter a key currently being built -- an input cycle would otherwise
    recurse forever (the ASP layer never emits such an influence, but a cyclic
    operator graph can emit nodes).
    """
    if operator_key in context.key_to_operator:
        return context.key_to_operator[operator_key]
    if operator_key in building:
        return None
    atom = atoms_by_key.get(operator_key)
    if atom is None:
        return None
    building.add(operator_key)
    try:
        operator = _get_or_make_operator(
            context,
            atom.type_,
            inputs_by_operator.get(operator_key, ()),
            operator_type_to_class,
            operator_input_class,
            atoms_by_key,
            inputs_by_operator,
            building,
        )
    finally:
        building.discard(operator_key)
    if operator is not None:
        context.key_to_operator[operator_key] = operator
    return operator


def _get_or_make_operator(
    context: pd2af.building.context.BuilderContext,
    operator_type: str,
    input_keys: collections.abc.Iterable[typing.Any],
    operator_type_to_class: dict[str, type],
    operator_input_class: type,
    atoms_by_key: dict,
    inputs_by_operator: dict,
    building: set,
) -> typing.Any:
    """Build (and intern) an operator of ``operator_type``.

    An activity-key input resolves through ``context.key_to_activity``; an
    operator-key input resolves recursively through
    :func:`_get_or_make_operator_recursive`, so a nested operator is built
    before its consumer. An unknown token yields ``None``.

    Each operator input and the operator itself are interned by content
    (``register_or_reuse``), so content-equal operators collapse to one
    canonical instance -- preserving the dedup-and-remap invariant for the
    influences that reference them.
    """
    operator_class = operator_type_to_class.get(operator_type)
    if operator_class is None:
        return None
    operator_inputs = []
    for input_key in input_keys:
        if isinstance(input_key, pd2af.asp.predicates.logicalOperatorKey):
            referred_element = _get_or_make_operator_recursive(
                context,
                input_key,
                atoms_by_key,
                inputs_by_operator,
                operator_type_to_class,
                operator_input_class,
                building,
            )
        else:
            referred_element = context.key_to_activity.get(input_key)
        if referred_element is None:
            continue
        operator_input = register_or_reuse(
            operator_input_class(referred_element=referred_element), context.cache
        )
        operator_inputs.append(operator_input)
    operator = operator_class(inputs=frozenset(operator_inputs))
    return register_or_reuse(operator, context.cache)
