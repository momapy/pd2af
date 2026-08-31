"""The parts of the model pass that are the same for either language.

:mod:`pd2af.celldesigner.building_model` and :mod:`pd2af.sbgn.building_model`
build their own compartments, species and influences, but read the clingo atoms,
resolve an influence source and walk the input subunit trees the same way. The logical-operator pass is
shared too: the CellDesigner builder makes ``BooleanLogicGate`` objects and the
SBGN-AF builder ``LogicalOperator`` objects, and otherwise the two do the
same thing -- read the operator atoms, resolve each operator's inputs through
``key_to_activity``, and add to the model those operators that actually source
an influence. Each builder passes in the classes it wants and the model
collection to add to.
"""

import pd2af.predicates
from pd2af.utils import add_model_element_if_new, register_or_reuse


def collect_atoms(context, clingo_model):
    """Sort the ``new(...)`` atoms into the context scratch lists by payload
    type: activities, influences, logical operators and operator inputs."""
    for atom in clingo_model.query(pd2af.predicates.new).all():
        payload = atom.object_
        if isinstance(payload, pd2af.predicates.activity):
            context.activity_atoms.append(payload)
        elif isinstance(payload, pd2af.predicates.INFLUENCE_PREDICATES):
            context.influence_atoms.append(payload)
        elif isinstance(payload, pd2af.predicates.logicalOperator):
            context.operator_atoms.append(payload)
        elif isinstance(payload, pd2af.predicates.logicalOperatorInput):
            context.operator_input_atoms.append(payload)


def resolve_influence_source(context, source_key):
    """Resolve an influence ``source`` key to its model element: a logical
    operator resolves through ``key_to_operator``, any activity key through
    ``key_to_activity``. ``None`` when the gate or the activity was not built,
    so the caller skips the edge."""
    if isinstance(source_key, pd2af.predicates.logicalOperatorKey):
        return context.key_to_operator.get(source_key)
    return context.key_to_activity.get(source_key)


def build_subunit_to_top_level(top_level_elements):
    """Map ``id(subunit)`` -> the top-level element it belongs to, for every
    subunit of ``top_level_elements`` at any depth.

    A subunit is a structural component, never an independent activity, so both
    builders attribute it to its outermost element: the CellDesigner builder
    takes the top-level species' compartment, the SBGN-AF builder the top-level
    entity pool's.
    """
    subunit_to_top_level = {}

    def walk(element, top_level_element):
        for subunit in getattr(element, "subunits", None) or ():
            subunit_to_top_level[id(subunit)] = top_level_element
            walk(subunit, top_level_element)

    for element in top_level_elements:
        walk(element, element)
    return subunit_to_top_level


def make_and_add_operators(
    context, operator_type_to_class, operator_input_class, model_operators
):
    """Build an operator for every authored logical operator and add those that
    actually source an influence to ``model_operators``.

    Every operator is built into ``context.key_to_operator``, so the influence
    pass can resolve an operator source. Only the operators that appear as an
    influence source are added to the model and recorded in
    ``context.operator_emissions`` for the layout pass: an operator whose target
    is not an activity yields no influence and would otherwise be a dangling
    node. The CellDesigner writer emits a gate only through its modulation, so
    dropping it also keeps the output round-trip-safe.
    """
    inputs_by_operator = {}
    for input_atom in context.operator_input_atoms:
        inputs_by_operator.setdefault(input_atom.operator, []).append(
            input_atom.input
        )
    used_operator_keys = {
        atom.source
        for atom in context.influence_atoms
        if isinstance(atom.source, pd2af.predicates.logicalOperatorKey)
    }
    seen_operator_identities = set()
    for atom in context.operator_atoms:
        operator = _get_or_make_operator(
            context,
            atom.type_,
            inputs_by_operator.get(atom.key, ()),
            operator_type_to_class,
            operator_input_class,
        )
        if operator is None:
            continue
        context.key_to_operator[atom.key] = operator
        if atom.key not in used_operator_keys:
            continue
        if add_model_element_if_new(
            model_operators, operator, seen_operator_identities
        ):
            input_operator = context.clingo_id_to_model_element[atom.key.gate]
            context.operator_emissions.append((operator, input_operator))


def _get_or_make_operator(
    context,
    operator_type,
    input_keys,
    operator_type_to_class,
    operator_input_class,
):
    """Build (and intern) an operator of ``operator_type`` whose inputs resolve
    through ``context.key_to_activity``. Returns ``None`` for an unknown token.

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
        activity = context.key_to_activity.get(input_key)
        if activity is None:
            continue
        operator_input = register_or_reuse(
            operator_input_class(referred_element=activity), context.cache
        )
        operator_inputs.append(operator_input)
    operator = operator_class(inputs=frozenset(operator_inputs))
    return register_or_reuse(operator, context.cache)
