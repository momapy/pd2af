"""Build the SBGN-AF model from clingo activity / influence atoms.

``make_and_add_model`` is the AF model pass: it walks the activity atoms
(``kept_species`` / ``promoted_subunit`` keys, each resolving to an input
SBGN-PD entity pool, phenotype or promoted subunit) and the influence atoms,
and populates ``context.model`` with canonical, content-deduped compartments,
activities and influences.

Each entity pool becomes a :class:`BiologicalActivity` carrying a typed
:class:`UnitOfInformation` (the entity class) and a label that is the canonical
serialization of the whole entity pool (:mod:`pd2af.sbgn.labels`). In the
merged modes (``normal``/``no-complex``) the label is built with state
variables stripped, so distinct proteoforms collapse to one merged activity;
otherwise they stay distinct under content-based model equality. A PD
:class:`Phenotype` process becomes an AF :class:`Phenotype` activity. Dedup is
honoured by interning every constructed element through the shared content
cache (``register_or_reuse``) and resolving influence endpoints through the
deduped activities, mirroring the model-element dedup invariant.
"""

import momapy.builder
import momapy.sbgn.af
import momapy.sbgn.pd

import pd2af.languages
import pd2af.predicates
import pd2af.sbgn.labels
from pd2af.utils import add_model_element_if_new, register_or_reuse


_ENTITY_CLASS_TO_UNIT_OF_INFORMATION_CLASS = {
    momapy.sbgn.pd.Macromolecule: momapy.sbgn.af.MacromoleculeUnitOfInformation,
    momapy.sbgn.pd.MacromoleculeMultimer: momapy.sbgn.af.MacromoleculeUnitOfInformation,
    momapy.sbgn.pd.NucleicAcidFeature: momapy.sbgn.af.NucleicAcidFeatureUnitOfInformation,
    momapy.sbgn.pd.NucleicAcidFeatureMultimer: momapy.sbgn.af.NucleicAcidFeatureUnitOfInformation,
    momapy.sbgn.pd.SimpleChemical: momapy.sbgn.af.SimpleChemicalUnitOfInformation,
    momapy.sbgn.pd.SimpleChemicalMultimer: momapy.sbgn.af.SimpleChemicalUnitOfInformation,
    momapy.sbgn.pd.Complex: momapy.sbgn.af.ComplexUnitOfInformation,
    momapy.sbgn.pd.Multimer: momapy.sbgn.af.ComplexUnitOfInformation,
    momapy.sbgn.pd.ComplexMultimer: momapy.sbgn.af.ComplexUnitOfInformation,
    momapy.sbgn.pd.UnspecifiedEntity: momapy.sbgn.af.UnspecifiedEntityUnitOfInformation,
    momapy.sbgn.pd.PerturbingAgent: momapy.sbgn.af.PerturbationUnitOfInformation,
    # Subunit classes (a promoted subunit becomes a top-level activity in the
    # no-complex modes) map to the same unit of information as their entity-pool
    # counterpart.
    momapy.sbgn.pd.MacromoleculeSubunit: momapy.sbgn.af.MacromoleculeUnitOfInformation,
    momapy.sbgn.pd.MacromoleculeMultimerSubunit: momapy.sbgn.af.MacromoleculeUnitOfInformation,
    momapy.sbgn.pd.NucleicAcidFeatureSubunit: momapy.sbgn.af.NucleicAcidFeatureUnitOfInformation,
    momapy.sbgn.pd.NucleicAcidFeatureMultimerSubunit: momapy.sbgn.af.NucleicAcidFeatureUnitOfInformation,
    momapy.sbgn.pd.SimpleChemicalSubunit: momapy.sbgn.af.SimpleChemicalUnitOfInformation,
    momapy.sbgn.pd.SimpleChemicalMultimerSubunit: momapy.sbgn.af.SimpleChemicalUnitOfInformation,
    momapy.sbgn.pd.ComplexSubunit: momapy.sbgn.af.ComplexUnitOfInformation,
    momapy.sbgn.pd.ComplexMultimerSubunit: momapy.sbgn.af.ComplexUnitOfInformation,
    momapy.sbgn.pd.UnspecifiedEntitySubunit: momapy.sbgn.af.UnspecifiedEntityUnitOfInformation,
}

# An unrecognised entity class falls back to an unspecified-entity unit of
# information rather than failing the whole transformation.
_FALLBACK_UNIT_OF_INFORMATION_CLASS = (
    momapy.sbgn.af.UnspecifiedEntityUnitOfInformation
)

# Influence predicate -> AF influence class. SBGN-PD only ever emits the four
# left-hand kinds (it has no "unknown" modulation twins), but the unknown
# variants are mapped too so the table is total over the predicate set.
_INFLUENCE_PREDICATE_TO_AF_CLASS = {
    pd2af.predicates.positivelyInfluences: momapy.sbgn.af.PositiveInfluence,
    pd2af.predicates.negativelyInfluences: momapy.sbgn.af.NegativeInfluence,
    pd2af.predicates.triggers: momapy.sbgn.af.NecessaryStimulation,
    pd2af.predicates.modulates: momapy.sbgn.af.UnknownInfluence,
    pd2af.predicates.unknownPositivelyInfluences: momapy.sbgn.af.PositiveInfluence,
    pd2af.predicates.unknownNegativelyInfluences: momapy.sbgn.af.NegativeInfluence,
    pd2af.predicates.unknownTriggers: momapy.sbgn.af.NecessaryStimulation,
    pd2af.predicates.unknownModulates: momapy.sbgn.af.UnknownInfluence,
}

_INFLUENCE_PREDICATE_CLASSES = tuple(_INFLUENCE_PREDICATE_TO_AF_CLASS)

# Operator-type token (``logicalOperator.type_``) -> SBGN-AF operator class.
# The NOT token is ``not_`` because bare ``not`` is a reserved clingo keyword.
# SBGN-PD authors only AND/OR/NOT operators (no unknown/delay), so the table
# covers exactly the tokens the ``sbgn_pd`` rule variant can emit.
_OPERATOR_TYPE_TO_OPERATOR_CLASS = {
    "and": momapy.sbgn.af.AndOperator,
    "or": momapy.sbgn.af.OrOperator,
    "not_": momapy.sbgn.af.NotOperator,
}


def make_and_add_model(context, clingo_model):
    context.model = momapy.builder.get_or_make_builder_cls(
        momapy.sbgn.af.SBGNAFModel
    )()
    _collect_atoms(context, clingo_model)
    _build_subunit_compartment_map(context)
    _make_and_add_compartments(context)
    _make_and_add_activities(context)
    _make_and_add_operators(context)
    _make_and_add_influences(context)


def _build_subunit_compartment_map(context):
    """Map ``id(subunit)`` -> the compartment of its enclosing complex.

    SBGN-PD subunit classes have no ``compartment`` field, so a promoted
    subunit activity would otherwise get ``compartment=None`` and never merge
    with a top-level twin. Mirror the CellDesigner builder: a subunit inherits
    its parent complex's compartment. Built once by walking the input model's
    complexes recursively (nested complexes pass their compartment down)."""
    context.subunit_id_to_parent_compartment = {}

    def walk(entity, inherited_compartment):
        compartment = getattr(entity, "compartment", None) or inherited_compartment
        for subunit in getattr(entity, "subunits", None) or ():
            context.subunit_id_to_parent_compartment[id(subunit)] = compartment
            walk(subunit, compartment)

    for entity in context.input_map.model.entity_pools:
        walk(entity, None)


def _collect_atoms(context, clingo_model):
    context.activity_atoms = []
    for atom in clingo_model.query(pd2af.predicates.new).all():
        payload = atom.object_
        if isinstance(payload, pd2af.predicates.activity):
            context.activity_atoms.append(payload)
        elif isinstance(payload, _INFLUENCE_PREDICATE_CLASSES):
            context.influence_atoms.append(payload)
        elif isinstance(payload, pd2af.predicates.logicalOperator):
            context.operator_atoms.append(payload)
        elif isinstance(payload, pd2af.predicates.logicalOperatorInput):
            context.operator_input_atoms.append(payload)


def _input_element_for_key(context, key):
    return context.clingo_id_to_model_element[key.species]


def _make_and_add_compartments(context):
    seen_compartment_identities = set()
    for atom in context.activity_atoms:
        input_element = _input_element_for_key(context, atom.key)
        input_compartment = getattr(input_element, "compartment", None)
        if input_compartment is None:
            input_compartment = context.subunit_id_to_parent_compartment.get(
                id(input_element)
            )
        if input_compartment is None:
            continue
        af_compartment = _get_or_make_compartment(context, input_compartment)
        add_model_element_if_new(
            context.model.compartments,
            af_compartment,
            seen_compartment_identities,
        )


def _get_or_make_compartment(context, input_compartment):
    canonical = context.input_compartment_to_af_compartment.get(
        id(input_compartment)
    )
    if canonical is None:
        candidate = momapy.sbgn.af.Compartment(label=input_compartment.label)
        canonical = register_or_reuse(candidate, context.cache)
        context.input_compartment_to_af_compartment[id(input_compartment)] = (
            canonical
        )
    # Reverse index so the layout pass can recover the input compartment (and
    # thus its glyph) from the canonical AF compartment. Keyed by id(canonical):
    # if several input compartments dedup to one AF compartment by label, the
    # last registered wins -- fine for plain-mode geometry.
    context.af_compartment_to_input_compartment[id(canonical)] = input_compartment
    return canonical


def _make_and_add_activities(context):
    strip = context.mode in pd2af.languages.MERGED_PROTEOFORM_MODES
    seen_activity_identities = set()
    for atom in context.activity_atoms:
        if atom.key in context.key_to_activity:
            continue
        input_element = _input_element_for_key(context, atom.key)
        activity = _make_activity(context, input_element, strip=strip)
        context.key_to_activity[atom.key] = activity
        if add_model_element_if_new(
            context.model.activities, activity, seen_activity_identities
        ):
            context.activity_emissions.append((activity, input_element))


def _make_activity(context, input_element, strip=False):
    """Build (and intern) the AF activity for ``input_element``.

    ``strip=True`` (the merged modes ``normal``/``no-complex``) drops state
    variables from the label so distinct proteoforms collapse into one merged
    activity under content-based model equality. ``strip=False`` keeps the full
    label, so distinct proteoforms stay distinct (keep-species behaviour)."""
    if isinstance(input_element, momapy.sbgn.pd.Phenotype):
        candidate = momapy.sbgn.af.Phenotype(label=input_element.label)
        return register_or_reuse(candidate, context.cache)
    unit_of_information_class = _ENTITY_CLASS_TO_UNIT_OF_INFORMATION_CLASS.get(
        type(input_element), _FALLBACK_UNIT_OF_INFORMATION_CLASS
    )
    unit_of_information = register_or_reuse(
        unit_of_information_class(), context.cache
    )
    input_compartment = getattr(input_element, "compartment", None)
    if input_compartment is None:
        # A subunit has no compartment field; inherit its parent complex's so a
        # promoted subunit can merge with a top-level twin (parity).
        input_compartment = context.subunit_id_to_parent_compartment.get(
            id(input_element)
        )
    compartment = None
    if input_compartment is not None:
        compartment = _get_or_make_compartment(context, input_compartment)
    candidate = momapy.sbgn.af.BiologicalActivity(
        label=pd2af.sbgn.labels.build_label(
            input_element, include_state_variables=not strip
        ),
        compartment=compartment,
        units_of_information=frozenset([unit_of_information]),
    )
    return register_or_reuse(candidate, context.cache)


def _make_and_add_operators(context):
    """Build a ``LogicalOperator`` for every authored operator and add those
    that actually source an influence to ``model.logical_operators``.

    Mirrors the CellDesigner gate builder: every operator is built into
    ``key_to_operator`` (so the influence pass can resolve an operator source),
    but only operators that source an influence are added to the model and
    recorded for the layout pass -- an operator whose target is not an activity
    yields no influence and would otherwise be a dangling node."""
    inputs_by_operator = {}
    for input_atom in context.operator_input_atoms:
        inputs_by_operator.setdefault(input_atom.operator, []).append(
            input_atom.input
        )
    used_operator_keys = {
        atom.source
        for atom in context.influence_atoms
        if isinstance(atom.source, pd2af.predicates.logical_operator_key)
    }
    seen_operator_identities = set()
    for atom in context.operator_atoms:
        operator = get_or_make_operator(
            atom.type_, inputs_by_operator.get(atom.key, ()), context
        )
        if operator is None:
            continue
        context.key_to_operator[atom.key] = operator
        if atom.key not in used_operator_keys:
            continue
        if add_model_element_if_new(
            context.model.logical_operators, operator, seen_operator_identities
        ):
            input_operator = context.clingo_id_to_model_element[atom.key.gate]
            context.operator_emissions.append((operator, input_operator))


def get_or_make_operator(operator_type, input_keys, context):
    """Build (and intern) an SBGN-AF ``LogicalOperator`` of ``operator_type``
    whose inputs resolve through ``key_to_activity``. Returns ``None`` for an
    unknown token. Each ``LogicalOperatorInput`` and the operator itself are
    interned by content (``register_or_reuse``), so content-equal operators
    collapse to one canonical instance."""
    operator_class = _OPERATOR_TYPE_TO_OPERATOR_CLASS.get(operator_type)
    if operator_class is None:
        return None
    operator_inputs = []
    for input_key in input_keys:
        activity = context.key_to_activity.get(input_key)
        if activity is None:
            continue
        operator_input = register_or_reuse(
            momapy.sbgn.af.LogicalOperatorInput(element=activity), context.cache
        )
        operator_inputs.append(operator_input)
    operator = operator_class(inputs=frozenset(operator_inputs))
    return register_or_reuse(operator, context.cache)


def _resolve_influence_source(context, source_key):
    """Resolve an influence ``source`` key to its model element: a logical
    operator resolves through ``key_to_operator`` (``None`` if not built); any
    activity key resolves through ``key_to_activity`` (``None`` if not emitted)."""
    if isinstance(source_key, pd2af.predicates.logical_operator_key):
        return context.key_to_operator.get(source_key)
    return context.key_to_activity.get(source_key)


def _make_and_add_influences(context):
    seen_influence_identities = set()
    for atom in context.influence_atoms:
        source = _resolve_influence_source(context, atom.source)
        target = context.key_to_activity.get(atom.target)
        if source is None or target is None:
            continue
        influence_class = _INFLUENCE_PREDICATE_TO_AF_CLASS[type(atom)]
        candidate = influence_class(source=source, target=target)
        canonical = register_or_reuse(candidate, context.cache)
        add_model_element_if_new(
            context.model.influences, canonical, seen_influence_identities
        )
