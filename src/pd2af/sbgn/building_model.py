"""Build the SBGN-AF model from clingo activity / influence atoms.

``make_and_add_model`` is the AF model pass: it walks the activity atoms
(``keptSpeciesKey`` / ``promotedSubunitKey`` keys, each resolving to an input
SBGN-PD entity pool, phenotype or promoted subunit) and the influence atoms,
and populates ``context.model`` with canonical, content-deduped compartments,
activities and influences.

Each entity pool becomes a :class:`BiologicalActivity` carrying a typed
:class:`UnitOfInformation` (the entity class) and a label that is the canonical
serialization of the whole entity pool (:mod:`pd2af.sbgn.building_labels`). In
the merged modes (``normal``/``normal-no-complex``) the label is built with state
variables stripped, so distinct proteoforms collapse to one merged activity,
and the entity's unit-of-information block is moved off the label onto the typed
:class:`UnitOfInformation` glyph (a curator's AF map carries ``ct:mRNA`` on the
glyph -- no brackets -- not ``[ct:mRNA]`` in the label); otherwise the label
keeps both blocks and
distinct proteoforms stay distinct under content-based model equality. A PD
:class:`Phenotype` process becomes an AF :class:`Phenotype` activity. Dedup is
honoured by interning every constructed element through the shared content
cache (``register_or_reuse``) and resolving influence endpoints through the
deduped activities, mirroring the model-element dedup invariant.
"""

import typing

import momapy.builder
import momapy.sbgn.af
import momapy.sbgn.pd

import pd2af.building_model
import pd2af.context
import pd2af.predicates
import pd2af.sbgn.building_labels
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
    # normal-no-complex modes) map to the same unit of information as their entity-pool
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
_FALLBACK_UNIT_OF_INFORMATION_CLASS = momapy.sbgn.af.UnspecifiedEntityUnitOfInformation

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

# Operator-type token (``logicalOperator.type_``) -> SBGN-AF operator class.
# The NOT token is ``not_`` because bare ``not`` is a reserved clingo keyword.
# SBGN-PD authors only AND/OR/NOT operators (no unknown/delay), so the table
# covers exactly the tokens the ``sbgn_pd`` rule variant can emit.
_OPERATOR_TYPE_TO_OPERATOR_CLASS = {
    "and": momapy.sbgn.af.AndOperator,
    "or": momapy.sbgn.af.OrOperator,
    "not_": momapy.sbgn.af.NotOperator,
}


def make_and_add_model(context: pd2af.context.BuilderContext, clingo_model: typing.Any):
    """Build ``context.model`` from the clingo atoms (pass 1)."""
    context.model = momapy.builder.get_or_make_builder_cls(momapy.sbgn.af.SBGNAFModel)()
    pd2af.building_model.collect_atoms(context, clingo_model)
    context.subunit_to_top_level = pd2af.building_model.build_subunit_to_top_level(
        context.input_map.model.entity_pools
    )
    _make_and_add_compartments(context)
    _make_and_add_activities(context)
    pd2af.building_model.make_and_add_operators(
        context,
        _OPERATOR_TYPE_TO_OPERATOR_CLASS,
        momapy.sbgn.af.LogicalOperatorInput,
        context.model.logical_operators,
    )
    _make_and_add_influences(context)


def _compartment_for_input_element(
    context: pd2af.context.BuilderContext, input_element: typing.Any
) -> typing.Any:
    """The compartment an element belongs to.

    SBGN-PD subunit classes have no ``compartment`` field, so a subunit inherits
    its top-level entity pool's: otherwise a promoted subunit activity would get
    ``compartment=None`` and never merge with a top-level twin.
    """
    compartment = getattr(input_element, "compartment", None)
    if compartment is not None:
        return compartment
    top_level_element = context.subunit_to_top_level.get(id(input_element))
    return getattr(top_level_element, "compartment", None)


def _input_element_for_key(
    context: pd2af.context.BuilderContext, key: typing.Any
) -> typing.Any:
    return context.clingo_id_to_model_element[key.species]


def _make_and_add_compartments(context: pd2af.context.BuilderContext):
    seen_compartment_identities = set()
    for atom in context.activity_atoms:
        input_element = _input_element_for_key(context, atom.key)
        input_compartment = _compartment_for_input_element(context, input_element)
        if input_compartment is None:
            continue
        af_compartment = _get_or_make_compartment(context, input_compartment)
        # A PD compartment becomes a distinct AF compartment; record the pair
        # for provenance so its annotations/notes carry (several PD compartments
        # may dedup to one AF compartment by label -- the bucket union handles it).
        context.compartment_emissions.append((input_compartment, af_compartment))
        add_model_element_if_new(
            context.model.compartments,
            af_compartment,
            seen_compartment_identities,
        )


def _get_or_make_compartment(
    context: pd2af.context.BuilderContext, input_compartment: typing.Any
) -> typing.Any:
    canonical = context.input_compartment_to_af_compartment.get(id(input_compartment))
    if canonical is None:
        candidate = momapy.sbgn.af.Compartment(label=input_compartment.label)
        canonical = register_or_reuse(candidate, context.cache)
        context.input_compartment_to_af_compartment[id(input_compartment)] = canonical
    # Reverse index so the layout pass can recover the input compartment (and
    # thus its glyph) from the canonical AF compartment. Keyed by id(canonical):
    # if several input compartments dedup to one AF compartment by label, the
    # last registered wins -- fine for plain-mode geometry.
    context.af_compartment_to_input_compartment[id(canonical)] = input_compartment
    return canonical


def _make_and_add_activities(context: pd2af.context.BuilderContext):
    strip = context.mode.merges_proteoforms
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


def _make_activity(
    context: pd2af.context.BuilderContext,
    input_element: typing.Any,
    strip: bool = False,
) -> typing.Any:
    """Build (and intern) the AF activity for ``input_element``.

    ``strip=True`` (the merged modes ``normal``/``normal-no-complex``) drops state
    variables from the label so distinct proteoforms collapse into one merged
    activity under content-based model equality, and relocates the entity's
    unit-of-information block off the label and onto the typed unit-of-information
    glyph (where a curator drawing AF from scratch would put it -- ``ct:mRNA`` on
    a nucleic-acid-feature glyph, no brackets, rather than inline in the label).
    Distinct
    units still keep activities distinct: ``UnitOfInformation.label`` is part of
    its content, so the dedup granularity is unchanged, only relocated.
    ``strip=False`` keeps the full label and a bare typed glyph, so distinct
    proteoforms stay distinct (keep-species behaviour).
    """
    if isinstance(input_element, momapy.sbgn.pd.Phenotype):
        candidate = momapy.sbgn.af.Phenotype(label=input_element.label)
        return register_or_reuse(candidate, context.cache)
    unit_of_information_class = _ENTITY_CLASS_TO_UNIT_OF_INFORMATION_CLASS.get(
        type(input_element), _FALLBACK_UNIT_OF_INFORMATION_CLASS
    )
    unit_of_information_label = (
        pd2af.sbgn.building_labels.make_units_of_information_label(input_element)
        if strip
        else None
    )
    unit_of_information = register_or_reuse(
        unit_of_information_class(label=unit_of_information_label), context.cache
    )
    input_compartment = _compartment_for_input_element(context, input_element)
    compartment = None
    if input_compartment is not None:
        compartment = _get_or_make_compartment(context, input_compartment)
    candidate = momapy.sbgn.af.BiologicalActivity(
        label=pd2af.sbgn.building_labels.make_label(
            input_element,
            include_state_variables=not strip,
            include_units_of_information=not strip,
        ),
        compartment=compartment,
        units_of_information=frozenset([unit_of_information]),
    )
    return register_or_reuse(candidate, context.cache)


def _make_and_add_influences(context: pd2af.context.BuilderContext):
    seen_influence_identities = set()
    for atom in context.influence_atoms:
        source = pd2af.building_model.resolve_influence_source(context, atom.source)
        target = context.key_to_activity.get(atom.target)
        if source is None or target is None:
            continue
        influence_class = _INFLUENCE_PREDICATE_TO_AF_CLASS[type(atom)]
        candidate = influence_class(source=source, target=target)
        canonical = register_or_reuse(candidate, context.cache)
        add_model_element_if_new(
            context.model.influences, canonical, seen_influence_identities
        )
