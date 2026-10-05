"""The clorm predicates the ASP program is written against."""

import clorm
import momapy.celldesigner


class keptSpeciesKey(clorm.Predicate):
    """Activity-key wrapper: the AF activity reuses an input PD species by identity.

    The single argument is the original CellDesigner species ID. The
    solver looks it up in ``id_to_model_element`` to recover the
    ``Species`` object.
    """

    species: clorm.ConstantStr


class promotedSubunitKey(clorm.Predicate):
    """Activity-key wrapper: a PD subunit of a dissolved complex is promoted.

    The subunit becomes a top-level activity, added to ``model.species``. The
    single argument is the synthetic ASP ID of the subunit species.
    Emitted by the ``no-complex`` mode.
    """

    species: clorm.ConstantStr


_ACTIVITY_KEY = keptSpeciesKey | promotedSubunitKey


class logicalOperatorKey(clorm.Predicate):
    """Activity-source wrapper: the AF influence is sourced by a logical operator.

    The operator (AND / OR / NOT / unknown) is authored in the input PD map. The
    single argument is the original gate ID -- a CellDesigner
    ``BooleanLogicGate`` or an SBGN-PD ``LogicalOperator``. The solver
    looks it up in ``clingo_id_to_model_element`` to recover the gate
    object. An operator is only ever an influence *source*, never a
    target, so it widens ``_INFLUENCE_SOURCE`` but not the activity key.
    """

    gate: clorm.ConstantStr


# The source of an influence edge: an activity (the two activity-key
# wrappers) or a logical operator. The target is always an activity, so
# only the influence ``source`` field is widened to this union.
_INFLUENCE_SOURCE = _ACTIVITY_KEY | logicalOperatorKey


class activity(clorm.Predicate):
    """An activity node in the new AF map.

    Always wrapped by ``new(...)`` in rule heads. The ``key`` field
    carries the activity's identity — one of the two activity-key
    wrappers.
    """

    key: _ACTIVITY_KEY


# Influence edges in the new AF map. There is one typed predicate per
# output influence class; the rules first converge on an internal
# ``influences(SOURCE, TARGET, KIND)`` relation (ASP-only, never wrapped
# in ``new``) and then fan it out to these typed heads. Each is always
# wrapped by ``new(...)`` in rule heads.


class positivelyInfluences(clorm.Predicate):
    """A positive-influence edge (source activates target)."""

    source: _INFLUENCE_SOURCE
    target: _ACTIVITY_KEY


class negativelyInfluences(clorm.Predicate):
    """A negative-influence edge (source inhibits target)."""

    source: _INFLUENCE_SOURCE
    target: _ACTIVITY_KEY


class modulates(clorm.Predicate):
    """A modulation edge: source influences target with an effect of unknown sign."""

    source: _INFLUENCE_SOURCE
    target: _ACTIVITY_KEY


class triggers(clorm.Predicate):
    """A triggering edge (necessary stimulation): source is required for target."""

    source: _INFLUENCE_SOURCE
    target: _ACTIVITY_KEY


class unknownPositivelyInfluences(clorm.Predicate):
    """A positive-influence edge whose existence is uncertain."""

    source: _INFLUENCE_SOURCE
    target: _ACTIVITY_KEY


class unknownNegativelyInfluences(clorm.Predicate):
    """A negative-influence edge whose existence is uncertain."""

    source: _INFLUENCE_SOURCE
    target: _ACTIVITY_KEY


class unknownModulates(clorm.Predicate):
    """A modulation edge whose existence is uncertain."""

    source: _INFLUENCE_SOURCE
    target: _ACTIVITY_KEY


class unknownTriggers(clorm.Predicate):
    """A triggering edge whose existence is uncertain."""

    source: _INFLUENCE_SOURCE
    target: _ACTIVITY_KEY


class logicalOperator(clorm.Predicate):
    """A logical operator node in the new AF map.

    Always wrapped by ``new(...)`` in rule heads. ``key`` carries the
    operator's identity; ``type_`` is the operator-type token
    (``and`` | ``or`` | ``not_`` | ``unknown``). The NOT token is spelt
    ``not_`` -- with a trailing underscore -- because bare ``not`` is a
    reserved clingo keyword (default negation) and cannot appear as a
    term. The builder maps the token to the momapy operator class per
    output language
    (``BooleanLogicGate`` subclass for CellDesigner, ``LogicalOperator``
    subclass for SBGN-AF). A single token field -- rather than one
    predicate per type -- mirrors the internal ``influences/3``
    KIND-token idiom and avoids same-name collisions with the input
    ontology's per-type operator functors.
    """

    key: logicalOperatorKey
    type_: clorm.ConstantStr


class logicalOperatorInput(clorm.Predicate):
    """An input edge of a logical operator: one activity feeding it.

    Always wrapped by ``new(...)`` in rule heads. The input is an activity
    key, or -- for SBGN-PD, where a logical operator's input may refer to
    another logical operator -- the key of that operator, mirroring how
    ``_INFLUENCE_SOURCE`` widens an influence source.
    """

    operator: logicalOperatorKey
    input: _ACTIVITY_KEY | logicalOperatorKey


class new(clorm.Predicate):
    """Top-level marker: this fact belongs to the *new* AF map being constructed.

    As opposed to facts about the input PD map.
    """

    object_: (
        activity
        | positivelyInfluences
        | negativelyInfluences
        | modulates
        | triggers
        | unknownPositivelyInfluences
        | unknownNegativelyInfluences
        | unknownModulates
        | unknownTriggers
        | logicalOperator
        | logicalOperatorInput
    )


# The typed influence predicates emitted into ``new(...)``.
INFLUENCE_PREDICATES = (
    positivelyInfluences,
    negativelyInfluences,
    modulates,
    triggers,
    unknownPositivelyInfluences,
    unknownNegativelyInfluences,
    unknownModulates,
    unknownTriggers,
)


predicate_to_model_element_class = {
    positivelyInfluences: momapy.celldesigner.PositiveInfluence,
    negativelyInfluences: momapy.celldesigner.NegativeInfluence,
    modulates: momapy.celldesigner.Modulation,
    triggers: momapy.celldesigner.Triggering,
    unknownPositivelyInfluences: momapy.celldesigner.UnknownPositiveInfluence,
    unknownNegativelyInfluences: momapy.celldesigner.UnknownNegativeInfluence,
    unknownModulates: momapy.celldesigner.UnknownModulation,
    unknownTriggers: momapy.celldesigner.UnknownTriggering,
}
