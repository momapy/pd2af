import clorm

import momapy.celldesigner


class kept_species(clorm.Predicate):
    """Activity-key wrapper: the AF activity reuses an input PD species
    by identity.

    The single argument is the original CellDesigner species ID. The
    solver looks it up in ``id_to_model_element`` to recover the
    ``Species`` object.
    """
    species: clorm.ConstantStr


class kept_subunit(clorm.Predicate):
    """Activity-key wrapper: a PD subunit of a kept complex contributes
    its own activity but is *not* added to ``model.species`` — it is
    already carried inside its parent complex's ``.subunits``.

    The single argument is the synthetic ASP ID of the subunit species.
    Emitted by ``normal`` and ``keep-species`` modes.
    """
    species: clorm.ConstantStr


class promoted_subunit(clorm.Predicate):
    """Activity-key wrapper: a PD subunit of a dissolved complex is
    promoted to a top-level activity (added to ``model.species``).

    The single argument is the synthetic ASP ID of the subunit species.
    Emitted by ``no-complex`` and ``keep-species-no-complex`` modes.
    """
    species: clorm.ConstantStr


class new_species_from_template(clorm.Predicate):
    """Activity-key wrapper: a PD species is replaced by a synthesized
    species whose template has been stripped of proteoform decorations
    (modification residues, regions) and whose compartment is the
    species' effective compartment.

    Used in ``normal`` and ``no-complex`` modes for templated species
    so content-equal proteoforms collapse via Python content interning
    at construction time.
    """
    species: clorm.ConstantStr


_ACTIVITY_KEY = (
    kept_species | kept_subunit | promoted_subunit | new_species_from_template
)


class activity(clorm.Predicate):
    """An activity node in the new AF map.

    Always wrapped by ``new(...)`` in rule heads. The ``key`` field
    carries the activity's identity — one of the three activity-key
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
    source: _ACTIVITY_KEY
    target: _ACTIVITY_KEY


class negativelyInfluences(clorm.Predicate):
    """A negative-influence edge (source inhibits target)."""
    source: _ACTIVITY_KEY
    target: _ACTIVITY_KEY


class modulates(clorm.Predicate):
    """A modulation edge: source influences target with an effect of
    unknown sign.
    """
    source: _ACTIVITY_KEY
    target: _ACTIVITY_KEY


class triggers(clorm.Predicate):
    """A triggering edge (necessary stimulation): source is required for
    target.
    """
    source: _ACTIVITY_KEY
    target: _ACTIVITY_KEY


class unknownPositivelyInfluences(clorm.Predicate):
    """A positive-influence edge whose existence is uncertain."""
    source: _ACTIVITY_KEY
    target: _ACTIVITY_KEY


class unknownNegativelyInfluences(clorm.Predicate):
    """A negative-influence edge whose existence is uncertain."""
    source: _ACTIVITY_KEY
    target: _ACTIVITY_KEY


class unknownModulates(clorm.Predicate):
    """A modulation edge whose existence is uncertain."""
    source: _ACTIVITY_KEY
    target: _ACTIVITY_KEY


class unknownTriggers(clorm.Predicate):
    """A triggering edge whose existence is uncertain."""
    source: _ACTIVITY_KEY
    target: _ACTIVITY_KEY


class new(clorm.Predicate):
    """Top-level marker: this fact belongs to the *new* AF map being
    constructed (as opposed to facts about the input PD map).
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
