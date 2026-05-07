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


_ACTIVITY_KEY = kept_species | new_species_from_template


class activity(clorm.Predicate):
    """An activity node in the new AF map.

    Always wrapped by ``new(...)`` in rule heads. The ``key`` field
    carries the activity's identity — one of the three activity-key
    wrappers.
    """
    key: _ACTIVITY_KEY


class positivelyInfluences(clorm.Predicate):
    """A positive-influence edge in the new AF map (source activates
    target). Always wrapped by ``new(...)`` in rule heads.
    """
    source: _ACTIVITY_KEY
    target: _ACTIVITY_KEY


class negativelyInfluences(clorm.Predicate):
    """A negative-influence (inhibition) edge in the new AF map
    (source inhibits target). Always wrapped by ``new(...)`` in rule
    heads.
    """
    source: _ACTIVITY_KEY
    target: _ACTIVITY_KEY


class new(clorm.Predicate):
    """Top-level marker: this fact belongs to the *new* AF map being
    constructed (as opposed to facts about the input PD map).
    """
    object_: activity | positivelyInfluences | negativelyInfluences


predicate_to_model_element_class = {
    positivelyInfluences: momapy.celldesigner.PositiveInfluence,
    negativelyInfluences: momapy.celldesigner.Inhibition,
}

model_element_class_to_layout_element_class = {
    momapy.celldesigner.PositiveInfluence: momapy.celldesigner.PositiveInfluenceLayout,
    momapy.celldesigner.Inhibition: momapy.celldesigner.InhibitionLayout,
}
