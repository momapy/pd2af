import clorm

import momapy.celldesigner


class kept_species(clorm.Predicate):
    """Identity wrapper for an AF entity that refers to a species
    carried over from the input PD map.

    Used for every activity in `keep-species` and `keep-species-no-complex`
    modes, and for templateless species (phenotypes, ions, simple
    molecules, surviving complexes) in `normal` and `no-complex` modes.
    The single argument is the original CellDesigner species ID — the
    solver looks it up in `id_to_model_element` to recover the `Species`
    object.

    Emitted by the `activity_key:kept` rule group (keep-species and
    keep-species-no-complex) and by `activity_key:merged` for
    templateless species; consumed by `solver._make_new_cd_model` and
    `solver._make_influences`.
    """
    species: clorm.ConstantStr


class derived_proteoform_class(clorm.Predicate):
    """Identity wrapper for an AF entity synthesized by merging
    proteoforms of the same template in the same compartment.

    Used in `normal` and `no-complex` modes. Two proteoforms of
    template `P` both in compartment `C` map to
    `derived_proteoform_class(P, C)` and collapse via clingo's set
    semantics. `compartment` may be the sentinel constant
    `no_compartment` for templated species without a compartment. The
    solver synthesizes a fresh `Species` for each such atom (see
    `solver._make_synthetic_species`).

    Emitted by the `activity_key:merged` rule group; consumed by
    `solver._make_new_cd_model` and `solver._make_influences`.
    """
    template: clorm.ConstantStr
    compartment: clorm.ConstantStr


class activity(clorm.Predicate):
    """An activity node in the new AF map.

    Always wrapped by `new(...)` in rule heads. The `key` field carries
    the activity's identity — either a `kept_species` (an input species
    promoted to an activity) or a `derived_proteoform_class` (a merged
    proteoform class, in `normal` and `no-complex` modes).

    Emitted by the shared `activity_derivation` rule group; consumed
    by `solver.make_new_cd_model` to populate `cd_model.species`.
    """
    key: kept_species | derived_proteoform_class


class positivelyInfluences(clorm.Predicate):
    """A positive-influence edge in the new AF map (source activates
    target).

    Always wrapped by `new(...)` in rule heads. `source` and `target`
    are activity keys (same wrapper types as `activity.key`). Emitted
    by `influences_from_paths` (positive paths between two activities)
    and `influences_consumption` (inhibitor-spares-reactant rule);
    consumed by `solver._make_influences` to build `PositiveInfluence`
    model elements.
    """
    source: kept_species | derived_proteoform_class
    target: kept_species | derived_proteoform_class


class negativelyInfluences(clorm.Predicate):
    """A negative-influence (inhibition) edge in the new AF map
    (source inhibits target).

    Always wrapped by `new(...)` in rule heads. `source` and `target`
    are activity keys (same wrapper types as `activity.key`). Emitted
    by `influences_from_paths` (negative paths between two activities)
    and by `influences_consumption` (catalyzer/physical-stimulator/
    trigger consume reactant); consumed by `solver._make_influences`
    to build `Inhibition` model elements.
    """
    source: kept_species | derived_proteoform_class
    target: kept_species | derived_proteoform_class


class new(clorm.Predicate):
    """Top-level marker: this fact belongs to the *new* AF map being
    constructed (as opposed to facts about the input PD map).

    The Python solver queries clingo for `new(_)` atoms and dispatches
    on the wrapped predicate type to build the new `CellDesignerModel`.
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
