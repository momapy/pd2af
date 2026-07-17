"""Clingo glue: build a control, solve, return the model + element registry."""

import clorm
import clorm.clingo
import clingo.ast

import momapy.celldesigner
import momapy.sbgn.pd
import momapy_kb.clingo.core

import pd2af.languages
import pd2af.ontology
import pd2af.predicates
import pd2af.rules


def _get_profile_from_mode(mode):
    return mode.replace("-", "_")


def _make_control(
    model, clingo_id_to_model_element, mode, language, active_ids, inactive_ids
):
    profile = _get_profile_from_mode(mode)
    control = clorm.clingo.Control(
        ["--warn=no-atom-undefined"],
        unifier=[pd2af.predicates.new],
    )
    with momapy_kb.clingo.core.Session() as session:
        ontology_rules = pd2af.ontology.make_rules(session, language)
        facts = session.make_facts_from_object(
            model, id_to_object=clingo_id_to_model_element
        )
    fact_base = clorm.FactBase(facts)
    with clingo.ast.ProgramBuilder(control) as control_builder:
        for ontology_rule in ontology_rules:
            clingo.ast.parse_string(ontology_rule, control_builder.add)
    control.add("base", [], pd2af.rules.build_program(profile, language))
    control.add_facts(fact_base)
    _add_activity_override_facts(
        control, clingo_id_to_model_element, active_ids, inactive_ids
    )
    return control


def _build_id_to_generated_constant(clingo_id_to_model_element):
    """Reverse `clingo_id_to_model_element` to `id_ -> generated ASP constant`.

    Only species / entity pools are kept, since those are the elements a user
    can mark active or inactive.
    """
    return {
        model_element.id_: generated_constant
        for generated_constant, model_element in clingo_id_to_model_element.items()
        if isinstance(
            model_element,
            (momapy.celldesigner.Species, momapy.sbgn.pd.EntityPool),
        )
    }


def _resolve_ids_to_generated_constants(
    ids, id_to_generated_constant, option_name
):
    """Validate `ids` name known species / entity pools and return their constants.

    Raises `ValueError` (naming `option_name`, e.g. `--active`) if any id is
    unknown. Returns the list of generated ASP constants, ready to interpolate
    bare into a fact.
    """
    unknown_ids = [
        identifier
        for identifier in ids
        if identifier not in id_to_generated_constant
    ]
    if unknown_ids:
        raise ValueError(
            f"{option_name} ids not found as a species or entity pool: "
            f"{sorted(unknown_ids)}"
        )
    return [id_to_generated_constant[identifier] for identifier in ids]


def _add_activity_override_facts(
    control, clingo_id_to_model_element, active_ids, inactive_ids
):
    """Inject the user's `--active` / `--inactive` activity overrides.

    `--active` ids become `hasActivityCandidate(ELEMENT, isInputParameter)`
    atoms (a candidate, so the `--inactive` veto can still suppress them).
    `--inactive` ids become `suppressActivity(ELEMENT)` atoms, which the
    bridging rule in `activity_base` reads to block *any* activity for that
    element. An id passed to both options is contradictory and raises.
    """
    active_ids = active_ids or []
    inactive_ids = inactive_ids or []
    conflicting_ids = set(active_ids) & set(inactive_ids)
    if conflicting_ids:
        raise ValueError(
            f"ids passed to both --active and --inactive: "
            f"{sorted(conflicting_ids)}"
        )
    if not active_ids and not inactive_ids:
        return
    id_to_generated_constant = _build_id_to_generated_constant(
        clingo_id_to_model_element
    )
    facts = []
    if active_ids:
        for generated_constant in _resolve_ids_to_generated_constants(
            active_ids, id_to_generated_constant, "--active"
        ):
            facts.append(
                f"hasActivityCandidate({generated_constant}, isInputParameter)."
            )
    if inactive_ids:
        for generated_constant in _resolve_ids_to_generated_constants(
            inactive_ids, id_to_generated_constant, "--inactive"
        ):
            facts.append(f"suppressActivity({generated_constant}).")
    control.add("base", [], "\n".join(facts))


def solve(map_, mode, active_ids=None, inactive_ids=None):
    clingo_id_to_model_element = {}
    language = pd2af.languages.language_from_map(map_)
    control = _make_control(
        map_.model,
        clingo_id_to_model_element,
        mode=mode,
        language=language,
        active_ids=active_ids,
        inactive_ids=inactive_ids,
    )
    control.ground([("base", [])])
    clingo_models = []
    control.solve(on_model=lambda model: clingo_models.append(model.facts(atoms=True)))
    return clingo_models[0], clingo_id_to_model_element
