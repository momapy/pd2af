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


def _make_control(model, clingo_id_to_model_element, mode, language, active_ids):
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
    _add_input_parameter_facts(control, clingo_id_to_model_element, active_ids)
    return control


def _add_input_parameter_facts(control, clingo_id_to_model_element, active_ids):
    """Inject `hasActivity(ELEMENT, isInputParameter)` for user-declared ids.

    `active_ids` are the original `id_` strings of species / entity pools the
    user marked active. We resolve each to the generated ASP constant by
    reversing `clingo_id_to_model_element` (`generated_constant -> object`),
    validate that every id names a species or entity pool, and inject the
    matching `hasActivity` atom. The generated constant is already an ASP
    constant token, so it is interpolated bare.
    """
    if not active_ids:
        return
    id_to_generated_constant = {
        model_element.id_: generated_constant
        for generated_constant, model_element in clingo_id_to_model_element.items()
        if isinstance(
            model_element,
            (momapy.celldesigner.Species, momapy.sbgn.pd.EntityPool),
        )
    }
    unknown_ids = [
        active_id
        for active_id in active_ids
        if active_id not in id_to_generated_constant
    ]
    if unknown_ids:
        raise ValueError(
            f"--active ids not found as a species or entity pool: "
            f"{sorted(unknown_ids)}"
        )
    input_parameter_facts = "\n".join(
        f"hasActivity({id_to_generated_constant[active_id]}, isInputParameter)."
        for active_id in active_ids
    )
    control.add("base", [], input_parameter_facts)


def solve(map_, mode, active_ids=None):
    clingo_id_to_model_element = {}
    language = pd2af.languages.language_from_map(map_)
    control = _make_control(
        map_.model,
        clingo_id_to_model_element,
        mode=mode,
        language=language,
        active_ids=active_ids,
    )
    control.ground([("base", [])])
    clingo_models = []
    control.solve(on_model=lambda model: clingo_models.append(model.facts(atoms=True)))
    return clingo_models[0], clingo_id_to_model_element
