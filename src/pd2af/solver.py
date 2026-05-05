"""Clingo glue: build a control, solve, return the model + element registry."""

import clorm
import clorm.clingo
import clingo.ast

import momapy_kb.clingo.core

import pd2af.ontology
import pd2af.predicates
import pd2af.rules


def _mode_to_profile(mode):
    return mode.replace("-", "_")


def _make_control(model, clingo_id_to_model_element, mode):
    profile = _mode_to_profile(mode)
    control = clorm.clingo.Control(
        ["--warn=no-atom-undefined"],
        unifier=[pd2af.predicates.new],
    )
    with momapy_kb.clingo.core.Session() as session:
        ontology_rules = pd2af.ontology.make_rules(session)
        facts = session.make_facts_from_object(
            model, id_to_object=clingo_id_to_model_element
        )
    fact_base = clorm.FactBase(facts)
    with clingo.ast.ProgramBuilder(control) as control_builder:
        for ontology_rule in ontology_rules:
            clingo.ast.parse_string(ontology_rule, control_builder.add)
    control.add("base", [], pd2af.rules.build_program(profile))
    control.add_facts(fact_base)
    return control


def solve(map_, mode):
    clingo_id_to_model_element = {}
    control = _make_control(
        map_.model, clingo_id_to_model_element, mode=mode
    )
    control.ground([("base", [])])
    clingo_models = []
    control.solve(
        on_model=lambda model: clingo_models.append(model.facts(atoms=True))
    )
    return clingo_models[0], clingo_id_to_model_element
