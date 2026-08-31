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


def _make_control(
    model,
    clingo_id_to_model_element,
    mode,
    language,
    set_active,
    set_inactive,
    set_all_active,
    set_all_inactive,
    exclude_groups=(),
    exclude_rules=(),
):
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
    control.add(
        "base",
        [],
        pd2af.rules.build_program(
            mode,
            language,
            exclude_groups=exclude_groups,
            exclude_rules=exclude_rules,
        ),
    )
    control.add_facts(fact_base)
    _add_activity_override_facts(
        control,
        clingo_id_to_model_element,
        set_active,
        set_inactive,
        set_all_active,
        set_all_inactive,
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


def _resolve_ids_to_generated_constants(ids, id_to_generated_constant, option_name):
    """Validate `ids` name known species / entity pools and return their constants.

    Raises `ValueError` (naming `option_name`, e.g. `--set-active`) if any id is
    unknown. Returns the list of generated ASP constants, ready to interpolate
    bare into a fact.
    """
    unknown_ids = [
        identifier for identifier in ids if identifier not in id_to_generated_constant
    ]
    if unknown_ids:
        raise ValueError(
            f"{option_name} ids not found as a species or entity pool: "
            f"{sorted(unknown_ids)}"
        )
    return [id_to_generated_constant[identifier] for identifier in ids]


def _add_activity_override_facts(
    control,
    clingo_id_to_model_element,
    set_active,
    set_inactive,
    set_all_active,
    set_all_inactive,
):
    """Inject the user's activity overrides as ASP facts.

    Four options, in strict precedence per-id > global > rules:

    * `--set-active` ids become `hasActivityCandidate(ELEMENT,
      isInputParameter)` *and* `forceActive(ELEMENT)`. The candidate makes
      the element active on its own; `forceActive` shields it from a
      concurrent `--set-all-inactive`.
    * `--set-inactive` ids become `suppressActivity(ELEMENT)`, which the
      bridging rule reads to block *any* activity for that element (it wins
      over `--set-all-active`, whose candidate is still suppressed).
    * `--set-all-active` emits the single fact `globalActivate.`, turning
      every top-level species / entity pool into a candidate.
    * `--set-all-inactive` emits the single fact `globalSuppress.`,
      suppressing every candidate not shielded by `forceActive`.

    An id passed to both `--set-active` and `--set-inactive` is
    contradictory and raises, as is asking for both global toggles at once.
    """
    set_active = set_active or []
    set_inactive = set_inactive or []
    if set_all_active and set_all_inactive:
        raise ValueError("cannot pass both --set-all-active and --set-all-inactive")
    conflicting_ids = set(set_active) & set(set_inactive)
    if conflicting_ids:
        raise ValueError(
            f"ids passed to both --set-active and --set-inactive: "
            f"{sorted(conflicting_ids)}"
        )
    facts = []
    if set_all_active:
        facts.append("globalActivate.")
    if set_all_inactive:
        facts.append("globalSuppress.")
    if set_active or set_inactive:
        id_to_generated_constant = _build_id_to_generated_constant(
            clingo_id_to_model_element
        )
        for generated_constant in _resolve_ids_to_generated_constants(
            set_active, id_to_generated_constant, "--set-active"
        ):
            facts.append(
                f"hasActivityCandidate({generated_constant}, isInputParameter)."
            )
            facts.append(f"forceActive({generated_constant}).")
        for generated_constant in _resolve_ids_to_generated_constants(
            set_inactive, id_to_generated_constant, "--set-inactive"
        ):
            facts.append(f"suppressActivity({generated_constant}).")
    if not facts:
        return
    control.add("base", [], "\n".join(facts))


def solve(
    map_,
    mode,
    set_active=None,
    set_inactive=None,
    set_all_active=False,
    set_all_inactive=False,
    exclude_groups=(),
    exclude_rules=(),
):
    clingo_id_to_model_element = {}
    language = pd2af.languages.get_language_from_map_or_model(map_)
    control = _make_control(
        map_.model,
        clingo_id_to_model_element,
        mode=mode,
        language=language,
        set_active=set_active,
        set_inactive=set_inactive,
        set_all_active=set_all_active,
        set_all_inactive=set_all_inactive,
        exclude_groups=exclude_groups,
        exclude_rules=exclude_rules,
    )
    control.ground([("base", [])])
    clingo_models = []
    control.solve(on_model=lambda model: clingo_models.append(model.facts(atoms=True)))
    return clingo_models[0], clingo_id_to_model_element
