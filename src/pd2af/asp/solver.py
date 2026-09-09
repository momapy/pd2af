"""Clingo glue: build a control, solve, return the model + element registry."""

import collections.abc
import functools
import types
import typing

import clorm
import clorm.clingo
import clingo.ast

import momapy.celldesigner
import momapy.core.elements
import momapy.core.map
import momapy.core.model
import momapy.sbgn.pd
import momapy_kb.clingo.core

import pd2af.asp.predicates
import pd2af.asp.rules
import pd2af.modes


_ONTOLOGY_BASES = (momapy.core.elements.ModelElement, momapy.core.model.Model)


def _iter_types(module: types.ModuleType) -> collections.abc.Iterator[type]:
    for attr_name in dir(module):
        if attr_name.startswith("_"):
            continue
        attr_value = getattr(module, attr_name)
        if isinstance(attr_value, type) and issubclass(attr_value, _ONTOLOGY_BASES):
            yield attr_value


def _make_ontology_rules(
    session: typing.Any, language: pd2af.modes.Language
) -> list[str]:
    """The sorted ontology rules for every model class of ``language``."""
    module = pd2af.modes.LANGUAGES[language]["momapy_module"]
    rules = set()
    for type_ in _iter_types(module):
        session.get_or_make_predicate_classes_from_type(
            type_, make_predicate_classes_recursively=True
        )
        rules.update(session.make_ontology_rules_from_type(type_))
    return sorted(rules)


@functools.lru_cache(maxsize=None)
def _get_ontology_rules_for_language(
    language: pd2af.modes.Language,
) -> tuple[str, ...]:
    """The ontology rules of ``language``, generated once and reused.

    The rules depend on the language's model classes only, so they are cached
    as an immutable tuple. They are generated in their own short-lived session:
    a session also holds the element-id counter and the predicate classes of
    the objects it converted, and neither may be shared between transforms.
    """
    with momapy_kb.clingo.core.Session() as session:
        return tuple(_make_ontology_rules(session, language))


def _make_control(
    model: momapy.core.model.Model,
    clingo_id_to_model_element: dict,
    mode: pd2af.modes.TransformationMode,
    language: pd2af.modes.Language,
    set_active: list[str] | None,
    set_inactive: list[str] | None,
    set_all_active: bool,
    set_all_inactive: bool,
    exclude_groups: tuple[str, ...] = (),
    exclude_rules: tuple[str, ...] = (),
) -> clorm.clingo.Control:
    control = clorm.clingo.Control(
        ["--warn=no-atom-undefined"],
        unifier=[pd2af.asp.predicates.new],
    )
    ontology_rules = _get_ontology_rules_for_language(language)
    with momapy_kb.clingo.core.Session() as session:
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
        pd2af.asp.rules.build_program(
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


def _build_id_to_generated_constant(clingo_id_to_model_element: dict) -> dict[str, str]:
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
    ids: collections.abc.Iterable[str],
    id_to_generated_constant: dict[str, str],
    option_name: str,
) -> list[str]:
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
    control: clorm.clingo.Control,
    clingo_id_to_model_element: dict,
    set_active: list[str] | None,
    set_inactive: list[str] | None,
    set_all_active: bool,
    set_all_inactive: bool,
) -> None:
    """Inject the user's activity overrides as ASP facts.

    Four options, in strict precedence per-id > global > rules:

    * `--set-active` ids become `hasActivityCandidate(ELEMENT,
      isInputParameter)` *and* `forceActive(ELEMENT)`. The candidate makes
      the element active on its own; `forceActive` shields it from a
      concurrent `--set-all-inactive`.
    * `--set-inactive` ids become `suppressActivity(ELEMENT)`, which the
      bridging rule reads to block *any* activity for that element: the veto
      is blanket, so it drops every candidate the element has, whichever rule
      proposed it (active flag or state, modulation source, reaction modifier,
      gate input, phenotype), and it wins over `--set-all-active`, whose
      candidate is suppressed like the others.
    * `--set-all-active` emits the single fact `globalActivate.`, turning
      every top-level species / entity pool into a candidate.
    * `--set-all-inactive` emits the single fact `globalSuppress.`,
      suppressing every candidate not shielded by `forceActive`.

    An id passed to both `--set-active` and `--set-inactive` is
    contradictory and raises, as is asking for both global toggles at once;
    both raise here, before the program is grounded or solved.
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
    input_map: momapy.core.map.Map,
    mode: pd2af.modes.TransformationMode,
    set_active: list[str] | None = None,
    set_inactive: list[str] | None = None,
    set_all_active: bool = False,
    set_all_inactive: bool = False,
    exclude_groups: tuple[str, ...] = (),
    exclude_rules: tuple[str, ...] = (),
) -> tuple[clorm.FactBase, dict]:
    """Solve the ASP program for ``input_map`` in ``mode``.

    Returns the single clingo model of derived atoms together with the
    ``clingo_id -> model_element`` map the build pass resolves keys through.
    A satisfiable program deriving no atom returns an empty fact base; a
    program with no answer set raises ``ValueError``.
    """
    clingo_id_to_model_element = {}
    language = pd2af.modes.get_language_from_map_or_model(input_map)
    control = _make_control(
        input_map.model,
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
    solve_result = control.solve(
        on_model=lambda model: clingo_models.append(model.facts(atoms=True))
    )
    if not clingo_models:
        raise ValueError(
            f"the ASP program has no answer set (solve result: {solve_result}); "
            f"an unsatisfiable program cannot be built into a map"
        )
    return clingo_models[0], clingo_id_to_model_element
