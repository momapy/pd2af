import clorm
import clorm.clingo
import clingo.ast

import momapy.builder
import momapy.celldesigner

import momapy_kb.clingo.core

import pd2af.ontology
import pd2af.predicates
import pd2af.rules


_MODE_TO_PROFILE = {
    "pd2af": "default",
    "pd2af-no-complex": "no_complex",
}


def _make_control(cd_model, id_to_model_element, mode):
    profile = _MODE_TO_PROFILE.get(mode)
    if profile is None:
        raise ValueError(f"mode {mode!r} is not supported")
    control = clorm.clingo.Control(
        ["--warn=no-atom-undefined"],
        unifier=[pd2af.predicates.new],
    )
    with momapy_kb.clingo.core.Session() as session:
        ontology_rules = pd2af.ontology.make_rules(session)
        facts = session.make_facts_from_object(
            cd_model, id_to_object=id_to_model_element
        )
    fact_base = clorm.FactBase(facts)
    with clingo.ast.ProgramBuilder(control) as control_builder:
        for ontology_rule in ontology_rules:
            clingo.ast.parse_string(ontology_rule, control_builder.add)
    control.add("base", [], pd2af.rules.build_program(profile))
    control.add_facts(fact_base)
    return control


def solve(cd_map, mode):
    id_to_model_element = {}
    control = _make_control(cd_map.model, id_to_model_element, mode=mode)
    control.ground([("base", [])])
    models = []
    control.solve(on_model=lambda model: models.append(model.facts(atoms=True)))
    return models[0], id_to_model_element


def _get_activity_atoms(clingo_model):
    return [
        atom.object_
        for atom in clingo_model.query(pd2af.predicates.new).all()
        if isinstance(atom.object_, pd2af.predicates.activity)
    ]


def _get_influence_atoms(clingo_model):
    return [
        atom.object_
        for atom in clingo_model.query(pd2af.predicates.new).all()
        if isinstance(
            atom.object_,
            (
                pd2af.predicates.positivelyInfluences,
                pd2af.predicates.negativelyInfluences,
            ),
        )
    ]


def _make_influences(influence_atoms, id_to_model_element):
    influences = {}
    for atom in influence_atoms:
        cls = pd2af.predicates.predicate_to_model_element_class[type(atom)]
        source = id_to_model_element[atom.source]
        target = id_to_model_element[atom.target]
        influence = cls(source=source, target=target)
        influences[influence.id_] = influence
    return influences


def make_new_cd_model(clingo_model, id_to_model_element):
    cd_model_builder_cls = momapy.builder.get_or_make_builder_cls(
        momapy.celldesigner.CellDesignerModel
    )
    cd_model_builder = cd_model_builder_cls()
    activity_atoms = _get_activity_atoms(clingo_model)
    species = [id_to_model_element[atom.name] for atom in activity_atoms]
    cd_model_builder.species = type(cd_model_builder.species)(species)
    compartments = set(
        s.compartment for s in species if s.compartment is not None
    )
    compartments_to_check = compartments
    while True:
        compartments_to_check = set(
            c.outside
            for c in compartments_to_check
            if c not in compartments_to_check
        )
        if not compartments_to_check:
            break
        compartments |= compartments_to_check
    cd_model_builder.compartments = type(cd_model_builder.compartments)(
        compartments
    )
    species_templates = set(
        s.template
        for s in species
        if hasattr(s, "template") and s.template is not None
    )
    cd_model_builder.species_templates = type(cd_model_builder.species_templates)(
        species_templates
    )
    influence_atoms = _get_influence_atoms(clingo_model)
    influences = _make_influences(influence_atoms, id_to_model_element)
    cd_model_builder.modulations = type(cd_model_builder.modulations)(
        influences.values()
    )
    return momapy.builder.object_from_builder(cd_model_builder)
