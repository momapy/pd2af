import dataclasses

import clorm
import clorm.clingo
import clingo.ast

import momapy.celldesigner

import momapy_kb.clingo.core

import pd2af.ontology
import pd2af.predicates
import pd2af.rules


_VALID_MODES = frozenset(
    {"normal", "no-complex", "keep-species", "keep-species-no-complex", "casq"}
)

_NO_COMPARTMENT_SENTINEL = "no_compartment"

_SYNTHESIZED_ID_PREFIX = "merged__"


def _mode_to_profile(mode):
    return mode.replace("-", "_")


_TEMPLATE_TO_SPECIES_CLASS = {
    momapy.celldesigner.GenericProteinTemplate: momapy.celldesigner.GenericProtein,
    momapy.celldesigner.TruncatedProteinTemplate: momapy.celldesigner.TruncatedProtein,
    momapy.celldesigner.ReceptorTemplate: momapy.celldesigner.Receptor,
    momapy.celldesigner.IonChannelTemplate: momapy.celldesigner.IonChannel,
    momapy.celldesigner.GeneTemplate: momapy.celldesigner.Gene,
    momapy.celldesigner.RNATemplate: momapy.celldesigner.RNA,
    momapy.celldesigner.AntisenseRNATemplate: momapy.celldesigner.AntisenseRNA,
}


def _make_control(cd_model, id_to_model_element, mode):
    if mode not in _VALID_MODES:
        raise ValueError(f"mode {mode!r} is not supported")
    profile = _mode_to_profile(mode)
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
    if mode not in _VALID_MODES:
        raise ValueError(f"mode {mode!r} is not supported")
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


def _stripped_template_for(template, template_id, cache):
    cached = cache.get(template_id)
    if cached is not None:
        return cached
    fields_to_clear = {}
    if hasattr(template, "modification_residues"):
        fields_to_clear["modification_residues"] = frozenset()
    if hasattr(template, "regions"):
        fields_to_clear["regions"] = frozenset()
    stripped = dataclasses.replace(
        template,
        id_=f"merged_template__{template_id}",
        **fields_to_clear,
    )
    cache[template_id] = stripped
    return stripped


def _make_synthetic_species(key, id_to_model_element, stripped_template_cache):
    template_id = key.template
    compartment_id = key.compartment
    template = id_to_model_element.get(template_id)
    if template is None:
        raise ValueError(
            f"merged-proteoform key references unknown template id {template_id!r}"
        )
    species_cls = _TEMPLATE_TO_SPECIES_CLASS.get(type(template))
    if species_cls is None:
        raise ValueError(
            f"cannot synthesize a merged-proteoform species for template "
            f"class {type(template).__name__}"
        )
    if compartment_id == _NO_COMPARTMENT_SENTINEL:
        compartment = None
    else:
        compartment = id_to_model_element.get(compartment_id)
        if compartment is None:
            raise ValueError(
                f"merged-proteoform key references unknown compartment id "
                f"{compartment_id!r}"
            )
    stripped_template = _stripped_template_for(
        template, template_id, stripped_template_cache
    )
    return species_cls(
        id_=f"{_SYNTHESIZED_ID_PREFIX}{template_id}__{compartment_id}",
        name=template.name,
        template=stripped_template,
        compartment=compartment,
    )
