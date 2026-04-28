import dataclasses

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
    "normal": "default",
    "no-complex": "no_complex",
    "pure-af": "pure_af",
}

_NO_COMPARTMENT_SENTINEL = "no_compartment"

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
    if mode not in _MODE_TO_PROFILE:
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
        id_=f"pure_af_template__{template_id}",
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
            f"pure-af synthesized key references unknown template id {template_id!r}"
        )
    species_cls = _TEMPLATE_TO_SPECIES_CLASS.get(type(template))
    if species_cls is None:
        raise ValueError(
            f"pure-af mode does not know how to synthesize a species for template "
            f"class {type(template).__name__}"
        )
    if compartment_id == _NO_COMPARTMENT_SENTINEL:
        compartment = None
    else:
        compartment = id_to_model_element.get(compartment_id)
        if compartment is None:
            raise ValueError(
                f"pure-af synthesized key references unknown compartment id "
                f"{compartment_id!r}"
            )
    stripped_template = _stripped_template_for(
        template, template_id, stripped_template_cache
    )
    return species_cls(
        id_=f"pure_af__{template_id}__{compartment_id}",
        name=template.name,
        template=stripped_template,
        compartment=compartment,
    )


def _resolve_activity_key(key, id_to_model_element, key_to_species, stripped_template_cache):
    cached = key_to_species.get(key)
    if cached is not None:
        return cached
    if isinstance(key, pd2af.predicates.kept_species):
        species = id_to_model_element[key.species]
    elif isinstance(key, pd2af.predicates.derived_proteoform_class):
        species = _make_synthetic_species(
            key, id_to_model_element, stripped_template_cache
        )
    else:
        raise ValueError(f"unknown activity key wrapper {type(key).__name__}")
    key_to_species[key] = species
    return species


def _make_influences(influence_atoms, id_to_model_element, key_to_species, stripped_template_cache):
    influences = {}
    for atom in influence_atoms:
        cls = pd2af.predicates.predicate_to_model_element_class[type(atom)]
        source = _resolve_activity_key(
            atom.source,
            id_to_model_element,
            key_to_species,
            stripped_template_cache,
        )
        target = _resolve_activity_key(
            atom.target,
            id_to_model_element,
            key_to_species,
            stripped_template_cache,
        )
        influence = cls(source=source, target=target)
        influences[influence.id_] = influence
    return influences


def make_new_cd_model(clingo_model, id_to_model_element):
    cd_model_builder_cls = momapy.builder.get_or_make_builder_cls(
        momapy.celldesigner.CellDesignerModel
    )
    cd_model_builder = cd_model_builder_cls()
    activity_atoms = _get_activity_atoms(clingo_model)
    key_to_species = {}
    stripped_template_cache = {}
    species = [
        _resolve_activity_key(
            atom.key,
            id_to_model_element,
            key_to_species,
            stripped_template_cache,
        )
        for atom in activity_atoms
    ]
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
    influences = _make_influences(
        influence_atoms,
        id_to_model_element,
        key_to_species,
        stripped_template_cache,
    )
    cd_model_builder.modulations = type(cd_model_builder.modulations)(
        influences.values()
    )
    return momapy.builder.object_from_builder(cd_model_builder)
