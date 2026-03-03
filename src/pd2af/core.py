import typing
import importlib.resources
import collections

import clorm
import clorm.clingo
import clingo.ast

import momapy.core
import momapy.geometry
import momapy.builder
import momapy.celldesigner.core

import momapy_kb.clingo.core
import momapy_kb.clingo.celldesigner

import pd2af.cd2asp
import pd2af.utils

_mode_to_file_name = {
    "pd2af": "pd2af.lp",
    "casq": "casq.lp",
    "pd2af-no-complex": "pd2af_no_complex.lp",
}


class activity(clorm.Predicate):
    name: clorm.ConstantStr


pos_path = clorm.simple_predicate("posPath", 2)
neg_path = clorm.simple_predicate("negPath", 2)


class positivelyInfluences(clorm.Predicate):
    source: clorm.ConstantStr
    target: clorm.ConstantStr


class negativelyInfluences(clorm.Predicate):
    source: clorm.ConstantStr
    target: clorm.ConstantStr


class triggers(clorm.Predicate):
    source: clorm.ConstantStr
    target: clorm.ConstantStr


class new(clorm.Predicate):
    object_: activity | positivelyInfluences | negativelyInfluences | triggers


predicate_to_model_element_class = {
    positivelyInfluences: momapy.celldesigner.core.PositiveInfluence,
    negativelyInfluences: momapy.celldesigner.core.Inhibition,
    triggers: momapy.celldesigner.core.Triggering,
}

model_element_class_to_layout_element_class = {
    momapy.celldesigner.core.PositiveInfluence: momapy.celldesigner.core.PositiveInfluenceLayout,
    momapy.celldesigner.core.NegativeInfluence: momapy.celldesigner.core.InhibitionLayout,
    momapy.celldesigner.core.Inhibition: momapy.celldesigner.core.InhibitionLayout,
    momapy.celldesigner.core.Triggering: momapy.celldesigner.core.TriggeringLayout,
}


def _get_layout_elements_from_model_element(
    map_, model_element
):  # TODO: delete when issue is solved in momapy
    layout_elements = map_.layout_model_mapping.get_mapping(model_element)
    if layout_elements is None:
        for key in map_.layout_model_mapping.inverse.keys():
            if isinstance(key, tuple) and key[0] == model_element:
                layout_elements = map_.layout_model_mapping.get_mapping(key)
                break
    return layout_elements


def _cd_model_to_facts(cd_model, id_to_model_element):
    return momapy_kb.clingo.core.make_facts_from_object(cd_model, id_to_model_element)


def _make_fact_base_from_cd_model(cd_model, id_to_model_element):
    facts = _cd_model_to_facts(cd_model, id_to_model_element)
    return clorm.FactBase(facts)


def _make_control_from_cd_model(cd_model, id_to_model_element, mode, active):
    control = clorm.clingo.Control(
        ["--warn=no-atom-undefined"],
        unifier=[new, neg_path, pos_path],
    )
    fact_base = _make_fact_base_from_cd_model(cd_model, id_to_model_element)
    with clingo.ast.ProgramBuilder(control) as control_builder:
        for ontology_rule in momapy_kb.clingo.celldesigner.ontology_rules:
            clingo.ast.parse_string(ontology_rule, control_builder.add)
    file_name = _mode_to_file_name.get(mode)
    if file_name is None:
        raise ValueError(f"mode {mode} is not supported")
    program_path = importlib.resources.files("pd2af.data") / file_name
    control.load(str(program_path))
    control.add_facts(fact_base)
    return control


def _solve_from_cd_map(cd_map, id_to_model_element, mode, active):
    cd_model = cd_map.model
    control = _make_control_from_cd_model(
        cd_model, id_to_model_element, mode=mode, active=active
    )
    control.ground([("base", [])])
    models = []
    control.solve(on_model=lambda model: models.append(model.facts(atoms=True)))
    model = models[0]
    return model


def _get_activity_atoms(model):
    activity_atoms = [
        atom.object_
        for atom in model.query(new).all()
        if isinstance(
            atom.object_,
            (activity,),
        )
    ]
    return activity_atoms


def _get_influence_atoms(model):
    influence_atoms = [
        atom.object_
        for atom in model.query(new).all()
        if isinstance(
            atom.object_,
            (positivelyInfluences, negativelyInfluences),
        )
    ]
    return influence_atoms


def _get_model_element_ids_to_layout_elements_from_model_element_ids(ids, map_):
    model_element_ids_to_layout_elements = collections.defaultdict(list)
    ids = set(ids)
    for model_element in map_.layout_model_mapping.inverse.keys():
        if isinstance(model_element, tuple):
            model_element = model_element[0]
        if model_element.id_ in ids:
            layout_elements = map_.layout_model_mapping.get_mapping(model_element)
            model_element_ids_to_layout_elements[model_element.id_] += layout_elements
    return model_element_ids_to_layout_elements


def _make_layout_elements_from_influence_atoms(
    influence_atoms, active_species_ids_to_layout_elements
):
    layout_elements = []
    for influence_atom in influence_atoms:
        source_species_id = influence_atom.source
        target_species_id = influence_atom.target
        for source_layout_element in active_species_ids_to_layout_elements[
            source_species_id
        ]:
            for target_layout_element in active_species_ids_to_layout_elements[
                target_species_id
            ]:
                model_element_class = predicate_to_model_element_class[
                    type(influence_atom)
                ]
                layout_element_class = model_element_class_to_layout_element_class[
                    model_element_class
                ]
                segment = momapy.geometry.Segment(
                    source_layout_element.border(target_layout_element.center()),
                    target_layout_element.border(source_layout_element.center()),
                )
                arc = layout_element_class(
                    source=source_layout_element,
                    target=target_layout_element,
                    segments=[segment],
                )
                layout_elements.append(arc)
    return layout_elements


def _get_id_to_model_element_from_ids(map_, ids):
    id_to_model_element = {}
    for model_element in map_.layout_model_mapping.inverse.keys():
        if isinstance(model_element, tuple):
            model_element = model_element[0]
        if model_element.id_ in ids:
            id_to_model_element[model_element.id_] = model_element
    return id_to_model_element


def _make_influence_id_to_model_element_from_atoms(
    influence_atoms, id_to_model_element
):
    influence_id_to_model_element = {}
    for influence_atom in influence_atoms:
        influence_model_element_class = predicate_to_model_element_class[
            type(influence_atom)
        ]
        source_model_element = id_to_model_element[influence_atom.source]
        target_model_element = id_to_model_element[influence_atom.target]
        influence_model_element = influence_model_element_class(
            source=source_model_element,
            target=target_model_element,
        )
        influence_id_to_model_element[influence_model_element.id_] = (
            influence_model_element
        )
    return influence_id_to_model_element


def _make_new_cd_model_from_clingo_model(cd_map, clingo_model, id_to_model_element):
    # for _ in clingo_model:
    #     print(_)
    cd_model_builder = momapy.celldesigner.core.CellDesignerModelBuilder()
    activity_atoms = _get_activity_atoms(clingo_model)
    active_species_ids = [activity_atom.name for activity_atom in activity_atoms]
    species = [
        id_to_model_element[active_species_id]
        for active_species_id in active_species_ids
    ]
    cd_model_builder.species = type(cd_model_builder.species)(species)
    compartments = set(
        [species.compartment for species in species if species.compartment is not None]
    )
    compartments_to_check = compartments
    while True:
        compartments_to_check = set(
            [
                compartment.outside
                for compartment in compartments_to_check
                if compartment not in compartments_to_check
            ]
        )
        if not compartments_to_check:
            break
        compartments |= compartments_to_check
    cd_model_builder.compartments = type(cd_model_builder.compartments)(compartments)
    species_templates = set(
        [
            species.template
            for species in species
            if hasattr(species, "template") and species.template is not None
        ]
    )
    cd_model_builder.species_templates = type(cd_model_builder.species_templates)(
        species_templates
    )
    influence_atoms = _get_influence_atoms(clingo_model)
    influence_id_to_model_element = _make_influence_id_to_model_element_from_atoms(
        influence_atoms, id_to_model_element
    )
    cd_model_builder.modulations = type(cd_model_builder.modulations)(
        influence_id_to_model_element.values()
    )
    cd_model = momapy.builder.object_from_builder(cd_model_builder)
    return cd_model


def _make_modulation_layout_element_tuples(cd_map, modulation):
    layout_element_tuples = []
    source_model_element = modulation.source
    target_model_element = modulation.target
    source_layout_elements = _get_layout_elements_from_model_element(
        cd_map, source_model_element
    )
    target_layout_elements = _get_layout_elements_from_model_element(
        cd_map, target_model_element
    )
    for source_layout_element in source_layout_elements:
        for target_layout_element in target_layout_elements:
            layout_element_tuple = _make_modulation_layout_element_tuple(
                modulation, source_layout_element, target_layout_element
            )
            layout_element_tuples.append(layout_element_tuple)
    return layout_element_tuples


def _make_modulation_layout_element_tuple(
    modulation, source_layout_element, target_layout_element
):
    layout_element_class = model_element_class_to_layout_element_class[type(modulation)]
    segment = momapy.geometry.Segment(
        source_layout_element.border(target_layout_element.center()),
        target_layout_element.border(source_layout_element.center()),
    )
    layout_element = layout_element_class(
        source=source_layout_element,
        target=target_layout_element,
        segments=tuple([segment]),
    )
    layout_element_tuple = (
        layout_element,
        source_layout_element,
        target_layout_element,
    )
    return layout_element_tuple


def _make_new_overlay_cd_map(cd_map, new_cd_model):
    new_layout = cd_map.layout
    new_layout_model_mapping_builder = momapy.core.LayoutModelMappingBuilder()
    new_layout_elements = []
    layout_elements_to_ignore = set()
    for compartment in new_cd_model.compartments:
        compartment_layouts = cd_map.layout_model_mapping.get_mapping(compartment)
        if compartment_layouts is not None:
            for compartment_layout in compartment_layouts:
                new_layout_model_mapping_builder.add_mapping(
                    compartment_layout, compartment
                )
                layout_elements_to_ignore.add(compartment_layout)
    for species in new_cd_model.species:
        species_layouts = _get_layout_elements_from_model_element(cd_map, species)
        for species_layout in species_layouts:
            new_layout_model_mapping_builder.add_mapping(species_layout, species)
            layout_elements_to_ignore.add(species_layout)
    for modulation in new_cd_model.modulations:
        modulation_layout_sets = cd_map.layout_model_mapping.get_mapping(modulation)
        if modulation_layout_sets is None:
            modulation_layout_tuples = _make_modulation_layout_element_tuples(
                cd_map, modulation
            )
            for modulation_layout_tuple in modulation_layout_tuples:
                new_layout_elements.append(modulation_layout_tuple[0])
            modulation_layout_sets = [
                frozenset(modulation_layout_tuple)
                for modulation_layout_tuple in modulation_layout_tuples
            ]
        for modulation_layout_set in modulation_layout_sets:
            new_layout_model_mapping_builder.add_mapping(
                modulation_layout_set, modulation
            )
            layout_elements_to_ignore.update(modulation_layout_set)
    new_layout_builder = momapy.builder.builder_from_object(new_layout)
    new_layout_builder.layout_elements += new_layout_elements
    builder_to_object = {}
    new_layout = momapy.builder.object_from_builder(
        new_layout_builder, builder_to_object=builder_to_object
    )
    new_layout = pd2af.utils.highlight_layout_elements(
        layout_elements_to_ignore, new_layout
    )
    new_layout_model_mapping = momapy.builder.object_from_builder(
        new_layout_model_mapping_builder, builder_to_object=builder_to_object
    )
    new_map = momapy.celldesigner.core.CellDesignerMap(
        model=new_cd_model,
        layout=new_layout,
        layout_model_mapping=new_layout_model_mapping,
    )
    return new_map


def _make_new_auto_cd_map(cd_map, new_cd_model):
    new_layout_builder = momapy.core.LayoutBuilder()
    new_layout_model_mapping_builder = momapy.core.LayoutModelMappingBuilder()
    for compartment in new_cd_model.compartments:
        compartment_layouts = cd_map.layout_model_mapping.get_mapping(compartment)
        if compartment_layouts is not None:
            for compartment_layout in compartment_layouts:
                new_layout_builder.layout_elements.append(compartment_layout)
                new_layout_model_mapping_builder.add_mapping(
                    compartment_layout, compartment
                )
                break  # we add just one
    species_to_layout_element = {}  # we need it to build modulation layout elements with correct source and target
    for species in new_cd_model.species:
        species_layouts = cd_map.layout_model_mapping.get_mapping(species)
        for species_layout in species_layouts:
            new_layout_builder.layout_elements.append(species_layout)
            new_layout_model_mapping_builder.add_mapping(species_layout, species)
            species_to_layout_element[species] = species_layout
            break  # we add just one
    for modulation in new_cd_model.modulations:
        source_layout_element = species_to_layout_element[modulation.source]
        target_layout_element = species_to_layout_element[modulation.target]
        modulation_layout_tuple = _make_modulation_layout_element_tuple(
            modulation, source_layout_element, target_layout_element
        )
        new_layout_builder.layout_elements.append(modulation_layout_tuple[0])
        modulation_layout_set = frozenset(modulation_layout_tuple)
        new_layout_model_mapping_builder.add_mapping(modulation_layout_set, modulation)
    builder_to_object = {}
    new_layout = momapy.builder.object_from_builder(
        new_layout_builder, builder_to_object=builder_to_object
    )
    new_layout_model_mapping = momapy.builder.object_from_builder(
        new_layout_model_mapping_builder, builder_to_object=builder_to_object
    )
    new_map = momapy.celldesigner.core.CellDesignerMap(
        model=new_cd_model,
        layout=new_layout,
        layout_model_mapping=new_layout_model_mapping,
    )
    new_map = pd2af.utils.auto_layout(new_map)
    return new_map


def _make_new_cd_map_no_layout(new_cd_model):
    new_map = momapy.celldesigner.core.CellDesignerMap(
        model=new_cd_model,
    )
    return new_map


def transform_map(
    cd_map,
    mode: typing.Literal["pd2af", "casq", "pd2af-no-complex"] = "pd2af",
    with_layout=True,
    layout_mode: typing.Literal["overlay", "auto", "all"] = "overlay",
    active: list[str] | None = None,
):
    if active is None:
        active = []
    new_maps = []
    id_to_model_element = {}
    clingo_model = _solve_from_cd_map(
        cd_map, id_to_model_element=id_to_model_element, mode=mode, active=active
    )
    new_cd_model = _make_new_cd_model_from_clingo_model(
        cd_map, clingo_model, id_to_model_element
    )
    if with_layout:
        if layout_mode == "overlay" or layout_mode == "all":
            new_overlay_cd_map = _make_new_overlay_cd_map(cd_map, new_cd_model)
            new_maps.append(new_overlay_cd_map)
        if layout_mode == "auto" or layout_mode == "all":
            new_auto_cd_map = _make_new_auto_cd_map(cd_map, new_cd_model)
            new_maps.append(new_auto_cd_map)
    else:
        new_map_no_layout = _make_new_cd_map_no_layout(new_cd_model)
        new_maps.append(new_map_no_layout)
    return new_maps


def _make_new_casq_cd_map(cd_map, new_cd_model):
    new_layout = cd_map.layout
    new_layout_model_mapping_builder = momapy.core.LayoutModelMappingBuilder()
    new_layout_elements = []
    layout_elements_to_ignore = set()
    for compartment in new_cd_model.compartments:
        compartment_layouts = cd_map.layout_model_mapping.get_mapping(compartment)
        if compartment_layouts is not None:
            for compartment_layout in compartment_layouts:
                new_layout_model_mapping_builder.add_mapping(
                    compartment_layout, compartment
                )
                layout_elements_to_ignore.add(compartment_layout)
    for species in new_cd_model.species:
        species_layouts = cd_map.layout_model_mapping.get_mapping(species)
        for species_layout in species_layouts:
            new_layout_model_mapping_builder.add_mapping(species_layout, species)
            layout_elements_to_ignore.add(species_layout)
    for modulation in new_cd_model.modulations:
        modulation_layout_sets = cd_map.layout_model_mapping.get_mapping(modulation)
        if modulation_layout_sets is None:
            modulation_layout_tuples = _make_modulation_layout_element_tuples(
                cd_map, modulation
            )
            for modulation_layout_tuple in modulation_layout_tuples:
                new_layout_elements.append(modulation_layout_tuple[0])
            modulation_layout_sets = [
                frozenset(modulation_layout_tuple)
                for modulation_layout_tuple in modulation_layout_tuples
            ]
        for modulation_layout_set in modulation_layout_sets:
            new_layout_model_mapping_builder.add_mapping(
                modulation_layout_set, modulation
            )
            layout_elements_to_ignore.update(modulation_layout_set)
    new_layout_builder = momapy.builder.builder_from_object(new_layout)
    new_layout_builder.layout_elements += new_layout_elements
    builder_to_object = {}
    new_layout = momapy.builder.object_from_builder(
        new_layout_builder, builder_to_object=builder_to_object
    )
    new_layout = pd2af.utils.highlight_layout_elements(
        layout_elements_to_ignore, new_layout
    )
    new_layout_model_mapping = momapy.builder.object_from_builder(
        new_layout_model_mapping_builder, builder_to_object=builder_to_object
    )
    new_map = momapy.celldesigner.core.CellDesignerMap(
        model=new_cd_model,
        layout=new_layout,
        layout_model_mapping=new_layout_model_mapping,
    )
    return new_map
