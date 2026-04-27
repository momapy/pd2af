import tempfile
import pathlib
import subprocess
import os.path
import typing

import momapy.io.core
import momapy.builder
import momapy.celldesigner

import pd2af.layouts


def _get_id_to_layout_element(cd_map):
    id_to_layout_element = {}
    for layout_element in cd_map.layout.layout_elements:
        id_to_layout_element[layout_element.id_] = layout_element
    return id_to_layout_element


sif_predicate_to_influence_class = {
    "POSITIVE": momapy.celldesigner.PositiveInfluence,
    "NEGATIVE": momapy.celldesigner.NegativeInfluence,
}


def _get_active_species(cd_map, sif_relation_tuples, id_to_layout_element):
    active_species_layout_ids = set([])
    for sif_relation_tuple in sif_relation_tuples:
        active_species_layout_ids.add(sif_relation_tuple[0])
        active_species_layout_ids.add(sif_relation_tuple[-1])
    active_species_layout_elements = [
        id_to_layout_element[active_species_layout_id]
        for active_species_layout_id in active_species_layout_ids
    ]
    active_species = [
        cd_map.layout_model_mapping.get_mapping(active_species_layout_element)
        for active_species_layout_element in active_species_layout_elements
    ]
    return active_species


def _make_influences(cd_map, sif_relation_tuples, id_to_layout_element):
    influences = []
    for sif_relation_tuple in sif_relation_tuples:
        source_layout_element = id_to_layout_element[sif_relation_tuple[0]]
        source = cd_map.layout_model_mapping.get_mapping(source_layout_element)
        target_layout_element = id_to_layout_element[sif_relation_tuple[-1]]
        target = cd_map.layout_model_mapping.get_mapping(target_layout_element)
        influence_class = sif_predicate_to_influence_class[sif_relation_tuple[1]]
        influence = influence_class(source=source, target=target)
        influences.append(influence)
    return influences


def _make_new_cd_model(cd_map, sif_relations):
    cd_model_builder_cls = momapy.builder.get_or_make_builder_cls(
        momapy.celldesigner.CellDesignerModel
    )
    cd_model_builder = cd_model_builder_cls()
    sif_relation_tuples = [
        tuple(sif_relation.split(" ")) for sif_relation in sif_relations
    ]
    id_to_layout_element = _get_id_to_layout_element(cd_map)
    active_species = _get_active_species(
        cd_map, sif_relation_tuples, id_to_layout_element
    )
    cd_model_builder.species = type(cd_model_builder.species)(active_species)
    compartments = set(
        [
            species.compartment
            for species in active_species
            if species.compartment is not None
        ]
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
            for species in active_species
            if hasattr(species, "template") and species.template is not None
        ]
    )
    cd_model_builder.species_templates = type(cd_model_builder.species_templates)(
        species_templates
    )
    influences = _make_influences(cd_map, sif_relation_tuples, id_to_layout_element)
    cd_model_builder.modulations = type(cd_model_builder.modulations)(influences)
    cd_model = momapy.builder.object_from_builder(cd_model_builder)
    return cd_model


def transform(
    input_file_path,
    layout_mode: typing.Literal["plain", "overlay", "auto"] | None = "plain",
):
    cd_map = momapy.io.core.read(input_file_path).obj
    _, output_file_path = tempfile.mkstemp(suffix=".sbml")
    output_file_path = pathlib.Path(output_file_path)
    output_file_parent = output_file_path.parent
    output_file_stem = output_file_path.stem
    output_raw_sif_file_stem = f"{output_file_stem}_raw"
    output_sif_file_name = pathlib.Path(output_raw_sif_file_stem).with_suffix(".sif")
    output_sif_file_path = os.path.join(output_file_parent, output_sif_file_name)
    command = ["casq", "--sif", input_file_path, output_file_path]
    subprocess.call(command)
    with open(output_sif_file_path) as f:
        sif_relations = [
            line.rstrip("\n") for line in f.readlines() if not line.startswith("#")
        ]
    new_cd_model = _make_new_cd_model(cd_map, sif_relations)
    return pd2af.layouts.build_map(cd_map, new_cd_model, layout_mode)
