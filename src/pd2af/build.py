import momapy.builder
import momapy.celldesigner

import pd2af.dedup
import pd2af.layouts
import pd2af.model


def build_map(map_, layout_mode, clingo_model, clingo_id_to_model_element):
    resolution = pd2af.model.resolve(
        clingo_model, clingo_id_to_model_element, map_
    )
    layout = pd2af.layouts.STRATEGIES[layout_mode]()

    model_builder = momapy.builder.get_or_make_builder_cls(
        momapy.celldesigner.CellDesignerModel
    )()
    layout_builder, mapping_builder = layout.start(map_)

    for compartment in resolution.compartments:
        model_builder.compartments.add(compartment)
        layout.on_compartment(map_, compartment, layout_builder, mapping_builder)

    for template in resolution.templates:
        model_builder.species_templates.add(template)

    for species, existing_species, is_subunit in resolution.species:
        if not is_subunit:
            model_builder.species.add(species)
        layout.on_species(
            map_,
            species,
            existing_species,
            is_subunit,
            layout_builder,
            mapping_builder,
        )

    layout.on_species_done(map_, layout_builder, mapping_builder)

    for modulation in resolution.modulations:
        model_builder.modulations.add(modulation)
        layout.on_modulation(modulation, layout_builder, mapping_builder)

    pd2af.dedup.dedup_and_remap_model(model_builder, mapping_builder)
    return layout.finish(map_, model_builder, layout_builder, mapping_builder)
