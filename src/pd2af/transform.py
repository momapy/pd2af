import typing

import momapy.builder
import momapy.celldesigner

import pd2af.layouts
import pd2af.model
import pd2af.solver
from pd2af.dedup import dedup_and_remap_model


_MERGED_PROTEOFORM_MODES = frozenset({"normal", "no-complex"})

_VALID_MODES = frozenset(
    {"normal", "no-complex", "keep-species", "keep-species-no-complex", "casq"}
)

_VALID_LAYOUT_MODES = frozenset({"plain", "overlay", "auto", "none", None})


def _normalize_layout_mode(layout_mode):
    if layout_mode == "none":
        return None
    return layout_mode


def _validate(mode, layout_mode):
    if mode not in _VALID_MODES:
        raise ValueError(f"mode {mode!r} is not supported")
    if layout_mode not in _VALID_LAYOUT_MODES:
        raise ValueError(f"layout_mode {layout_mode!r} is not supported")
    if mode in _MERGED_PROTEOFORM_MODES and layout_mode not in (None, "auto"):
        raise ValueError(
            f"mode {mode!r} requires layout_mode='auto' or 'none' (no PD "
            "layout to reuse for synthesized merged-proteoform activities)"
        )


def _new_model_builder():
    cls = momapy.builder.get_or_make_builder_cls(
        momapy.celldesigner.CellDesignerModel
    )
    return cls()


def transform(
    map_,
    mode: typing.Literal[
        "normal",
        "no-complex",
        "keep-species",
        "keep-species-no-complex",
        "casq",
    ] = "normal",
    layout_mode: typing.Literal[
        "plain", "overlay", "auto", "none"
    ] | None = "auto",
):
    layout_mode = _normalize_layout_mode(layout_mode)
    _validate(mode, layout_mode)

    clingo_model, clingo_id_to_model_element = pd2af.solver.solve(map_, mode)
    resolution = pd2af.model.resolve(
        clingo_model, clingo_id_to_model_element, mode
    )
    layout = pd2af.layouts.STRATEGIES[layout_mode]()

    model_builder = _new_model_builder()
    layout_builder, mapping_builder = layout.start(map_)

    for compartment in resolution.compartments:
        model_builder.compartments.add(compartment)
        layout.on_compartment(
            map_, compartment, layout_builder, mapping_builder
        )

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

    dedup_and_remap_model(model_builder, mapping_builder)
    return layout.finish(map_, model_builder, layout_builder, mapping_builder)
