import typing

import momapy.builder
import momapy.celldesigner

from pd2af.placers import PlaceContext, placer_for
from pd2af.walkers import BuildStep, BuildStepKind, walker_for

# Side-effect imports to register walkers/placers in their factories.
import pd2af.walkers.keep_species  # noqa: F401
import pd2af.walkers.normal  # noqa: F401
import pd2af.placers.auto  # noqa: F401
import pd2af.placers.no_layout  # noqa: F401
import pd2af.placers.overlay  # noqa: F401
import pd2af.placers.plain  # noqa: F401


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


def _add_to_model(model_builder, build_step: BuildStep):
    element = build_step.new_element
    if build_step.kind is BuildStepKind.COMPARTMENT:
        model_builder.compartments.add(element)
    elif build_step.kind is BuildStepKind.TEMPLATE:
        model_builder.species_templates.add(element)
    elif build_step.kind is BuildStepKind.SPECIES:
        if build_step.parent is None:
            model_builder.species.add(element)
    elif build_step.kind is BuildStepKind.MODULATION:
        model_builder.modulations.add(element)


def _new_model_builder():
    cls = momapy.builder.get_or_make_builder_cls(
        momapy.celldesigner.CellDesignerModel
    )
    return cls()


def _finalize_map(model_builder, layout_builder, mapping_builder):
    new_model = momapy.builder.object_from_builder(model_builder)
    if layout_builder is None:
        return momapy.celldesigner.CellDesignerMap(model=new_model)
    builder_to_object = {}
    new_layout = momapy.builder.object_from_builder(
        layout_builder, builder_to_object=builder_to_object
    )
    new_mapping = momapy.builder.object_from_builder(
        mapping_builder, builder_to_object=builder_to_object
    )
    return momapy.celldesigner.CellDesignerMap(
        model=new_model,
        layout=new_layout,
        layout_model_mapping=new_mapping,
    )


def transform(
    cd_map,
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
    walker = walker_for(mode)
    placer = placer_for(layout_mode)
    if walker is None:
        raise ValueError(f"no walker registered for mode {mode!r}")
    if placer is None:
        raise ValueError(
            f"no placer registered for layout_mode {layout_mode!r}"
        )
    model_builder = _new_model_builder()
    layout_builder, mapping_builder = placer.make_layout_builders(cd_map)
    place_context = PlaceContext(
        cd_map=cd_map,
        layout_mode=layout_mode,
        new_layout_builder=layout_builder,
        new_mapping_builder=mapping_builder,
    )
    for build_step in walker.walk(cd_map):
        _add_to_model(model_builder, build_step)
        placer.place(build_step, place_context)
    placer.finalize(place_context, model_builder)
    new_map = _finalize_map(
        model_builder,
        place_context.new_layout_builder,
        place_context.new_mapping_builder,
    )
    return placer.post_finalize(place_context, new_map)
