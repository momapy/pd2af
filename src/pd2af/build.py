"""Coordinator: drive the two-phase BuilderContext pipeline.

Pass 1 (:func:`pd2af._building_model.make_and_add_model`) walks clingo
activity / influence atoms and populates ``context.model`` with canonical,
content-deduped compartments, templates, species and modulations.
References between elements are wired to canonical instances at
construction time.

Pass 2 (:func:`pd2af._building_layout.make_and_add_layout`) -- skipped
entirely when ``layout_mode is None`` -- populates ``context.layout`` and
``context.layout_model_mapping``, branching on ``layout_mode``.

``build_map`` creates the ``BuilderContext``, runs the two passes, then
assembles the final CellDesignerMap from the three context slots and
returns the frozen map.
"""

import dataclasses

import momapy.builder
import momapy.celldesigner

import pd2af._building_layout
import pd2af._building_model
import pd2af.utils


@dataclasses.dataclass
class BuilderContext:
    # --- inputs ---
    input_map: object
    layout_mode: str | None
    clingo_id_to_model_element: dict
    influence_pairing: str = "cross"

    # --- outputs being built ---
    model: object = None
    layout: object = None
    layout_model_mapping: object = None

    # --- Pass-1 -> Pass-2 handoff ---
    species_emissions: list = dataclasses.field(default_factory=list)
    input_model_element_to_canonical_model_element: dict = dataclasses.field(
        default_factory=dict
    )

    # --- Pass-1 scratch ---
    cache: dict = dataclasses.field(default_factory=dict)
    subunit_to_top_level: dict = None
    activity_atoms_by_key_class: dict = dataclasses.field(default_factory=dict)
    influence_atoms: list = dataclasses.field(default_factory=list)
    key_to_species: dict = dataclasses.field(default_factory=dict)

    # --- Pass-2 scratch ---
    model_element_to_layout_elements: dict = dataclasses.field(default_factory=dict)
    object_to_builder: dict = dataclasses.field(default_factory=dict)
    synthetic_index: int = 0


def build_map(
    map_,
    layout_mode,
    clingo_model,
    clingo_id_to_model_element,
    influence_pairing="cross",
):
    context = BuilderContext(
        input_map=map_,
        layout_mode=layout_mode,
        clingo_id_to_model_element=clingo_id_to_model_element,
        influence_pairing=influence_pairing,
    )
    pd2af._building_model.make_and_add_model(context, clingo_model)
    if layout_mode is not None:
        pd2af._building_layout.make_and_add_layout(context)

    map_builder_class = momapy.builder.get_or_make_builder_cls(
        momapy.celldesigner.CellDesignerMap
    )
    map_builder = map_builder_class(
        model=context.model,
        layout=context.layout,
        layout_model_mapping=context.layout_model_mapping,
    )
    new_map = momapy.builder.object_from_builder(map_builder)

    if layout_mode == "auto":
        new_map = pd2af.utils.auto_layout(new_map)
    return new_map
