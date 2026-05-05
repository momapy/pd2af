import momapy.builder
import momapy.core.mapping

import pd2af.utils
from pd2af.placers import PlaceContext, Placer, register
from pd2af.placers._shared import (
    add_modulation_mapping,
    get_input_layouts_for,
    make_modulation_arc,
)
from pd2af.walkers import BuildStep, BuildStepKind


@register("overlay")
class OverlayPlacer(Placer):
    """Clones the full input layout up front, then dims layout elements
    no build step consumed. This preserves the original "highlight what
    survives" rendering behaviour."""

    def __init__(self):
        self._object_to_builder = {}
        self._kept_builders = set()
        # Maps id(new_element) → cloned layout builder, used by modulation
        # build steps to find arc endpoints.
        self._species_layout_for_new = {}

    def make_layout_builders(self, cd_map):
        layout_builder = momapy.builder.builder_from_object(
            cd_map.layout, object_to_builder=self._object_to_builder
        )
        mapping_builder = (
            momapy.core.mapping.LayoutModelMappingBuilder.from_object(
                cd_map.layout_model_mapping,
                object_to_builder=self._object_to_builder,
            )
        )
        return layout_builder, mapping_builder

    def _cloned_layouts_for(self, cd_map, input_element):
        layouts = get_input_layouts_for(cd_map, input_element)
        if not layouts:
            return ()
        cloned = []
        for layout in layouts:
            builder = self._object_to_builder.get(id(layout))
            if builder is not None:
                cloned.append(builder)
        return tuple(cloned)

    def place(self, build_step: BuildStep, context: PlaceContext) -> None:
        if build_step.kind is BuildStepKind.COMPARTMENT:
            self._place_compartment(build_step, context)
        elif build_step.kind is BuildStepKind.TEMPLATE:
            return None
        elif build_step.kind is BuildStepKind.SPECIES:
            self._place_species(build_step, context)
        elif build_step.kind is BuildStepKind.MODULATION:
            self._place_modulation(build_step, context)

    def _place_compartment(self, build_step, context):
        for input_compartment in build_step.provenances:
            for builder in self._cloned_layouts_for(context.cd_map, input_compartment):
                self._kept_builders.add(builder)

    def _place_species(self, build_step, context):
        cloned = ()
        for input_species in build_step.provenances:
            cloned = self._cloned_layouts_for(context.cd_map, input_species)
            if cloned:
                break
        if not cloned:
            return
        species_layout = cloned[0]
        self._species_layout_for_new[id(build_step.new_element)] = species_layout
        for builder in cloned:
            self._kept_builders.add(builder)

    def _place_modulation(self, build_step, context):
        modulation = build_step.new_element
        source_layout = self._species_layout_for_new.get(id(modulation.source))
        target_layout = self._species_layout_for_new.get(id(modulation.target))
        if source_layout is None or target_layout is None:
            return
        # Build the arc against cloned source/target builders. The arc's
        # internal Segment/Point entries get throwaway addresses; isolate
        # them in a per-arc cache so they don't poison object_to_builder.
        per_arc_cache = dict(self._object_to_builder)
        arc = make_modulation_arc(modulation, source_layout, target_layout)
        arc_builder = momapy.builder.builder_from_object(
            arc, object_to_builder=per_arc_cache
        )
        source_builder = momapy.builder.builder_from_object(
            source_layout, object_to_builder=self._object_to_builder
        )
        target_builder = momapy.builder.builder_from_object(
            target_layout, object_to_builder=self._object_to_builder
        )
        context.new_layout_builder.layout_elements.append(arc_builder)
        add_modulation_mapping(
            context.new_mapping_builder,
            arc_builder,
            source_builder,
            target_builder,
            modulation,
        )
        self._kept_builders.update(
            [arc_builder, source_builder, target_builder]
        )

    def finalize(self, context: PlaceContext, new_map_builder) -> None:
        if context.new_layout_builder is None:
            return
        replaced = pd2af.utils.highlight_layout_elements(
            self._kept_builders, context.new_layout_builder
        )
        context.new_layout_builder = replaced
        pd2af.utils.harmonize_root_layout(context.new_layout_builder)
