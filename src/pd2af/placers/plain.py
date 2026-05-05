import pd2af.utils
from pd2af.placers import PlaceContext, Placer, register
from pd2af.placers._shared import (
    add_modulation_mapping,
    copy_subtree_singleton_mappings,
    get_input_layouts_for,
    make_modulation_arc,
)
from pd2af.walkers import BuildStep, BuildStepKind


@register("plain")
class PlainPlacer(Placer):
    """Reuses input layouts wholesale for compartments and top-level
    species; synthesises fresh modulation arcs.

    Species build steps are collected during the walk; layouts are
    appended in a single partitioned pass so that nested layouts (a
    subunit's layout that is already a descendant of its parent's
    layout) aren't duplicated at the top level.
    """

    def __init__(self):
        # Pending top-level species: (build_step, input_species_layout) pairs.
        self._pending_top_level = []
        # Pending subunit species (parent is a kept complex): same shape.
        self._pending_subunits = []
        # Maps id(new_element) → layout element (post-flush) used as
        # arc endpoint.
        self._species_layout_for_new = {}
        self._flushed = False

    def place(self, build_step: BuildStep, context: PlaceContext) -> None:
        if build_step.kind is BuildStepKind.COMPARTMENT:
            self._place_compartment(build_step, context)
        elif build_step.kind is BuildStepKind.TEMPLATE:
            return None
        elif build_step.kind is BuildStepKind.SPECIES:
            self._enqueue_species(build_step, context)
        elif build_step.kind is BuildStepKind.MODULATION:
            self._ensure_species_flushed(context)
            self._place_modulation(build_step, context)

    def _place_compartment(self, build_step, context):
        for input_compartment in build_step.provenances:
            layouts = get_input_layouts_for(context.cd_map, input_compartment)
            if not layouts:
                continue
            compartment_layout = layouts[0]
            context.new_layout_builder.layout_elements.append(compartment_layout)
            copy_subtree_singleton_mappings(
                context.cd_map.layout_model_mapping,
                compartment_layout,
                context.new_mapping_builder,
            )
            break

    def _enqueue_species(self, build_step, context):
        layout = self._first_input_layout(context, build_step)
        if build_step.parent is None:
            self._pending_top_level.append((build_step, layout))
        else:
            self._pending_subunits.append((build_step, layout))

    def _first_input_layout(self, context, build_step):
        for input_species in build_step.provenances:
            layouts = get_input_layouts_for(context.cd_map, input_species)
            if layouts:
                return layouts[0]
        return None

    def _ensure_species_flushed(self, context):
        if self._flushed:
            return
        self._flushed = True
        self._flush_species(context)

    def _flush_species(self, context):
        # Build the descendant id set of every candidate top-level species
        # layout (collected up front so order is irrelevant). A candidate
        # whose layout is in someone else's descendant set is "nested".
        candidate_layouts = [
            layout for _, layout in self._pending_top_level if layout is not None
        ]
        descendant_ids = set()
        for layout in candidate_layouts:
            for descendant in layout.descendants():
                descendant_ids.add(id(descendant))
        for build_step, layout in self._pending_top_level:
            if layout is None:
                self._on_species_without_layout(build_step, context)
                continue
            if id(layout) in descendant_ids:
                # Nested: its parent's layout (also a candidate) will bring
                # it along as a descendant. Don't append here, but record
                # the mapping for modulation lookup.
                self._species_layout_for_new[id(build_step.new_element)] = layout
                continue
            context.new_layout_builder.layout_elements.append(layout)
            copy_subtree_singleton_mappings(
                context.cd_map.layout_model_mapping,
                layout,
                context.new_mapping_builder,
            )
            self._species_layout_for_new[id(build_step.new_element)] = layout
        for build_step, layout in self._pending_subunits:
            if layout is not None:
                self._species_layout_for_new[id(build_step.new_element)] = layout

    def _on_species_without_layout(self, build_step, context):
        # Plain mode does not synthesise layouts. Subclass hook for auto.
        return None

    def _place_modulation(self, build_step, context):
        modulation = build_step.new_element
        source_layout = self._species_layout_for_new.get(id(modulation.source))
        target_layout = self._species_layout_for_new.get(id(modulation.target))
        if source_layout is None or target_layout is None:
            return
        arc = make_modulation_arc(modulation, source_layout, target_layout)
        context.new_layout_builder.layout_elements.append(arc)
        add_modulation_mapping(
            context.new_mapping_builder,
            arc,
            source_layout,
            target_layout,
            modulation,
        )

    def finalize(self, context: PlaceContext, new_map_builder) -> None:
        self._ensure_species_flushed(context)
        if context.new_layout_builder is not None:
            pd2af.utils.harmonize_root_layout(context.new_layout_builder)
