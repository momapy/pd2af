"""Layout strategies for the four ``layout_mode`` values.

Duck-typed event handlers — no ABC. Each strategy exposes the same
methods invoked by :func:`pd2af.transform.transform`:

* ``start(map_)`` returns ``(layout_builder, mapping_builder)``.
* ``on_compartment``, ``on_species``, ``on_modulation`` place each
  resolved element.
* ``on_species_done`` is the explicit flush between species and
  modulations.
* ``finish`` produces the final :class:`CellDesignerMap`.

The ``STRATEGIES`` dict at the bottom maps ``layout_mode`` → class.
"""

import momapy.builder
import momapy.celldesigner
import momapy.core.layout
import momapy.core.mapping
import momapy.geometry

import pd2af.predicates
import pd2af.utils


def _get_input_layouts_for(map_, input_element):
    layouts = map_.layout_model_mapping.get_mapping(input_element)
    if layouts is not None:
        return layouts
    for key in map_.layout_model_mapping.inverse.keys():
        if isinstance(key, tuple) and key[0] == input_element:
            return map_.layout_model_mapping.get_mapping(key)
    return None


def _make_modulation_arc(modulation, source_layout, target_layout):
    arc_class = pd2af.predicates.model_element_class_to_layout_element_class[
        type(modulation)
    ]
    if source_layout is target_layout:
        start_point = source_layout.anchor_point("north_north_west")
        end_point = source_layout.anchor_point("north_north_east")
    else:
        start_point = source_layout.border(target_layout.center())
        end_point = target_layout.border(source_layout.center())
        if start_point is None:
            start_point = source_layout.north_west()
        if end_point is None:
            end_point = target_layout.north_east()
    segment = momapy.geometry.Segment(start_point, end_point)
    return arc_class(
        source=source_layout,
        target=target_layout,
        segments=(segment,),
    )


def _add_modulation_mapping(
    mapping_builder, arc, source_layout, target_layout, modulation
):
    source_key = mapping_builder._singleton_to_key.get(source_layout)
    source_cluster = (
        source_key if source_key is not None else frozenset([source_layout])
    )
    target_key = mapping_builder._singleton_to_key.get(target_layout)
    target_cluster = (
        target_key if target_key is not None else frozenset([target_layout])
    )
    mapping_builder.add_mapping(
        frozenset([arc]) | source_cluster | target_cluster,
        modulation,
        anchor=arc,
    )


def _copy_subtree_singleton_mappings(
    source_mapping, layout_element, target_builder
):
    elements = [layout_element] + list(layout_element.descendants())
    for element in elements:
        if element in source_mapping:
            target_builder.add_mapping(element, source_mapping[element])


def _new_layout_builders():
    layout_builder_class = momapy.builder.get_or_make_builder_cls(
        momapy.core.layout.Layout
    )
    return (
        layout_builder_class(),
        momapy.core.mapping.LayoutModelMappingBuilder(),
    )


_SPECIES_CLASS_TO_LAYOUT_CLASS = {
    momapy.celldesigner.GenericProtein: momapy.celldesigner.GenericProteinLayout,
    momapy.celldesigner.TruncatedProtein: momapy.celldesigner.TruncatedProteinLayout,
    momapy.celldesigner.Receptor: momapy.celldesigner.ReceptorLayout,
    momapy.celldesigner.IonChannel: momapy.celldesigner.IonChannelLayout,
    momapy.celldesigner.Gene: momapy.celldesigner.GeneLayout,
    momapy.celldesigner.RNA: momapy.celldesigner.RNALayout,
    momapy.celldesigner.AntisenseRNA: momapy.celldesigner.AntisenseRNALayout,
    momapy.celldesigner.Phenotype: momapy.celldesigner.PhenotypeLayout,
    momapy.celldesigner.Ion: momapy.celldesigner.IonLayout,
    momapy.celldesigner.SimpleMolecule: momapy.celldesigner.SimpleMoleculeLayout,
    momapy.celldesigner.Drug: momapy.celldesigner.DrugLayout,
    momapy.celldesigner.Unknown: momapy.celldesigner.UnknownLayout,
    momapy.celldesigner.Complex: momapy.celldesigner.ComplexLayout,
}


def _make_synthetic_species_layout(species, index):
    layout_class = _SPECIES_CLASS_TO_LAYOUT_CLASS.get(type(species))
    if layout_class is None:
        raise ValueError(
            f"no default layout class registered for species type "
            f"{type(species).__name__}"
        )
    position = momapy.geometry.Point(float(index), 0.0)
    label = momapy.core.layout.TextLayout(
        text=getattr(species, "name", "") or "",
        position=position,
    )
    return layout_class(position=position, label=label)


class NoLayout:
    """No layout output — only the model is built."""

    def start(self, map_):
        return None, None

    def on_compartment(self, *_args, **_kwargs):
        return None

    def on_species(self, *_args, **_kwargs):
        return None

    def on_species_done(self, *_args, **_kwargs):
        return None

    def on_modulation(self, *_args, **_kwargs):
        return None

    def finish(self, map_, model_builder, layout_builder, mapping_builder):
        new_model = momapy.builder.object_from_builder(model_builder)
        return momapy.celldesigner.CellDesignerMap(model=new_model)


class PlainLayout:
    """Reuses input layouts wholesale for compartments and top-level
    species; synthesises fresh modulation arcs.

    A model species can have several alias layouts in the input map
    (e.g. one Complex drawn at two locations). All of them are emitted
    into the output: the writer needs every alias declared so that
    modulation arcs referencing any of them resolve on read-back.

    Species are buffered through :meth:`on_species` so that
    :meth:`on_species_done` can partition top-level vs. nested layouts
    in a single pass: an alias whose layout is a descendant of another
    candidate's layout is nested and brought along by its parent rather
    than appended at the top level.
    """

    def __init__(self):
        self._pending_top_level = []
        self._pending_subunits = []
        self._species_layout_for_new = {}

    def start(self, map_):
        return _new_layout_builders()

    def on_compartment(
        self, map_, compartment, layout_builder, mapping_builder
    ):
        layouts = _get_input_layouts_for(map_, compartment)
        if not layouts:
            return
        compartment_layout = layouts[0]
        layout_builder.layout_elements.append(compartment_layout)
        _copy_subtree_singleton_mappings(
            map_.layout_model_mapping, compartment_layout, mapping_builder
        )

    def on_species(
        self,
        map_,
        species,
        existing_species,
        is_subunit,
        layout_builder,
        mapping_builder,
    ):
        layouts = (
            tuple(_get_input_layouts_for(map_, existing_species) or ())
            if existing_species is not None
            else ()
        )
        if is_subunit:
            self._pending_subunits.append((species, layouts))
        else:
            self._pending_top_level.append((species, layouts))

    def on_species_done(self, map_, layout_builder, mapping_builder):
        # Build the descendant id set of every candidate top-level
        # alias layout (collected up front so order is irrelevant). An
        # alias whose layout is in someone else's descendant set is
        # "nested" — its parent (also a candidate) brings it along.
        descendant_ids = set()
        for _, layouts in self._pending_top_level:
            for candidate_layout in layouts:
                for descendant in candidate_layout.descendants():
                    descendant_ids.add(id(descendant))
        for species, layouts in self._pending_top_level:
            if not layouts:
                self._on_species_without_layout(
                    species, layout_builder, mapping_builder
                )
                continue
            chosen_layout = layouts[0]
            for alias_layout in layouts:
                if id(alias_layout) in descendant_ids:
                    continue
                layout_builder.layout_elements.append(alias_layout)
                _copy_subtree_singleton_mappings(
                    map_.layout_model_mapping, alias_layout, mapping_builder
                )
            self._species_layout_for_new[id(species)] = chosen_layout
        for species, layouts in self._pending_subunits:
            if layouts:
                self._species_layout_for_new[id(species)] = layouts[0]

    def on_modulation(
        self, modulation, layout_builder, mapping_builder
    ):
        source_layout = self._species_layout_for_new.get(id(modulation.source))
        target_layout = self._species_layout_for_new.get(id(modulation.target))
        if source_layout is None or target_layout is None:
            return
        arc = _make_modulation_arc(modulation, source_layout, target_layout)
        layout_builder.layout_elements.append(arc)
        _add_modulation_mapping(
            mapping_builder, arc, source_layout, target_layout, modulation
        )

    def finish(self, map_, model_builder, layout_builder, mapping_builder):
        pd2af.utils.harmonize_root_layout(layout_builder)
        return _finalize_map(model_builder, layout_builder, mapping_builder)

    def _on_species_without_layout(
        self, species, layout_builder, mapping_builder
    ):
        # Plain mode does not synthesise layouts. Hook for AutoLayout.
        return None


class AutoLayout(PlainLayout):
    """Like :class:`PlainLayout` but synthesises fallback layouts for
    species without input layouts and runs the auto-layout solver
    after finalisation."""

    def __init__(self):
        super().__init__()
        self._synthetic_index = 0

    def _on_species_without_layout(
        self, species, layout_builder, mapping_builder
    ):
        synthetic = _make_synthetic_species_layout(
            species, self._synthetic_index
        )
        self._synthetic_index += 1
        layout_builder.layout_elements.append(synthetic)
        mapping_builder.add_mapping(synthetic, species)
        self._species_layout_for_new[id(species)] = synthetic

    def finish(self, map_, model_builder, layout_builder, mapping_builder):
        pd2af.utils.harmonize_root_layout(layout_builder)
        new_map = _finalize_map(model_builder, layout_builder, mapping_builder)
        return pd2af.utils.auto_layout(new_map)


class OverlayLayout:
    """Clones the full input layout up front, then dims layout elements
    no event consumed. Preserves the original "highlight what survives"
    rendering."""

    def __init__(self):
        self._object_to_builder = {}
        self._kept_builders = set()
        # Maps id(species) → ORIGINAL (frozen) input layout. We keep
        # the frozen one so :func:`_make_modulation_arc` builds an arc
        # whose source/target are frozen objects;
        # :func:`builder_from_object` then reuses the cloned builders
        # via the per-arc cache. Storing the builder here would produce
        # a frozen arc with builder children and break the conversion
        # to a builder tree.
        self._species_layout_for_new = {}

    def start(self, map_):
        layout_builder = momapy.builder.builder_from_object(
            map_.layout, object_to_builder=self._object_to_builder
        )
        # We deliberately avoid `LayoutModelMappingBuilder.from_object`:
        # it round-trips model-side values through `builder_from_object`
        # / `object_from_builder`, producing fresh frozen clones with
        # different `id()`. In overlay-compatible modes resolved species
        # reuse input species objects by identity for the AF model, so
        # cloned mapping values become invisible to the writer's
        # identity lookup and `<listOfSpeciesAliases>` ends up empty.
        # Build the mapping by hand: convert layout-side keys to
        # builders, leave model-side values as the input frozen objects.
        mapping_builder = momapy.core.mapping.LayoutModelMappingBuilder()
        input_mapping = map_.layout_model_mapping
        anchor_for_key_id = {}
        for input_anchor, input_key in input_mapping._singleton_to_key.items():
            anchor_for_key_id[id(input_key)] = input_anchor
        for input_layout_key, input_model_value in input_mapping.items():
            new_layout_key = self._convert_layout_side(input_layout_key)
            input_anchor = anchor_for_key_id.get(id(input_layout_key))
            if input_anchor is not None:
                new_anchor = self._object_to_builder.get(
                    id(input_anchor), input_anchor
                )
                mapping_builder.add_mapping(
                    new_layout_key, input_model_value, anchor=new_anchor
                )
            else:
                mapping_builder[new_layout_key] = input_model_value
        return layout_builder, mapping_builder

    def on_compartment(
        self, map_, compartment, layout_builder, mapping_builder
    ):
        for builder in self._cloned_layouts_for(map_, compartment):
            self._kept_builders.add(builder)

    def on_species(
        self,
        map_,
        species,
        existing_species,
        is_subunit,
        layout_builder,
        mapping_builder,
    ):
        if existing_species is None:
            return
        original_layouts = (
            _get_input_layouts_for(map_, existing_species) or ()
        )
        cloned_layouts = self._cloned_layouts_for(map_, existing_species)
        if not original_layouts:
            return
        # Store the FROZEN original layout for modulation arc building;
        # cloned builders go into kept_builders for the dim pass.
        self._species_layout_for_new[id(species)] = original_layouts[0]
        for builder in cloned_layouts:
            self._kept_builders.add(builder)

    def on_species_done(self, map_, layout_builder, mapping_builder):
        return None

    def on_modulation(
        self, modulation, layout_builder, mapping_builder
    ):
        source_layout = self._species_layout_for_new.get(id(modulation.source))
        target_layout = self._species_layout_for_new.get(id(modulation.target))
        if source_layout is None or target_layout is None:
            return
        # Build the arc against cloned source/target builders. The arc's
        # internal Segment/Point entries get throwaway addresses; isolate
        # them in a per-arc cache so they don't poison object_to_builder.
        per_arc_cache = dict(self._object_to_builder)
        arc = _make_modulation_arc(modulation, source_layout, target_layout)
        arc_builder = momapy.builder.builder_from_object(
            arc, object_to_builder=per_arc_cache
        )
        source_builder = momapy.builder.builder_from_object(
            source_layout, object_to_builder=self._object_to_builder
        )
        target_builder = momapy.builder.builder_from_object(
            target_layout, object_to_builder=self._object_to_builder
        )
        layout_builder.layout_elements.append(arc_builder)
        _add_modulation_mapping(
            mapping_builder,
            arc_builder,
            source_builder,
            target_builder,
            modulation,
        )
        self._kept_builders.update(
            [arc_builder, source_builder, target_builder]
        )

    def finish(self, map_, model_builder, layout_builder, mapping_builder):
        layout_builder = pd2af.utils.highlight_layout_elements(
            self._kept_builders, layout_builder
        )
        pd2af.utils.harmonize_root_layout(layout_builder)
        return _finalize_map(model_builder, layout_builder, mapping_builder)

    def _convert_layout_side(self, layout_key):
        if isinstance(layout_key, frozenset):
            return frozenset(
                self._object_to_builder.get(id(element), element)
                for element in layout_key
            )
        return self._object_to_builder.get(id(layout_key), layout_key)

    def _cloned_layouts_for(self, map_, input_element):
        layouts = _get_input_layouts_for(map_, input_element)
        if not layouts:
            return ()
        cloned = []
        for layout in layouts:
            builder = self._object_to_builder.get(id(layout))
            if builder is not None:
                cloned.append(builder)
        return tuple(cloned)


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


STRATEGIES = {
    None: NoLayout,
    "plain": PlainLayout,
    "auto": AutoLayout,
    "overlay": OverlayLayout,
}
