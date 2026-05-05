import momapy.geometry

import pd2af.utils
from pd2af.placers import PlaceContext, register
from pd2af.placers.plain import PlainPlacer
from pd2af.walkers import BuildStep, BuildStepKind


_SPECIES_CLASS_TO_LAYOUT_CLASS = {}


def _layout_class_for(species):
    import momapy.celldesigner

    global _SPECIES_CLASS_TO_LAYOUT_CLASS
    if not _SPECIES_CLASS_TO_LAYOUT_CLASS:
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
    return _SPECIES_CLASS_TO_LAYOUT_CLASS.get(type(species))


def _make_synthetic_species_layout(species, index):
    import momapy.core.layout

    layout_cls = _layout_class_for(species)
    if layout_cls is None:
        raise ValueError(
            f"no default layout class registered for species type "
            f"{type(species).__name__}"
        )
    position = momapy.geometry.Point(float(index), 0.0)
    label = momapy.core.layout.TextLayout(
        text=getattr(species, "name", "") or "",
        position=position,
    )
    return layout_cls(position=position, label=label)


@register("auto")
class AutoPlacer(PlainPlacer):
    """Like `PlainPlacer` but synthesises fallback layouts for species
    without input layouts and runs the auto-layout solver after
    finalisation."""

    def __init__(self):
        super().__init__()
        self._synthetic_index = 0

    def _on_species_without_layout(self, build_step, context):
        synthetic = _make_synthetic_species_layout(
            build_step.new_element, self._synthetic_index
        )
        self._synthetic_index += 1
        context.new_layout_builder.layout_elements.append(synthetic)
        context.new_mapping_builder.add_mapping(synthetic, build_step.new_element)
        self._species_layout_for_new[id(build_step.new_element)] = synthetic

    def post_finalize(self, context: PlaceContext, new_map):
        return pd2af.utils.auto_layout(new_map)
