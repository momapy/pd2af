"""Input-language tokens and inference from the input map type.

The input language drives both the ontology vocabulary
(:mod:`pd2af.ontology`) and which rule variant is resolved
(:mod:`pd2af.rules`). The output language is *deduced* from the input:
``celldesigner`` input -> CellDesigner output, ``sbgn_pd`` input ->
SBGN-AF output.
"""

import momapy.celldesigner
import momapy.sbgn.pd


CELLDESIGNER = "celldesigner"
SBGN_PD = "sbgn_pd"


# The "merged" (true-AF) transformation modes: proteoforms collapse and all
# post-translational decorations are stripped. The build stage strips PTMs iff
# ``context.mode`` is one of these; the complementary modes keep decorations.
# (Tokens are the hyphenated transform-mode names, not the underscored ASP
# profile names.)
MERGED_PROTEOFORM_MODES = frozenset({"normal", "normal-no-complex"})


def language_from_map(map_):
    """Infer the input language token from the input map's type."""
    if isinstance(map_, momapy.sbgn.pd.SBGNPDMap):
        return SBGN_PD
    if isinstance(map_, momapy.celldesigner.CellDesignerMap):
        return CELLDESIGNER
    raise ValueError(
        f"unsupported input map type {type(map_).__name__!r}; "
        f"expected SBGNPDMap or CellDesignerMap"
    )


def language_from_model(model):
    """Infer the input language token from a bare input model's type."""
    if isinstance(model, momapy.sbgn.pd.SBGNPDModel):
        return SBGN_PD
    if isinstance(model, momapy.celldesigner.CellDesignerModel):
        return CELLDESIGNER
    raise ValueError(
        f"unsupported input model type {type(model).__name__!r}; "
        f"expected SBGNPDModel or CellDesignerModel"
    )
