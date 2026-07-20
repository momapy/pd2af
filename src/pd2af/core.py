import dataclasses
import typing

import momapy.celldesigner
import momapy.core.model
import momapy.sbgn.pd
import momapy.utils

import pd2af.build
import pd2af.languages
import pd2af.solver


@dataclasses.dataclass
class TransformerResult:
    """Result of :func:`transform`: the built AF map plus provenance.

    Attributes:
        obj: The transformed map (``CellDesignerMap`` or ``SBGNAFMap``), or
            the transformed model (``CellDesignerModel`` or ``SBGNAFModel``)
            when the input was a bare model rather than a map.
        provenance: Maps each input PD/CD model element to the
            ``frozenset`` of output AF model elements derived from it.
            Several input elements may map to one output element because
            the merged modes content-dedup their results (many-to-one);
            use ``.inverse`` (``id(output_element) -> frozenset`` of input
            elements) for the result-to-source direction. Input elements
            are valid keys because momapy's reader enforces the model-dedup
            invariant, so no two content-equal-but-distinct elements coexist
            in a single input map.
    """

    obj: typing.Any = None
    provenance: momapy.utils.FrozenIdentityMultiDict | None = None


_TRANSFORMATION_MODES = frozenset(
    {"normal", "normal-no-complex", "keep-species", "keep-species-no-complex", "casq"}
)
_MERGED_PROTEOFORM_MODES = pd2af.languages.MERGED_PROTEOFORM_MODES
# Accepted `layout_mode` input values: the concrete modes plus the `"auto"`
# meta value (resolved from the input type) and the `None` sentinel.
_LAYOUT_MODES = frozenset({"plain", "overlay", "dot", "auto", None})
_INFLUENCE_PAIRINGS = frozenset({"cross", "nearest"})

# The concrete layout modes (excluding the `None` sentinel and the `"auto"`
# meta value), in display order.
LAYOUT_MODES = ("plain", "overlay", "dot")

# SBGN-AF output supports the curated-geometry `plain` mode and the graphviz
# `dot` mode; the `overlay` dimming is CellDesigner-only.
SBGN_AF_LAYOUT_MODES = frozenset({"plain", "dot", None})
_SBGN_AF_LAYOUT_MODES = SBGN_AF_LAYOUT_MODES


def get_compatible_layout_modes_for_transformation_mode(mode):
    """Return the concrete layout modes valid for a transformation mode.

    The merged-proteoform modes only accept `dot`, because their synthesized
    merged activities have no original geometry to reuse; every other mode
    accepts all three layout modes.
    """
    if mode in _MERGED_PROTEOFORM_MODES:
        return ("dot",)
    return LAYOUT_MODES


def _normalize_layout_mode(layout_mode):
    if layout_mode == "none":
        return None
    return layout_mode


def _validate_layout_mode(layout_mode, mode):
    if mode in _MERGED_PROTEOFORM_MODES and layout_mode not in (None, "dot"):
        raise ValueError(f"mode {mode} requires layout_mode 'dot' or None")


def _validate_layout_mode_for_language(layout_mode, language):
    if (
        language == pd2af.languages.SBGN_PD
        and layout_mode not in _SBGN_AF_LAYOUT_MODES
    ):
        raise ValueError(
            f"SBGN-AF output supports layout_mode 'plain', 'dot' or None, "
            f"got {layout_mode!r} ('overlay' is unsupported)"
        )


def _wrap_model_in_map(model, language):
    """Wrap a bare input model in a layout-less map of the matching language.

    ``Map`` is a keyword-only frozen dataclass whose ``layout`` and
    ``layout_model_mapping`` default to ``None``, so a model-only map is a
    direct construction.
    """
    if language == pd2af.languages.SBGN_PD:
        return momapy.sbgn.pd.SBGNPDMap(model=model)
    return momapy.celldesigner.CellDesignerMap(model=model)


def transform(
    map_or_model,
    mode: typing.Literal[
        "normal",
        "normal-no-complex",
        "keep-species",
        "keep-species-no-complex",
        "casq",
    ] = "normal",
    layout_mode: typing.Literal["auto", "dot", "plain", "overlay"] | None = "auto",
    influence_pairing: typing.Literal["cross", "nearest"] = "cross",
    set_active: list[str] | None = None,
    set_inactive: list[str] | None = None,
    set_all_active: bool = False,
    set_all_inactive: bool = False,
):
    layout_mode = _normalize_layout_mode(layout_mode)
    is_model_input = isinstance(map_or_model, momapy.core.model.Model)
    if is_model_input:
        if layout_mode not in (None, "auto"):
            raise ValueError(
                f"model input supports only layout_mode 'auto' or None, "
                f"got {layout_mode!r}"
            )
        layout_mode = None
        map_ = _wrap_model_in_map(
            map_or_model, pd2af.languages.language_from_model(map_or_model)
        )
    else:
        map_ = map_or_model
        if layout_mode == "auto":
            layout_mode = "dot"
    _validate_layout_mode(layout_mode, mode)
    _validate_layout_mode_for_language(
        layout_mode, pd2af.languages.language_from_map(map_)
    )
    if influence_pairing not in _INFLUENCE_PAIRINGS:
        raise ValueError(
            f"influence_pairing must be one of {sorted(_INFLUENCE_PAIRINGS)}, "
            f"got {influence_pairing!r}"
        )
    clingo_model, clingo_id_to_model_element = pd2af.solver.solve(
        map_,
        mode,
        set_active=set_active,
        set_inactive=set_inactive,
        set_all_active=set_all_active,
        set_all_inactive=set_all_inactive,
    )
    result = pd2af.build.build_map(
        map_,
        layout_mode,
        clingo_model,
        clingo_id_to_model_element,
        influence_pairing,
        mode=mode,
    )
    if is_model_input:
        return TransformerResult(
            obj=result.obj.model, provenance=result.provenance
        )
    return result
