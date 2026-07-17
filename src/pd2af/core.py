import dataclasses
import typing

import momapy.utils

import pd2af.build
import pd2af.languages
import pd2af.solver


@dataclasses.dataclass
class TransformerResult:
    """Result of :func:`transform`: the built AF map plus provenance.

    Attributes:
        obj: The transformed map (``CellDesignerMap`` or ``SBGNAFMap``).
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
_LAYOUT_MODES = frozenset({"plain", "overlay", "auto", None})
_INFLUENCE_PAIRINGS = frozenset({"cross", "nearest"})

# The concrete layout modes (excluding the `None` sentinel), in display order.
LAYOUT_MODES = ("plain", "overlay", "auto")

# SBGN-AF output supports the curated-geometry `plain` mode and the graphviz
# `auto` mode; the `overlay` dimming is CellDesigner-only.
SBGN_AF_LAYOUT_MODES = frozenset({"plain", "auto", None})
_SBGN_AF_LAYOUT_MODES = SBGN_AF_LAYOUT_MODES


def get_compatible_layout_modes_for_transformation_mode(mode):
    """Return the concrete layout modes valid for a transformation mode.

    The merged-proteoform modes only accept `auto`, because their synthesized
    merged activities have no original geometry to reuse; every other mode
    accepts all three layout modes.
    """
    if mode in _MERGED_PROTEOFORM_MODES:
        return ("auto",)
    return LAYOUT_MODES


def _normalize_layout_mode(layout_mode):
    if layout_mode == "none":
        return None
    return layout_mode


def _validate_layout_mode(layout_mode, mode):
    if mode in _MERGED_PROTEOFORM_MODES and layout_mode not in (None, "auto"):
        raise ValueError(f"mode {mode} requires layout_mode 'auto' or None")


def _validate_layout_mode_for_language(layout_mode, language):
    if (
        language == pd2af.languages.SBGN_PD
        and layout_mode not in _SBGN_AF_LAYOUT_MODES
    ):
        raise ValueError(
            f"SBGN-AF output supports layout_mode 'plain', 'auto' or None, "
            f"got {layout_mode!r} ('overlay' is unsupported)"
        )


def transform(
    map_,
    mode: typing.Literal[
        "normal",
        "normal-no-complex",
        "keep-species",
        "keep-species-no-complex",
        "casq",
    ] = "normal",
    layout_mode: typing.Literal["auto", "plain", "overlay"] | None = "auto",
    influence_pairing: typing.Literal["cross", "nearest"] = "cross",
    active_ids: list[str] | None = None,
    inactive_ids: list[str] | None = None,
):
    layout_mode = _normalize_layout_mode(layout_mode)
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
        map_, mode, active_ids=active_ids, inactive_ids=inactive_ids
    )
    return pd2af.build.build_map(
        map_,
        layout_mode,
        clingo_model,
        clingo_id_to_model_element,
        influence_pairing,
        mode=mode,
    )
