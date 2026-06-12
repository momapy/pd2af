import typing

import pd2af.build
import pd2af.languages
import pd2af.solver


_TRANSFORMATION_MODES = frozenset(
    {"normal", "normal-no-complex", "keep-species", "keep-species-no-complex", "casq"}
)
_MERGED_PROTEOFORM_MODES = pd2af.languages.MERGED_PROTEOFORM_MODES
_LAYOUT_MODES = frozenset({"plain", "overlay", "auto", None})
_INFLUENCE_PAIRINGS = frozenset({"cross", "nearest"})

# SBGN-AF output supports the curated-geometry `plain` mode and the graphviz
# `auto` mode; the `overlay` dimming is CellDesigner-only.
_SBGN_AF_LAYOUT_MODES = frozenset({"plain", "auto", None})


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
    clingo_model, clingo_id_to_model_element = pd2af.solver.solve(map_, mode)
    return pd2af.build.build_map(
        map_,
        layout_mode,
        clingo_model,
        clingo_id_to_model_element,
        influence_pairing,
        mode=mode,
    )
