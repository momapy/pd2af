import typing

import pd2af.build
import pd2af.languages
import pd2af.solver


_TRANSFORMATION_MODES = frozenset(
    {"normal", "no-complex", "keep-species", "keep-species-no-complex", "casq"}
)
_MERGED_PROTEOFORM_MODES = frozenset({"normal", "no-complex"})
_LAYOUT_MODES = frozenset({"plain", "overlay", "auto", None})
_INFLUENCE_PAIRINGS = frozenset({"cross", "nearest"})

# SBGN-AF output currently only reuses the curated input geometry (plain);
# graphviz `auto` and the `overlay` dimming are CellDesigner-only for now.
_SBGN_AF_LAYOUT_MODES = frozenset({"plain", None})


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
            f"SBGN-AF output currently supports layout_mode 'plain' or None, "
            f"got {layout_mode!r} ('auto'/'overlay' are future work)"
        )


def transform(
    map_,
    mode: typing.Literal[
        "normal",
        "no-complex",
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
    )
