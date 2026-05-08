import typing

import pd2af.build
import pd2af.solver


_TRANSFORMATION_MODES = frozenset(
    {"normal", "no-complex", "keep-species", "keep-species-no-complex", "casq"}
)
_MERGED_PROTEOFORM_MODES = frozenset({"normal", "no-complex"})
_LAYOUT_MODES = frozenset({"plain", "overlay", "auto", None})


def _normalize_layout_mode(layout_mode):
    if layout_mode == "none":
        return None
    return layout_mode


def _validate_layout_mode(layout_mode, mode):
    if mode in _MERGED_PROTEOFORM_MODES and layout_mode not in (None, "auto"):
        raise ValueError(f"mode {mode} requires layout_mode 'auto' or None")


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
):
    layout_mode = _normalize_layout_mode(layout_mode)
    _validate_layout_mode(layout_mode, mode)
    clingo_model, clingo_id_to_model_element = pd2af.solver.solve(map_, mode)
    return pd2af.build.build_map(
        map_, layout_mode, clingo_model, clingo_id_to_model_element
    )
