import typing

import pd2af.layouts
import pd2af.solver


_MERGED_PROTEOFORM_MODES = frozenset({"normal", "no-complex"})


def transform(
    cd_map,
    mode: typing.Literal[
        "normal",
        "no-complex",
        "keep-species",
        "keep-species-no-complex",
        "casq",
    ] = "normal",
    layout_mode: typing.Literal["plain", "overlay", "auto"] | None = "auto",
):
    if mode in _MERGED_PROTEOFORM_MODES and layout_mode != "auto":
        raise ValueError(
            f"mode {mode!r} requires layout_mode='auto' (no PD layout to "
            "reuse for synthesized merged-proteoform activities)"
        )
    clingo_model, id_to_model_element = pd2af.solver.solve(cd_map, mode=mode)
    new_cd_model = pd2af.solver.make_new_cd_model(
        clingo_model, id_to_model_element
    )
    return pd2af.layouts.build_map(cd_map, new_cd_model, layout_mode)
