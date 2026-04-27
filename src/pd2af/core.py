import typing

import pd2af.layouts
import pd2af.solver


def transform(
    cd_map,
    mode: typing.Literal["pd2af", "pd2af-no-complex"] = "pd2af",
    layout_mode: typing.Literal["plain", "overlay", "auto"] | None = "plain",
):
    clingo_model, id_to_model_element = pd2af.solver.solve(cd_map, mode=mode)
    new_cd_model = pd2af.solver.make_new_cd_model(
        clingo_model, id_to_model_element
    )
    return pd2af.layouts.build_map(cd_map, new_cd_model, layout_mode)
