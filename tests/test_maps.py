"""Smoke tests over the bundled CellDesigner maps.

For every XML map under tests/maps/{pd_dm,covid_dm}, the `normal` and
`no-complex` transformation modes are exercised in the no-layout mode.
Cases are parametrized so each map and mode is reported as its own
pytest node id, e.g.::

    tests/test_maps.py::test_smoke[pd_dm/Glycolysis.xml-normal]
"""

import os

import pytest

import momapy.celldesigner

import pd2af

from tests._helpers import (
    COVID_DM_DIR,
    MAPS_DIR,
    PD_DM_DIR,
    list_xml_files,
    read_cd_map,
)


def _collect_cases():
    cases = []
    for directory in (PD_DM_DIR, COVID_DM_DIR):
        for path in list_xml_files(directory):
            for mode in ("normal", "no-complex"):
                rel = os.path.relpath(path, MAPS_DIR)
                cases.append(pytest.param(path, mode, id=f"{rel}-{mode}"))
    return cases


_CASES = _collect_cases()


@pytest.mark.skipif(not _CASES, reason="no XML maps under tests/maps/")
@pytest.mark.parametrize("path,mode", _CASES)
def test_smoke(path, mode):
    cd_map = read_cd_map(path)
    out = pd2af.transform(cd_map, mode=mode, layout_mode=None)
    assert isinstance(out, momapy.celldesigner.CellDesignerMap)
    assert isinstance(out.model, momapy.celldesigner.CellDesignerModel)
    for mod in out.model.modulations:
        assert isinstance(mod.source, momapy.celldesigner.Species)
        assert isinstance(mod.target, momapy.celldesigner.Species)
        assert isinstance(
            mod,
            (
                momapy.celldesigner.PositiveInfluence,
                momapy.celldesigner.Inhibition,
            ),
        )
