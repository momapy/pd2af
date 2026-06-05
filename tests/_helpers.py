import os
import shutil

import momapy.io.core


_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
MAPS_DIR = os.path.join(_TESTS_DIR, "maps", "celldesigner")
EXAMPLE_MAP_PATH = os.path.join(MAPS_DIR, "example.xml")


def read_cd_map(path):
    return momapy.io.core.read(path).obj


def has_dot_binary():
    return shutil.which("dot") is not None


def modulation_set(model):
    return {
        (type(m).__name__, m.source.name, m.target.name)
        for m in model.modulations
    }


def species_names(model):
    return sorted(s.name for s in model.species)
