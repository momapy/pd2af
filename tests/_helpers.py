import os
import shutil

import momapy.io.core


_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
MAPS_DIR = os.path.join(_TESTS_DIR, "maps", "celldesigner")
EXAMPLE_MAP_PATH = os.path.join(MAPS_DIR, "example.xml")

SBGN_MAPS_DIR = os.path.join(_TESTS_DIR, "maps", "sbgn")
SBGN_EXAMPLE_MAP_PATH = os.path.join(
    SBGN_MAPS_DIR, "insulin-like_growth_factor_signaling.sbgn"
)
# A small SBGN-PD map that carries explicit `compartmentRef` attributes (the
# curated maps do not), so the compartment-handling paths can be exercised.
SBGN_WITH_COMPARTMENTS_MAP_PATH = os.path.join(
    SBGN_MAPS_DIR, "with_compartments.sbgn"
)


def read_cd_map(path):
    return momapy.io.core.read(path).obj


def read_sbgn_map(path):
    return momapy.io.core.read(path, reader="sbgnml").obj


def has_dot_binary():
    return shutil.which("dot") is not None


def modulation_set(model):
    return {
        (type(modulation).__name__, modulation.source.name, modulation.target.name)
        for modulation in model.modulations
    }


def species_names(model):
    return sorted(species.name for species in model.species)
