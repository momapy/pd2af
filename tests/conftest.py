import pytest

from tests._helpers import (
    EXAMPLE_MAP_PATH,
    SBGN_EXAMPLE_MAP_PATH,
    read_cd_map,
    read_sbgn_map,
)


@pytest.fixture(scope="session")
def example_map_path():
    return EXAMPLE_MAP_PATH


@pytest.fixture(scope="session")
def example_cd_map():
    return read_cd_map(EXAMPLE_MAP_PATH)


@pytest.fixture(scope="session")
def sbgn_example_map():
    return read_sbgn_map(SBGN_EXAMPLE_MAP_PATH)
