import pytest
import tomlkit
from naunet.network import Network
from naunet.configuration import BaseConfiguration, NetworkConfiguration


@pytest.fixture
def example_base_config():
    return BaseConfiguration("test_project")


@pytest.fixture
def example_network_config():
    return NetworkConfiguration("test_project", Network())


def test_content(example_base_config):
    content = tomlkit.loads(example_base_config.content)

    general = content["general"]
    chemistry = content["chemistry"]
    chem_species = chemistry["species"]
    assert general["name"] == "test_project"
    assert chem_species["elements"] == []
    assert chem_species["allowed_species"] == []
