import tomlkit
from naunet.configuration import Configuration


def test_init_configuration():
    config = Configuration("newproject")


def test_content():
    config = Configuration("newproject")
    content = config.content

    loaded = tomlkit.loads(content)
    general = loaded["general"]
    chemistry = loaded["chemistry"]
    assert general["name"] == "newproject"
    assert chemistry["elements"] == []
    assert chemistry["species"] == []