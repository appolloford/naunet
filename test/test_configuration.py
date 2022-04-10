from naunet.configuration import Configuration


def test_init_configuration():
    config = Configuration("newproject")


def test_content():
    config = Configuration("newporject")
    content = config.content