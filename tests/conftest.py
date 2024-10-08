import pytest
from pathlib import Path
from cleo.application import Application
from naunet.console.commands import (
    ExampleCommand,
    ExtendCommand,
    InitCommand,
    NewCommand,
    RenderCommand,
)


def pytest_addoption(parser):
    parser.addoption(
        "--runslow", action="store_true", default=False, help="run slow tests"
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: mark test as slow to run")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runslow"):
        # --runslow given in cli: do not skip slow tests
        return
    skip_slow = pytest.mark.skip(reason="need --runslow option to run")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)


@pytest.fixture
def datadir():
    return Path(__file__).absolute().parent / "data"


@pytest.fixture
def application():
    app = Application()
    commands = [
        ExampleCommand(),
        ExtendCommand(),
        InitCommand(),
        NewCommand(),
        RenderCommand(),
    ]

    for command in commands:
        app.add(command)

    return app
