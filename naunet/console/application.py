import sys

from cleo import Application as BaseApplication

from ..__version__ import __version__

from .commands.init import InitCommand
from .commands.render import RenderCommand


class Application(BaseApplication):
    def __init__(self):
        super(Application, self).__init__("naunet", __version__)

        for command in self.get_default_commands():
            self.add(command)

    def get_default_commands(self):  # type: () -> list
        commands = [
            InitCommand(),
            RenderCommand(),
        ]

        return commands


if __name__ == "__main__":
    Application().run()