import sys
from cleo import Application as BaseApplication
from importlib.metadata import version

from .commands.init import InitCommand
from .commands.render import RenderCommand
from .commands.example import ExampleCommand


class Application(BaseApplication):
    def __init__(self):
        super(Application, self).__init__("naunet", version("naunet"))

        for command in self.get_default_commands():
            self.add(command)

    def get_default_commands(self):  # type: () -> list
        commands = [
            InitCommand(),
            RenderCommand(),
            ExampleCommand(),
        ]

        return commands


if __name__ == "__main__":
    Application().run()