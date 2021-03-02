# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import os
import re
import sys
import urllib.parse

from cleo import argument, option
from tomlkit import dumps
from tomlkit.toml_file import TOMLFile

from .command import Command


class RenderCommand(Command):
    """
    Render source codes according to the project setting

    render
    """

    def __init__(self):
        super(RenderCommand, self).__init__()

    def handle(self):

        config = TOMLFile("naunet_config.toml")

        content = config.read()

        chemistry = content["chemistry"]
        network = chemistry["network"]
        database = chemistry["database"]
        element = chemistry["elements"]
        pseudo_element = chemistry["pseudo_elements"]

        from pathlib import Path
        import naunet

        naunet.settings.initialize(element=element, pseudo_element=pseudo_element)

        from naunet.network import Network

        net = Network(network, database)

        net.check_duplicate_reaction()
        net.info.to_ccode(
            to_file=True,
            prefix=os.path.join(Path.cwd()),
            file_name="naunet_constants.h",
        )
        net.ode_expression.to_ccode(
            function="fex",
            to_file=True,
            prefix=os.path.join(Path.cwd()),
            file_name=f"fex.cpp",
            header=True,
            header_file=f"fex.h",
        )

        # progress = self.progress_bar()
        # progress.finish()
