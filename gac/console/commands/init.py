# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import os
import re
import sys
import urllib.parse
import shutil

from cleo import option

from .command import Command


class InitCommand(Command):
    name = "init"

    def __init__(self):
        super(InitCommand, self).__init__()

    def handle(self):
        from pathlib import Path
        import gac

        src_parent_path = Path(gac.__file__).parent
        csrc_path = str(src_parent_path) + "/cxx_src"

        dest_path = str(Path.cwd())

        for file in os.listdir(csrc_path):
            src = "/".join([csrc_path, file])
            dest = "/".join([dest_path, file])
            if os.path.isdir(src):
                shutil.copytree(src, dest)
            elif os.path.isfile(src):
                shutil.copyfile(src, dest)
            # else:
            #     raise TypeError
