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

        target_path = str(Path.cwd())

        for file in os.listdir(csrc_path):
            shutil.copyfile("/".join([csrc_path, file]), "/".join([target_path, file]))
