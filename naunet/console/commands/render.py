# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import errno
import os
import re
import sys
import urllib.parse
import shutil

from cleo import argument, option
from tomlkit import dumps
from tomlkit.toml_file import TOMLFile

from .command import Command


class RenderCommand(Command):
    """
    Render source codes according to the project setting

    render
        {--f|force : forced to override the existing files}
        {--update-species=? : allow to update the species list in toml configure file}
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
        species = chemistry["species"]

        odesolver = content["ODEsolver"]
        required = odesolver["required"]

        from pathlib import Path
        import naunet

        if len(species) == 0:
            species = None

        naunet.settings.initialize(element=element, pseudo_element=pseudo_element)

        from naunet.network import Network

        header_prefix = os.path.join(Path.cwd(), "include")
        source_prefix = os.path.join(Path.cwd(), "src")

        if os.path.exists(header_prefix):
            if not os.path.isdir(header_prefix):
                raise FileNotFoundError(
                    errno.ENOENT, os.strerror(errno.ENOENT), header_prefix
                )

            elif os.listdir(header_prefix):
                overwrite = self.option("force") or self.confirm(
                    "Non-empty include directory. Overwrite?", False
                )

                if not overwrite:
                    sys.exit()

        else:
            os.mkdir(header_prefix)

        if os.path.exists(source_prefix):
            if not os.path.isdir(source_prefix):
                raise FileNotFoundError(
                    errno.ENOENT, os.strerror(errno.ENOENT), source_prefix
                )

            elif os.listdir(source_prefix):
                overwrite = self.option("force") or self.confirm(
                    "Non-empty source directory. Overwrite?", False
                )

                if not overwrite:
                    sys.exit()

        else:
            os.mkdir(source_prefix)

        net = Network(network, database, species=species)

        net.check_duplicate_reaction()
        net.info.to_ccode(
            to_file=True,
            prefix=header_prefix,
            file_name="naunet_constants.h",
        )
        net.userdata.to_ccode(
            to_file=True,
            prefix=header_prefix,
            file_name="naunet_userdata.h",
        )
        for func in required:
            net.ode_expression.to_ccode(
                function=func,
                to_file=True,
                prefix=source_prefix,
                file_name=f"{func}.cpp",
                header=True,
                header_prefix=header_prefix,
                header_file=f"{func}.h",
            )

        src_parent_path = Path(naunet.__file__).parent
        csrc_path = os.path.join(src_parent_path, "cxx_src", "cvode_example", "src")
        dest_path = os.path.join(Path.cwd(), "src")

        incfile = os.path.join(
            src_parent_path, "cxx_src", "cvode_example", "include", "naunet.h"
        )
        dest = os.path.join(Path.cwd(), "include", "naunet.h")
        shutil.copyfile(incfile, dest)

        for src in ["main.cpp", "naunet.cpp"]:
            srcfile = os.path.join(csrc_path, src)
            dest = os.path.join(dest_path, src)
            shutil.copyfile(srcfile, dest)

        cmakefile = os.path.join(
            src_parent_path, "cxx_src", "cvode_example", "CMakeLists.txt"
        )
        dest = os.path.join(Path.cwd(), "CMakeLists.txt")
        shutil.copyfile(cmakefile, dest)

        print(self.option("update-species"))
        update = self.option("update-species")
        if not update:
            update = self.confirm("Update species in configure file?", False)

        if update or self.option("update-species").lower != "false":
            chemistry["species"] = [x.name for x in net.info.net_species]

            content["chemistry"] = chemistry

            config_file = os.path.join(Path.cwd(), "naunet_config.toml")
            with open(config_file, "w", encoding="utf-8") as f:
                f.write(dumps(content))

        # csrc_path = str(src_parent_path) + "/cxx_src"
        # dest_path = str(Path.cwd())
        # for file in os.listdir(csrc_path):
        #     src = "/".join([csrc_path, file])
        #     dest = "/".join([dest_path, file])
        #     if os.path.isdir(src):
        #         shutil.copytree(src, dest)
        #     elif os.path.isfile(src):
        #         shutil.copyfile(src, dest)
        #     # else:
        #     #     raise TypeError

        # progress = self.progress_bar()
        # progress.finish()
