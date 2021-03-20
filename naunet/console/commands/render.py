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

        for subdir in ["include", "src", "test"]:
            prefix = os.path.join(Path.cwd(), subdir)

            if os.path.exists(prefix):
                if not os.path.isdir(prefix):
                    raise FileNotFoundError(
                        errno.ENOENT, os.strerror(errno.ENOENT), prefix
                    )

                elif os.listdir(prefix):
                    overwrite = self.option("force") or self.confirm(
                        f"Non-empty {subdir} directory. Overwrite?", False
                    )

                    if not overwrite:
                        sys.exit()

            else:
                os.mkdir(prefix)

        net = Network(network, database, species=species)

        header_prefix = os.path.join(Path.cwd(), "include")
        source_prefix = os.path.join(Path.cwd(), "src")

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

        net.ode_expression.to_ccode(
            to_file=True,
            prefix=source_prefix,
            file_name=f"naunet_ode.cpp",
            header=True,
            header_prefix=header_prefix,
            header_file=f"naunet_ode.h",
        )

        src_parent_path = Path(naunet.__file__).parent
        template_path = os.path.join(src_parent_path, "cxx_src", "cvode_example")
        csrc_path = os.path.join(template_path, "src")

        for src in ["naunet.cpp"]:
            srcfile = os.path.join(csrc_path, src)
            dest = os.path.join(Path.cwd(), "src", src)
            shutil.copyfile(srcfile, dest)

        inc_path = os.path.join(template_path, "include")
        incfile = os.path.join(inc_path, "naunet.h")
        dest = os.path.join(Path.cwd(), "include", "naunet.h")
        shutil.copyfile(incfile, dest)

        testfile = os.path.join(template_path, "test", "main.cpp")
        dest = os.path.join(Path.cwd(), "test", "main.cpp")
        shutil.copyfile(testfile, dest)

        parfile = os.path.join(template_path, "test", "timeres.dat")
        dest = os.path.join(Path.cwd(), "test", "timeres.dat")
        shutil.copyfile(parfile, dest)

        for cmakesrc in ["CMakeLists.txt", "src/CMakeLists.txt", "test/CMakeLists.txt"]:
            cmakefile = os.path.join(template_path, cmakesrc)
            dest = os.path.join(Path.cwd(), cmakesrc)
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
