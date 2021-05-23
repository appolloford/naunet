# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import errno
import os
import re
import sys
import urllib.parse
import shutil

from cleo import argument, option
from importlib.metadata import version
from tomlkit import dumps
from tomlkit.toml_file import TOMLFile

from .command import Command


class RenderCommand(Command):
    """
    Render source codes according to the project setting

    render
        {--f|force : forced to override the existing files}
        {--update-species=? : allow to update the species list in toml configure file}
        {--patch= : create patch for target code}
    """

    def __init__(self):
        super(RenderCommand, self).__init__()

    def handle(self):

        ver = version("naunet")

        config = TOMLFile("naunet_config.toml")

        content = config.read()

        chemistry = content["chemistry"]
        network = chemistry["network"]
        database = chemistry["database"]
        element = chemistry["elements"]
        pseudo_element = chemistry["pseudo_elements"]
        species = chemistry["species"]
        dust = chemistry["dust"]
        binding = chemistry["binding_energy"]
        yields = chemistry["photon_yield"]
        rate_modifier = chemistry["rate_modifier"]
        ode_modifier = chemistry["ode_modifier"]

        odesolver = content["ODEsolver"]
        solver = odesolver["solver"]
        method = odesolver["method"]
        device = odesolver["device"]
        # required = odesolver["required"]

        from pathlib import Path
        import naunet

        if len(species) == 0:
            species = None

        from naunet.network import Network, supported_reaction_class
        from naunet.species import Species

        Species.set_known_elements(element)
        Species.set_known_pseudoelements(pseudo_element)

        from naunet.chemistry import update_binding_energy, update_photon_yield

        update_binding_energy(binding)
        update_photon_yield(yields)

        if not supported_reaction_class.get(database):
            from importlib import util

            spec = util.spec_from_file_location(database, f"{database}.py")
            module = util.module_from_spec(spec)
            spec.loader.exec_module(module)

        net = Network(species=species, dusttype=dust["type"])
        net.add_reaction_from_file(network, database)
        net.rate_modifier = rate_modifier
        net.ode_modifier = ode_modifier

        # patches are rendered independently
        patch = self.option("patch")
        if patch:
            pm = net.patchmaker(patch)
            pm.render(os.path.join(Path.cwd(), patch))
            return

        # Check whether include and src folders exist, test folder is checked in example.py
        for subdir in ["include", "src"]:
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

        header_prefix = os.path.join(Path.cwd(), "include")
        source_prefix = os.path.join(Path.cwd(), "src")

        net.check_duplicate_reaction()
        tl = net.templateloader(solver=solver, method=method, device=device)
        tl.render_constants(prefix=source_prefix, headerprefix=header_prefix)
        tl.render_macros(prefix=header_prefix)
        tl.render_naunet(prefix=source_prefix, headerprefix=header_prefix)
        tl.render_ode(prefix=source_prefix, headerprefix=header_prefix)
        tl.render_physics(prefix=source_prefix, headerprefix=header_prefix)
        tl.render_userdata(prefix=header_prefix)
        tl.render_cmake(prefix=Path.cwd(), version=ver)

        src_parent_path = Path(naunet.__file__).parent
        template_path = os.path.join(src_parent_path, "templates", solver)
        csrc_path = os.path.join(template_path, "src")

        # for src in ["naunet.cpp"]:
        #     srcfile = os.path.join(csrc_path, src)
        #     dest = os.path.join(Path.cwd(), "src", src)
        #     shutil.copyfile(srcfile, dest)

        inc_path = os.path.join(template_path, "include")
        for inc in ["naunet_timer.h"]:
            incfile = os.path.join(inc_path, inc)
            dest = os.path.join(Path.cwd(), "include", inc)
            shutil.copyfile(incfile, dest)

        # testfile = os.path.join(template_path, "test", "main.cpp")
        # dest = os.path.join(Path.cwd(), "test", "main.cpp")
        # shutil.copyfile(testfile, dest)

        # parfile = os.path.join(template_path, "test", "timeres.dat")
        # dest = os.path.join(Path.cwd(), "test", "timeres.dat")
        # shutil.copyfile(parfile, dest)

        # for cmakesrc in ["CMakeLists.txt", "src/CMakeLists.txt", "test/CMakeLists.txt"]:
        #     cmakefile = os.path.join(template_path, cmakesrc)
        #     dest = os.path.join(Path.cwd(), cmakesrc)
        #     shutil.copyfile(cmakefile, dest)

        update = self.option("update-species")
        if not update:
            update = self.confirm("Update species in configure file?", False)

        if update and not (type(update) is str and update.lower() in ["false", "no"]):

            chemistry["species"] = [x.name for x in net.info.species]

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
