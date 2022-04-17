# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import errno
import os
import re
import sys
import urllib.parse
import shutil

import tomlkit

from importlib.metadata import version
from tomlkit.toml_file import TOMLFile

from .command import Command


class RenderCommand(Command):
    """
    Render source codes according to the project setting

    render
        {--f|force : forced to override the existing files}
        {--patch= : create patch for target code}
        {--with-pattern : render the jacobian pattern}
    """

    def __init__(self):
        super(RenderCommand, self).__init__()

    def handle(self):

        ver = version("naunet")

        config = TOMLFile("naunet_config.toml")

        content = config.read()

        chemistry = content["chemistry"]
        network = chemistry["network"]
        format = chemistry["format"]
        element = chemistry["elements"]
        pseudo_element = chemistry["pseudo_elements"]
        species = chemistry["species"]
        extra_species = chemistry["extra_species"]
        heating = chemistry["heating"]
        cooling = chemistry["cooling"]
        binding = chemistry["binding_energy"]
        yields = chemistry["photon_yield"]
        shielding = chemistry["shielding"]
        rate_modifier = chemistry["rate_modifier"]
        ode_modifier = chemistry["ode_modifier"]

        dust = chemistry["dust"]
        dustmodel = dust["model"]
        dustspecies = dust["species"]

        odesolver = content["ODEsolver"]
        solver = odesolver["solver"]
        method = odesolver["method"]
        device = odesolver["device"]
        # required = odesolver["required"]

        from pathlib import Path
        import naunet

        from naunet.species import Species

        Species.set_known_elements(element)
        Species.set_known_pseudoelements(pseudo_element)
        Species.set_dust_species(dustspecies)

        from naunet.network import Network, supported_reaction_class
        from naunet.chemistry import update_binding_energy, update_photon_yield

        update_binding_energy(binding)
        update_photon_yield(yields)

        for fmt in format:
            if not supported_reaction_class.get(fmt):
                from importlib import util

                spec = util.spec_from_file_location(fmt, f"{fmt}.py")
                module = util.module_from_spec(spec)
                spec.loader.exec_module(module)

        net = Network(
            filelist=network,
            fileformats=format,
            allowed_species=species + extra_species,
            required_species=extra_species,
            dustmodel=dustmodel,
            dustparams=dust,
            heating=heating,
            cooling=cooling,
            shielding=shielding,
        )

        net.check_duplicate_reaction()

        # Don't change the include/src/test when rendering patches
        patch = self.option("patch")
        if patch:
            pm = net.patchmaker(patch, device)
            pm.render(Path.cwd() / patch)
            return

        # If not creating patch, check whether include and src folders exist
        header_prefix = Path.cwd() / "include"
        source_prefix = Path.cwd() / "src"
        for prefix in [header_prefix, source_prefix]:

            if prefix.exists():
                if not os.path.isdir(prefix):
                    raise FileNotFoundError(
                        errno.ENOENT, os.strerror(errno.ENOENT), prefix
                    )

                elif os.listdir(prefix):
                    overwrite = self.option("force") or self.confirm(
                        f"Non-empty {prefix.name} directory. Overwrite?", False
                    )

                    if not overwrite:
                        sys.exit()

            else:
                os.mkdir(prefix)

        rate_modifier = {int(key): value for key, value in rate_modifier.items()}
        tl = net.templateloader(
            solver=solver,
            method=method,
            device=device,
            ratemodifier=rate_modifier,
            odemodifier=ode_modifier,
        )

        tl.render_constants(prefix=source_prefix, headerprefix=header_prefix)
        tl.render_macros(prefix=header_prefix)
        tl.render_naunet(prefix=source_prefix, headerprefix=header_prefix)
        tl.render_ode(prefix=source_prefix, headerprefix=header_prefix)
        tl.render_physics(prefix=source_prefix, headerprefix=header_prefix)
        tl.render_utilities(prefix=source_prefix, headerprefix=header_prefix)
        tl.render_data(prefix=header_prefix)
        tl.render_cmake(prefix=Path.cwd(), version=ver)

        pkgpath = Path(naunet.__file__).parent
        timerfile = pkgpath / "templates/common/include/naunet_timer.h"
        shutil.copyfile(timerfile, header_prefix / "naunet_timer.h")

        demo = Path.cwd() / "demo.ipynb"
        if not demo.exists():
            shutil.copyfile(pkgpath / "templates/common/demo.ipynb", demo)

        pattern = self.option("with-pattern")
        if pattern:
            tl.render_sparsity(prefix=Path.cwd())

        summary = tomlkit.table()
        dust = net.info.dust
        gas_species = [s.name for s in net.info.species if not s.is_surface]
        ice_species = [g.name for g in net.info.species if g.is_surface]
        dust_species = [d.name for d in dust.species] if dust else []
        summary["num_of_species"] = len(net.info.species)
        summary["num_of_gas_species"] = len(gas_species)
        summary["num_of_ice_species"] = len(ice_species)
        summary["num_of_dust_species"] = len(dust_species)
        summary["num_of_reactions"] = len(net.info.reactions)
        summary["list_of_species"] = [x.name for x in net.info.species]
        summary["list_of_species_alias"] = [x.alias for x in net.info.species]
        summary["list_of_gas_species"] = gas_species
        summary["list_of_ice_species"] = ice_species
        summary["list_of_dust_species"] = dust_species

        content["summary"] = summary

        config_file = Path.cwd() / "naunet_config.toml"
        with open(config_file, "w", encoding="utf-8") as f:
            f.write(tomlkit.dumps(content))

        # progress = self.progress_bar()
        # progress.finish()
