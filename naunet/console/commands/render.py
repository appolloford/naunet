# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import errno
import os
import re
import sys
import urllib.parse
import shutil

import tomlkit

from tomlkit.toml_file import TOMLFile

from .command import Command


class RenderCommand(Command):
    """
    Render source codes according to the project setting

    render
        {--f|force : forced to override the existing files}
        {--patch= : create patch for target code}
        {--patch-source= : source directory of patch}
        {--with-pattern : render the jacobian pattern}
        {--with-summary-py : render with a pymodule having network summary information}
    """

    def __init__(self):
        super(RenderCommand, self).__init__()

    def handle(self):

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

        dupes = net.find_duplicate_reaction()
        print(f"The following {len(dupes)} reactions are duplicate:\n")
        print("\n".join([str(dup[1]) for dup in dupes]))

        # if elements are all capital, checking the repeat species in capital
        cap_species_name = all([e == e.upper() for e in element])
        patchname = self.option("patch")
        source = self.option("patch-source")
        # Don't change the include/src/test when rendering patches
        if patchname:
            patch = net.patch(patchname, device, cap_species_name, source)
            patch.render(path=Path.cwd() / patchname)
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

        tl.render(path=Path.cwd())

        pkgpath = Path(naunet.__file__).parent

        demo = Path.cwd() / "demo.ipynb"
        if not demo.exists():
            shutil.copyfile(pkgpath / "templates/base/demo.ipynb", demo)

        pattern = self.option("with-pattern")
        if pattern:
            tl.render_jac_pattern()

        summary = tomlkit.table()
        dust = net.info.dust
        all_species = [x.name for x in net.info.species]
        all_alias = [x.alias for x in net.info.species]
        gas_species = [s.name for s in net.info.species if not s.is_surface]
        ice_species = [g.name for g in net.info.species if g.is_surface]
        dust_species = [d.name for d in dust.species] if dust else []
        summary["num_of_species"] = len(net.info.species)
        summary["num_of_gas_species"] = len(gas_species)
        summary["num_of_ice_species"] = len(ice_species)
        summary["num_of_dust_species"] = len(dust_species)
        summary["num_of_reactions"] = len(net.info.reactions)
        summary["list_of_species"] = all_species
        summary["list_of_species_alias"] = all_alias
        summary["list_of_gas_species"] = gas_species
        summary["list_of_ice_species"] = ice_species
        summary["list_of_dust_species"] = dust_species

        content["summary"] = summary

        config_file = Path.cwd() / "naunet_config.toml"
        with open(config_file, "w", encoding="utf-8") as f:
            f.write(tomlkit.dumps(content))

        if self.option("with-summary-py"):
            with open("summary.py", "w") as outf:

                outf.write(f"nspec = {len(all_species)}\n")
                outf.write(f"ngas = {len(gas_species)}\n")
                outf.write(f"nice = {len(ice_species)}\n")
                outf.write(f"ndust = {len(dust_species)}\n")
                outf.write(f"nreac = {len(net.info.reactions)}\n")
                outf.write("\n\n")

                speclistlist = [
                    all_species,
                    all_alias,
                    gas_species,
                    ice_species,
                    dust_species,
                ]
                namelist = [
                    "all_species",
                    "all_alias",
                    "gas_species",
                    "ice_species",
                    "dust_species",
                ]

                for speclist, name in zip(speclistlist, namelist):

                    specstr = ",\n    ".join(f"'{x}'" for x in speclist)
                    if specstr:
                        outf.write("".join([f"{name} = [\n    ", specstr, ",\n]"]))
                        outf.write("\n\n")

                for ele in Species.known_elements():
                    specnatom = [s.element_count.get(ele, 0) for s in net.info.species]
                    eledictstr = f",\n    ".join(
                        f"'{s}': {natom}"
                        for s, natom in zip(all_species, specnatom)
                        if natom
                    )
                    if eledictstr:
                        outf.write(
                            "".join(
                                [
                                    f"{ele}_dict = {{\n    ",
                                    eledictstr,
                                    f",\n}}",
                                ]
                            )
                        )
                        outf.write("\n\n")

        # progress = self.progress_bar()
        # progress.finish()
