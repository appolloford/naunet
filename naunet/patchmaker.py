import os
import shutil
from pathlib import Path
from jinja2 import Environment, PackageLoader

from .species import Species
from .templateloader import NetworkInfo
from .utilities import _stmwrap


class PatchMaker:
    def __init__(
        self,
        netinfo: NetworkInfo,
        couple: str,
        device: str,
        *args,
        **kwargs,
    ) -> None:

        if couple.lower() == "enzo":
            self.patchtemplateloader = EnzoPatch(netinfo, device, *args, **kwargs)
        else:
            raise ValueError("Not supported target!")

    def render(self, prefix: str = "./"):
        self.patchtemplateloader.render(prefix)


class EnzoPatch:

    # define both e- and E- to exclude electron
    enzo_defined_species_name = [
        "e-",
        "E-",
        "H",
        "H+",
        "He",
        "He+",
        "He++",
        "H-",
        "H2",
        "H2+",
        "D",
        "D+",
        "HD",
        "C",
        "C+",
        "O",
        "O+",
        "Si",
        "Si+",
        "Si++",
        "CH",
        "CH2",
        "CH3+",
        "C2",
        "CO",
        "HCO+",
        "OH",
        "H2O",
        "O2",
    ]

    enzo_defined_alias = [
        "Electron",
        "HI",
        "HII",
        "HeI",
        "HeII",
        "HeIII",
        "HM",
        "H2I",
        "H2II",
        "DI",
        "DII",
        "HDI",
        "CI",
        "CII",
        "OI",
        "OII",
        "SiI",
        "SiII",
        "SiIII",
        "CHI",
        "CH2I",
        "CH3II",
        "C2I",
        "COI",
        "HCOII",
        "OHI",
        "H2OI",
        "O2I",
    ]

    # repeated define e- and E- to exclude electron
    grackle_species_name = [
        "e-",
        "E-",
        "H",
        "H+",
        "He",
        "He+",
        "He++",
        "H-",
        "H2",
        "H2+",
        "D",
        "D+",
        "HD",
    ]

    grackle_defined_alias = [
        "De",
        "HI",
        "HII",
        "HeI",
        "HeII",
        "HeIII",
        "HM",
        "H2I",
        "H2II",
        "DI",
        "DII",
        "HDI",
    ]

    def __init__(
        self,
        netinfo: object,
        device: str,
        defined_species: list = None,
        startidx: int = 104,
    ) -> None:

        loader = PackageLoader("naunet", "templates/patches/enzo")
        self._env = Environment(loader=loader)
        self._env.filters["stmwrap"] = _stmwrap
        self._env.globals.update(zip=zip)
        self._env.trim_blocks = True
        self._env.rstrip_blocks = True

        self.netinfo = netinfo
        self.device = device

        # filter the enzo-defined species if species are not in current network
        specname_in_network = [s.name for s in self.netinfo.species]
        self.defined_species_name = (
            defined_species if defined_species else EnzoPatch.enzo_defined_species_name
        )
        self.defined_species_name = [
            s for s in self.defined_species_name if s in specname_in_network
        ]
        self.grackle_defined_species_name = [
            s for s in EnzoPatch.grackle_species_name if s in specname_in_network
        ]

        self.defined_species = [Species(s) for s in self.defined_species_name]
        self.grackle_species = [Species(s) for s in self.grackle_defined_species_name]
        self.startidx = startidx

    def _copy_static_patches(self, prefix="./"):
        import naunet

        src_parent_path = Path(naunet.__file__).parent
        static_patch_path = os.path.join(
            src_parent_path, "templates", "patches", "enzo"
        )

        for src in [
            "Make.config.assemble",
            "Make.config.objects",
            "Make.config.settings",
            "Make.mach.linux-gnu",
            "global_data.h",
            "InitializeNew.C",
            "ReadParameterFile.C",
            "WriteParameterFile.C",
            "SetDefaultGlobalValues.C",
            "Grid_MultiSpeciesHandler.C",
            "hydro_rk/Grid_UpdatePrim.C",
            "hydro_rk/Grid_UpdateMHDPrim.C",
        ]:
            srcfile = os.path.join(static_patch_path, src)
            dstfile = os.path.join(prefix, src)
            shutil.copyfile(srcfile, dstfile)

    def _render_derived_field(self, prefix: str = "./") -> None:
        """
        Generate the derived fields of number density of species for yt. Saved
        in name of "derived_fields_of_network.py"

        Args:
            prefix (str, optional): Path to save output file. Defaults to "./".

        """

        derived_species_field = []

        for s in self.netinfo.species:
            alias = "Electron" if s.iselectron else s.alias
            derived_species_field.append(
                "\n".join(
                    [
                        f"@derived_field(name='{alias}_ndensity', sampling_type='cell')",
                        f"def {alias}_ndensity(field, data):",
                        f"    if 'enzo' not in data.ds.dataset_type:",
                        f"        return",
                        f"    if data.ds.parameters['MultiSpecies'] < 4:",
                        f"        return",
                        f"    dunit = data.ds.mass_unit/data.ds.length_unit**3",
                        f"    num_unit = dunit / mh_cgs / {1.0 if s.iselectron else s.massnumber}",
                        f"    arr = (num_unit*data['{alias}_Density']).to_ndarray()",
                        f"    return arr",
                    ]
                )
            )

        with open(os.path.join(prefix, "derived_fields_of_network.py"), "w") as outf:
            outf.write("import yt\n")
            outf.write("from yt import derived_field\n\n")
            outf.write("mh_cgs = float(yt.units.mh_cgs)\n\n")
            outf.write("\n\n".join(derived_species_field))

    def _render_enzoheader(self, prefix: str = "./"):
        template = self._env.get_template("naunet_enzo.h.j2")

        nspec = len(self.netinfo.species)
        has_thermal = True if self.netinfo.heating or self.netinfo.cooling else False
        neqns = nspec + has_thermal

        repeatspec = [
            s
            for s in self.netinfo.species
            if s.alias in self.grackle_defined_alias or s.iselectron
        ]
        nspec_in_enzo = nspec + len(self.grackle_defined_alias) - len(repeatspec)
        nspec_in_enzo -= 1  # don't count electron, enzo treat it specially

        result = template.render(
            nspec=nspec_in_enzo, neqns=neqns, species=self.netinfo.species
        )
        with open(os.path.join(prefix, "naunet_enzo.h"), "w") as out:
            out.write(result)

    def _render_grid(self, prefix: str = "./"):
        template = self._env.get_template("Grid.h.j2")

        specnum = [
            "DeNum" if s.iselectron else f"{s.alias}Num" for s in self.netinfo.species
        ]
        arguments = ", ".join(f"int &{n}" for n in specnum)

        result = template.render(arguments=arguments)
        with open(os.path.join(prefix, "Grid.h"), "w") as out:
            out.write(result)

    def _render_identify(self, prefix: str = "./"):

        template = self._env.get_template("Grid_IdentifyNaunetSpeciesFields.C.j2")

        specnum = [
            "DeNum" if s.iselectron else f"{s.alias}Num" for s in self.netinfo.species
        ]
        typedefidx = [
            "ElectronDensity" if s.iselectron else f"{s.alias}Density"
            for s in self.netinfo.species
        ]
        arguments = ", ".join(f"int &{n}" for n in specnum)
        initial = " = ".join(specnum)
        initial = " = ".join([initial, "0"])
        findfield = [
            f"{n} = FindField({t}, FieldType, NumberOfBaryonFields);"
            for n, t in zip(specnum, typedefidx)
        ]

        check = [f'if ({n} < 0) ENZO_VFAIL("{n}=%" ISYM "\\n", {n})' for n in specnum]

        result = template.render(
            arguments=arguments,
            initial=initial,
            findfield=findfield,
            check=check,
        )
        with open(
            os.path.join(prefix, "Grid_IdentifyNaunetSpeciesFields.C"), "w"
        ) as out:
            out.write(result)

    def _render_rkptr(self, prefix: str = "./"):

        addspec = [
            s
            for s in self.netinfo.species
            if s.alias not in self.grackle_defined_alias and not s.iselectron
        ]
        addspecnum = [f"{s.alias}Num" for s in addspec]

        declare = ", ".join(addspecnum)
        prim = [f"Prim[nfield++] = BaryonField[{n}]" for n in addspecnum]
        oldprim = [f"Prim[nfield++] = OldBaryonField[{n}]" for n in addspecnum]

        specnum = [
            "DeNum" if s.iselectron else f"{s.alias}Num" for s in self.netinfo.species
        ]
        identify = ", ".join(specnum)

        template = self._env.get_template("hydro_rk/Grid_ReturnHydroRKPointers.C.j2")
        result = template.render(declare=declare, identify=identify, prim=prim)

        with open(
            os.path.join(prefix, "hydro_rk/Grid_ReturnHydroRKPointers.C"), "w"
        ) as out:
            out.write(result)

        template = self._env.get_template("hydro_rk/Grid_ReturnOldHydroRKPointers.C.j2")
        result = template.render(declare=declare, identify=identify, oldprim=oldprim)

        with open(
            os.path.join(prefix, "hydro_rk/Grid_ReturnOldHydroRKPointers.C"), "w"
        ) as out:
            out.write(result)

    def _render_test(self, prefix: str = "./"):

        addspec = [
            s
            for s in self.netinfo.species
            if s.alias not in self.grackle_defined_alias and not s.iselectron
        ]
        addspecnum = [f"{s.alias}Num" for s in addspec]

        addname = [f"{f'{s.alias}Name':<20}   = \"{s.alias}_Density\"" for s in addspec]
        addlabel = [f"DataLabel[count++] = (char*) {s.alias}Name" for s in addspec]

        declare = ", ".join(addspecnum)

        addfieldtype = [
            f"FieldType[{s.alias}Num  = NumberOfBaryonFields++] = {s.alias}Density"
            for s in addspec
        ]

        template = self._env.get_template("hydro_rk/CollapseMHD3DInitialize.C.j2")
        result = template.render(addname=addname, addlabel=addlabel)
        with open(
            os.path.join(prefix, "hydro_rk/CollapseMHD3DInitialize.C"), "w"
        ) as out:
            out.write(result)

        template = self._env.get_template(
            "hydro_rk/Grid_CollapseMHD3DInitializeGrid.C.j2"
        )
        result = template.render(declare=declare, addfieldtype=addfieldtype)
        with open(
            os.path.join(prefix, "hydro_rk/Grid_CollapseMHD3DInitializeGrid.C"), "w"
        ) as out:
            out.write(result)

        template = self._env.get_template("hydro_rk/TurbulenceInitialize.C.j2")
        result = template.render(addname=addname, addlabel=addlabel)
        with open(os.path.join(prefix, "hydro_rk/TurbulenceInitialize.C"), "w") as out:
            out.write(result)

        template = self._env.get_template("hydro_rk/Grid_TurbulenceInitializeGrid.C.j2")
        result = template.render(declare=declare, addfieldtype=addfieldtype)
        with open(
            os.path.join(prefix, "hydro_rk/Grid_TurbulenceInitializeGrid.C"), "w"
        ) as out:
            out.write(result)

    def _render_typedef(self, prefix: str = "./"):

        template = self._env.get_template("typedefs.h.j2")

        start = self.startidx
        addspec = [
            s
            for s in self.netinfo.species
            if s.alias not in self.enzo_defined_alias and not s.iselectron
        ]
        idxdefine = [f"{s.alias}Density = {si+start}," for si, s in enumerate(addspec)]
        idxdefine.append(f"FieldUndefined = {start+len(addspec)};")

        result = template.render(idxdefine=idxdefine)
        with open(os.path.join(prefix, "typedefs.h"), "w") as out:
            out.write(result)

    def _render_updateelectron(self, prefix: str = "./") -> None:
        """
        Render the Grid_UpdateElectronDensity.C file for Enzo

        Args:
            prefix (str, optional): Path to save output file. Defaults to "./".
        """

        template = self._env.get_template("hydro_rk/Grid_UpdateElectronDensity.C.j2")

        addspec = [
            s
            for s in self.netinfo.species
            if s.alias not in self.grackle_defined_alias and not s.iselectron
        ]
        addspecnum = [f"{s.alias}Num" for s in addspec]

        declare = ", ".join(addspecnum)

        specnum = [
            "DeNum" if s.iselectron else f"{s.alias}Num" for s in self.netinfo.species
        ]
        identify = ", ".join(specnum)

        ionchargesum = " + ".join(
            f"{s.charge:.1f} * BaryonField[{n}][i] / {s.massnumber}"
            for s, n in zip(addspec, addspecnum)
            if s.charge > 0
        )

        result = template.render(
            declare=declare, identify=identify, ionchargesum=ionchargesum
        )

        with open(
            os.path.join(prefix, "hydro_rk/Grid_UpdateElectronDensity.C"), "w"
        ) as out:
            out.write(result)

    def _render_wrapper(self, prefix: str = "./"):
        """
        Render Grid_NaunetWrapper.C for Enzo

        Args:
            prefix (str, optional): Path to save output file. Defaults to "./".
        """
        template = self._env.get_template("Grid_NaunetWrapper.C.j2")

        specnum = [
            "DeNum" if s.iselectron else f"{s.alias}Num" for s in self.netinfo.species
        ]
        declare = ", ".join(specnum)
        initial = " = ".join([*specnum, "0"])

        abund = []
        invabund = []

        for s, n in zip(self.netinfo.species, specnum):

            # In Enzo, electron has the same mass number as hydrogen
            massnum = 1.0 if s.iselectron else s.massnumber

            if self.device == "cpu":
                abund.append(
                    f"y[IDX_{s.alias}] = max(BaryonField[{n}][igrid], 1e-40) * NumberDensityUnits / {massnum}"
                )
                invabund.append(
                    f"BaryonField[{n}][igrid] = max(y[IDX_{s.alias}] * {massnum} / NumberDensityUnits, 1e-40)"
                )

            elif self.device == "gpu":
                abund.append(
                    f"y[sidx + IDX_{s.alias}] = max(BaryonField[{n}][igrid], 1e-40) * NumberDensityUnits / {massnum}"
                )
                invabund.append(
                    f"BaryonField[{n}][igrid] = max(y[sidx + IDX_{s.alias}] * {massnum} / NumberDensityUnits, 1e-40)"
                )

        result = template.render(
            device=self.device,
            declare=declare,
            initial=initial,
            abund=abund,
            invabund=invabund,
        )
        with open(os.path.join(prefix, "Grid_NaunetWrapper.C"), "w") as out:
            out.write(result)

    def render(self, prefix: str = "./"):
        os.makedirs(os.path.join(prefix))
        os.makedirs(os.path.join(prefix, "hydro_rk"))
        self._copy_static_patches(prefix)
        self._render_enzoheader(prefix)
        self._render_grid(prefix)
        self._render_identify(prefix)
        self._render_rkptr(prefix)
        self._render_updateelectron(prefix)
        self._render_test(prefix)
        self._render_typedef(prefix)
        self._render_wrapper(prefix)
        self._render_derived_field(prefix)
