import os
import shutil
from pathlib import Path
from jinja2 import Environment, PackageLoader

from .species import Species


class PatchMaker:
    def __init__(
        self,
        netinfo: object,
        couple: str,
        *args,
        **kwargs,
    ) -> None:

        if couple.lower() == "enzo":
            self.patchtemplateloader = EnzoPatch(netinfo, *args, **kwargs)
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

    def __init__(
        self,
        netinfo: object,
        defined_species: list = None,
        startidx: int = 104,
    ) -> None:

        loader = PackageLoader("naunet", "templates/patches/enzo")
        self._env = Environment(loader=loader)
        self._env.globals.update(zip=zip)
        self._env.trim_blocks = True
        self._env.rstrip_blocks = True

        self.netinfo = netinfo

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
            "Grid_MultiSpeciesHandler.C",
            "hydro_rk/Grid_UpdatePrim.C",
            "hydro_rk/Grid_UpdateMHDPrim.C",
        ]:
            srcfile = os.path.join(static_patch_path, src)
            dstfile = os.path.join(prefix, src)
            shutil.copyfile(srcfile, dstfile)

    def _render_enzoheader(self, prefix: str = "./"):
        template = self._env.get_template("naunet_enzo.h.j2")

        nspec = len(self.netinfo.species)
        neq = len(self.netinfo.species)

        result = template.render(nspec=nspec, neq=neq, species=self.netinfo.species)
        with open(os.path.join(prefix, "naunet_enzo.h"), "w") as out:
            out.write(result)

    def _render_grid(self, prefix: str = "./"):
        template = self._env.get_template("Grid.h.j2")

        specnum = [
            "DeNum" if s.name in ["e-", "E-"] else f"{s.alias}Num"
            for s in self.netinfo.species
        ]
        arguments = ", ".join(f"int &{n}" for n in specnum)

        result = template.render(arguments=arguments)
        with open(os.path.join(prefix, "Grid.h"), "w") as out:
            out.write(result)

    def _render_identify(self, prefix: str = "./"):

        template = self._env.get_template("Grid_IdentifyNaunetSpeciesFields.C.j2")

        specnum = [
            "DeNum" if s.name in ["e-", "E-"] else f"{s.alias}Num"
            for s in self.netinfo.species
        ]
        typedefidx = [
            "ElectronDensity" if s.name in ["e-", "E-"] else f"{s.alias}Density"
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

        addspec = [s for s in self.netinfo.species if s not in self.grackle_species]
        addspecnum = [f"{s.alias}Num" for s in addspec]

        declare = ", ".join(addspecnum)
        prim = [f"Prim[nfield++] = BaryonField[{n}]" for n in addspecnum]
        oldprim = [f"Prim[nfield++] = OldBaryonField[{n}]" for n in addspecnum]

        specnum = [
            "DeNum" if s.name in ["e-", "E-"] else f"{s.alias}Num"
            for s in self.netinfo.species
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

        addspec = [s for s in self.netinfo.species if s not in self.grackle_species]
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
        addspec = [s for s in self.netinfo.species if s not in self.defined_species]
        idxdefine = [f"{s.alias}Density = {si+start}," for si, s in enumerate(addspec)]
        idxdefine.append(f"FieldUndefined = {start+len(addspec)};")

        result = template.render(idxdefine=idxdefine)
        with open(os.path.join(prefix, "typedefs.h"), "w") as out:
            out.write(result)

    def _render_wrapper(self, prefix: str = "./"):
        template = self._env.get_template("Grid_NaunetWrapper.C.j2")

        specnum = [
            "DeNum" if s.name in ["e-", "E-"] else f"{s.alias}Num"
            for s in self.netinfo.species
        ]
        declare = ", ".join(specnum)
        initial = " = ".join([*specnum, "0"])

        # In Enzo, electron has the same mass number as hydrogen
        abund = [
            f"y[IDX_{s.alias}] =  BaryonField[{n}][i] / 1.0"
            if s.name in ["e-", "E-"]
            else f"y[IDX_{s.alias}] =  BaryonField[{n}][i] / {s.massnumber}"
            for s, n in zip(self.netinfo.species, specnum)
        ]
        invabund = [
            f"BaryonField[{n}][i] = y[IDX_{s.alias}] * 1.0"
            if s.name in ["e-", "E-"]
            else f"BaryonField[{n}][i] = y[IDX_{s.alias}] * {s.massnumber}"
            for s, n in zip(self.netinfo.species, specnum)
        ]

        result = template.render(
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
        self._render_test(prefix)
        self._render_typedef(prefix)
        self._render_wrapper(prefix)
