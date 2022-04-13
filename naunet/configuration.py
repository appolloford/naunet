from __future__ import annotations

import tomlkit

from datetime import datetime

from .species import Species
from .templateloader import NetworkInfo

NAUNET_CONFIG_DEFAULT = """\
# Naunet config document

[general]
creation_time = ""
name = ""
description = ""

[chemistry]
elements = []
pseudo_elements = []
species = []
extra_species = []
network = []
format = []
heating = []
cooling = []
ode_modifier = []

[chemistry.dust]
model = "none"
species = []

[chemistry.binding_energy]

[chemistry.photon_yield]

[chemistry.shielding]

[chemistry.rate_modifier]

[ODEsolver]
solver = ""
device = ""
method = ""

[summary]
num_of_species = -1
num_of_gas_species = -1
num_of_ice_species = -1
num_of_dust_species = -1
num_of_reactions = -1
list_of_species = []
list_of_species_alias = []
list_of_gas_species = []
list_of_ice_species = []
list_of_dust_species = []
"""


class Configuration:
    def __init__(
        self,
        project: str,
        description: str = "",
        element: list[str] = None,
        pseudo_element: list[str] = None,
        allowed_species: list[str] = None,
        required_species: list[str] = None,
        binding_energy: dict[str, float] = None,
        photon_yield: dict[str, float] = None,
        network: list[str] = None,
        format: list[str] = None,
        heating: list[str] = None,
        cooling: list[str] = None,
        shielding: dict[str, str] = None,
        dustmodel: str = "none",
        dustspecies: list[str] = None,
        rate_modifier: dict[int, str] = None,
        ode_modifier: list[str] = None,
        solver: str = "cvode",
        device: str = "cpu",
        method: str = "dense",
        networkinfo: NetworkInfo = None,
    ) -> None:

        self._name = project
        self._description = description
        self._element = element if element else Species.known_elements()
        self._pseudoelement = (
            pseudo_element if pseudo_element else Species.known_pseudoelements()
        )
        self._species = allowed_species if allowed_species else []
        self._extraspecies = required_species if required_species else []
        self._bindingenergy = binding_energy if binding_energy else {}
        self._photonyield = photon_yield if photon_yield else {}
        self._network = network if network else []
        self._format = format if format else []
        self._heating = heating if heating else []
        self._cooling = cooling if cooling else []
        self._shielding = shielding if shielding else {}
        self._dustmodel = dustmodel
        self._dustspecies = dustspecies if dustspecies else []
        self._ratemodifier = rate_modifier if rate_modifier else {}
        self._odemodifier = ode_modifier if ode_modifier else []
        self._solver = solver
        self._device = device
        self._method = method
        self._networkinfo = networkinfo

    @property
    def content(self) -> str:
        content = tomlkit.loads(NAUNET_CONFIG_DEFAULT)

        general = content["general"]

        general["creation_time"] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        general["name"] = self._name
        general["description"] = self._description

        chemistry = content["chemistry"]
        chemistry["elements"] = self._element
        chemistry["pseudo_elements"] = self._pseudoelement
        chemistry["species"] = self._species
        chemistry["extra_species"] = self._extraspecies
        chemistry["binding_energy"] = self._bindingenergy
        chemistry["photon_yield"] = self._photonyield
        chemistry["network"] = self._network
        chemistry["format"] = self._format
        chemistry["heating"] = self._heating
        chemistry["cooling"] = self._cooling
        chemistry["shielding"] = self._shielding
        chemistry["dust"]["model"] = self._dustmodel
        chemistry["dust"]["species"] = self._dustspecies
        chemistry["rate_modifier"] = self._ratemodifier
        chemistry["ode_modifier"] = self._odemodifier

        odesolver = content["ODEsolver"]
        odesolver["solver"] = self._solver
        odesolver["device"] = self._device
        odesolver["method"] = self._method

        summary = content["summary"]
        info = self._networkinfo
        if info:
            dust = info.dust
            gas_species = [s.name for s in info.species if not s.is_surface]
            ice_species = [g.name for g in info.species if g.is_surface]
            dust_species = [d.name for d in dust.species] if dust else []
            summary["num_of_species"] = len(info.species)
            summary["num_of_gas_species"] = len(gas_species)
            summary["num_of_ice_species"] = len(ice_species)
            summary["num_of_dust_species"] = len(dust_species)
            summary["num_of_reactions"] = len(info.reactions)
            summary["list_of_species"] = [x.name for x in info.species]
            summary["list_of_species_alias"] = [x.alias for x in info.species]
            summary["list_of_gas_species"] = gas_species
            summary["list_of_ice_species"] = ice_species
            summary["list_of_dust_species"] = dust_species

        return tomlkit.dumps(content)
