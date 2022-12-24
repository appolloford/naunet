from __future__ import annotations

import tomlkit

from datetime import datetime
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .network import Network

NAUNET_CONFIG_DEFAULT = """\
# Naunet config document

[general]
creation_time = ""
name = ""
description = ""
loads = []

[chemistry]
[chemistry.species]
elements = []
pseudo_elements = []
allowed_species = []
extra_species = []

[chemistry.species.symbol]
surface = ""
bulk = ""

[chemistry.species.binding_energy]

[chemistry.species.photon_yield]

[chemistry.grain]
symbol = ""
model = ""

[chemistry.network]
files = []
formats = []

[chemistry.thermal]
heating = []
cooling = []

[chemistry.shielding]

[chemistry.rate_modifier]

[chemistry.ode_modifier]

[ODEsolver]
solver = ""
device = ""
method = ""

[summary]
num_of_elements = 0
num_of_species = 0
num_of_grains = 0
num_of_gas_species = 0
num_of_ice_species = 0
num_of_grain_species = 0
num_of_reactions = 0
list_of_elements = []
list_of_species = []
list_of_species_alias = []
list_of_gas_species = []
list_of_ice_species = []
list_of_grain_species = []
"""


class BaseConfiguration:
    def __init__(
        self,
        project: str,
        description: str = "",
        load: list[str] = None,
        element: list[str] = None,
        pseudo_element: list[str] = None,
        allowed_species: list[str] = None,
        required_species: list[str] = None,
        species_kwargs: dict[str, str] = None,
        binding_energy: dict[str, float] = None,
        photon_yield: dict[str, float] = None,
        filenames: list[str] = None,
        formats: list[str] = None,
        heating: list[str] = None,
        cooling: list[str] = None,
        shielding: dict[str, str] = None,
        grain_model: str = "",
        rate_modifier: dict[int, str] = None,
        ode_modifier: dict[str, dict[str, list[str | list[str]]]] = None,
        solver: str = "cvode",
        device: str = "cpu",
        method: str = "dense",
    ) -> None:

        self._name = project
        self._description = description
        self._load = load or []
        self._element = element.copy() if element else []
        self._pseudoelement = pseudo_element.copy() if pseudo_element else []
        self._allowedspecies = allowed_species.copy() if allowed_species else []
        self._extraspecies = required_species.copy() if required_species else []
        self._species_kwargs = species_kwargs.copy() if species_kwargs else {}
        self._bindingenergy = binding_energy.copy() if binding_energy else {}
        self._photonyield = photon_yield.copy() if photon_yield else {}
        self._filenames = filenames.copy() if filenames else []
        self._formats = formats.copy() if formats else []
        self._heating = heating.copy() if heating else []
        self._cooling = cooling.copy() if cooling else []
        self._shielding = shielding.copy() if shielding else {}
        self._grain_model = grain_model
        self._ratemodifier = rate_modifier.copy() if rate_modifier else {}
        self._odemodifier = ode_modifier.copy() if ode_modifier else {}
        self._solver = solver
        self._device = device
        self._method = method

        self._num_of_reactions = 0
        self._network_elements = []
        self._network_species = []
        self._network_alias = []
        self._network_gas_species = []
        self._network_ice_species = []
        self._network_grain_species = []
        self._network_grains = []

    @property
    def content(self) -> str:
        content = tomlkit.loads(NAUNET_CONFIG_DEFAULT)

        general = content["general"]

        general["creation_time"] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        general["name"] = self._name
        general["description"] = self._description
        general["loads"] = self._load

        chemistry = content["chemistry"]

        chem_species = chemistry["species"]
        chem_species["elements"] = self._element
        chem_species["pseudo_elements"] = self._pseudoelement
        chem_species["allowed_species"] = self._allowedspecies
        chem_species["extra_species"] = self._extraspecies
        chem_species["symbol"]["surface"] = self._species_kwargs.get(
            "surface_prefix", "#"
        )
        chem_species["symbol"]["bulk"] = self._species_kwargs.get("builk_prefix", "@")
        chem_species["binding_energy"] = self._bindingenergy
        chem_species["photon_yield"] = self._photonyield

        chem_grain = chemistry["grain"]
        chem_grain["symbol"] = self._species_kwargs.get("grain_symbol", "GRAIN")
        chem_grain["model"] = self._grain_model

        chem_network = chemistry["network"]
        chem_network["files"] = self._filenames
        chem_network["formats"] = self._formats

        chem_thermal = chemistry["thermal"]
        chem_thermal["heating"] = self._heating
        chem_thermal["cooling"] = self._cooling

        chemistry["shielding"] = self._shielding
        chemistry["rate_modifier"] = self._ratemodifier
        chemistry["ode_modifier"] = self._odemodifier

        odesolver = content["ODEsolver"]
        odesolver["solver"] = self._solver
        odesolver["device"] = self._device
        odesolver["method"] = self._method

        summary = content["summary"]
        summary["num_of_elements"] = len(self._network_elements)
        summary["num_of_species"] = len(self._network_species)
        summary["num_of_grains"] = len(self._network_grains)
        summary["num_of_gas_species"] = len(self._network_gas_species)
        summary["num_of_ice_species"] = len(self._network_ice_species)
        summary["num_of_grain_species"] = len(self._network_grain_species)
        summary["num_of_reactions"] = self._num_of_reactions
        summary["list_of_elements"] = self._network_elements
        summary["list_of_species"] = self._network_species
        summary["list_of_species_alias"] = self._network_alias
        summary["list_of_gas_species"] = self._network_gas_species
        summary["list_of_ice_species"] = self._network_ice_species
        summary["list_of_grain_species"] = self._network_grain_species

        return tomlkit.dumps(content)


class NetworkConfiguration(BaseConfiguration):
    def __init__(
        self,
        project: str,
        network: Network,
        description: str = "",
        solver: str = "cvode",
        device: str = "cpu",
        method: str = "dense",
    ) -> None:
        super().__init__(
            project,
            description=description,
            solver=solver,
            device=device,
            method=method,
        )

        binding = {s.name: s.eb for s in network.species if s.is_surface}
        yields = {s.name: s.photon_yield for s in network.species if s.is_surface}

        self._element = network._known_elements.copy()
        self._pseudoelement = network._known_pseudo_elements.copy()
        self._allowedspecies = network.allowed_species.copy()
        self._extraspecies = network.required_species.copy()
        self._species_kwargs = network._species_kwargs.copy()
        self._bindingenergy = binding
        self._photonyield = yields
        self._filenames = ["reactions.naunet"]
        self._format = ["naunet"]
        self._heating = network._heating_names.copy()
        self._cooling = network._cooling_names.copy()
        self._shielding = network.shielding.copy()
        self._grain_model = network.grain_model
        self._ratemodifier = network.rate_modifier.copy()
        self._odemodifier = network.ode_modifier.copy()

        species = network.species
        grains = network.grains
        self._num_of_reactions = len(network.reactions)
        self._network_elements = [x.name for x in network.elements]
        self._network_species = [x.name for x in species]
        self._network_alias = [x.alias for x in species]
        self._network_gas_species = [s.name for s in species if not s.is_surface]
        self._network_ice_species = [s.name for s in species if s.is_surface]
        self._network_grain_species = [s.name for g in grains for s in g.species]
