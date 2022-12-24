description = "Example: primordial chemical network from KROME"
files = "primordial.krome"
formats = "krome"
grain_model = ""
grain_symbol = "GRAIN"
surface_prefix = "#"
bulk_prefix = "@"
elements = ["e", "H", "D", "He"]
pseudo_elements = ["Photon"]
allowed_species = [
    "e-",
    "H",
    "H+",
    "H-",
    "D",
    "D+",
    "He",
    "He+",
    "He++",
    "H2",
    "H2+",
    "HD",
    "Photon",
]
extra_species = []
heating = []
cooling = [
    "CIC_HI",
    "CIC_HeI",
    "CIC_HeII",
    "CIC_He_2S",
    "RC_HII",
    "RC_HeI",
    "RC_HeII",
    "RC_HeIII",
    "CEC_HI",
    "CEC_HeI",
    "CEC_HeII",
]
shielding = {}
binding_energy = {}
photon_yield = {}
rate_modifier = {}
ode_modifier = {}
