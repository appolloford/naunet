description = "Example: deuterium network used in Hsu et al. 2021"
files = "deuterium.krome"
formats = "krome"
grain_model = ""
grain_symbol = "GRAIN"
surface_prefix = "#"
bulk_prefix = "@"

elements = ["e", "H", "D", "He", "C", "N", "O"]

pseudo_elements = ["o", "p", "m"]

element_replacement = {}

allowed_species = [
    "C",
    "C+",
    "C-",
    "C2",
    "C2+",
    "C2D",
    "C2D+",
    "C2H",
    "C2H+",
    "C2N",
    "C2N+",
    "C2O+",
    "C3",
    "C3+",
    "CCO",
    "CD",
    "CD+",
    "CD2",
    "CD2+",
    "CH",
    "CH+",
    "CH2",
    "CH2+",
    "CHD",
    "CHD+",
    "CN",
    "CN+",
    "CN-",
    "CNC+",
    "CO",
    "CO+",
    "CO2",
    "CO2+",
    "D",
    "D+",
    "D-",
    "D2O",
    "D2O+",
    "D3O+",
    "DCN",
    "DCN+",
    "DCO",
    "DCO+",
    "DNC",
    "DNC+",
    "DNO",
    "DNO+",
    "DOC+",
    "GRAIN-",
    "GRAIN0",
    "H",
    "H+",
    "H-",
    "H2DO+",
    "H2O",
    "H2O+",
    "H3O+",
    "HCN",
    "HCN+",
    "HCO",
    "HCO+",
    "HD",
    "HD+",
    "HD2O+",
    "HDO",
    "HDO+",
    "HNC",
    "HNC+",
    "HNO",
    "HNO+",
    "HOC+",
    "He",
    "He+",
    "HeD+",
    "HeH+",
    "N",
    "N+",
    "N2",
    "N2+",
    "N2D+",
    "N2H+",
    "N2O",
    "NCO+",
    "ND",
    "ND+",
    "ND2",
    "ND2+",
    "NH",
    "NH+",
    "NH2",
    "NH2+",
    "NHD",
    "NHD+",
    "NO",
    "NO+",
    "NO2",
    "NO2+",
    "O",
    "O+",
    "O-",
    "O2",
    "O2+",
    "O2D",
    "O2D+",
    "O2H",
    "O2H+",
    "OCN",
    "OD",
    "OD+",
    "OD-",
    "OH",
    "OH+",
    "OH-",
    "e-",
    "mD3+",
    "oD2",
    "oD2+",
    "oD2H+",
    "oD3+",
    "oH2",
    "oH2+",
    "oH2D+",
    "oH3+",
    "pD2",
    "pD2+",
    "pD2H+",
    "pD3+",
    "pH2",
    "pH2+",
    "pH2D+",
    "pH3+",
]

extra_species = []

heating = []

cooling = []

shielding = {}

binding_energy = {}

photon_yield = {}

rate_modifier = {}

ode_modifier = {}
