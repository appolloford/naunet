import sympy as sym

element_list = None

pseudo_element_list = None

species_list = None

surface_symbol = "#"

charge_symbols = ["+", "-"]

user_symbols = {
    "CRIR": sym.symbols("zeta"),
    "Temperature": sym.symbols("Tgas"),
    "VisualExtinction": sym.symbols("Av"),
    "UVPHOT": sym.symbols("uv"),
}
