from enum import IntEnum

# If the reaction types have the same formalism, they share the same value in the enum class
class ReactionType(IntEnum):
    """
    The definition of reaction types. Types are categorized by phases
    (gas, grain, surface) and then divided to sub-groups
    """

    # 10x: common types of gas phase reactions
    GAS_TWOBODY = 100
    GAS_COSMICRAY = 101
    GAS_PHOTON = 102
    GAS_THREEBODY = 103

    # 11x: special types in KIDA
    GAS_KIDA_IP1 = 110
    GAS_KIDA_IP2 = 111

    # 12x: special types in UMIST
    GAS_UMIST_CRPHOT = 120

    # 13x: special types from Walsh+2015
    GAS_LEEDS_XRAY = 130
    GAS_LEEDS_RECOM = 131
    GAS_LEEDS_ECAPTURE = 132

    # 20x: common types of gas-grain interaction
    GRAIN_FREEZE = 200
    GRAIN_DESORPT_THERMAL = 201
    GRAIN_DESORPT_COSMICRAY = 202
    GRAIN_DESORPT_PHOTON = 203

    # 21x: special desorption from Walsh+2015
    GRAIN_DESORPT_REACTIVE = 210

    # 22x: special desoprtion from UCLCHEM
    GRAIN_DESORPT_H2 = 220

    # 30x: common types of surface reactions
    SURFACE_TWOBODY = 300
    SURFACE_COSMICRAY = 301
    SURFACE_PHOTON = 302

    # 31x: surface diffusion
    SURFACE_DIFFUSION = 310

    UNKNOWN = 999