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

    # 13x: other special types
    GAS_XRAY = 130

    # 20x: common types of gas-grain interaction
    GRAIN_FREEZE = 200
    GRAIN_DESORB_THERMAL = 201
    GRAIN_DESORB_COSMICRAY = 202
    GRAIN_DESORB_PHOTON = 203
    GRAIN_DESORB_REACTIVE = 204

    # 21x: special desoprtion
    GRAIN_DESORB_H2 = 210  # H2 formation desorption (Roberts et al. 2007)

    # 22x: special grain involved reactions
    GRAIN_RECOMINE = 220
    GRAIN_ECAPTURE = 221

    # 30x: common types of surface reactions
    SURFACE_TWOBODY = 300
    SURFACE_COSMICRAY = 301
    SURFACE_PHOTON = 302

    # 31x: surface diffusion
    SURFACE_DIFFUSION = 310

    UNKNOWN = 999
    DUMMY = 1000
