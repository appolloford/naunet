import pytest
from pathlib import Path
from naunet.species import Species
from naunet.dusts.hh93dust import HH93Dust
from naunet.reactions.reaction import Reaction
from naunet.reactions.reactiontype import ReactionType


@pytest.fixture
def example_depletion_reaction():
    return Reaction(
        ["C"],
        ["#C"],
        5.0,
        41000.0,
        1.0,
        0.0,
        0.0,
        reaction_type=ReactionType.GRAIN_FREEZE,
    )


@pytest.fixture
def example_thermal_desorption_reaction():
    return Reaction(
        ["#C"],
        ["C"],
        5.0,
        41000.0,
        1.0,
        0.0,
        0.0,
        reaction_type=ReactionType.GRAIN_DESORB_THERMAL,
    )


@pytest.fixture
def hh93dust():
    return HH93Dust([Species("GRAIN0"), Species("GRAIN-")])


def test_init_dust(hh93dust):
    assert "gdens" not in hh93dust.varis
    assert "double gdens = y[IDX_GRAIN0I] + y[IDX_GRAINM]" in hh93dust.locvars


def test_depletion_rate(hh93dust, example_depletion_reaction):
    with pytest.raises(ValueError, match="Gas temperature symbol is not set"):
        example_depletion_reaction.rateexpr(dust=hh93dust)

    hh93dust.sym_tgas = "Tgas"
    assert (
        example_depletion_reaction.rateexpr(dust=hh93dust)
        == "opt_frz * 1.0 * pi * rG * rG * gdens * sqrt(8.0 * kerg * Tgas/ (pi*amu*12.0))"
    )


def test_thermal_desorption_rate(hh93dust, example_thermal_desorption_reaction):
    hh93dust.sym_tdust = "Tdust"
    assert (
        example_thermal_desorption_reaction.rateexpr(dust=hh93dust)
        == "opt_thd * cov * nMono * densites * "
        "sqrt(2.0*sites*kerg*eb_GCI/(pi*pi*amu*12.0)) * exp(-eb_GCI/(Tdust))"
    )
