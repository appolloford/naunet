import pytest
from pathlib import Path
from naunet.species import Species
from naunet.grains.hh93grain import HH93Grain
from naunet.reactions.reaction import Reaction
from naunet.reactiontype import ReactionType


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
def hh93grain():
    return HH93Grain([Species("GRAIN0"), Species("GRAIN-")])


def test_init_grain(hh93grain):
    assert "gdens" not in hh93grain.params.keys()
    assert "gdens" in hh93grain.deriveds.keys()
    assert "y[IDX_GRAIN0I] + y[IDX_GRAINM]" in hh93grain.deriveds.values()


def test_depletion_rate(hh93grain, example_depletion_reaction):
    assert (
        example_depletion_reaction.rateexpr(grain=hh93grain)
        == "opt_frz * 1.0 * pi * rG * rG * gdens * sqrt(8.0 * kerg * Tgas/ (pi*amu*12.0))"
    )


def test_thermal_desorption_rate(hh93grain, example_thermal_desorption_reaction):
    assert (
        example_thermal_desorption_reaction.rateexpr(grain=hh93grain)
        == "opt_thd * cov * nMono * densites * "
        "sqrt(2.0*sites*kerg*eb_GCI/(pi*pi*amu*12.0)) * exp(-eb_GCI/(Tgas))"
    )
