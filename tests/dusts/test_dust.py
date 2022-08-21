import pytest
from pathlib import Path
from naunet.species import Species
from naunet.grains.grain import Grain
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
def grain():
    return Grain()


def test_init_dust(grain):
    assert "gdens" in grain.varis
    # assert "gdens" in grain.locvars


def test_depletion_rate(grain, example_depletion_reaction):
    assert (
        example_depletion_reaction.rateexpr(grain=grain)
        == "1.0 * pi * rG * rG * gdens * sqrt(8.0 * kerg * Tgas/ (pi*amu*12.0))"
    )


def test_thermal_desorption_rate(grain, example_thermal_desorption_reaction):
    with pytest.raises(NotImplementedError, match="not implemented in base"):
        example_thermal_desorption_reaction.rateexpr(grain=grain)
