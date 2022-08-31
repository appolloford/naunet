import pytest
from pathlib import Path
from naunet.species import Species
from naunet.grains.rr07grain import RR07Grain, RR07XGrain
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
def rr07grain():
    return RR07Grain()


@pytest.fixture
def rr07xgrain():
    return RR07XGrain()


def test_init_grain(rr07grain, rr07xgrain):
    assert "gdens" in rr07grain.params.keys()
    assert "gdens" in rr07xgrain.params.keys()
    assert "opt_thd" in rr07xgrain.params.keys()


def test_depletion_rate(rr07grain, example_depletion_reaction):
    assert (
        example_depletion_reaction.rateexpr(grain=rr07grain)
        == "4.57e4 * 1.0 * gxsec * fr * sqrt(Tgas / 12.0)"
    )


def test_thermal_desorption_rate(
    rr07grain, rr07xgrain, example_thermal_desorption_reaction
):
    with pytest.raises(NotImplementedError, match="not implemented in rr07"):
        example_thermal_desorption_reaction.rateexpr(grain=rr07grain)

    assert (
        example_thermal_desorption_reaction.rateexpr(grain=rr07xgrain)
        == "mantabund > 1e-30 ? (opt_thd * sqrt(2.0*sites*kerg*eb_GCI/"
        "(pi*pi*amu*12.0)) * 2.0 * densites * exp(-eb_GCI/Tgas)) : 0.0"
    )
