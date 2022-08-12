import pytest
from pathlib import Path
from naunet.species import Species
from naunet.dusts.rr07dust import RR07Dust, RR07DustX
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
def rr07dust():
    return RR07Dust([])


@pytest.fixture
def rr07dustx():
    return RR07DustX([])


def test_init_dust(rr07dust, rr07dustx):
    assert "gdens" in rr07dust.varis
    assert "gdens" in rr07dustx.varis
    assert "opt_thd" in rr07dustx.varis


def test_depletion_rate(rr07dust, example_depletion_reaction):
    with pytest.raises(ValueError, match="Gas temperature symbol is not set"):
        example_depletion_reaction.rateexpr(dust=rr07dust)

    rr07dust.sym_tgas = "Tgas"
    assert (
        example_depletion_reaction.rateexpr(dust=rr07dust)
        == "4.57e4 * 1.0 * gxsec * fr * sqrt(Tgas / 12.0)"
    )


def test_thermal_desorption_rate(
    rr07dust, rr07dustx, example_thermal_desorption_reaction
):
    rr07dust.sym_tdust = "Tdust"
    with pytest.raises(NotImplementedError, match="not implemented in rr07"):
        example_thermal_desorption_reaction.rateexpr(dust=rr07dust)

    rr07dustx.sym_tdust = "Tdust"
    assert (
        example_thermal_desorption_reaction.rateexpr(dust=rr07dustx)
        == "mantabund > 1e-30 ? (opt_thd * sqrt(2.0*sites*kerg*eb_GCI/"
        "(pi*pi*amu*12.0)) * 2.0 * densites * exp(-eb_GCI/Tdust)) : 0.0"
    )
