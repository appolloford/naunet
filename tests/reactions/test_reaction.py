import pytest
from naunet.species import Species
from naunet.reactions.reaction import Reaction
from naunet.reactiontype import ReactionType


@pytest.fixture
def reference():
    return Reaction(
        ["He", "CR"],
        ["He+", "e-"],
        -9999.0,
        9999.0,
        reaction_type=ReactionType.GAS_COSMICRAY,
    )


# A fixture wrapper to create reaction instance
@pytest.fixture
def reaction(request):
    return Reaction(*request.param)


@pytest.mark.parametrize(
    "reaction, reacts, prods, tmin, tmax, a, b, c, rtype, idx",
    [
        ((), [], [], -1.0, -1.0, 0.0, 0.0, 0.0, ReactionType.UNKNOWN, -1),
        (
            (
                ["He"],
                ["He+", "e-"],
                5.0,
                20.0,
                0.5,
                0.0,
                0.0,
                ReactionType.GAS_COSMICRAY,
                3,
            ),
            [Species("He")],
            [Species("He+"), Species("e-")],
            5.0,
            20.0,
            0.5,
            0.0,
            0.0,
            ReactionType.GAS_COSMICRAY,
            3,
        ),
        (
            (
                [Species("He")],
                [Species("He+"), Species("e-")],
                5.0,
                20.0,
                0.5,
                0.0,
                0.0,
                ReactionType.GAS_COSMICRAY,
                3,
            ),
            [Species("He")],
            [Species("He+"), Species("e-")],
            5.0,
            20.0,
            0.5,
            0.0,
            0.0,
            101,
            3,
        ),
    ],
    indirect=["reaction"],
)
def test_init_reaction(reaction, reacts, prods, tmin, tmax, a, b, c, rtype, idx):
    assert reaction.reactants == reacts
    assert reaction.products == prods
    assert reaction.temp_min == tmin
    assert reaction.temp_max == tmax
    assert reaction.alpha == a
    assert reaction.beta == b
    assert reaction.gamma == c
    assert reaction.reaction_type == rtype
    assert reaction.idxfromfile == idx


@pytest.mark.parametrize(
    "reaction, is_equal",
    [
        (
            (
                ["He", "CR"],
                ["He+", "e-"],
                -9999.0,
                9999.0,
                1.0,
                0.0,
                0.0,
                ReactionType.GAS_COSMICRAY,
            ),
            True,
        ),
        (
            (
                ["He", "CR"],
                ["He+", "e-"],
                -9999.0,
                9999.0,
                1.0,
                0.0,
                0.0,
                ReactionType.UNKNOWN,
            ),
            True,
        ),
        (
            (
                ["He", "Photon"],
                ["He+", "e-"],
                -9999.0,
                9999.0,
                1.0,
                0.0,
                0.0,
                ReactionType.GAS_COSMICRAY,
            ),
            True,  # pseudo-element is not in the reactants, not compared
        ),
        (
            (
                ["He", "CR"],
                ["He+", "e-"],
                0.0,
                9999.0,
                1.0,
                0.0,
                0.0,
                ReactionType.GAS_COSMICRAY,
            ),
            False,  # different temperature
        ),
        (
            (
                ["He", "CR"],
                ["He+", "e-"],
                -9999.0,
                9999.0,
                1.0,
                0.0,
                0.0,
                ReactionType.GAS_PHOTON,
            ),
            False,  # different reaction type
        ),
    ],
    indirect=["reaction"],
)
def test_eq_reaction(reference, reaction, is_equal):

    assert (reference == reaction) == is_equal
    assert (set([reference, reaction]) == set([reference])) == is_equal


def test_contains(reference):

    assert Species("He") in reference
    assert Species("He+") in reference
    assert Species("e") in reference
    assert Species("H") not in reference


def test_reaction_str():
    reac = Reaction(
        ["He", "CR"],
        ["He+", "e-"],
        0.0,
        9999.0,
        reaction_type=ReactionType.GAS_COSMICRAY,
    )

    longstr = "".join(
        [
            "He               -> He+ + e-                        ",
            ",     0.0 < T <  9999.0",
            ", Type: GAS_COSMICRAY            ",
            ", Source: unknown",
            ", Index: -1",
        ]
    )

    minimalstr = "He -> He+ + e-"

    assert str(reac) == longstr
    assert f"{reac:minimal}" == minimalstr
