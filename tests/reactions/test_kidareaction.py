import pytest
from naunet.species import Species
from naunet.reactions.reaction import Reaction
from naunet.reactiontype import ReactionType
from naunet.reactions.kidareaction import KIDAReaction


def test_init_kidareaction():

    # initialize from reaction string
    react_string = "".join(
        [
            "He         CR                     ",
            "He+        e-                                            ",
            "5.000e-01  0.000e+00  0.000e+00 ",
            "2.00e+00 0.00e+00 ",
            "logn  1  -9999   9999  1     3 1  1 ",
        ]
    )
    kidareac = KIDAReaction(react_string)
    assert kidareac.reactants == [Species("He")]
    assert kidareac.products == [Species("He+"), Species("e-")]
    assert kidareac.temp_min == -9999.0
    assert kidareac.temp_max == 9999.0
    assert kidareac.alpha == 0.5
    assert kidareac.beta == 0.0
    assert kidareac.gamma == 0.0
    assert kidareac.reaction_type == ReactionType.GAS_COSMICRAY
    assert kidareac.source == "kida"
    assert kidareac.idxfromfile == 3


def test_eq_reaction_kidareaction():

    react_string = "".join(
        [
            "He         CR                     ",
            "He+        e-                                            ",
            "5.000e-01  0.000e+00  0.000e+00 ",
            "2.00e+00 0.00e+00 ",
            "logn  1  -9999   9999  1     3 1  1 ",
        ]
    )
    kidareac = KIDAReaction(react_string)

    reac = Reaction(
        ["He", "CR"],
        ["He+", "e-"],
        -9999.0,
        9999.0,
        reaction_type=ReactionType.GAS_COSMICRAY,
    )
    assert reac == kidareac

    reac = Reaction(
        ["He", "CR"],
        ["He+", "e-"],
        -9999.0,
        9999.0,
        reaction_type=ReactionType.GAS_COSMICRAY,
        idxfromfile=999,
    )
    assert reac == kidareac
    assert set([reac, kidareac]) == set([reac])

    reac = Reaction(
        ["He", "CR"],
        ["He+", "e-"],
        0.0,
        9999.0,
        reaction_type=ReactionType.GAS_COSMICRAY,
    )
    assert reac != kidareac


@pytest.mark.slow
def test_init_kidareaction_from_file(datadir):
    with open(datadir / "deuspin.kida", "r") as react_file:
        for line in react_file.readlines():
            reaction = KIDAReaction(line)
            # check all non-empty reactions produce a reaction rate
            assert reaction.is_empty or reaction.rateexpr()
