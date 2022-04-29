from naunet.reactions.umistreaction import UMISTReaction


def test_init_umistreaction_from_file(datadir):
    with open(datadir / "rate12.umist", "r") as react_file:
        for line in react_file.readlines():
            reaction = UMISTReaction(line)
            # check all non-empty reactions produce a reaction rate
            assert reaction.is_empty or reaction.rateexpr()
