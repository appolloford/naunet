from naunet.reactions.umistreaction import UMISTReaction


def test_init_umistreaction_from_file(datadir):
    with open(datadir / "rate12.umist", "r") as react_file:
        for line in react_file.readlines():
            reaction = UMISTReaction(line)
            rate = reaction.rateexpr() if not reaction.is_empty else "0.0"
            assert rate  # check all reactions can produce a reaction rate
