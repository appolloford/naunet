from naunet.reactions.umistreaction import UMISTReaction


def test_init_umistreaction_from_file():
    umist_reaction = []
    with open("test/test_input/rate12.rates", "r") as react_file:
        for line in react_file.readlines():
            umist_reaction.append(UMISTReaction(line))
