# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from .command import Command
from ...species import Species
from ...reactions.reaction import Reaction
from ...reactions.reactiontype import ReactionType


class ExtendCommand(Command):
    """
    Extend an existing chemical network

    extend
        {input : The input filename}
        {output : The output filename}
        {--input-format=naunet : The format of the input file}
        {--output-format=naunet : The format of the output file}
        {--append-depletion : Append neutral species depletion reactions}
        {--append-thermal-desorption : Append thermal desorption reaction for surface species}
        {--append-photon-desorption : Append photon desorption reaction for surface species}
        {--append-cosmic-ray-desorption : Append cosmic-ray desorption reaction for surface species}
        {--limit-species= : Limit the reactions only involving the species}
        {--remove-duplicate : Remove the duplicate reactions from the network}
        {--remove-species= : Remove all reactions involving the species}
        {--surface-prefix= : Prefix of surface species}
    """

    def __init__(self):
        super(ExtendCommand, self).__init__()

    def handle(self):

        inp = self.argument("input")
        out = self.argument("output")
        informat = self.option("input-format")
        outformat = self.option("output-format")

        if self.option("surface-prefix"):
            Species.surface_prefix = self.option("surface-prefix")
        sprefix = Species.surface_prefix

        from naunet.network import Network, supported_reaction_class

        for fmt in [informat, outformat]:
            if not supported_reaction_class.get(fmt):
                from importlib import util

                spec = util.spec_from_file_location(fmt, f"{fmt}.py")
                module = util.module_from_spec(spec)
                spec.loader.exec_module(module)

        net = Network(filelist=inp, fileformats=informat)

        allowed_species = self.option("limit-species")
        if allowed_species:
            reactions = net.reaction_list
            allowed_species = [s.strip() for s in allowed_species.split(",") if s]
            newlist = [
                reac
                for reac in reactions
                if all(
                    rp.name in allowed_species for rp in reac.reactants + reac.products
                )
            ]
            net = Network(newlist)

        remove_species = self.option("remove-species")
        if remove_species:
            remove_species = [s.strip() for s in remove_species.split(",") if s]
            remove_index = set()
            for spec in remove_species:
                remove_index.update(net.where_species(spec))
            net.remove_reaction(list(remove_index))

        if self.option("remove-duplicate"):
            dupes = net.find_duplicate_reaction()
            net.remove_reaction([dup[0] for dup in dupes])

        initializer = supported_reaction_class.get(informat, Reaction)
        kept_reactions = [reac.react_string for reac in net.reaction_list]
        with open(f"dropped_reactions.{informat}", "w") as outf:
            with open(f"{inp}") as inpf:
                for line in inpf.readlines():
                    if initializer.preprocessing(line) not in kept_reactions:
                        outf.write(line)

        if self.option("append-depletion"):
            species = net.reactants | net.products
            for spec in species:
                # check the species is neutral gas species
                if spec.name == spec.gasname and spec.charge == 0:
                    prod = f"{sprefix}{spec}"
                    reaction = Reaction(
                        [spec],
                        [prod],
                        alpha=1.0,
                        reaction_type=ReactionType.GRAIN_FREEZE,
                    )
                    net.add_reaction(reaction)
                # TODO: electron?
                # elif spec.iselectron:
                #     pass

        options = ["thermal", "photon", "cosmic-ray"]
        rtypes = [
            ReactionType.GRAIN_DESORB_THERMAL,
            ReactionType.GRAIN_DESORB_PHOTON,
            ReactionType.GRAIN_DESORB_COSMICRAY,
        ]
        for option, rtype in zip(options, rtypes):
            if self.option(f"append-{option}-desorption"):
                species = net.reactants | net.products
                for spec in species:
                    if spec.is_surface:
                        prod = f"{spec.gasname}"
                        reaction = Reaction(
                            [spec],
                            [prod],
                            alpha=1.0,
                            reaction_type=rtype,
                        )
                        net.add_reaction(reaction)

        source, sink = net.find_source_sink()
        if source:
            print(f"Found sources: {source}")
        if sink:
            print(f"Found sinks: {sink}")

        net.reindex()
        net.write(out, outformat)