import pandas as pd
import itertools
import numpy as np
from collections import defaultdict
import re


def parseExportPhospho(phosphoCataSolo, phosphoCataComp, stateOfGeneReg, phosphoGeneIn, phosphoGeneOut, outputFile):
    '''
    parse exported files from neo4j requests to get genes controlling phospholrylation 
    output is a sif file with 'source relationType target reactionName ModifiedPositionOfPhosphoSites'
    relationType = 'controls-phospho-of', 'controls-dephospho-of','regulates-phospho-of', 'regulates-dephospho-of'
    '''
    
    # ===========================================
    # READIND EXPORTS ===========================
    # ===========================================
    regulatorsAndCataInReaction = {}

    with open(stateOfGeneReg, "r") as regGene:
        lines = regGene.readlines()

    for line in lines[1::]:
        parts = line.strip().split("\t")
        reaction_reg = parts[0]
        regulators = parts[1][1:-1].split(",")
        for regulator in regulators:
            regulator = regulator.strip(' "').replace("'", '')  # remove double quote from the metabo in list
            if reaction_reg not in regulatorsAndCataInReaction:
                regulatorsAndCataInReaction[reaction_reg] = [regulator]
            else:
                regulatorsAndCataInReaction[reaction_reg].append((regulator))

    # ------------------------- catalyseur part 

    catalystOnly = {}  # to differentiate with regulators

    with open(phosphoCataSolo, "r") as cataSolo:
        lines = cataSolo.readlines()

    for line in lines[1::]:
        parts = line.strip().split("\t")
        reaction_cata = parts[0]
        catalysers = parts[1][1:-1].split(",")
        for catalyser in catalysers:
            catalyser = catalyser.strip(' "').replace("'", '')
            if reaction_cata not in regulatorsAndCataInReaction:
                regulatorsAndCataInReaction[reaction_cata] = [catalyser]
            else:
                regulatorsAndCataInReaction[reaction_cata].append((catalyser))
            if reaction_cata not in catalystOnly:
                catalystOnly[reaction_cata] = [catalyser]
            else:
                catalystOnly[reaction_cata].append((catalyser))

    with open(phosphoCataComp, "r") as cataComp:
        lines = cataComp.readlines()

    for line in lines[1::]:
        parts = line.strip().split("\t")
        reaction_cata = parts[0]
        catalysers = parts[1][1:-1].split(",")
        for catalyser in catalysers:
            catalyser = catalyser.strip(' "').replace("'", '')
            if reaction_cata not in regulatorsAndCataInReaction:
                regulatorsAndCataInReaction[reaction_cata] = [catalyser]
            else:
                regulatorsAndCataInReaction[reaction_cata].append((catalyser))
            if reaction_cata not in catalystOnly:
                catalystOnly[reaction_cata] = [catalyser]
            else:
                catalystOnly[reaction_cata].append((catalyser))

    # ---------------- parse position of modified residue(s)

    dictPositionsOut = defaultdict(dict)

    with open(phosphoGeneOut, "r") as pos:
        lines = pos.readlines()

    for line in lines[1::]:
        parts = line.strip().split("\t")
        reaction = parts[0]
        geneP = parts[1]
        positions = parts[2][1:-1].split(",")
        for position in positions:
            position = position.strip(' "').replace("'", '')
            if position == "":
                position = "NA"
            if reaction in dictPositionsOut and geneP in dictPositionsOut[reaction]:
                dictPositionsOut[reaction][geneP].append(position)
            else:
                dictPositionsOut[reaction][geneP] = [position]


    # ---------------- parse position of modified residu in input

    dictPositionsIn = defaultdict(dict)

    with open(phosphoGeneIn, "r") as posin:
        lines = posin.readlines()

    for line in lines[1::]:
        parts = line.strip().split("\t")
        reaction = parts[0]
        geneInP = parts[1]
        positions = parts[2][1:-1].split(",")
        for position in positions:
            position = position.strip(' "').replace("'", '')
            if position == "":
                position = "NA"
            if reaction in dictPositionsIn and geneInP in dictPositionsIn[reaction]:
                dictPositionsIn[reaction][geneInP].append(position)
            else:
                dictPositionsIn[reaction][geneInP] = [position]


    # ===========================================
    # PARSING DICTS =============================
    # ===========================================

    # ------------ Create dictionary of new phosphorylated position

    newPhosphoSite = defaultdict(dict)

    for (
        reac,
        nested,
    ) in (
        dictPositionsOut.items()
    ):  # comparer les react, les genes puis les positions en output par rapport à ceux en input
        for gene, pos in nested.items():
            if (
                reac in dictPositionsIn and gene in dictPositionsIn[reac]
            ):  # soustraire les sites deja phosphorylés den input a ceux phosphorylé(s) dans cette reaction reac
                posTokeep = [p for p in pos if p not in dictPositionsIn[reac][gene]]
                if (
                    len(posTokeep) > 0
                ):  # pour eviter ce genre d'out {'P54646': [], 'Q13131': [], 'O75143': ['NA']} ou on conserve les genes en key une fois toute les positions supprimées
                    newPhosphoSite[reac][gene] = posTokeep
            if (
                reac not in dictPositionsIn or gene not in dictPositionsIn[reac]
            ):  # si la reaction n'a pas de gene phospho en input ou si le gene n'est pas phospho en input alors le(s) site(s) phosph en out le sont pendant cette reaction
                newPhosphoSite[reac][gene] = pos

    # ------- Create dictionary of dephosphorylate site

    dephosphoSite = defaultdict(dict)

    for reac, nested in dictPositionsIn.items():
        for gene, pos in nested.items():
            if reac in dictPositionsOut and gene in dictPositionsOut[reac]:
                posTokeep = [p for p in pos if p not in dictPositionsOut[reac][gene]]
                if len(posTokeep) > 0:
                    dephosphoSite[reac][gene] = posTokeep
            if reac not in dictPositionsOut or gene not in dictPositionsOut[reac]:
                dephosphoSite[reac][gene] = pos

    # ------------------------ parsing part -> Create Output for phosporylation

    coupleTokeep = []

    for react, nested in newPhosphoSite.items():
        targets = nested.keys()
        if react in regulatorsAndCataInReaction:
            sources = regulatorsAndCataInReaction.get(react, [])
            for target in targets:
                for source in sources:
                    if (
                        react in catalystOnly
                    ):  # some reaction have no activeUnit which catalyst = R-HSA-69227
                        if source in catalystOnly[react]:
                            relationType = "controls-phospho-of"
                        else:
                            relationType = "regulates-phospho-of"
                    else:
                        relationType = "regulates-phospho-of"
                    if target != source:
                        coupleTokeep.append(
                            (
                                source,
                                relationType,
                                target,
                                react,
                                newPhosphoSite[react][target],
                            )
                        )

    # ---------------- parsing part : Create Output for dephosporylation

    for react, nested in dephosphoSite.items():
        targets = nested.keys()
        if react in regulatorsAndCataInReaction:
            sources = regulatorsAndCataInReaction.get(react, [])
            for target in targets:
                for source in sources:
                    if (
                        react in catalystOnly
                    ):  # some reaction has no activeUnit which catalyst = R-HSA-69227
                        if source in catalystOnly[react]:
                            relationType = "controls-dephospho-of"
                        else:
                            relationType = "regulates-dephospho-of"
                    else:
                        relationType = "regulates-dephospho-of"
                    if target != source:
                        coupleTokeep.append(
                            (
                                source,
                                relationType,
                                target,
                                react,
                                dephosphoSite[react][target],
                            )
                        )

    # ---------------- writing outputs

    with open(outputFile, "w") as sifOut:
        sifOut.write("source\trelationType\ttarget\treaction\tpositionResiduP\n")
        for couple in coupleTokeep:
            sifOut.write(f"{couple[0]}\t{couple[1]}\t{couple[2]}\t{couple[3]}\t{couple[4]}\n")



