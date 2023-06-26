#!/usr/bin/env python3
import sys
from rdkit import Chem


def removeBadGroups(inpSmi, badGroupsPath, outFilePath):
    """
    Detects bad groups in a ref file, removes them and writes a new file without any bad groups-containing molecules

    Parameters
    ----------
    inpSmi : string
        path to the input file containing the candidate molecules. One line per SMILES
    badGroupsPath : string
        path to a text file that can be evaluated as a list of SMARTS of bad groups to search and eliminate
        from the original file.
    outFilePath: string
        path to the output file to be written.

    Returns
    -------
    None
        Will write a csv file without any bad groups according to the given list

    """

    fh = open(inpSmi, "r")
    inpSmiText = fh.read().split("\n")[:-1]
    fh.close()

    fh = open(badGroupsPath, "r")
    badGroups = eval(fh.read())
    fh.close()
    outText = ""

    for line in inpSmiText:
        bad = False
        smi = line.split(" ")[0]
        mol = Chem.MolFromSmiles(smi)
        for badGroup in badGroups:
            badMol = Chem.MolFromSmarts(badGroup)
            matches = mol.GetSubstructMatches(badMol)
            if len(matches) != 0:
                bad = True
                break
        if bad is False:
            outText = outText + line + "\n"
    

    fh = open(outFilePath, "w")
    fh.write(outText)
    fh.close()

removeBadGroups(sys.argv[1], sys.argv[2], sys.argv[3])