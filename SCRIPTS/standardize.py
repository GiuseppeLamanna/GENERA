#!/usr/bin/env python3
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.info')
import rdkit
import sys

uncharger = rdMolStandardize.Uncharger()


def molStandardize(inpSmiPath, outPath):
    """
    Standardardize the molecules in the input file. The input file must be a smi-like file.

    Parameters
    ----------
    inpSmiPath : string
        path to the .csv input file containing the smiles as the first column and nothing else.
    outPath : string
        path to the output file where the output will be written

    Returns
    -------
    None
        but will write a file in the specified location and name.

    """

    fh = open(inpSmiPath, "r")
    smis = fh.read().split("\n")
    fh.close()

    outFile = ""
    for line in smis:
        smi = line.split(",")[0].split(" ")[0]
        m = Chem.MolFromSmiles(smi,sanitize=False)
        m.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(m,sanitizeOps=(Chem.SANITIZE_ALL^Chem.SANITIZE_CLEANUP^Chem.SANITIZE_PROPERTIES))
        cm = rdMolStandardize.Normalize(m)
        cm = uncharger.uncharge(cm)
        cm = rdMolStandardize.Reionize(cm)
        outFile += Chem.MolToSmiles(cm) + "\n"

    fh = open(outPath, "w")
    fh.write(outFile)
    fh.close()
    
molStandardize(sys.argv[1], sys.argv[2])

