"""
Module for finding solvent accessible surface area for residues/atoms in proteins,
and other burial related analysis
"""
import sys
import os

from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.core.pose import remove_nonprotein_residues
init('-mute all')


def residue_sasas(pdb):
    """
    Calculates the solvent accessible surface areas for all residues in pdb. Returns
    dictionary of corresponding sasas

    Parameters:
        pdb (str, required): path to the pdb to be analyzed

    Returns:
        dict: keys are residue number, values are solven accesible surface area
    """

    # create pose, rename sasa calculator, score pose
    pose = pose_from_pdb(pdb)
    remove_nonprotein_residues(pose)
    scorefxn = get_fa_scorefxn()
    scorefxn(pose)

    sasa_calc = protocols.vardist_solaccess.VarSolDistSasaCalculator()
    # create sasas map
    sasa_map = sasa_calc.calculate(pose)

    # dictionary to store the resulting sasas
    residue_sasas = {}

    for i in range(1, sasa_map.size()+1):
        atom_sasas = sasa_map[i]

        for j in range(1, len(atom_sasas)+1):
            if j == 1:
                residue_sasas[i] = atom_sasas[j]
            else:
                residue_sasas[i] += atom_sasas[j]

    return residue_sasas
