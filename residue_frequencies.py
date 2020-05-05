from pyrosetta import *
import os
init()

import numpy as np
import json
import argparse


def get_pdb_list(path):
    """
    List of all pdbs in a folder.
    """
    print(os.listdir(path))
    return [path + '/' + p for p in os.listdir(path) if 'pdb' in p]


def append_protein_residues(pdb, count_dict):
    """
    Given the path to a pdb, count all the amino acids and return a dictionary of
    the counts.
    """
    pose = pose_from_pdb(pdb)
    nres = pyrosetta.rosetta.core.pose.nres_protein(pose) # number of residues
    for i in range(1, nres+1):
        r = pose.residue(i)     # residue obj
        rname = r.name()[:3]    # residue name

        count_dict[rname] += 1

    return count_dict


def normalize_frequencies(count_dict):
    """
    Given the total counts of all residues in a set, return a dictionary of their
    normalized frequencies summing to 1.
    """
    total = 0.
    frequencies = {}
    for key in count_dict.keys():
        total += count_dict[key]
        frequencies[key] = 0

    for key in frequencies.keys():
        frequencies[key] = count_dict[key] / total

    return frequencies


def make_json_frequencies(input_folder, output_name):
    # add a command line argument parser here
    # path = args.path

    # some conversion dictionaries
    alpha_1 = list("ARNDCQEGHILKMFPSTWYV-")
    states = len(alpha_1)
    alpha_3 = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
               'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','GAP']

    aa_1_N = {a:n for n,a in enumerate(alpha_1)}
    aa_N_1 = {n:a for n,a in enumerate(alpha_1)}
    aa_1_3 = {a:b for a,b in zip(alpha_1,alpha_3)}
    aa_3_1 = {b:a for a,b in zip(alpha_1,alpha_3)}

    # starting counts for all amino acids in set
    total_counts = {aa:0 for aa in alpha_3}

    # get the pdb_list
    pdb_list = get_pdb_list(input_folder)

    # add up the occurence of every AA in a pose
    for i, pdb in enumerate(pdb_list):
        print('On protein {}.'.format(i))
        total_counts = append_protein_residues(pdb, total_counts)


    frequencies = normalize_frequencies(total_counts)

    with open(output_name, 'w') as f:
        json.dump(frequencies, f)


def main():
    parser = argparse.ArgumentParser(description="Amino acid frequency calculator for pdbs in folder",
                                     epilog="Have fun")

    parser.add_argument('--path',
                        '-p',
                        action='store',
                        type=str,
                        help='Path to source pdbs from.')
    
    parser.add_argument('--out',
                        '-o',
                        action='store',
                        type=str,
                        default='freqs.json',
                        help='Path to output json file.')
    
    args = parser.parse_args()
    
    input_folder = args.path
    output_name = args.out

    make_json_frequencies(input_folder, output_name)

main()
