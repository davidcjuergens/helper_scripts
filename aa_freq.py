from pyrosetta import *
import os
init()
import numpy as np
import json

# some conversion dictionaries 
aas = 'ARNDCEQGHILKMFPSTWYV'
aas = [char for char in aas]
# dict for turing letter into a single number
aa2idx = {char:i for i, char in enumerate(aas)}
idx2aa = {i:char for i, char in enumerate(aas)}




def get_pdb_list(path):
    """
    List of all pdbs in a folder.
    """
    print(os.listdir(path))
    return [path + '/' + p for p in os.listdir(path) if 'npz' in p]

def count_aa(npz):
    """
    Counts the amino acids in a single pose
    """
    data = np.load(npz)

    seq = str(data['seq'])

    start_count = {char:i for i, char in enumerate(aas)}

    for aa in seq:
        start_count[aa] += 1

    return start_count

def total_aa_freq(path):
    """
    Count all the amino acids in a folder of pdbs. Return relative frequencies 
    """
    master_count = {char:i for i, char in enumerate(aas)}
    npzs = get_pdb_list(path)

    for npz in npzs:

        aa_dict = count_aa(npz)
        for aa in aa_dict.keys():
            master_count[aa] += aa_dict[aa]

    # at this point master_count should have a count of all amino acids

    total_residues = np.sum(master_count[key] for key in master_count.keys())

    master_count = {key:master_count[key]/total_residues for key in master_count.keys()}

    return master_count

def main():
    path = '/home/davidcj/input_features/900_centroid_features_train'

    freqs = total_aa_freq(path)
    print('THIS IS FREQS')
    print(freqs)

    with open('freqs_900.json', 'w') as f:
        json.dump(freqs, f)
    

main()    
