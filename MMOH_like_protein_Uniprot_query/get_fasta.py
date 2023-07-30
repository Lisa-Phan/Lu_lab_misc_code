#!/bin/python3

""" 
Create a script that takes a directory of .pdb
and return a fasta file
Last mod: 2023/06/12
"""

import os
import sys
from Bio.PDB import PDBParser


OUTFILE = 'sMMO_seq.fasta'
AA_DICT = {
    'ALA': 'A', 
    'ARG': 'R', 
    'ASN': 'N', 
    'ASP': 'D', 
    'CYS': 'C',
    'GLU': 'E',
    'GLN': 'Q',
    'GLY': 'G',
    'HIS': 'H',
    'ILE': 'I',
    'LEU': 'L',
    'LYS': 'K',
    'MET': 'M',
    'PHE': 'F',
    'PRO': 'P',
    'SER': 'S',
    'THR': 'T',
    'TRP': 'W',
    'TYR': 'Y',
    'VAL': 'V'
}


def get_name(pdb_file):
    """
    extract Uniprot ID from alphafold pdb file
    """
    return pdb_file.split('-')[1]

def pdb2fasta(pdb_dir):
    with open(OUTFILE, 'w') as outFile:
        file_list = os.listdir(pdb_dir)
        for file in file_list:
            if file.endswith('.pdb'):
                pdb_file = os.path.join(pdb_dir, file)
                name = get_name(file)
                seq = get_sequence(pdb_file)s
                outFile.write(f'>{name}\n{seq}\n')

def get_sequence(pdb_file, aa_alphabet = AA_DICT):
    """
    Return sequence from pdb file
    Assuming that the pdb is single model single chain, with chain id A
    """
    parser = PDBParser()
    structure = parser.get_structure('X', pdb_file)
    res_list = [res.get_resname() for res in structure[0]['A'].get_residues()]
    seq = ''.join([aa_alphabet[res] for res in res_list])
    return seq

def main():
    pdb_dir = sys.argv[1]
    pdb2fasta(pdb_dir)

if __name__ == '__main__':
    main()

