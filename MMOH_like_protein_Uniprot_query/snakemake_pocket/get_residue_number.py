"""
Selecting residues of interest from a set of 
Ca- superimposed structures

input: residues numbering, query pdb, and a directory of pdbs
output: a table with one column per queried residue and one row per pdb
assumption: all single chain pdbs, Ca backbone only

"""

import os, sys
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, NeighborSearch

QUERY_PDB = sys.argv[1]
TARGET_DATABASE = sys.argv[2]
RESLIST = sys.argv[3]

def get_min_distance(coord, structure):
    """
    param res1: numpy array of coordinates (x,y,z)
    param res2: pdb.object
    return: tuple of (min_dist, resname, resnumber)
    where min_dist is the minimum distance between the queried coord and one atoms in structure
    resname is the residue name of the closest atom
    resnumber is the residue number of the closest atom
    """
    min_dist = 1000
    for atom in structure.get_atoms():
        dist = atom.coord - coord
        if dist < min_dist:
            min_dist = dist
            #three letter name
            resname = atom.get_parent().get_resname()
            resnumber = atom.get_parent().get_id()[1]
    return (min_dist, resname, resnumber)

def euclidean_distance(coord1, coord2):
    return np.linalg.norm(coord1 - coord2)

def get_Ca_coord_from_number(pdb, res_number):
    """
    #assuming that this is a single-chain pdb object
    pdb: biopython structure object
    res_number: int
    return: numpy array of coordinates (x,y,z)
    """
    for atom in pdb.get_atoms():
        if atom.get_parent().get_id()[1] == res_number:
            return atom.coord

def get_atom_attributes(atom):
    """
    Return an atom's parent's res number, res type
    """
    return atom.get_parent().get_id()[1], atom.get_parent().get_resname()

def res_number_from_query_to_target(query_pdb, target_pdb, res_number):
    """
    target_pdb: biopython structure object
    res_list: list of residue numbers

    query pdb: biopython structure object
    target  pdb: biopython structure object
    res_number: int

    return the number of the nearest residue number in target pdb
    """
    ns = NeighborSearch(list(target_pdb.get_atoms()))    
    coord = get_Ca_coord_from_number(query_pdb, res_number)
    nearby_atoms = ns.search(coord, 3.0, level='A')
    if len(nearby_atoms):
        nearby_atoms.sort(key=lambda atom: euclidean_distance(atom.coord, coord))
        return get_atom_attributes(nearby_atoms[0])
    else:
        print('no nearby atoms found')
        return ('NA', 'NA')

def make_col(query_pdb, target_pdb_list, query_res_number):
    """
    Generate one col per query
    query_pdb: biopython structure object
    target_pdb_list: list of biopython structure object
    res_list: list of residue numbers for each target pdb
    """
    res_num_record = []
    res_type_record = []
    for target_structure in target_pdb_list:
        res_num, res_type = res_number_from_query_to_target(query_pdb, target_structure, query_res_number)
        res_num_record.append(res_num)
        res_type_record.append(res_type)

    return res_num_record, res_type_record

def file_to_pdb_obj(pdb_file):
    """
    read a file and only select smmoh chain
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_file, pdb_file)
    return structure

def main():
    pdbs = os.listdir(TARGET_DATABASE)
    full_path_pdbs = [os.path.join(TARGET_DATABASE, pdb) for pdb in pdbs]
    res_list = [int(res) for res in RESLIST.split(',')]
    #the reference structure 
    query_struct = file_to_pdb_obj(QUERY_PDB)

    #list of searched structure
    target_pdb_list = [file_to_pdb_obj(pdb) for pdb in full_path_pdbs]
    
    table = {}
    for res in res_list:
        res_num_record, res_type_record = make_col(query_struct, target_pdb_list, res)
        table[f'{res}_num'] = res_num_record
        table[f'{res}_type'] = res_type_record
    table['pdb_name'] = pdbs
    pd.DataFrame.from_dict(table).to_csv('residue_table.csv', sep = '\t', index=False)

if __name__ == '__main__':
    main()
