#!/bin/python3

import os
import sys
from Bio.PDB import PDBParser, PDBIO

QUERY_PDB = sys.argv[1]
TARGET_DATABASE = sys.argv[2]
FILENAME = 'sMMO_align.tsv'

#code chain mapping is dictionary of pdb:chain corresponding to sMMO
from fpocket_lots import code_chain_mapping

def get_smmoh_chain(alignment_file):
    """
    Using the dictionary, search in foldseek generated alignment file
    param: alignment_file
    tsv file with col1: querypdb_chain, 
                  col2: target_chain,
                  col3: rotation matrix in csv format, xxxyyyzzz
                  col4: translation matrix
    return the translation vector and rotation matrix in a dictionary
    """
    file_mapping = {}
    with open(alignment_file, 'r') as f:
        for line in f:
            #check if the target and query chains are sMMO chains
            query, target, rotation, translation = line.split('\t')
            #get pdb codes
            query_code = query.split('.')[0]
            target_code = target.split('.')[0]
            if code_chain_mapping[query_code] == query[-1]:
                if code_chain_mapping[target_code] == target[-1]:
                    file_mapping[target_code] = create_ro_trans_dict(rotation, translation)
    
    return file_mapping

#https://github.com/joannalange/pdb-transform/blob/master/pdb-transform.py
def apply_transformation(structure, matrix):
    """
    Transforms structure based on the transformation matrix
    :param structure: biopython's structure object
    :param matrix: transformation matrix dict
    :return: transformed structure
    """
    # rotation = np.asmatrix(np.array(matrix['rotation']))
    rotation = matrix['rotation']
    # translation = np.array(matrix['translation'])
    translation = matrix['translation']

    # apply transformation to each atom
    map(lambda atom: atom.transform(rotation, translation), structure.get_atoms())

    return structure

def create_ro_trans_dict(rotation_str, translation_str):
    """
    create a dictionary matrix of 
    {
    'rotation': [[x1,x2,x3], [y1,y2,y3], [z1,z2,z3]],
    'translation: [x,y,z]
    } from strings
    """
    rotation = convert_rotation_string(rotation_str)
    translation = [float(value) for value in translation_str.split(',')]
    matrix = {
            'rotation' : rotation, 
            'translation': translation
            }

    return matrix

def convert_rotation_string(coord_string):
    coordinates = coord_string.split(',')
    result = []
    
    for i in range(0, 9, 3):
        sublist = [float(coordinates[i]), float(coordinates[i+1]), float(coordinates[i+2])]
        result.append(sublist)
    
    return result

def get_pdb_obj(pdb_file):
    """
    read a file and only select smmoh chain
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_file, pdb_file)
    chain = code_chain_mapping[pdb_file.split('.')[0]]
    return structure[0][chain]

def write_pdb_file(structure, pdb_file):
    """
    write pdb file
    params structure: biopython structure object
    params pdb_file: string
    """
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_file)


def main():
    #Calculate translation and rotation matrix for structure of interest

    os.system(f'foldseek easy-search {QUERY_PDB} {TARGET_DATABASE} {FILENAME} tmp --format-output query,target,u,t')
    os.system('mkdir -p foldseek_superimpose')
    os.system(f'mv {FILENAME} tmp foldseek_superimpose')
    
    #get the rotation and translation matrix for the sMMO chains
    file_mapping = get_smmoh_chain(f'foldseek_superimpose/{FILENAME}')

    #read pdb files, select necessary chains, apply transformation and write to file
    for pdb_files in os.listdir(TARGET_DATABASE):
        structure = get_pdb_obj(pdb_files)
        
        #apply transformation
        pdb_code = pdb_files.split('.')[0]
        structure = apply_transformation(structure, file_mapping[pdb_code])
        
        #write to file
        write_pdb_file(structure, f'foldseek_superimpose/{pdb_code}_chain{code_chain_mapping[pdb_code]}.pdb')





