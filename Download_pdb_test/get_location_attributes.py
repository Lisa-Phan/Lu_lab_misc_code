import Bio.PDB
from Bio.PDB import PDBParser
from Bio.PDB import NeighborSearch
from Bio.PDB import Selection
from Bio.PDB.Selection import unfold_entities

from itertools import chain
import pandas as pd
import pickle
import sys

""" 
Take in a file of protein structure annotation
Output a binary object with annotation info on the heteroatoms of each residue

"""


def read_pdb(PDB_code, path):
    """
    param: PDB code: the PDB code string
    param: path: file path
    get_structure takes two args, the id for the structure and the file
    """
    name = PDB_code
    pdbparser = Bio.PDB.PDBParser(QUIET=True)   # suppress PDBConstructionWarning
    try: 
        struct = pdbparser.get_structure(name, path)
        return struct
    except Exception as err:
        print(err)
        return None

exclude_list = ['SER', 'ASP', 'GLU', 'LEU',
                'PHE', 'THR', 'ILE', 'ASN',
                'VAL', 'CYS', 'PRO', 'ALA', 
                'GLN', 'GLY', 'LYS', 'TYR', 
                'ARG', 'TRP', 'HIS', 'MET', 'HOH']
    
def find_hetatom(pdb, exclude_list = exclude_list):
    """
    Find all residues in PDB files that are not in the exclude_list
    param: pdb: pdb structure object, named by the 4-letter pdb code
    param: exclude_list: list of strings that are residue names to exclude from the search
    return: a hetatom list
    """
    hetatlist = []
    if pdb is None:
        return None
    for model in pdb:
        for chain in model:
            for residue in chain:
                if residue.resname not in exclude_list:
                    hetatlist.append(residue)
    return hetatlist

def get_neighboring_coordinates(res, search_radius = 3.5):
    """Get neighboring residues within search radius of res 
    param: res: an BioPDB.Residue object, can be an amino acid or hetero atom
    return list of residue object
    """
    structure = unfold_entities(res, target_level = 'M')[0]
    ne_se = NeighborSearch(list(structure.get_atoms()))
    atoms_coord = [atoms.coord for atoms in res]
    neighbors = [ne_se.search(coord, search_radius,'R') for coord in atoms_coord ]
    unlist = list(set(chain(*neighbors)))
    return unlist

def get_resi_location_attribute(res):
    """ Take a list residue object
        Return a set of attributes for each residues
        param: res: Bio.PDB.Structure.Residue object
    """
    if res is None:
        return None
    many_resi_dict = {}
    attribute = {'Chain': unfold_entities(res, target_level = 'C')[0].id,
                                'Model': unfold_entities(res, target_level = 'M')[0].id,
                                'Neighbor': get_neighboring_coordinates(res)}
    return attribute

def apply_read_pdb(tab):
    return read_pdb(tab['Single_pdb'], tab['PDB_paths'])

def main():
    #### Read in file
    tab = pd.read_csv(path, sep = '\t')
    
    #### Create structure column
    tab['PDB_Structure_Object'] = tab.apply(apply_read_pdb, axis = 1)
    
    #### Create heteroatom column
    tab['Hetero_atom'] = tab['PDB_Structure_Object'].apply(find_hetatom)
    
    #### Make each residue a row, and get location attribute of each residue 
    longtab = tab.explode('Hetero_atom').reset_index(drop = True)
    longtab['Hetero_location_attribute'] = longtab['Hetero_atom'].apply(get_resi_location_attribute)
    
    with open('Pickle_structure_with_residue_annotation', 'wb') as outfile:
        pickle.dump(longtab, outfile)
        outfile.close()

############################################

path = sys.argv[1]

if __name__ == '__main__':
    main()