"""
Converting xyz to series of pdb

Module that takes original pdb, a set of element-name-only xyz, 
return a set of pdb files

"""


from itertools import islice
import numpy as np
from dataclasses import dataclass

#CA are C-alphas, connected to prior N and forward C
#C are C=O, needs an O at corresponding position
#N is connected to prior C, one H an

TYPE_MAPPING = {'C': ['CA', 'C', 'NON'],
                'N': ['N', 'NON'], 
                'O': ['O', 'NON'], 
                'H': ['H', 'NON']}

BACKBONE_ATOM_CONNECTIVITY_OPTION = {'CA': ['N', 'C', 'H', 'NON'], 
                     'N': ['CA', 'H', 'C'], 
                     'C': ['CA', 'O', 'N'], 
                     'O': ['C'], 
                     'H': ['N', 'CA'], 
                     'NON': ['CA', 'NON']}

ANGLE_CONSTRAINT = {
    'CA': {"N_CA_C":(100, 120)},
    'CB': {"CA_CB_N": (110, 125)},
    'C': {"CA_C_O": (114, 131)},
}

ATOM_CONNECTIVITY_NUMBER = {
    'CA': 4,
    'C': 3,
    'O': 2,
    'N': 3,
    'H': 1, 
    'NON': 1
}

@dataclass(repr=True)
class XYZ_object: 
    """
    One xyz object contains the coordinates that makes up a molecule
    """
    name: str
    atoms: list

    def __str__(self):
        return f'{self.name} {self.coordinates}'
    
    def get_xyz(self):
        return self.atoms
    
@dataclass
class XYZ_atoms:
    """
    Atoms in XYZ file, contains name and coordinate of atoms
    """
    name: str
    coordinates: tuple
    domain: list = TYPE_MAPPING[name]
    connected_atoms: list = BACKBONE_ATOM_CONNECTIVITY_OPTION[name]
    
    def __str__(self):
        return f'{self.name} {self.coordinates}'

    def __repr__(self):
        return f'{self.name} {self.coordinates}'
    
    
#each atom is a dictionary that contains: 
#name, coordinates, connectivity, angle constraint

#Constraints to satisfy
# + Connectivity
# + Angle check

#how to break the problem down into functions
"""
Pseudo code
#if initial statting point is not specified
for every atom: 
    if atom type is C:
        assign atom to CA:
            infer corresponding CB N
            check constraint satisfaction for atoms in assignment
            if constraint is satisfied:
                find next ca based on assignment
            if not satisfied: 
                backtrack

def check_constrait(assignment):
    for atom in assignment:
        if atom is CA: 
        
"""




def get_distance(atom1: XYZ_atoms, atom2: XYZ_atoms):
    """
    Calculate distance between two atoms
    """
    return np.lin.alg.norm(atom1 - atom2)

def find_neighbor_atoms(xyz_object, distance, atom_name):
    """
    """
    pass


def create_coordinate(coord_list):
    """
    Create a coordinate object from a string
    """
    return tuple([float(value) for value in coord_list])

def create_XYZ_objs(file):
    """
    Read a file
    Get number of atoms
    Create XYZ objects
    return list of xyz objects
    """
    with open(file, 'r') as f:
        atom_number = int(f.readline().strip())
        #read in chunks of atom_number
        while True:
            lines = list(islice(f, atom_number + 2))

            for line in lines:
                print(line[0], line[1], line[-1])
            XYZ_atom_list = [XYZ_atoms(line.split()[0], create_coordinate(line.split()[1:])) for line in lines]
            print(len(XYZ_atom_list), atom_number)
            assert len(XYZ_atom_list) == atom_number, "Number of atoms in XYZ file does not match the number of atoms in the file"
            print(XYZ_atom_list)
            if not lines:
                break






    