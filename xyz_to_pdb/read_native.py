"""
Given an intact input pdb, read through the file and 

check for
Ca-Cb-N bond angle distribution
Cb-N-Ca bond angle distribution
Ca-Cb-N-Ca dihedral angle distribution
"""
from warnings import warn
from Bio.PDB import PDBParser
import numpy as np
import matplotlib.pyplot as plt



pdb_link  = r"C:\Users\lisap\Lab_code\xyz_to_pdb\MCPB_OUT_solv.pdb"

def get_bond_angle(atom1, atom2, atom3):
    """
    Calculate angle created by vectors constructed by <atom1, atom2> and <atom2, atom3>
    
    https://stackoverflow.com/questions/14066933/direct-way-of-computing-the-clockwise-angle-between-two-vectors
    dot = x1*x2 + y1*y2 + z1*z2    # Between [x1, y1, z1] and [x2, y2, z2]
    lenSq1 = x1*x1 + y1*y1 + z1*z1
    lenSq2 = x2*x2 + y2*y2 + z2*z2
    angle = acos(dot/sqrt(lenSq1 * lenSq2))
    """

    vector_a = atom1.coord - atom2.coord
    vector_b = atom3.coord - atom2.coord
    dot = np.dot(vector_a, vector_b)
    lenSq1 = np.dot(vector_a, vector_a)
    lenSq2 = np.dot(vector_b, vector_b)
    angle = np.arccos(dot/np.sqrt(lenSq1 * lenSq2))

    return angle*180/np.pi

def get_angle_constraint_CA(link):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('test', link)
    #get all CA atoms
    ca_atoms = [atom for atom in structure.get_atoms() if atom.get_name() == 'CA']
    return ca_atoms


def get_connected_atoms_to_ca(atom):
    """
    get atoms connecting to given ca
    """
    connected_atoms = []
    for atom in atom.get_parent().get_atoms():
        if atom.get_name() in ['C', 'N']:
            connected_atoms.append(atom)
    return connected_atoms

def check_CB_connection(atomC):
    """
    C is double bonded to O, to one Ca, and N from new atom
    
    """
    # get numeric assignment for residue
    ca = [atom for atom in atomC.get_parent().get_atoms() if atom.get_name() == 'CA'][0]
    new_atom_num = atomC.get_parent().id[1] + 1
    neighbor_res = atomC.get_parent().get_parent()[new_atom_num]
    try: 
        neighbor_res_N = [atom for atom in neighbor_res.get_atoms() if atom.get_name() == 'N'][0]
    except IndexError:
        warn("IndexError: {} does not have a N atom".format(neighbor_res))
        return None
    #check angle
    angle = get_bond_angle(neighbor_res_N, atomC, ca)
    return angle

def check_N_connection(atomN):
    """
    Check N-Ca-N bond angle
    
    """
    if atomN.get_parent().id[1] > 1:
        prior_neighbor = atomN.get_parent().get_parent()[atomN.get_parent().id[1] - 1]
        
        prior_neighbor_CB = [atom for atom in prior_neighbor.get_atoms() if atom.get_name() == 'C'][0]
        ca = [atom for atom in atomN.get_parent().get_atoms() if atom.get_name() == 'CA'][0]
        angle = get_bond_angle(prior_neighbor_CB, atomN, ca)
        return angle
    else:
        return None

cas = get_angle_constraint_CA(pdb_link)
angles = []
for ca in cas:
    N, C = get_connected_atoms_to_ca(ca)
    ang = check_N_connection(N)
    if ang:
        angles.append(ang)

    # ang = check_CB_connection(C)
    # if ang:
    #     angles.append(ang)
    
print(min(angles), max(angles))
#Since this number is likely to fluctuate, there needs to be a buffer range
# value for osprey is "106.29983231717422" "115.49650444967797"
#optimistic range values: 100 to 120

#Cbeta angles 
#115.46985305153659 121.54734171463559
#optimistic range values: 110 to 125