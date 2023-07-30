"""
Test pyuul with bounding box conditions
last mod: 2023/6/12
"""    
from pyuul import utils, VolumeMaker
import os,urllib,torch
import numpy as np

import matplotlib
import matplotlib.pyplot as plt


######## GLOBAL VARIABLES ########
#Example center coordinate is the Fe of heme in 5YCE
CENTER = (0.615, 13.576, -15.285)
PDB = "5YCE.pdb"
BOX = (5,5,5)

#on Lonestar6, can run on gpu/cuda via gpu-a100 node
#set as cpu for now
DEVICE = 'cpu'


######### FUNCTIONS ##########
def get_atoms(pdb, center, measurements):
    """
    pdb is a pdb inputfile
    center is the center of the bounding box
    measurement is a set of (length, width, height) tuples
    
    Return a tensor of shape (1, n_atoms, 3) for coords
    and a list of list for atname
    """
    #get pdb via pyuul utils
    #coords is a tensor of shape (1, n_atoms, 3)
    #atname is a list of atom names
    coords, atname = utils.parsePDB(pdb)
    #iterate over the coords and check if they are part of the binding box
    
    atomList = []
    atomName = []
    #for every atom
    for i in range(coords[0].shape[0]):

        #individual atom coordinates
        coord = coords[0][i]

        #iterate over x, y, z to check if atom is in the box
        if is_atom_in_box(coord[0], center[0], measurements[0]):
            if is_atom_in_box(coord[1], center[1], measurements[1]): 
                if is_atom_in_box(coord[2], center[2], measurements[2]):
                    atomList.append(coord.tolist())
                    atomName.append(atname[0][i])
    
    return torch.Tensor([atomList]), [atomName]
                    
def is_atom_in_box(atom_coord, center, measurements):
    """
    atom_coord is either the x, y, or z value 
    center is the corresponding x, y or z center of the bounding box
    measurement is either the length(x), width(y), or height(z) of the bounding box
    return True if the atom is in the box
    """
    if atom_coord > center - measurements/2 and atom_coord < center + measurements/2:
        return True
    else:
        return False
    

def main(): 
    
    #get atoms in binding box param
    coord_tensor, atom_name = get_atoms(PDB, CENTER, BOX)
    device = DEVICE


    #standard syntax based on tutorial
    atom_channel = utils.atomlistToChannels(atom_name).to(device)
    radius = utils.atomlistToRadius(atom_name).to(device)
    
    volmaker = VolumeMaker.Voxels(device=device)
    voxellizedVolume = volmaker(coord_tensor, radius, atom_channel, resolution=0.5).to_dense()
    
    ### TODO: find out what to use as labeled trainig data

if __name__ == "__main__":
    main()
