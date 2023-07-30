#!/bin/python3

# a script that performs many pocket characterization on a set of pdb files

#1. collect all the sMMO pdb, make note of what the hydrolase chain is 
#2. Determine the important params, -m and -M combo
#3. Run sequentially, rename each output and move to designated folder

## Command for solvent radius stuff


#Diferrous XFEL, Diferric XFEL, ___, ___, OB3b with carboxylate 
#Bath sMMO with reg subunit, 


import os, sys
import numpy as np

code_chain_mapping = {"6YD0":'D', 
               "6YDI":'D', 
               "1MTY":'D', 
               "6VK4":'A',
               "6VK8":'A', 
               "4GAM":'A'}

UPPER = (4.5, 5)
LOWER = (3, 3.5)


#via testing with 6YD0, the best param is determined to be the default -m 3 -M 5
# combo = [ (i,j) for i in np.linspace(LOWER[0], LOWER[1], 3) \
#          for j in np.linspace(UPPER[0], UPPER[1], 3) ]
#run fpocket search

MIN = sys.argv[1]
MAX = sys.argv[2]


def test_many_param_fpocket(combo, code_chain_mapping, output_dir):
    """
    Combo is a list of (cartesian product) of the lower and upper bound search radii
    code_chain_mapping is a dictionary of pdb code and chain to use for fpocket 
    -m is the min radius, -M is the max radius

    """
    for pdb_code in list(code_chain_mapping.keys()):
        for param in combo:
            os.system(f"fpocket -f {pdb_code}.pdb " + 
                      f"-k {code_chain_mapping[pdb_code]} " +
                      f"-m {param[0]} -M {param[1]} ") 
            #rename output to reflect param used
            os.system(f"mv {pdb_code}_out {pdb_code}_m{param[0]}_M{param[1]}.pdb_out")
            #move to designated folder
            os.system('mv *_out ' + output_dir)

def main():
    os.system('mkdir fpocket_output')
    test_many_param_fpocket([(3,5)], code_chain_mapping, "fpocket_output")



if __name__ == "__main__":
    main()
