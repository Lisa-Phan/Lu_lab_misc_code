import pickle
import sys
import os
import pyrosetta
from pyrosetta import *
from pyrosetta.teaching import *
init()

def main():
        if len(sys.argv) != 2:
                print('Usage: python3 script <pdbfile>')

        sfxn = get_score_function()
        pose = pose_from_pdb(sys.argv[1])
        relax = pyrosetta.rosetta.protocols.relax.FastRelax()
        relax.set_scorefxn(sfxn)
	relax.apply(pose)
	pose.pdb_dump('relaxed.pdb')


if __name__ == '__main__':
        main()

