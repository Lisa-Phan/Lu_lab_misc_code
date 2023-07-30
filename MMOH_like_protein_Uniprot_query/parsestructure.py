"""
Test code to parse PDB structure and chain, omitting non-amino acid residues
Last mod 6/6/23
"""
import Bio.PDB
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO

import sys
import os

#read pdb file and get chain
def get_chain(pdb, chain):
    """
    pdb is the path to the pdb file
    chain is the single letter string for chain id
    """
    parser = PDBParser()
    structure = parser.get_structure(file=pdb, id=pdb.strip('.pdb'))
    try: 
        spec_chain = structure[0][chain]
    except KeyError:
        print(f'Chain {chain} not found in model 1 of pdb')
        sys.exit()
    #get amino acid residues
    reslist = [res for res in spec_chain.get_residues() if res.id[0] == ' ']
    return create_chain(reslist)

def create_chain(reslist, chain_id='A'):
    """
    Take residue list, a chain name string
    return a chain object
    """
    chain = Bio.PDB.Chain.Chain(chain_id)
    for count, res in enumerate(reslist):
        chain.insert(count, res)
    return chain

#write chain to pdb file
def write_chain(chain,filename):
    """
    chain is a chain object
    pdb is the path to the pdb file
    """
    io = PDBIO()
    io.set_structure(chain)
    io.save(filename)
