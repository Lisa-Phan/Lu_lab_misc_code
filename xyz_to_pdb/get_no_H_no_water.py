import Bio
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB import *


complete_pdb = r"C:\Users\lisap\Lab_code\xyz_to_pdb\MCPB_OUT_solv.pdb"

parser = PDBParser(QUIET=True)
structure = parser.get_structure('test', complete_pdb)
#remove all hetam residues
print([residue.get_full_id() for residue in structure.get_residues()])

not_hetatm = [residue for residue in structure.get_residues() if residue.get_id()[0] == ' ']

def create_struct_from_reslist(reslist):
    """
    create structure from list of atoms
    """
    new_structure = Bio.PDB.Structure.Structure('test')
    new_model = Bio.PDB.Model.Model(0)
    new_chain = Bio.PDB.Chain.Chain('A')
    new_structure.add(new_model)
    new_model.add(new_chain)
    for res in reslist:
        try:
            new_chain.add(res)
        except Bio.PDB.PDBExceptions.PDBConstructionException:
            continue
    return new_structure

#write to new pdb file
io = PDBIO()
io.set_structure(create_struct_from_reslist(not_hetatm))
io.save(r"C:\Users\lisap\Lab_code\xyz_to_pdb\MCPB_nowater1.pdb")

