"""
Silly implementation
Does not work since bioPDB can't properly read the converted xyz file or original xyz file
"""

from Bio.PDB import PDBParser, PDBIO, Superimposer

intact_PDB = r"C:\Users\lisap\Lab_code\xyz_to_pdb\MCPB_OUT_solv.pdb"
frag_xyz = r"C:\Users\lisap\Downloads\NPT-pos-1-trimmed.xyz"
frag_converted = r"C:\Users\lisap\Lab_code\xyz_to_pdb\62_Osprey.pdb"


parser = PDBParser()
struct1 = parser.get_structure('intact', intact_PDB)
struct2 = parser.get_structure('frag', frag_converted)

C_NAME = ['CB', 'C', 'CA']

def align_on_cs(struct1, struct2):
    c1 = [atoms for atoms in struct1.get_atoms() if atoms.get_name() in C_NAME]
    c2 = [atoms for atoms in struct2.get_atoms() if atoms.get_name() in C_NAME]
    print(len(c1), len(c2))
    assert len(c1) == len(c2), 'Number of C atoms do not match'
    super_imposer = Superimposer()
    super_imposer.set_atoms(c1, c2)
    super_imposer.apply(struct2.get_atoms())

align_on_cs(struct1, struct2)