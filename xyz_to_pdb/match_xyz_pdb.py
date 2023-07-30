"""
read two files consecutively

"""
pdb_file = r"C:\Users\lisap\Lab_code\xyz_to_pdb\MCPB_OUT_solv_noH_nowater.pdb"
xyz_file = r"C:\Users\lisap\Lab_code\xyz_to_pdb\NPT_xyz_noH_nowater.xyz"

#solving atom connectivity in order

def read_one_pdb_residue(pdb_file: str, number: int):
    """
    Read one pdb file and retrieve atoms for one residue of specified number
    """
    lines = []
    for line in pdb_file: 
        all_fields = line.split()
        if all_fields[0] == 'ATOM' and all_fields[4] == str(number):
            lines.append(line)
    return lines

#line count both files
def line_count(file):
    with open(file, 'r') as f:
        return sum(1 for line in f)

#infer residue naming by C O N pattern
#get atom name
def get_atom_pattern(file, field_num):
    pattern = ''
    with open(file, 'r') as f:
        for line in f:
            atom_name = line.split()[field_num]
            pattern += atom_name
    return pattern


def extract_substring(match_string, substring):
    """
    count instances of substring occuring within match_string
    """
    counter = 0
    for count, letter in enumerate(match_string):
        lower = count - 1
        upper = count + len(substring) - 1
        #print(match_string[lower : upper])
        if match_string[lower:upper] == substring:
            counter += 1
    return counter

xyz_pattern = get_atom_pattern(xyz_file, 0)
print(extract_substring(xyz_pattern, 'CCONC'))

pdb_pattern = get_atom_pattern(pdb_file, 2)
print(extract_substring(pdb_pattern, 'CCONC'))