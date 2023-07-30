"""
Takes in a directory of pdbs to search and a mapping file
for which chain to search
run pocket finder with given param
three and five was tested to be relatively ok 
"""
import sys
import os

MAPPING_FILE = sys.argv[1]
PDB_DIRECTORY = sys.argv[2]
MIN = sys.argv[3]
MAX = sys.argv[4]

assert len(sys.argv) == 5, 'Usage: python3 fpocket_runner.py mapping_file pdb_directory min max'

def create_mapping(file):
    """
    file is a csv file with first column for pdb code
    second column as chain identifier

    return a dictionary of code: chain mapping
    """
    mapping = {}
    with open(file, 'r') as f:
        for line in f:
            key, value = line.strip().split(',')
            mapping[key] = value
    return mapping

def run_fpocket(file, min, max, chain_mapping):
    """
    run fpocket on a file provided. Assuming file format is
    CODE.pdb
    """
    code = os.path.basename(file).strip().split('.')[0]
    chain = chain_mapping[code]
    os.system(f'fpocket -f {file} -m {min} -M {max} -k {chain}')

def run_on_many_files(directory, min, max, chain_mapping):
    """
    run fpocket on all files in directory
    """
    os.system('mkdir -p fpocket_output')
    for file in os.listdir(directory):
        if file.endswith('.pdb'):
            filepath = os.path.join(directory, file)
            run_fpocket(filepath, min, max, chain_mapping)
    
    os.system('mv *_out fpocket_output')

def main():
    chain_mapping = create_mapping(MAPPING_FILE)
    run_on_many_files(PDB_DIRECTORY, MIN, MAX, chain_mapping)

if __name__ == '__main__':
    main()


