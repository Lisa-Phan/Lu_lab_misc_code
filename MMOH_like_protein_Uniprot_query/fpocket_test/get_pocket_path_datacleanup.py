"""
2023/07/12
PATCHY DATA CLEAN UP LOG from ipython
"""
import os, sys
import pandas as pd

dir = sys.argv[1]

def list_files_and_directories(directory):
    everything = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            everything.append(os.path.join(root, file))
        for dir in dirs:
            everything.append(os.path.join(root, dir))
    return everything
files = list_files_and_directories(dir)

raw_out_pdb = [f for f in files if f.endswith('_out.pdb')]
raw_out_pdb = [f for f in raw_out_pdb if 'chain_only' in f]
raw_out_pdb = [f for f in raw_out_pdb if 'iterate_many_pocket' in f]

def get_params(path_string):
   from operator import itemgetter
   dir_name_with_param = os.path.split(path_string)[0]
   just_dir_name = os.path.split(dir_name_with_param)[1]
   m, M = itemgetter(1,2)(just_dir_name.split('_'))
   return m[1:], M.strip('.pdb')[1:]

min_sphere, max_sphere = [], []

for file in raw_out_pdb:
    min_par, max_par = get_params(file)
    min_sphere.append(min_par)
    max_sphere.append(max_par)

table = pd.DataFrame({'file_path': raw_out_pdb, 'min': min_sphere, 'max': max_sphere})
table.to_csv('/scratch/09069/dhp563/pathfinding_snakemake/all_pocket_output.tsv', sep = '\t', index = False)
