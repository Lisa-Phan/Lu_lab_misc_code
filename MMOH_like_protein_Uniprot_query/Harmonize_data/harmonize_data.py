#!/bin/python

"""
Script that creates a pipeline to harmonize data into one file
Takes into account the following:
1. Raw output from foldseek
2. Cluster file from foldseek
3. Uniprot query
Last mod 23/6/16
"""


import pandas as pd
import sys

DATA1, DATA2, DATA3 = sys.argv[1], sys.argv[2], sys.argv[3]
assert len(sys.argv) == 4, 'Usage: python harmonize_data.py \
    <foldseek_raw> <cluster_file> <uniprot_query>'

#Target name is AF-{name}-F1-model_v4
DATA1_HEADER = ['query', 'target', 
                'fident', 'alnlen',
                'mismatch', 'gapopen',
                'qstart', 'qend',
                'tstart', 'tend', 
                'evalue', 'bits']

#accession is uniprot id, like Q3YA75 
DATA3_HEADER = ['id', 'accession', 'length', 
                'protein_name', 'reviewed', 
                'gene_names', 'go_id', 'xref_pdb', 
                'xref_interpro', 'structure_3d', 
                'organism_name']

def make_join_table(data1, data2, data3, outname='Harmonized_foldseek_uniprot_file.tsv'):
    """
    Create a combined dataframe from all sources 
    param: data1 is the link to the raw output from foldseek
    param: data2 is the link to the cluster file from foldseek
    param: data3 is the link to the uniprot query
    return: a dataframe with all the information
    """
    tab_foldseek = pd.read_csv(data1, sep='\t', names=DATA1_HEADER)
    tab_foldseek['target'] = tab_foldseek['target'].str.split('-').str[1]
    
    cluster_map = get_dictionary(data2)
    
    tab_uniprot = pd.read_csv(data3, sep='\t', names=DATA3_HEADER)
    tab_uniprot['cluster_membership'] = tab_uniprot['accession'].map(cluster_map)
    
    merged_file = pd.merge(tab_foldseek, tab_uniprot, left_on='target', right_on='accession')
    merged_file.to_csv(outname, sep = '\t', index=False)

def get_dictionary(data2):
    """
    Create a dictionary from data2
    """ 
    cluster_map = {}
    with open(data2, 'r') as mapFile:
        for line in mapFile:
            value, key = line.strip().split('\t')
            #second column is the member 
            key = key.split('-')[1]
            #first column is the cluster
            value = value.split('-')[1]
            cluster_map[key] = value
    return cluster_map        

if __name__ == '__main__':
    make_join_table(DATA1, DATA2, DATA3)
    








