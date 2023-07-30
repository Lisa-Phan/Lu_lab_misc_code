import json
from Bio import Entrez
from ete3 import NCBITaxa, Tree, TreeStyle
import pandas as pd

ID_FILE = ''
TREE_FILE = ''

def accession2taxid(acc: str, db="protein") -> str:
    """Get refseq taxid from refseq sequence accession number"""
    handle = Entrez.esearch(db=db, term=acc)
    record = Entrez.read(handle)
    gi = record["IdList"][0]
    handle = Entrez.esummary(db=db, id=gi, retmode="json")
    result = json.load(handle)["result"]
    taxid = result[gi]["taxid"]
    return str(taxid)

#### Get mapping from refseq id with function 
mapping = {}
with open(ID_FILE, 'r') as f:
    for line in f:
        refseqid = line.strip()
        mapping[refseqid] = accession2taxid(refseqid)
    f.close()

ncbi = NCBITaxa()
tree_string = TREE_FILE.open().read()
tree = Tree(tree_string, format=1)

#get a sense of the organisms
for key in mapping:
    print(key, ncbi.get_taxid_translator([mapping[key]]))

def find_replace_tree_string(tab, str):
    for index, row in tab.iterrows():
        str = str.replace('"' + row['refseq_id'] + '"', row['name'])
    return str

def plot_tree(tree_string_file, annotated_file):
    import pandas as pd
    tree_string = open(tree_string_file).read()
    annotated_file = pd.read_csv(annotated_file, sep = '\t')
    named_string = find_replace_tree_string(annotated_file, tree_string)
    tree = Tree(named_string, format=1)
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.branch_vertical_margin = 10
    tree.show(tree_style=ts)

desired_ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

###### Converion to genera for comparison with published analyses ####### 

def get_desired_ranks(taxid, desired_ranks):
    """
    taxid is an int
    desired_ranks is a list of strings
    return a dictionary 
    """
    global ncbi
    lineage = ncbi.get_lineage(taxid)   
    names = ncbi.get_taxid_translator(lineage)
    lineage2ranks = ncbi.get_rank(names)
    ranks2lineage = dict((rank,taxid) for (taxid, rank) in lineage2ranks.items())
    return{'{}_id'.format(rank): ranks2lineage.get(rank, '<not present>') for rank in desired_ranks}


def do_mapping(taxids, desired_ranks):
    results = list()
    for taxid in taxids:
        results.append(list())
        results[-1].append(str(taxid))
        ranks = get_desired_ranks(taxid, desired_ranks)
        for key, rank in ranks.items():
            if rank != '<not present>':
                results[-1].append(list(ncbi.get_taxid_translator([rank]).values())[0])
            else:
                results[-1].append(rank)

    #generate the header
    header = ['Original_query_taxid']
    header.extend(desired_ranks)
    taxonomic_df = pd.DataFrame(results, columns = header)
    return taxonomic_df

def change_tree_annotation(desired_rank, taxonomic_df, tree_string_file):
    tree_string = open(tree_string_file).read()
    for index, row in taxonomic_df.iterrows():
        tree_string = tree_string.replace('"' + row['refseq_id'] + '"', row[desired_rank])
    tree = Tree(tree_string, format=1)
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.branch_vertical_margin = 10
    tree.show(tree_style=ts)