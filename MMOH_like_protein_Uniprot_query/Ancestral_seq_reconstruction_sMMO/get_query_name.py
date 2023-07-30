link = r"C:\Users\lisap\Lab_code\MMOH_like_protein_Uniprot_query\Ancestral_seq_reconstruction_sMMO\FireProt_catalytic_constr_di_iron_site_20230620\tree.tre"


outfilename = r"C:\Users\lisap\Lab_code\MMOH_like_protein_Uniprot_query\Ancestral_seq_reconstruction_sMMO\FireProt_catalytic_constr_di_iron_site_20230620\leaf_id.txt"
import re

def extract_words_between_quotes(string):
    pattern = r'"([^"]*)"'
    matches = re.findall(pattern, string)
    return matches

with open(link, "r") as f:
    tree = f.read()
    #extract all the node names
    leaf = extract_words_between_quotes(tree)
    print(leaf)
    print(len(leaf))

    with open(outfilename, 'w') as outfile:
        for item in leaf:
            outfile.write(f"{item}\n")
    