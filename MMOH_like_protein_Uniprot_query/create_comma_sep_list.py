"""
Takes a file and create a csv string
"""
file = r"C:\Users\lisap\Lab_code\MMOH_like_protein_Uniprot_query\Foldseek_easysearch_result\full_id.txt"
outfile = r"C:\Users\lisap\Lab_code\MMOH_like_protein_Uniprot_query\Foldseek_easysearch_result\full_id.csv"

with open(file, 'r') as f:
    with open(outfile, 'w') as of:
        for line in f:
            of.write(line.strip() + ",") 
