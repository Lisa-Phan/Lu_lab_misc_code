import pandas as pd

link = r"C:\Users\lisap\Lab_code\MMOH_like_protein_Uniprot_query\fpocket_test\path_tracing.txt"

with open(link, 'r') as f:
    for line in f: 
        print(line.split(' '))

tab = pd.read_csv(link, sep=' ')
print(tab.head())