import pandas as pd

link=r'C:\Users\lisap\Lab_code\MMOH_like_protein_Uniprot_query\MMOH_Uniprot_query_20230623.tsv'

def main():
    df=pd.read_csv(link, sep='\t')
    df.to_excel(link.replace('.tsv','.xlsx'), index=False)

if __name__ == '__main__':
    main()
