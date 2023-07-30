import pandas as pd 
import matplotlib.pyplot as plt

link1=r"C:\Users\lisap\Lab_code\MMOH_like_protein_Uniprot_query\Foldseek_easysearch_result\aln_AF_P27353"
link2=r"C:\Users\lisap\Lab_code\MMOH_like_protein_Uniprot_query\Foldseek_easysearch_result\aln_chainD_6YD0.txt"
link3 = r"C:\Users\lisap\Lab_code\MMOH_like_protein_Uniprot_query\Foldseek_easysearch_result\aln_6YD0_exhaustive"
colnames = ['Query', 
            'Subject', 
            'Identity',
            'Alignment Length',
            'Mismatches',
            'Gap Openings',
            'Query Start',
            'Query End',
            'Subject Start',
            'Subject End',
            'E-Value',
            'Bit Score']

def dataformat(link):
    """
    Take a link to a tab delimited file
    Return a dataframe with the column names
    """
    tab = pd.read_csv(link, sep='\t', header=None)
    tab.columns = colnames
    return tab


def get_histograms(link, color):
    tab1=dataformat(link)
    #Plot histogram of column 2    
    tab1['Identity'].hist(color=color, bins=30)
    plt.show()
    tab1['Alignment Length'].hist(color=color, bins=30)
    plt.show()
    tab1['E-Value'].hist(color=color, bins=30)
    plt.show()

#get_histograms(link1)
#get_histograms(link2)

#How many entries of the two query matches? 
def check_match(link1, link2):
    tab1=dataformat(link1)
    tab2=dataformat(link2)

    biglist = tab1['Subject'].tolist()
    biglist.extend(tab2['Subject'].tolist())
    print(biglist)
    print(len(set(biglist)))

#check_match(link1, link2)
get_histograms(link3, 'green')