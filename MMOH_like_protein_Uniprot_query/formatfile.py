"""
2023/6/3
#data formatting analysis
#retrieve all entries marked with reviewed
#rectangularize the data
"""

import sys
import pandas as pd

link=sys.argv[1]
assert len(sys.argv)==2, 'Usage: python script.py link_to_file'

#ref
#https://stackoverflow.com/questions/50731229/split-cell-into-multiple-rows-in-pandas-dataframe


def main():
    tab = pd.read_excel(link)
    tab = tab[tab['Reviewed']=='reviewed']
    tab['PubMed ID'] = tab['PubMed ID'].str.split(';')

    #explode the PubMed ID column
    tab = tab.explode('PubMed ID')
    tab.to_excel(link.replace('.xlsx','_reviewed.xlsx'), index=False)

if __name__ == '__main__':
    main()


