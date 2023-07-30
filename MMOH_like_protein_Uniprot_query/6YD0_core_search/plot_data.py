"""
Data cleaning
Take the dataframe and obtain the name of the protein using field 
after field 'Molecule: '
"""

import pandas as pd
import sys, os
import matplotlib.pyplot as plt


name = sys.argv[1]

def plot_namecount(tab):
    tab.NAME.value_counts().sort_values().plot(kind='bar')
    plt.show()
    

def main():
    tab = pd.read_csv(name, sep='\t')
    plot_namecount(tab)

if __name__ == '__main__':
    main()