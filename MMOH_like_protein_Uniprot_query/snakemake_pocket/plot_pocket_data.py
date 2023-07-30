"""
Plotting code for pocket analyis
"""

#dictionary for selecting pockets belonging to main channel
#subjective, since the solvent accessible area varies a lot
#also, this is handpicked based on literature opinion, 
#not subjective

import pandas as pd
import os,sys
import matplotlib.pyplot as plt


main_channel_pockets = {'6VK8': [1], 
                        '6VK4': [1,3,4],
                        '6YDI': [2,24,25,20],
                        '6YD0': [1,10,24],
                        '4GAM': [1,2],
                        '1MTY': [1,2,4] }

file = sys.argv[1]


#for each pdb, select pocket
table = pd.read_csv(file, sep='\t')

main_pocket = pd.DataFrame()
for index, rows in table.iterrows():
    code = rows['pdb_code']
    if rows['pocket_number'] in main_channel_pockets[code]:
        print(code, rows['pocket_number'])
        main_pocket = main_pocket.append(rows)

main_pocket = main_pocket.drop_duplicates(subset=['pdb_code', 'pocket_number'])

#since some structure have multiple pockets 
#a rough estimate of volume would be to sum them
vol = main_pocket.groupby(['pdb_code']).sum(numeric_only=True)
print(vol[['Number of Alpha Spheres', 'Total SASA', 'Volume']])
print(vol.columns)