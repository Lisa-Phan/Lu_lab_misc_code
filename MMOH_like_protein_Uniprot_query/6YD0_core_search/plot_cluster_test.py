import pandas as pd
import matplotlib.pyplot as plt

link = r"C:\Users\lisap\Lab_code\MMOH_like_protein_Uniprot_query\6YD0_core_search\6YD0_core_clust_cluster.tsv"

tab = pd.read_csv(link, sep='\t')
#plot cluster count distribution
tab.columns = ['cluster', 'member']

tab['cluster'].value_counts().plot(kind='bar')

#most populated cluster
print(tab['cluster'].value_counts().head(10))



#plt.show()