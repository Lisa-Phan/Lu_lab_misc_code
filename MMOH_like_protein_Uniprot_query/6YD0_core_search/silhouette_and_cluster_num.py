"""
Plotting silhouette score as a function of cluster number
"""

import pandas as pd
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt

# Sample text data 
data_link = r"C:\Users\lisap\Lab_code\MMOH_like_protein_Uniprot_query\6YD0_core_search\header_info_name.tsv"

# Convert the text data into a DataFrame
df = pd.read_csv(data_link, sep='\t')

# TF-IDF vectorizer
vectorizer = TfidfVectorizer(stop_words='english')
tfidf_matrix = vectorizer.fit_transform(df['NAME'])

# List to store silhouette scores
sse_values = []

# Range of cluster numbers to try
cluster_range = range(10, 30)  

# Calculate silhouette score for each cluster number
for num_clusters in cluster_range:
    kmeans = KMeans(n_clusters=num_clusters, random_state=42)
    kmeans.fit(tfidf_matrix)
    sse_values.append(kmeans.inertia_)

# Plot the silhouette score as a function of the cluster number
plt.plot(cluster_range, sse_values, marker='o')
plt.xlabel('Number of Clusters')
plt.ylabel('Squared sum error')
plt.title('Squared sum error as a Function of Cluster Number')
plt.grid(True)
plt.show()