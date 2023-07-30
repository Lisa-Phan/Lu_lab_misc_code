import pandas as pd
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

# Sample text data
data_link = r"C:\Users\lisap\Lab_code\MMOH_like_protein_Uniprot_query\6YD0_core_search\header_info_name.tsv"

# Convert the text data into a DataFrame
df = pd.read_csv(data_link, sep='\t')

# TF-IDF vectorizer
vectorizer = TfidfVectorizer(stop_words='english')
tfidf_matrix = vectorizer.fit_transform(df['NAME'])

# KMeans clustering
num_clusters = 18  # Modifiable
kmeans = KMeans(n_clusters=num_clusters, random_state=42)
kmeans.fit(tfidf_matrix)

# Assign cluster labels to the DataFrame
df['cluster_label'] = kmeans.labels_

# Evaluate the clustering using silhouette score
silhouette_avg = silhouette_score(tfidf_matrix, kmeans.labels_)
print(f"Silhouette Score: {silhouette_avg}")

# Display the clusters
for cluster_num in range(num_clusters):
    print(f"\nCluster {cluster_num}:")
    cluster_data = df[df['cluster_label'] == cluster_num]
    print(cluster_data['NAME'].to_string(index=False))

print(df['cluster_label'].value_counts())
df.to_csv(r"C:\Users\lisap\Lab_code\MMOH_like_protein_Uniprot_query\6YD0_core_search\header_info_name_clustered.tsv", sep='\t', index=False)