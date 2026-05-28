"""The logic here is to understand which are the genes involved in most number of BPs at differnet cut off
Currently the top 50"""

import pandas as pd
from collections import Counter
from scipy.stats import binom_test
import os
import networkx as nx
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests

path = "C:/projects/Final Rat Trans Project Result2023/blindmolerat_rnaseqanalysis-2024/DMBA/Results_10092025_f2.5_p0.01_SELECTED/entrezID/"
os.chdir(path)

# Load CSV
df = pd.read_csv("S.control-S.wound/GeneOntologyBP.csv") 
# df = pd.read_csv("MM.control-M.wound/GeneOntologyBP.csv") 
df = df.iloc[0:30, ]
df["gene_list"] = df["geneID"].str.split("/")

# Count gene occurrences across all descriptions
all_genes = [gene for sublist in df["gene_list"] for gene in sublist]
gene_counts = Counter(all_genes)

# Convert to DataFrame
gene_count_df = pd.DataFrame.from_dict(gene_counts, orient="index", columns=["count"])
gene_count_df.index.name = "GeneSymbol"
print(gene_count_df)

# Total appearances
total_genes = len(all_genes)
unique_genes = len(gene_count_df)

# Expected frequency under null
expected_prob = 1 / unique_genes

# Binomial test
gene_count_df["p_value"] = gene_count_df["count"].apply(
    lambda x: binom_test(x, n=total_genes, p=expected_prob, alternative="greater")
)

# FDR correction
gene_count_df["p_adj"] = multipletests(gene_count_df["p_value"], method="fdr_bh")[1]
gene_count_df = gene_count_df.sort_values("p_adj")
print(gene_count_df.head(20))

# --- Top 10 most frequent genes (excluding Acod1) ---
top_genes = [g for g, _ in gene_counts.most_common(11) if g != "Acod1"][:10]
print("Top 10 genes:", top_genes)

# Filter rows where Acod1 is present
acod1_rows = df[df["gene_list"].apply(lambda genes: "Acod1" in genes)]

# Keep only rows where Acod1 + at least one top gene is present
acod1_rows = acod1_rows[acod1_rows["gene_list"].apply(lambda genes: any(g in genes for g in top_genes))]

# --- Build network ---
G = nx.Graph()
G.add_node("Acod1", color='red', size=800)

for _, row in acod1_rows.iterrows():
    bp = row['Description']
    G.add_node(bp, color='skyblue', size=1500)
    G.add_edge("Acod1", bp)
    # Connect top genes if present
    for g in row['gene_list']:
        if g in top_genes:
            if g not in G:
                G.add_node(g, color='orange', size=0)  # temp size 0
            G.add_edge(g, bp)

# Update gene node sizes proportional to number of edges (BP connections)
for g in top_genes:
    if g in G:
        num_edges = len([n for n in G.neighbors(g)])
        G.nodes[g]['size'] = 300 + 200 * num_edges  # base size + multiplier

# Node attributes
colors = [G.nodes[n]['color'] for n in G.nodes]
sizes = [G.nodes[n]['size'] for n in G.nodes]

# Plot
plt.figure(figsize=(14, 10))
pos = nx.spring_layout(G, k=0.5, seed=42)
nx.draw(G, pos, with_labels=True, node_color=colors, node_size=sizes,
        font_size=10, font_weight='bold', edge_color='gray')
plt.title("Top 10 Frequent Genes and Biological Processes involving Acod1", fontsize=14)
plt.show()
