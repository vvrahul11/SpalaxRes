"""The logic here is to plot a network of the spalax core immune activation and regulatory pathways genes"""

import pandas as pd
from collections import Counter
import os
import networkx as nx
import matplotlib.pyplot as plt

path = "C:/projects/Final Rat Trans Project Result2023/blindmolerat_rnaseqanalysis-2024/DMBA/Results_10092025_f2.5_p0.01_SELECTED/entrezID/"
os.chdir(path)

# Load CSV
df = pd.read_csv("S.control-S.wound/GeneOntologyBP.csv") 
#df = pd.read_csv("MM.control-M.wound/GeneOntologyBP.csv") 
df["gene_list"] = df["geneID"].str.split("/")

common_genes = ["Zfp36", "Socs3", "Nlrp3", "Nlrp12", "Anxa1", "Acod1"] #a & b
# Filter rows containing at least one common gene
df_filtered = df[df["gene_list"].apply(lambda genes: any(g in genes for g in common_genes))]
#import pdb; pdb.set_trace()
df_filtered = df_filtered.iloc[0:] #20

# Count occurrences for all common genes
all_genes = [gene for sublist in df_filtered["gene_list"] for gene in sublist if gene in common_genes]
gene_counts = Counter(all_genes)


# --- Build network ---
G = nx.Graph()

# Add BP nodes
for _, row in df_filtered.iterrows():
    bp = row['Description']
    if bp.startswith("negative"):
                G.add_node(bp, color='grey', size=1500)
    for g in row['gene_list']:
        if g in common_genes:
            if g not in G:
                G.add_node(g, color='red', size=0)  # temporary size
            if bp.startswith("negative"):
                import pdb
                #pdb.set_trace()
                G.add_edge(g, bp)

# Update gene node sizes proportional to number of BP connections
for g in common_genes:
    if g in G:
        num_edges = len([n for n in G.neighbors(g)])
        G.nodes[g]['size'] = 300 + 200 * num_edges  # base + multiplier

# Node attributes
colors = [G.nodes[n]['color'] for n in G.nodes]
sizes = [G.nodes[n]['size'] for n in G.nodes]

# Plot network
plt.figure(figsize=(30, 18))
pos = nx.spring_layout(G, k=0.9, seed=4)
pos = nx.kamada_kawai_layout(G)
nx.draw(G, pos, with_labels=True, node_color=colors, node_size=sizes,
        font_size=10, font_weight='bold', edge_color='gray')
plt.title("Biological Processes involving Selected Common Genes", fontsize=14)
#plt.show()
plt.savefig("C:/projects/Final Rat Trans Project Result2023/blindmolerat_rnaseqanalysis-2024/nlrp3_acod1.png", dpi=600, bbox_inches='tight')