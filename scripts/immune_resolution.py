"""The logic here is to understand which genes invloved in inflammatory BP that is common between Spalax and mouse and what are the other BPs
Currently the top 50"""

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

# Your common genes
a = "Acod1/Adam8/Adora2b/Anxa1/Apoa1/Apoe/C1qtnf3/Ccr5/Clcf1/Cnr2/Daglb/Duoxa2/Ggt1/Il17ra/Il1r1/Il6/Mefv/Mmp8/Nfkbiz/Nlrp12/Nlrp3/Nod2/Nr1d1/Osm/Pglyrp1/Ptgis/Ptgs2/S100a9/Serpine1/Snca/Socs3/Tlr2/Tnc/Tnf/Tnfaip3/Tslp/Zfp36".split("/")
# Mouse
b = "Acod1/Ada/Adam8/Adipoq/Adora1/Adora3/Agt/Aoc3/Ccl24/Clcf1/Esr1/Fcer1a/Ffar4/Ggt1/Gpx2/Il1b/Il22ra2/Il33/Il6/Isl1/Lpl/Mmp8/Nlrp3/Npy/Nr1d2/Osm/Pbk/Pglyrp1/Pla2g5/Ptger3/Ptgs2/S100a9/Scgb1a1/Serpine1/Sucnr1/Tnc/Tnfsf11/Tnfsf4/Wfdc1/Ccn4".split("/")
a = set(a)
b = set(b)

common_genes = ["Acod1", "Nlrp3", "Il6"] #a & b
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

plt.savefig("nlrp3_acod1.png", dpi=600, bbox_inches='tight')