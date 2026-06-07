import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from scipy.cluster.hierarchy import linkage, leaves_list

path = "C:/projects/Final Rat Trans Project Result2023/blindmolerat_rnaseqanalysis-2024/DMBA/Results_10092025_f2.5_p0.01_SELECTED/genesymbol/"
os.chdir(path)
##### All DEG BP comparison
# regulation of inflammatory response
# Spalax
a = "Acod1/Adam8/Adora2b/Anxa1/Apoa1/Apoe/C1qtnf3/Ccr5/Clcf1/Cnr2/Daglb/Duoxa2/Ggt1/Il17ra/Il1r1/Il6/Mefv/Mmp8/Nfkbiz/Nlrp12/Nlrp3/Nod2/Nr1d1/Osm/Pglyrp1/Ptgis/Ptgs2/S100a9/Serpine1/Snca/Socs3/Tlr2/Tnc/Tnf/Tnfaip3/Tslp/Zfp36".split("/")
# Mouse
b = "Acod1/Ada/Adam8/Adipoq/Adora1/Adora3/Agt/Aoc3/Ccl24/Clcf1/Esr1/Fcer1a/Ffar4/Ggt1/Gpx2/Il1b/Il22ra2/Il33/Il6/Isl1/Lpl/Mmp8/Nlrp3/Npy/Nr1d2/Osm/Pbk/Pglyrp1/Pla2g5/Ptger3/Ptgs2/S100a9/Scgb1a1/Serpine1/Sucnr1/Tnc/Tnfsf11/Tnfsf4/Wfdc1/Ccn4".split("/")


# Find intersection
a = set(a)
b = set(b)
common_genes = a & b
unique_a = a - b
unique_b = b - a

print(f"Number of common genes: {len(common_genes)}")

# Load data
df1 = pd.read_csv("S.control-S.wound/S.control-S.wound_genes.csv") 
df2 = pd.read_csv("MM.control-M.wound/MM.control-M.wound_genes.csv")

# Keep only common genes
df1 = df1[df1['genename'].isin(common_genes)][['genename','log2FoldChange']]
df2 = df2[df2['genename'].isin(common_genes)][['genename','log2FoldChange']]

df1.rename(columns={'log2FoldChange': 'Spalax wound log2FC'}, inplace=True)
df2.rename(columns={'log2FoldChange': 'Mouse wound log2FC'}, inplace=True)

df_merge = pd.merge(df1, df2, on='genename', how='outer').set_index('genename')

# Cluster rows (genes) using hierarchical clustering
link = linkage(df_merge.fillna(0), method='average', metric='euclidean')
gene_order = leaves_list(link)
df_merge = df_merge.iloc[gene_order]

# Plot composite figure
fig, axes = plt.subplots(1, 2, figsize=(16, 8))
#fig.subplots_adjust(left=0.08, right=0.95, top=0.95, bottom=0.08, wspace=0.8, hspace=0.4)

# --- Left: Venn diagram ---
common_genes_ordered = [gene for gene in df_merge.index if gene in common_genes]
venn = venn2([a, b], set_labels=("Spalax wound", "Mouse wound"), ax=axes[0])

# Inside labels: only counts
if venn.get_label_by_id("10"):
    venn.get_label_by_id("10").set_text(f"(n={len(unique_a)})")
if venn.get_label_by_id("01"):
    venn.get_label_by_id("01").set_text(f"(n={len(unique_b)})")
if venn.get_label_by_id("11"):
    venn.get_label_by_id("11").set_text("\n".join(sorted(common_genes)) + f"\n\n(n={len(common_genes)})")
if venn.get_label_by_id("11"):
    venn.get_label_by_id("11").set_text(
        "\n".join(common_genes_ordered) + f"\n\n(n={len(common_genes_ordered)})"
    )
# Place unique gene names OUTSIDE circles
# Place unique gene names OUTSIDE circles
c0 = venn.get_circle_center(0)
c1 = venn.get_circle_center(1)

axes[0].annotate("\n".join(sorted(unique_a)),
                 xy=(c0.x, c0.y),
                 xytext=(-0.7, 0), ha='right', va='center',
                 fontsize=9)

axes[0].annotate("\n".join(sorted(unique_b)),
                 xy=(c1.x, c1.y),
                 xytext=(0.7, 0), ha='left', va='center',
                 fontsize=9, )

#axes[0].set_title("Venn Diagram of Gene Sets", fontsize=14)

# --- Right: Heatmap with clustering ---
sns.heatmap(
    df_merge,
    cmap="vlag",
    center=0,
    annot=True,
    linewidths=0.5,
    ax=axes[1]
)
axes[1].set_ylabel("Gene names", fontsize=12, labelpad=30)  # increases distance from axis
axes[1].yaxis.set_label_coords(-0.15, 0.9)  # x-position, y-position in axes fraction

#fig.savefig("gene_venn_heatmap.png", dpi=600, bbox_inches='tight')
#fig.savefig("gene_venn_heatmap.pdf", bbox_inches='tight')

plt.tight_layout()
plt.show()
