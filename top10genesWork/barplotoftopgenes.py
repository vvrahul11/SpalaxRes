import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.patches import Patch

# Load data
csv_path = "C:/projects/Final Rat Trans Project Result2023/blindmolerat_rnaseqanalysis-2024/top10genesWork/Top20genesUpandDownAll.csv"
outdir = os.path.dirname(csv_path)
topgenes = pd.read_csv(csv_path)

# Extract subsets
s_wound = topgenes.iloc[:, [0, 1, 2]]
m_wound = topgenes.iloc[:, [3, 4, 5]]
s_scab  = topgenes.iloc[:, [6, 7, 8]]
m_scab  = topgenes.iloc[:, [9, 10, 11]]

datasets = [
    (s_wound, "Spalax wound up- and down-regulated genes"),
    (m_wound, "Mouse wound up- and down-regulated genes"),
    (s_scab,  "Spalax scab up- and down-regulated genes"),
    (m_scab,  "Mouse scab up- and down-regulated genes"),
]

# Global x-axis limits
all_vals = pd.concat([df.iloc[:, 1].astype(float) for df, _ in datasets])
x_min, x_max = all_vals.min(), all_vals.max()

# Set up figure
fig, axes = plt.subplots(1, 4, figsize=(22, 10), sharex=True)

# Define custom legend handles
legend_handles = [
    Patch(facecolor="red", label="Up"),
    Patch(facecolor="darkgreen", label="Down"),
]

for i, (ax, (df, title)) in enumerate(zip(axes, datasets)):
    genes = df.iloc[:, 0]
    logfc = df.iloc[:, 1].astype(float)
    y_pos = np.arange(len(genes))

    colors = ["red" if val > 0 else "darkgreen" for val in logfc]

    ax.barh(y_pos, logfc, color=colors)
    ax.axvline(0, color="black", linewidth=0.8)

    ax.set_xlim(x_min - 1, x_max + 1)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.set_xlabel("log2Foldchange", fontsize=12, fontweight="bold")

    # Y-axis labels only for first subplot    
    ax.set_yticks(y_pos)
    ax.set_yticklabels(genes, fontsize=10)

    # Add legend slightly lower than top-right
    
    ax.legend(
        handles=legend_handles,
        title="Direction of regulation",
        loc="upper right",
        bbox_to_anchor=(1, 0.85),  # push legend a bit down
        frameon=False,
        fontsize=11,
        title_fontsize=12,
    )

plt.tight_layout()

# Save high-quality outputs
outfile_base = os.path.join(outdir, "combined_genes")
plt.savefig(outfile_base + ".pdf")  # vector for Nature submission
plt.savefig(outfile_base + ".png", dpi=600)  # hi-res raster fallback
plt.close()
print(f"Saved figures: {outfile_base}.pdf and .png")
