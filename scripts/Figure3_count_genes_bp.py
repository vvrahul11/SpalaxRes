# """The logic here is to understand which are the genes involved in top 20 BP terms in both species"""

import os
import textwrap
from collections import Counter

import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests


# ============================================================
# SETTINGS
# ============================================================

base_path = "C:/projects/Final Rat Trans Project Result2023/blindmolerat_rnaseqanalysis-2024/DMBA/Results_10092025_f2.5_p0.01_SELECTED/entrezID/"
os.chdir(base_path)

species_files = {
    "Mouse": "MM.control-M.wound/GeneOntologyBP.csv",
    "Spalax": "S.control-S.wound/GeneOntologyBP.csv",
}

top_n_bp = 20
top_n_genes = 10
fdr_cutoff = 0.05
min_bp_degree = 1

output_file = "Mouse_vs_Spalax_clean_GO_BP_network_full_terms_600dpi.png"


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def wrap_label(label, width=28):
    return "\n".join(textwrap.wrap(str(label), width=width))


def load_go_table(file_path, top_n_bp=20):
    df = pd.read_csv(file_path)
    df = df.iloc[:top_n_bp].copy()

    df["gene_list"] = df["geneID"].astype(str).str.split("/")
    df["label"] = df["Description"].apply(lambda x: wrap_label(x, 28))

    return df


def recurrent_gene_test(df):
    all_genes = [g for genes in df["gene_list"] for g in genes]
    gene_counts = Counter(all_genes)

    gene_count_df = pd.DataFrame.from_dict(
        gene_counts,
        orient="index",
        columns=["count"]
    )

    gene_count_df.index.name = "GeneSymbol"
    gene_count_df = gene_count_df.sort_values("count", ascending=False)

    total_gene_appearances = len(all_genes)
    unique_genes = len(gene_count_df)
    expected_prob = 1 / unique_genes

    gene_count_df["p_value"] = gene_count_df["count"].apply(
        lambda x: binomtest(
            k=x,
            n=total_gene_appearances,
            p=expected_prob,
            alternative="greater"
        ).pvalue
    )

    gene_count_df["p_adj"] = multipletests(
        gene_count_df["p_value"],
        method="fdr_bh"
    )[1]

    gene_count_df = gene_count_df.sort_values(
        ["p_adj", "count"],
        ascending=[True, False]
    )

    return gene_count_df


def build_network(
    df,
    gene_stats,
    species,
    top_n_genes=7,
    fdr_cutoff=0.05,
    min_bp_degree=2
):

    sig_genes = gene_stats[gene_stats["p_adj"] < fdr_cutoff].copy()

    if sig_genes.empty:
        sig_genes = gene_stats.head(top_n_genes).copy()
    else:
        sig_genes = sig_genes.head(top_n_genes).copy()

    selected_genes = sig_genes.index.tolist()

    G = nx.Graph()

    for gene in selected_genes:
        count = gene_stats.loc[gene, "count"]
        padj = gene_stats.loc[gene, "p_adj"]

        G.add_node(
            f"{species}_GENE_{gene}",
            label=gene,
            node_type="gene",
            species=species,
            count=count,
            p_adj=padj,
            size=600 + count * 160,
        )

    for _, row in df.iterrows():
        bp = row["Description"]
        bp_label = row["label"]
        genes = row["gene_list"]

        matching_genes = [g for g in genes if g in selected_genes]

        if len(matching_genes) == 0:
            continue

        bp_node = f"{species}_BP_{bp}"

        G.add_node(
            bp_node,
            label=bp_label,
            node_type="bp",
            species=species,
            size=850,
        )

        for gene in matching_genes:
            gene_node = f"{species}_GENE_{gene}"
            G.add_edge(gene_node, bp_node)

    bp_nodes_to_remove = [
        n for n, d in G.nodes(data=True)
        if d["node_type"] == "bp" and G.degree(n) < min_bp_degree
    ]

    G.remove_nodes_from(bp_nodes_to_remove)

    isolated_genes = [
        n for n, d in G.nodes(data=True)
        if d["node_type"] == "gene" and G.degree(n) == 0
    ]

    G.remove_nodes_from(isolated_genes)

    selected_genes_final = [
        G.nodes[n]["label"]
        for n, d in G.nodes(data=True)
        if d["node_type"] == "gene"
    ]

    return G, selected_genes_final


def plot_clean_networks(graphs, selected_gene_dict, output_file=None):

    fig, axes = plt.subplots(1, 2, figsize=(28, 14))

    species_colors = {
        "Mouse": {
            "gene": "#F4A3A3",
            "gene_edge": "#B22222",
            "bp": "#D7EAFB",
            "bp_edge": "#4A90C2",
            "title": "#1F5A99",
        },
        "Spalax": {
            "gene": "#F4A3A3",
            "gene_edge": "#B22222",
            "bp": "#DFF0D8",
            "bp_edge": "#4F8A3D",
            "title": "#3C7A35",
        },
    }

    for ax, (species, G) in zip(axes, graphs.items()):

        colors = species_colors[species]

        pos = nx.spring_layout(
            G,
            seed=42,
            k=1.6,
            iterations=500
        )

        gene_nodes = [
            n for n, d in G.nodes(data=True)
            if d["node_type"] == "gene"
        ]

        bp_nodes = [
            n for n, d in G.nodes(data=True)
            if d["node_type"] == "bp"
        ]

        gene_sizes = [G.nodes[n]["size"] for n in gene_nodes]
        bp_sizes = [G.nodes[n]["size"] for n in bp_nodes]

        nx.draw_networkx_edges(
            G,
            pos,
            ax=ax,
            edge_color="gray",
            width=1.0,
            alpha=0.28
        )

        nx.draw_networkx_nodes(
            G,
            pos,
            nodelist=bp_nodes,
            node_size=bp_sizes,
            node_color=colors["bp"],
            edgecolors=colors["bp_edge"],
            linewidths=1.0,
            ax=ax
        )

        nx.draw_networkx_nodes(
            G,
            pos,
            nodelist=gene_nodes,
            node_size=gene_sizes,
            node_color=colors["gene"],
            edgecolors=colors["gene_edge"],
            linewidths=1.3,
            ax=ax
        )

        bp_labels = {
            n: G.nodes[n]["label"]
            for n in bp_nodes
        }

        gene_labels = {
            n: G.nodes[n]["label"]
            for n in gene_nodes
        }

        nx.draw_networkx_labels(
            G,
            pos,
            labels=bp_labels,
            font_size=14,
            font_weight="normal",
            ax=ax
        )

        nx.draw_networkx_labels(
            G,
            pos,
            labels=gene_labels,
            font_size=14,
            font_weight="bold",
            ax=ax
        )

        # ax.set_title(
        #     ""
        #     # f"{species}\nTop recurrent genes: {', '.join(selected_gene_dict[species])}",
        #     fontsize=15,
        #     fontweight="bold",
        #     color=colors["title"],
        #     pad=18
        # )

        ax.axis("off")

    # fig.suptitle(
    #     # f"Recurrent genes across top {top_n_bp} GO biological processes in wound response",
    #     fontsize=18,
    #     fontweight="bold",
    #     y=0.98
    # )

    plt.tight_layout(rect=[0, 0, 1, 0.94])

    if output_file:
        plt.savefig(
            output_file,
            dpi=600,
            bbox_inches="tight",
            facecolor="white"
        )

    plt.show()


# ============================================================
# RUN ANALYSIS
# ============================================================

graphs = {}
selected_gene_dict = {}
gene_stats_dict = {}

for species, file_path in species_files.items():

    df = load_go_table(file_path, top_n_bp=top_n_bp)

    gene_stats = recurrent_gene_test(df)
    gene_stats_dict[species] = gene_stats

    print("\n" + "=" * 80)
    print(species)
    print("=" * 80)
    print(gene_stats.head(20))

    G, selected_genes = build_network(
        df=df,
        gene_stats=gene_stats,
        species=species,
        top_n_genes=top_n_genes,
        fdr_cutoff=fdr_cutoff,
        min_bp_degree=min_bp_degree
    )

    graphs[species] = G
    selected_gene_dict[species] = selected_genes


plot_clean_networks(
    graphs=graphs,
    selected_gene_dict=selected_gene_dict,
    output_file=output_file
)

print(f"Saved figure to: {os.path.abspath(output_file)}")