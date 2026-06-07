#!/usr/bin/env python3
"""
Figure 5 – Coupled inflammatory activation and feedback regulation
             in Spalax wound tissue

Run from the same directory as this file; data lives in ./data_for_figure/
Outputs: Figure5_conceptual_spalax_wound.{png,pdf,svg}
"""

import warnings
warnings.filterwarnings("ignore")

import re
import sys
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
import networkx as nx

# ─── Paths ────────────────────────────────────────────────────────────────────
DATA_DIR = Path(__file__).parent / "data_for_figure"
OUT_STEM = DATA_DIR / "Figure5_conceptual_spalax_wound"

# ─── Global color palette (consistent across all panels) ─────────────────────
COLORS = {
    "immune_activation":   "#B03A2E",   # deep red
    "feedback_regulation": "#1F618D",   # steel blue
    "tissue_repair":       "#1E8449",   # forest green
    "other":               "#7F8C8D",   # grey
}
CAT_LABELS = {
    "immune_activation":   "Immune activation",
    "feedback_regulation": "Feedback regulation",
    "tissue_repair":       "Tissue repair / remodeling",
    "other":               "Other",
}

# ─── Supplementary gene → Entrez-ID map for target genes not in KEGG file ───
# Primarily standard NCBI mouse (Mus musculus) Entrez IDs.
# Where the KEGG file provides an ID it takes priority (see load_expression_files).
SUPPLEMENTARY_GENE_MAP = {
    "Tnf":      "21926",   # Tumor necrosis factor
    "Nlrp3":    "216799",  # NLR family pyrin domain-containing 3
    "Nod2":     "257629",  # Nucleotide-binding oligomerization domain 2
    "Cxcl10":   "15945",   # C-X-C motif chemokine ligand 10 (IP-10)
    "Tlr2":     "24088",   # Toll-like receptor 2 (also Trem1 in pipeline – check)
    "Tlr8":     "170743",  # Toll-like receptor 8
    "Il1b":     "15891",   # Interleukin 1 beta (also from KEGG)
    # Feedback regulators
    "Nfkbia":   "18035",   # NF-kappaB inhibitor alpha
    "Socs3":    "12702",   # Suppressor of cytokine signaling 3
    "Zfp36":    "22695",   # ZFP36 ring finger protein (Tristetraprolin)
    "Acod1":    "107869",  # Aconitate decarboxylase 1 (formerly Irg1)
    "Anxa1":    "16952",   # Annexin A1
    "Cnr2":     "12802",   # Cannabinoid receptor 2
    "Adora2b":  "11540",   # Adenosine A2b receptor
    "Nr1d1":    "217166",  # Nuclear receptor subfamily 1 group D member 1
    "Nlrp12":   "268749",  # NLR family pyrin domain-containing 12
    # Repair / remodeling
    "Mmp8":     "17394",   # Matrix metallopeptidase 8
    "Mmp12":    "17381",   # Matrix metallopeptidase 12
    "Serpine1": "18787",   # Serpin family E member 1 (PAI-1)
    "Thbs1":    "21825",   # Thrombospondin 1
    "Vtn":      "22371",   # Vitronectin (22371 per NCBI; 22370=Vcan in pipeline)
}

# ═════════════════════════════════════════════════════════════════════════════
#  1. load_enrichment
# ═════════════════════════════════════════════════════════════════════════════
def load_enrichment(filepath: Path) -> pd.DataFrame:
    """Load a clusterProfiler enrichment CSV with auto-detected columns."""
    df = pd.read_csv(filepath, index_col=0)
    print(f"\n[load_enrichment] {filepath.name}  shape={df.shape}")
    print(f"  columns: {list(df.columns)}")

    # ── q-value ──────────────────────────────────────────────────────────────
    for col in ("qvalue", "q_value", "q.value", "p.adjust", "padj", "FDR"):
        if col in df.columns:
            df["_qval"] = pd.to_numeric(df[col], errors="coerce")
            print(f"  q-value col  → '{col}'")
            break
    else:
        raise ValueError(
            f"No q-value column in {filepath.name}. "
            f"Found: {list(df.columns)}"
        )

    # ── description ──────────────────────────────────────────────────────────
    for col in ("Description", "description", "term", "pathway", "Term"):
        if col in df.columns:
            df["_desc"] = df[col].astype(str)
            break
    else:
        df["_desc"] = df.index.astype(str)

    # Strip " - Mus musculus (house mouse)" from KEGG descriptions
    df["_desc"] = df["_desc"].str.replace(
        r"\s*-\s*Mus musculus.*$", "", regex=True
    ).str.strip()

    # ── gene list ─────────────────────────────────────────────────────────────
    for col in ("geneID", "gene_names", "genes", "Genes", "geneList"):
        if col in df.columns:
            df["_genes"] = df[col].astype(str)
            print(f"  gene-list col → '{col}'")
            break
    else:
        df["_genes"] = ""

    # ── count & ratio ─────────────────────────────────────────────────────────
    for col in ("Count", "count", "GeneCount"):
        if col in df.columns:
            df["_count"] = pd.to_numeric(df[col], errors="coerce")
            break
    else:
        df["_count"] = np.nan

    for col in ("GeneRatio", "generatio", "gene_ratio"):
        if col in df.columns:
            df["_ratio"] = df[col].apply(
                lambda x: eval(str(x)) if "/" in str(x) else float(x)
            )
            break
    else:
        df["_ratio"] = np.nan

    df["-log10q"] = -np.log10(df["_qval"].clip(lower=1e-300))
    return df


# ═════════════════════════════════════════════════════════════════════════════
#  2. parse_gene_list
# ═════════════════════════════════════════════════════════════════════════════
def parse_gene_list(gene_str: str) -> list:
    """Split a gene string on '/', ';', or ',' and return clean symbols."""
    if pd.isna(gene_str) or str(gene_str).strip() in ("", "nan"):
        return []
    return [
        g.strip()
        for g in re.split(r"[/;,]", str(gene_str))
        if g.strip() and g.strip().lower() != "nan"
    ]


# ═════════════════════════════════════════════════════════════════════════════
#  3. assign_go_category
# ═════════════════════════════════════════════════════════════════════════════
def assign_go_category(description: str) -> str:
    """Map a GO / KEGG description string to a biological category key."""
    d = description.lower()

    feedback_kw = [
        "negative regulation", "neg. regulation", "inhibition",
        "suppression", "feedback",
    ]
    repair_kw = [
        "wound healing", "tissue remodel", "extracellular matrix",
        "matrix organization", "proteolysis", "cell migration",
        "angiogenesis", "ecm-receptor", "coagulation", "arachidonic",
        "lipid and athero", "repair", "peptidase activity",
        "response to external stimulus",
    ]
    immune_kw = [
        "inflammatory", "immune", "cytokine", "leukocyte", "bacterium",
        "innate", "toll-like", "tlr", "nlrp", "nf-kb", "il-17", "tnf",
        "jak-stat", "chemokine", "interferon", "neutrophil", "monocyte",
        "macrophage", "activation", "defense", "pathogen",
        "cytokine receptor", "viral", "infection", "malaria",
        "staphylococcus", "complement", "neuroactive", "allograft",
        "rheumatoid", "diabetes", "hepatitis", "trypanosomiasis",
    ]

    if any(k in d for k in feedback_kw):
        return "feedback_regulation"
    if any(k in d for k in repair_kw):
        return "tissue_repair"
    if any(k in d for k in immune_kw):
        return "immune_activation"
    return "other"


# ═════════════════════════════════════════════════════════════════════════════
#  4. assign_gene_module
# ═════════════════════════════════════════════════════════════════════════════
def assign_gene_module(gene: str) -> str:
    """Return the functional module key for a gene symbol."""
    g = gene.lower()
    inflammatory = {
        "il6", "tnf", "ptgs2", "osm", "nlrp3", "s100a9", "lcn2", "csf3",
        "cxcl2", "cxcl3", "cxcl10", "cxcl1", "il1b", "il1a", "il1f10",
        "tlr2", "tlr8", "nod2", "nlrc4", "nlrc5",
    }
    feedback = {
        "tnfaip3", "nfkbia", "nfkbiz", "socs3", "bcl3", "traf1", "zfp36",
        "acod1", "nlrp12", "anxa1", "cnr2", "adora2b", "nr1d1",
    }
    repair = {
        "arg1", "mmp3", "mmp8", "mmp9", "mmp12", "mmp13", "serpine1",
        "thbs1", "apoe", "plaur", "vtn", "tgfb1",
    }
    if g in inflammatory:
        return "immune_activation"
    if g in feedback:
        return "feedback_regulation"
    if g in repair:
        return "tissue_repair"
    return "other"


# ═════════════════════════════════════════════════════════════════════════════
#  5. plot_go_enrichment  (Panel A)
# ═════════════════════════════════════════════════════════════════════════════
def plot_go_enrichment(ax: plt.Axes, go_df: pd.DataFrame) -> None:
    """
    Horizontal dot plot with three explicit section groups:
      Immune activation | Feedback regulation | Tissue repair / remodeling
    Point size = gene count; x-axis = –log10(q-value).
    """
    # ── Explicit term selections per section ─────────────────────────────────
    # Each list = partial patterns; first match wins (case-insensitive search)
    SECTIONS = {
        "immune_activation": [
            "inflammatory response",
            "leukocyte activation",
            "cytokine production",
            "innate immune response",
            "response to bacterium",
            "leukocyte migration",
            "neutrophil activation",
        ],
        "feedback_regulation": [
            "negative regulation of cytokine production",
            "negative regulation of inflammatory response",
            "negative regulation of defense response",
            "negative regulation of nf-kb",
            "negative regulation of interleukin-1 production",
            "negative regulation of leukocyte activation",
            "negative regulation of immune response",
            "negative regulation of toll-like receptor",
        ],
        "tissue_repair": [
            "wound healing",
            "extracellular matrix organization",
            "cell migration",
            "blood coagulation",
            "tissue remodeling",
            "regulation of peptidase activity",
        ],
    }

    desc_lower = go_df["_desc"].str.lower()

    def pick_terms(patterns, n_max=8):
        picked, seen = [], set()
        for pat in patterns:
            hits = desc_lower[desc_lower.str.contains(pat, regex=False)]
            # Among hits take the one with the best q-value
            hits_sorted = go_df.loc[hits.index].sort_values("_qval")
            for idx in hits_sorted.index:
                if idx not in seen and len(picked) < n_max:
                    picked.append(idx)
                    seen.add(idx)
                    break
        return picked

    blocks = []
    for cat_key, patterns in SECTIONS.items():
        idxs = pick_terms(patterns)
        sub = go_df.loc[idxs].copy()
        sub["category"] = cat_key
        sub["_section_order"] = list(SECTIONS.keys()).index(cat_key)
        blocks.append(sub)

    df_plot = pd.concat(blocks)
    # Sort: section_order DESCENDING so immune_activation (0) gets the highest
    # y positions and appears at the TOP of the plot (repair at bottom).
    df_plot = df_plot.sort_values(
        ["_section_order", "_qval"], ascending=[False, True]
    ).reset_index(drop=False)

    n = len(df_plot)
    y_pos = np.arange(n)
    node_colors = [COLORS[c] for c in df_plot["category"]]
    cnt = df_plot["_count"].fillna(10).astype(float)
    cnt_range = max(cnt.max() - cnt.min(), 1)
    sizes = 40 + 130 * (cnt - cnt.min()) / cnt_range

    ax.scatter(df_plot["-log10q"], y_pos, c=node_colors, s=sizes,
               zorder=3, edgecolors="white", linewidths=0.6, alpha=0.92)

    labels = [textwrap.fill(d, width=40) for d in df_plot["_desc"]]
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=7.2)
    ax.set_xlabel("–log₁₀(q-value)", fontsize=9)
    ax.set_title(
        "GO Biological Process enrichment\nin Spalax wound tissue",
        fontsize=9.5, fontweight="bold", pad=5,
    )

    # Section order in plot (bottom → top): repair, feedback, immune
    plot_section_order = ["tissue_repair", "feedback_regulation", "immune_activation"]
    section_sizes = [
        len(df_plot[df_plot["category"] == k]) for k in plot_section_order
    ]

    # Draw subtle horizontal section dividers
    boundary = 0
    for sz in section_sizes[:-1]:
        boundary += sz
        ax.axhline(boundary - 0.5, color="#D5D8DC", lw=1.0, ls="--", zorder=0)

    # Colored section bracket on the left
    sec_labels = {
        "immune_activation":   "Immune activation",
        "feedback_regulation": "Negative regulation /\nfeedback",
        "tissue_repair":       "Tissue repair /\nremodeling",
    }
    boundary = 0
    for cat_key, sz in zip(plot_section_order, section_sizes):
        mid = boundary + sz / 2 - 0.5
        ax.text(
            -0.28, mid, sec_labels[cat_key],
            transform=ax.get_yaxis_transform(),
            ha="center", va="center", fontsize=6.5,
            color=COLORS[cat_key], fontweight="bold",
            clip_on=False,
        )
        boundary += sz

    # Legend (size encodes gene count)
    patches = [
        mpatches.Patch(color=COLORS[k], label=CAT_LABELS[k])
        for k in ("immune_activation", "feedback_regulation", "tissue_repair")
    ]
    ax.legend(handles=patches, fontsize=6.5, loc="lower right",
              framealpha=0.80, handlelength=1.1, borderpad=0.5)

    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", lw=0.3, color="lightgrey", zorder=0)
    ax.set_xlim(left=0)


# ═════════════════════════════════════════════════════════════════════════════
#  6. plot_gene_go_network  (Panel B)
# ═════════════════════════════════════════════════════════════════════════════
def plot_gene_go_network(ax: plt.Axes, go_df: pd.DataFrame) -> None:
    """
    Tripartite Gene–GO network matching conceptual_figure2.jpg:
      LEFT  → immune-activation GO terms (6 terms, red dots)
      CENTER → gene nodes (colored capsules, readable external labels)
      RIGHT  → negative regulation GO terms (blue, 8) + repair GO terms (green, 5)
    """
    # ── Gene sets ─────────────────────────────────────────────────────────────
    TARGET_GENES = [
        # Inflammatory effectors (top group, red)
        "Tnf", "Il6", "Tlr2", "Tlr8", "Nlrp3", "Nod2", "S100a9",
        # Feedback regulators (middle group, blue)
        "Tnfaip3", "Zfp36", "Socs3", "Nlrp12", "Nr1d1",
        "Anxa1", "Adora2b", "Cnr2", "Acod1",
        # Repair / remodeling (bottom group, green)
        "Arg1", "Mmp12", "Serpine1", "Thbs1", "Apoe", "Plaur",
    ]
    target_lower = {g.lower(): g for g in TARGET_GENES}

    # ── Explicit GO term selection ─────────────────────────────────────────────
    # Left: immune activation GO terms
    LEFT_PATTERNS = [
        "inflammatory response",
        "leukocyte activation",
        "cytokine production",
        "innate immune response",
        "response to bacterium",
        "leukocyte migration",
    ]
    # Right-top: negative regulation (feedback) GO terms
    RIGHT_FEEDBACK_PATTERNS = [
        "negative regulation of cytokine production",
        "negative regulation of nf-kb",
        "negative regulation of interleukin-1 production",
        "negative regulation of inflammatory response",
        "negative regulation of defense response",
        "negative regulation of leukocyte activation",
        "negative regulation of immune response",
        "negative regulation of toll-like receptor",
    ]
    # Right-bottom: tissue repair GO terms
    RIGHT_REPAIR_PATTERNS = [
        "wound healing",
        "extracellular matrix organization",
        "cell migration",
        "blood coagulation",
        "tissue remodeling",
    ]

    go_df = go_df.copy()
    go_df["category"] = go_df["_desc"].apply(assign_go_category)
    desc_lower = go_df["_desc"].str.lower()

    def pick_one(patterns):
        """Return one GO term index per pattern (best q-value match)."""
        result, seen = [], set()
        for pat in patterns:
            hits = go_df.loc[desc_lower.str.contains(pat, regex=False)].sort_values("_qval")
            for idx in hits.index:
                if idx not in seen:
                    result.append(idx)
                    seen.add(idx)
                    break
        return result

    left_idxs         = pick_one(LEFT_PATTERNS)
    right_fb_idxs     = pick_one(RIGHT_FEEDBACK_PATTERNS)
    right_repair_idxs = pick_one(RIGHT_REPAIR_PATTERNS)

    all_idxs = left_idxs + right_fb_idxs + right_repair_idxs
    selected_go = go_df.loc[all_idxs]

    # ── Build bipartite graph ──────────────────────────────────────────────────
    G = nx.Graph()
    go_label_cat = {}   # short label → category string

    for idx in all_idxs:
        row = go_df.loc[idx]
        full_desc = row["_desc"]
        # Shorten label to ≤34 chars for display
        label = full_desc if len(full_desc) <= 34 else full_desc[:32] + "…"
        cat = row["category"]
        if idx in right_fb_idxs:
            cat = "feedback_regulation"   # enforce even if pattern mis-classified
        elif idx in right_repair_idxs:
            cat = "tissue_repair"
        elif idx in left_idxs:
            cat = "immune_activation"
        G.add_node(label, ntype="go", category=cat, full=full_desc)
        go_label_cat[label] = cat

        for g in parse_gene_list(str(row["_genes"])):
            gl = g.lower()
            if gl in target_lower:
                canonical = target_lower[gl]
                mod = assign_gene_module(canonical)
                G.add_node(canonical, ntype="gene", module=mod)
                G.add_edge(label, canonical)

    go_nodes   = [n for n in G.nodes if G.nodes[n]["ntype"] == "go"]
    gene_nodes = [n for n in G.nodes if G.nodes[n]["ntype"] == "gene"]

    left_go         = [n for n in go_nodes if G.nodes[n]["category"] == "immune_activation"]
    right_fb_go     = [n for n in go_nodes if G.nodes[n]["category"] == "feedback_regulation"]
    right_repair_go = [n for n in go_nodes if G.nodes[n]["category"] == "tissue_repair"]

    # ── Layout constants ───────────────────────────────────────────────────────
    X_LEFT, X_GENE, X_RIGHT = 0.0, 0.5, 1.0
    Y_TOP, Y_BOT = 1.0, 0.0

    pos = {}

    def place_col(nodes, x, yhi=Y_TOP, ylo=Y_BOT):
        n = len(nodes)
        if n == 0:
            return
        ys = np.linspace(yhi, ylo, n) if n > 1 else [(yhi + ylo) / 2]
        for nd, y in zip(nodes, ys):
            pos[nd] = (x, float(y))

    place_col(left_go, X_LEFT)

    # Right column: feedback in upper half, repair in lower half
    n_fb = len(right_fb_go)
    n_rp = len(right_repair_go)
    split_y = 0.38 if n_rp > 0 else Y_BOT
    place_col(right_fb_go,     X_RIGHT, yhi=Y_TOP,   ylo=split_y + 0.06)
    place_col(right_repair_go, X_RIGHT, yhi=split_y, ylo=Y_BOT)

    # Genes sorted by module
    mod_order = {"immune_activation": 0, "feedback_regulation": 1,
                 "tissue_repair": 2, "other": 3}
    gene_sorted = sorted(
        gene_nodes,
        key=lambda g: (mod_order.get(G.nodes[g]["module"], 3), g)
    )
    place_col(gene_sorted, X_GENE)

    # ── Draw edges ─────────────────────────────────────────────────────────────
    for u, v in G.edges():
        go_nd = u if G.nodes[u]["ntype"] == "go" else v
        cat = G.nodes[go_nd]["category"]
        c = COLORS.get(cat, COLORS["other"])
        if go_nd in pos and (v if u == go_nd else u) in pos:
            x0, y0 = pos[go_nd]
            x1, y1 = pos[v if u == go_nd else u]
            ax.plot([x0, x1], [y0, y1], color=c, lw=0.7, alpha=0.30, zorder=1)

    # ── Draw GO term nodes (filled circles, larger) ───────────────────────────
    for nd in go_nodes:
        if nd not in pos:
            continue
        x, y = pos[nd]
        cat = G.nodes[nd]["category"]
        ax.scatter(x, y, s=180, marker="o", color=COLORS.get(cat, COLORS["other"]),
                   zorder=3, edgecolors="white", linewidths=0.7)
        # Label: text to the outside of the dot
        ha  = "right" if x < X_GENE else "left"
        off = -0.045 if x < X_GENE else 0.045
        ax.text(x + off, y, nd, ha=ha, va="center",
                fontsize=6.2, color="#1C2833", clip_on=False)

    # ── Draw gene nodes with EXTERNAL labels ──────────────────────────────────
    for nd in gene_sorted:
        if nd not in pos:
            continue
        x, y = pos[nd]
        mod  = G.nodes[nd].get("module", "other")
        col  = COLORS.get(mod, COLORS["other"])
        # Pill-shaped node via FancyBboxPatch in data coordinates
        w, h = 0.085, 0.028
        rect = mpatches.FancyBboxPatch(
            (x - w / 2, y - h / 2), w, h,
            boxstyle="round,pad=0.005",
            facecolor=col, edgecolor="white", linewidth=0.5,
            transform=ax.transData, clip_on=False, zorder=4,
        )
        ax.add_patch(rect)
        ax.text(x, y, nd, ha="center", va="center",
                fontsize=6.5, fontweight="bold", color="white",
                fontstyle="italic", zorder=5, clip_on=False)

    # ── Section labels on right column (horizontal, placed above each group) ──
    if right_fb_go:
        top_fb = max(pos[n][1] for n in right_fb_go)
        ax.text(X_RIGHT, top_fb + 0.06,
                "",
                #"Negative regulation /\nfeedback (GO terms)",
                ha="center", va="bottom", fontsize=6.2,
                fontweight="bold", color=COLORS["feedback_regulation"],
                clip_on=False)
    if right_repair_go:
        top_rp = max(pos[n][1] for n in right_repair_go)
        ax.text(X_RIGHT, top_rp + 0.06,
                "",
                #"Tissue remodeling /\nprotection (GO terms)",
                ha="center", va="bottom", fontsize=6.2,
                fontweight="bold", color=COLORS["tissue_repair"],
                clip_on=False)

    # ── Column headers ─────────────────────────────────────────────────────────
    ax.text(X_LEFT,  Y_TOP + 0.10, 
            "",
            #"Immune activation\n(GO terms)",
            ha="center", va="bottom", fontsize=7, fontweight="bold",
            color=COLORS["immune_activation"], clip_on=False)
    ax.text(X_GENE,  Y_TOP + 0.10, "Core driver genes",
            ha="center", va="bottom", fontsize=7, fontweight="bold",
            color="#2C3E50", clip_on=False)

    ax.set_xlim(-0.28, 1.28)
    ax.set_ylim(-0.06, 1.15)
    ax.axis("off")
    ax.set_title(
        "Gene–GO network: core driver genes connect\n"
        "activation, feedback and remodeling programs",
        fontsize=9.5, fontweight="bold", pad=2,
    )

    patches = [
        mpatches.Patch(color=COLORS[k], label=CAT_LABELS[k])
        for k in ("immune_activation", "feedback_regulation", "tissue_repair")
    ]
    ax.legend(handles=patches, fontsize=6.5, loc="upper right",
              framealpha=0.80, handlelength=1.1, borderpad=0.5)


# ═════════════════════════════════════════════════════════════════════════════
#  7. load_expression_files
# ═════════════════════════════════════════════════════════════════════════════
def load_expression_files(
    data_dir: Path, kegg_df: pd.DataFrame
):
    """
    Load four DEG files and return:
      expr_dfs  – {condition_label: DataFrame with '_lfc' column}
      sym_to_id – {gene_symbol: entrez_id_str} merged from KEGG + supplement
    """
    FILE_MAP = {
        "Mouse scab":   "MM.control-M.scar.csv",
        "Mouse wound":  "MM.control-M.wound.csv",
        "Spalax scab":  "S.control-S.scar.csv",
        "Spalax wound": "S.control-S.wound.csv",
    }

    expr_dfs = {}
    for label, fname in FILE_MAP.items():
        fp = data_dir / fname
        df = pd.read_csv(fp, index_col=0)
        print(f"\n[load_expression] {fname}  shape={df.shape}  cols={list(df.columns)}")
        # Detect LFC column
        for col in ("log2FoldChange", "logFC", "avg_log2FC", "Log2FC",
                    "log2FC", "LFC", "FoldChange"):
            if col in df.columns:
                df["_lfc"] = pd.to_numeric(df[col], errors="coerce")
                print(f"  LFC col → '{col}'")
                break
        else:
            df["_lfc"] = pd.to_numeric(df.iloc[:, 0], errors="coerce")
            print("  LFC col → first column (fallback)")
        # Ensure string index for consistent look-up
        df.index = df.index.astype(str)
        expr_dfs[label] = df

    # ── Build gene-symbol → Entrez-ID map ────────────────────────────────────
    # Step 1: from KEGG file (geneID numeric + gene_names symbol)
    sym_to_id: dict[str, str] = {}
    if "gene_names" in kegg_df.columns and "geneID" in kegg_df.columns:
        for _, row in kegg_df.iterrows():
            ids   = parse_gene_list(str(row.get("geneID", "")))
            names = parse_gene_list(str(row.get("gene_names", "")))
            for gid, gname in zip(ids, names):
                if "unknown" not in gname.lower():
                    sym_to_id[gname] = str(gid)

    # Step 2: supplement with hardcoded NCBI IDs for missing targets
    for sym, eid in SUPPLEMENTARY_GENE_MAP.items():
        if sym not in sym_to_id:
            sym_to_id[sym] = eid

    print(f"\n  Gene symbol map: {len(sym_to_id)} entries "
          f"({sum(1 for v in sym_to_id.values() if v)} non-empty)")
    return expr_dfs, sym_to_id


# ═════════════════════════════════════════════════════════════════════════════
#  8. plot_gene_heatmap  (Panel C)
# ═════════════════════════════════════════════════════════════════════════════
def plot_gene_heatmap(
    ax: plt.Axes, expr_dfs: dict, sym_to_id: dict
) -> None:
    """
    Z-scored log2FC heatmap of three gene modules across four conditions.
    Conditions ordered: Mouse scab | Mouse wound | Spalax scab | Spalax wound.
    """
    MODULES = {
        "Inflammatory\neffectors": [
            "Il6", "Tnf", "Ptgs2", "Osm", "Nlrp3", "S100a9", "Lcn2",
            "Csf3", "Cxcl2", "Cxcl3", "Cxcl10", "Tlr2", "Tlr8", "Nod2",
        ],
        "Feedback\nregulators": [
            "Tnfaip3", "Nfkbia", "Socs3", "Bcl3", "Traf1", "Zfp36",
            "Acod1", "Nlrp12", "Anxa1", "Cnr2", "Adora2b", "Nr1d1",
        ],
        "Repair /\nremodeling": [
            "Arg1", "Mmp3", "Mmp8", "Mmp9", "Mmp12", "Mmp13",
            "Serpine1", "Thbs1", "Apoe", "Plaur", "Vtn",
        ],
    }
    COL_ORDER = ["Mouse scab", "Mouse wound", "Spalax scab", "Spalax wound"]
    MOD_CAT = {
        "Inflammatory\neffectors": "immune_activation",
        "Feedback\nregulators":    "feedback_regulation",
        "Repair /\nremodeling":    "tissue_repair",
    }

    rows, row_labels, row_modules = [], [], []

    for mod_name, genes in MODULES.items():
        for gene in genes:
            gene_id = sym_to_id.get(gene)
            vals = []
            for cond in COL_ORDER:
                df = expr_dfs[cond]
                v = np.nan
                if gene_id and gene_id in df.index:
                    v = df.at[gene_id, "_lfc"]
                vals.append(v)

            n_found = sum(1 for v in vals if not np.isnan(v))
            if n_found >= 1:           # include if found in ≥1 condition
                rows.append(vals)
                row_labels.append(gene)
                row_modules.append(mod_name)
            else:
                print(f"  WARNING: {gene} (ID={gene_id}) – not found in any "
                      f"expression file, skipping")

    if not rows:
        ax.text(0.5, 0.5, "No genes mapped to expression files.\n"
                "Check gene-ID mapping.", ha="center", va="center",
                transform=ax.transAxes, fontsize=9, color="red")
        return

    mat = pd.DataFrame(rows, index=row_labels, columns=COL_ORDER)
    print(f"\n  Heatmap: {mat.shape[0]} genes × {mat.shape[1]} conditions")
    print(f"  Missing values: {mat.isna().sum().sum()} / {mat.size}")

    # Z-score each row (gene) across conditions; fill NaN with 0 after z-score
    def zscore_row(r):
        valid = r.dropna()
        if len(valid) < 2:
            return r.fillna(0)
        return (r - valid.mean()) / (valid.std() + 1e-9)

    mat_z = mat.apply(zscore_row, axis=1).fillna(0)

    # Draw heatmap
    cmap = sns.diverging_palette(240, 10, s=85, l=45, as_cmap=True)
    hm = sns.heatmap(
        mat_z, ax=ax, cmap=cmap, center=0, vmin=-2.5, vmax=2.5,
        cbar_kws={"shrink": 0.55, "label": "Z-score (log₂FC)",
                  "ticks": [-2, -1, 0, 1, 2]},
        xticklabels=True, yticklabels=True,
        linewidths=0.25, linecolor="white",
    )

    ax.set_xticklabels(
        ax.get_xticklabels(), rotation=30, ha="right", fontsize=8
    )
    ax.set_yticklabels(
        ax.get_yticklabels(), fontstyle="italic", fontsize=6.2, va="center"
    )
    ax.set_title(
        "Heatmap of key gene modules across conditions",
        fontsize=9, fontweight="bold", pad=5,
    )

    # ── Colored module side bars ───────────────────────────────────────────────
    bar_x  = -1.5    # axes-data units to the left of y-axis
    bar_w  = 0.6
    y_start = 0
    for mod_name in MODULES:
        n = sum(1 for m in row_modules if m == mod_name)
        if n == 0:
            continue
        cat = MOD_CAT[mod_name]
        rect = mpatches.FancyBboxPatch(
            (bar_x, y_start), bar_w, n,
            boxstyle="round,pad=0.05",
            facecolor=COLORS[cat], edgecolor="none",
            transform=ax.transData, clip_on=False, zorder=5,
        )
        ax.add_patch(rect)
        ax.text(
            bar_x + bar_w / 2, y_start + n / 2, mod_name,
            ha="center", va="center", fontsize=5.8,
            color="white", fontweight="bold", rotation=90,
            transform=ax.transData, clip_on=False, zorder=6,
        )
        y_start += n


# ═════════════════════════════════════════════════════════════════════════════
#  9. plot_kegg_enrichment  (Panel D)
# ═════════════════════════════════════════════════════════════════════════════
def plot_kegg_enrichment(ax: plt.Axes, kegg_df: pd.DataFrame) -> None:
    """
    Horizontal bar chart of top KEGG pathways.
    x = –log10(q-value), bars colored by biological category, labeled by n.
    """
    PRIORITY = [
        "cytokine-cytokine receptor",
        "il-17 signaling",
        "tnf signaling",
        "jak-stat",
        "pi3k-akt",
        "leukocyte transendothelial",
        "neutrophil extracellular",
        "complement and coagulation",
        "arachidonic acid",
        "ecm-receptor",
        "viral protein interaction",
        "neuroactive ligand",
        "staphylococcus",
        "malaria",
        "lipid and athero",
    ]

    desc_lower = kegg_df["_desc"].str.lower()
    selected, used = [], set()
    for pat in PRIORITY:
        for idx in desc_lower[desc_lower.str.contains(pat, regex=False)].index:
            if idx not in used:
                selected.append(idx)
                used.add(idx)
    for idx in kegg_df.sort_values("_qval").index:
        if len(selected) >= 15:
            break
        if idx not in used:
            selected.append(idx)
            used.add(idx)

    sub = kegg_df.loc[selected].copy()
    sub["category"] = sub["_desc"].apply(assign_go_category)
    sub = sub.sort_values("-log10q", ascending=True)   # bottom→top

    y_pos = np.arange(len(sub))
    bar_colors = [COLORS.get(c, COLORS["other"]) for c in sub["category"]]

    bars = ax.barh(y_pos, sub["-log10q"], color=bar_colors,
                   edgecolor="white", height=0.72, alpha=0.88)

    x_max = sub["-log10q"].max()
    for bar, cnt in zip(bars, sub["_count"].fillna(0).astype(int)):
        ax.text(
            bar.get_width() + x_max * 0.015,
            bar.get_y() + bar.get_height() / 2,
            f"n={cnt}", va="center", fontsize=6.5,
        )

    # Wrap long labels
    labels = [textwrap.fill(d, width=40) for d in sub["_desc"]]
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("–log₁₀(q-value)", fontsize=8.5)
    ax.set_title(
        "KEGG pathway enrichment\n(top pathways, Spalax wound)",
        fontsize=9, fontweight="bold", pad=5,
    )

    patches = [
        mpatches.Patch(color=COLORS[k], label=CAT_LABELS[k])
        for k in ("immune_activation", "feedback_regulation", "tissue_repair")
    ]
    ax.legend(handles=patches, fontsize=6.5, loc="lower right",
              framealpha=0.75, handlelength=1)

    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", lw=0.3, color="lightgrey", zorder=0)
    ax.set_xlim(left=0)


# ═════════════════════════════════════════════════════════════════════════════
#  10. main
# ═════════════════════════════════════════════════════════════════════════════
def main() -> None:
    matplotlib.rcParams.update({
        "font.family":  "DejaVu Sans",
        "font.size":    8,
        "axes.linewidth": 0.7,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
    })

    # ── Load data ─────────────────────────────────────────────────────────────
    go_df   = load_enrichment(DATA_DIR / "GeneOntologyBP_Spalax_wound.csv")
    kegg_df = load_enrichment(DATA_DIR / "Kegg_pathway_Spalax_wound.csv")
    expr_dfs, sym_to_id = load_expression_files(DATA_DIR, kegg_df)

    # ── Figure skeleton ───────────────────────────────────────────────────────
    fig = plt.figure(figsize=(18, 14))
    fig.suptitle(
        "Coupled inflammatory activation and feedback regulation\n"
        "in Spalax wound tissue",
        fontsize=13, fontweight="bold", y=0.99,
    )

    gs = gridspec.GridSpec(
        2, 2, figure=fig,
        hspace=0.40, wspace=0.38,
        left=0.10, right=0.97, top=0.93, bottom=0.06,
    )
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, 0])
    axD = fig.add_subplot(gs[1, 1])

    # ── Draw panels ───────────────────────────────────────────────────────────
    print("\n── Panel A: GO enrichment ──────────────────────────────────────")
    plot_go_enrichment(axA, go_df)

    print("\n── Panel B: Gene–GO network ────────────────────────────────────")
    plot_gene_go_network(axB, go_df)

    print("\n── Panel C: Expression heatmap ─────────────────────────────────")
    plot_gene_heatmap(axC, expr_dfs, sym_to_id)

    print("\n── Panel D: KEGG enrichment ─────────────────────────────────────")
    plot_kegg_enrichment(axD, kegg_df)

    # ── Bold panel labels ─────────────────────────────────────────────────────
    for ax, label in [(axA, "A"), (axB, "B"), (axC, "C"), (axD, "D")]:
        ax.text(
            -0.13, 1.05, label,
            transform=ax.transAxes,
            fontsize=18, fontweight="bold", va="top",
        )

    # ── Save ──────────────────────────────────────────────────────────────────
    for ext, dpi in [("png", 600), ("pdf", None), ("svg", None)]:
        out = Path(str(OUT_STEM) + f".{ext}")
        kwargs = dict(bbox_inches="tight", format=ext)
        if dpi:
            kwargs["dpi"] = dpi
        fig.savefig(out, **kwargs)
        print(f"\nSaved → {out}")

    plt.show()
    print("\nDone.")


if __name__ == "__main__":
    main()
