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
DATA_DIR = Path('C:/projects/Final Rat Trans Project Result2023/blindmolerat_rnaseqanalysis-2024/data_for_figure')
OUT_STEM = DATA_DIR / "Figure5_ABC"

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

# NOTE: the Entrez-ID-based SUPPLEMENTARY_GENE_MAP that used to live here
# has been removed. The expression files are indexed by gene symbol
# ('genename'), not Entrez ID, so load_expression_files() now looks genes
# up directly by (lowercased) symbol. See load_expression_files() below.


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
#  4. Gene → functional-module assignment
# ═════════════════════════════════════════════════════════════════════════════
# NOT a keyword heuristic. These 52 genes are the statistically-significant
# hub genes from an exact Poisson-binomial recurrence test (genes appearing
# in significantly more enriched GO Biological Process terms than expected
# by chance, given each term's size), run on a curated subset of the 602
# significant Spalax-wound GO BP terms restricted to immune-activation,
# feedback-regulation, and tissue-repair-relevant biology; genes were then
# assigned to whichever of those 3 categories their own significant-term
# memberships were proportionally most enriched in (Benjamini–Hochberg
# FDR < 0.05 throughout).
#
# GENE_SECONDARY flags genuine dual/"moonlighting" roles: either the gene's
# own data showed strong representation (>=15% of that category's term
# pool) in a second category, or the dual role is well established in the
# literature but under-represented by GO annotation wording in this dataset
# (flagged below as such). A manual override was applied to Tnf, whose raw
# category proportions were a near-three-way tie (49.0% repair / 42.3%
# immune / 39.0% feedback) — reassigned to its canonical immune-activation
# role, with tissue_repair kept as a documented secondary role (TNF drives
# MMP induction/ECM remodeling).
GENE_PRIMARY = {
    # feedback regulation
    "Apoa1": "feedback_regulation", "Arg1": "feedback_regulation",
    "Clec4g": "feedback_regulation", "Fas": "feedback_regulation",
    "Mefv": "feedback_regulation", "Pglyrp1": "feedback_regulation",
    "Ripor2": "feedback_regulation", "Tarm1": "feedback_regulation",
    "Tnfaip3": "feedback_regulation",
    # immune activation
    "Acod1": "immune_activation", "Adam8": "immune_activation",
    "Bcl3": "immune_activation", "Ccr5": "immune_activation",
    "Cd86": "immune_activation", "Clcf1": "immune_activation",
    "Clec7a": "immune_activation", "Coro1a": "immune_activation",
    "Il17ra": "immune_activation", "Il18r1": "immune_activation",
    "Il1a": "immune_activation", "Il1r1": "immune_activation",
    "Il6": "immune_activation", "Lep": "immune_activation",
    "Ltf": "immune_activation", "Nfkbid": "immune_activation",
    "Nfkbiz": "immune_activation", "Nlrp3": "immune_activation",
    "Nod2": "immune_activation", "Nr4a3": "immune_activation",
    "Ptafr": "immune_activation", "Rsad2": "immune_activation",
    "Runx1": "immune_activation", "Slc11a1": "immune_activation",
    "Tlr2": "immune_activation", "Tlr8": "immune_activation",
    "Tnf": "immune_activation",  # manual override, see note above
    "Tslp": "immune_activation", "Vav1": "immune_activation",
    # tissue repair / remodeling
    "Adora2b": "tissue_repair", "Adtrp": "tissue_repair",
    "Anxa1": "tissue_repair", "Apoe": "tissue_repair",
    "Grem1": "tissue_repair", "Hif1a": "tissue_repair",
    "Itgb2": "tissue_repair", "Mmp9": "tissue_repair",
    "Nos3": "tissue_repair", "Ptgs2": "tissue_repair",
    "S100a9": "tissue_repair", "Serpine1": "tissue_repair",
    "Thbs1": "tissue_repair", "Vtn": "tissue_repair",
}

GENE_SECONDARY = {
    "Adam8": "tissue_repair", "Adtrp": "feedback_regulation",
    "Anxa1": "immune_activation", "Apoe": "feedback_regulation",
    "Arg1": "immune_activation",
    "Bcl3": "feedback_regulation",       # literature override (atypical IkB)
    "Grem1": "feedback_regulation", "Il1a": "tissue_repair",
    "Il6": "tissue_repair", "Itgb2": "immune_activation",
    "Nfkbid": "feedback_regulation",
    "Nfkbiz": "feedback_regulation",     # literature override (IkBzeta)
    "Ptgs2": "feedback_regulation", "Serpine1": "feedback_regulation",
    "Thbs1": "feedback_regulation",
    "Tnf": "tissue_repair",              # manual override, see note above
    "Tnfaip3": "immune_activation", "Vtn": "feedback_regulation",
}


def assign_gene_module(gene: str) -> str:
    """Primary functional module for a gene symbol (falls back to 'other')."""
    return GENE_PRIMARY.get(gene, "other")


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
    ROW_SPACING = 1.7   # extra vertical room so wrapped multi-line labels don't collide
    y_pos = np.arange(n) * ROW_SPACING
    node_colors = [COLORS[c] for c in df_plot["category"]]
    cnt = df_plot["_count"].fillna(10).astype(float)
    cnt_range = max(cnt.max() - cnt.min(), 1)
    sizes = 40 + 130 * (cnt - cnt.min()) / cnt_range

    ax.scatter(df_plot["-log10q"], y_pos, c=node_colors, s=sizes,
               zorder=3, edgecolors="white", linewidths=0.6, alpha=0.92)

    labels = [textwrap.fill(d, width=30) for d in df_plot["_desc"]]
    ax.set_yticks(y_pos)
    # ax.set_yticklabels(labels, fontsize=7.2)
    yticklabels = ax.set_yticklabels(
    labels,
    fontsize=6.8
    )

    for tick, cat in zip(yticklabels, df_plot["category"]):
        tick.set_color(COLORS[cat])
        tick.set_fontweight("bold")

    ax.set_xlabel("–log₁₀(q-value)", fontsize=9)

    # Section order in plot (bottom → top): repair, feedback, immune
    plot_section_order = ["tissue_repair", "feedback_regulation", "immune_activation"]
    section_sizes = [
        len(df_plot[df_plot["category"] == k]) for k in plot_section_order
    ]

    # Draw subtle horizontal section dividers
    boundary = 0
    for sz in section_sizes[:-1]:
        boundary += sz
        ax.axhline(boundary * ROW_SPACING - ROW_SPACING / 2,
                   color="#D5D8DC", lw=1.0, ls="--", zorder=0)

    # Colored section bracket on the left
    # sec_labels = {
    #     "immune_activation":   "Immune activation",
    #     "feedback_regulation": "Negative regulation /\nfeedback",
    #     "tissue_repair":       "Tissue repair /\nremodeling",
    # }
    # boundary = 0
    # for cat_key, sz in zip(plot_section_order, section_sizes):
    #     mid = boundary + sz / 2 - 0.5
    #     ax.text(
    #         -0.28, mid, sec_labels[cat_key],
    #         transform=ax.get_yaxis_transform(),
    #         ha="center", va="center", fontsize=6.5,
    #         color=COLORS[cat_key], fontweight="bold",
    #         clip_on=False,
    #     )
    #     boundary += sz

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
    # ── Gene set (Panel B) ──────────────────────────────────────────────────
    # A readable subset of the full 52-gene, statistically-significant hub-gene
    # set (see GENE_PRIMARY/GENE_SECONDARY above) — all 9 feedback_regulation
    # genes are kept (small category, all informative), plus the top
    # immune_activation and tissue_repair genes by recurrence significance.
    # The full 52-gene set is shown in Panel C.
    TARGET_GENES = [
        # immune activation
        "Il6", "Tnf", "Nod2", "Nlrp3", "Tlr2", "Il1a", "Bcl3", "Nfkbiz", "Acod1",
        # feedback regulation (all)
        "Tnfaip3", "Arg1", "Ripor2", "Clec4g", "Apoa1", "Fas", "Mefv", "Pglyrp1", "Tarm1",
        # tissue repair / remodeling
        "Anxa1", "S100a9", "Serpine1", "Thbs1", "Mmp9", "Vtn",
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
        # Panel B is now much wider (2x the width of Panel A), so the full
        # GO term description fits without truncation for all but the most
        # extreme outliers; keep a generous safety cap only for those.
        label = full_desc if len(full_desc) <= 60 else full_desc[:58] + "…"
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
        off = -0.05 if x < X_GENE else 0.05
        ax.text(x + off, y, nd, ha=ha, va="center",
                fontsize=6.8, color="#1C2833", clip_on=False)

    # ── Draw gene nodes with EXTERNAL labels ──────────────────────────────────
    for nd in gene_sorted:
        if nd not in pos:
            continue
        x, y = pos[nd]
        mod  = G.nodes[nd].get("module", "other")
        col  = COLORS.get(mod, COLORS["other"])
        sec  = GENE_SECONDARY.get(nd)
        ring_col = COLORS.get(sec, "white") if sec else "white"
        ring_lw  = 1.8 if sec else 0.5
        # Pill-shaped node via FancyBboxPatch in data coordinates.
        # Width scales with the gene name length so longer symbols
        # (e.g. "Serpine1", "Tnfaip3", "Adora2b") aren't clipped.
        gene_fontsize = 6.3
        w = max(0.095, 0.017 * len(nd) + 0.045)
        h = 0.030
        rect = mpatches.FancyBboxPatch(
            (x - w / 2, y - h / 2), w, h,
            boxstyle="round,pad=0.005",
            facecolor=col, edgecolor=ring_col, linewidth=ring_lw,
            transform=ax.transData, clip_on=False, zorder=4,
        )
        ax.add_patch(rect)
        ax.text(x, y, nd, ha="center", va="center",
                fontsize=gene_fontsize, fontweight="bold", color="white",
                fontstyle="italic", zorder=5, clip_on=False)

    # Note explaining the dual-role ring, placed under "Core driver genes"
    ax.text(X_GENE, Y_TOP + 0.055, "(ring color = secondary/dual role)",
            ha="center", va="bottom", fontsize=5.6, color="#5D6D7E",
            style="italic", clip_on=False)

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

    ax.set_xlim(-0.55, 2.05)
    ax.set_ylim(-0.06, 1.15)
    ax.axis("off")

    patches = [
        mpatches.Patch(color=COLORS[k], label=CAT_LABELS[k])
        for k in ("immune_activation", "feedback_regulation", "tissue_repair")
    ]
    ax.legend(handles=patches, fontsize=6.5, loc="upper right",
              framealpha=0.80, handlelength=1.1, borderpad=0.5)


# ═════════════════════════════════════════════════════════════════════════════
#  7. load_expression_files
# ═════════════════════════════════════════════════════════════════════════════
# ═════════════════════════════════════════════════════════════════════════════
#  7. load_expression_files
# ═════════════════════════════════════════════════════════════════════════════
def load_expression_files(data_dir: Path):
    """
    Load four DEG files and return:
      expr_dfs – {condition_label: DataFrame with '_lfc' column, indexed by
                  lowercase gene symbol for case-insensitive lookup}

    NOTE: earlier versions of this script looked genes up by Entrez ID via a
    hand-maintained symbol->ID map (SUPPLEMENTARY_GENE_MAP). That assumed the
    expression files were indexed by Entrez ID. Checking the user's actual
    file confirmed they are indexed by gene SYMBOL ('genename' column), with
    'id' (Entrez) as a separate, unused column. Looking genes up directly by
    symbol removes an entire class of silent-drop bugs from incomplete or
    incorrect ID mappings.
    """
    FILE_MAP = {
        "Mouse scab":   "MM.control-M.scar_genes.csv",
        "Mouse wound":  "MM.control-M.wound_genes.csv",
        "Spalax scab":  "S.control-S.scar_genes.csv",
        "Spalax wound": "S.control-S.wound_genes.csv",
    }

    expr_dfs = {}
    for label, fname in FILE_MAP.items():
        fp = data_dir / fname
        df = pd.read_csv(fp)
        print(f"\n[load_expression] {fname}  shape={df.shape}  cols={list(df.columns)}")

        # Identify the gene-symbol column
        sym_col = None
        for col in ("genename", "gene_name", "Gene", "gene_symbol", "SYMBOL", "gene"):
            if col in df.columns:
                sym_col = col
                break
        if sym_col is None:
            sym_col = df.columns[0]
            print(f"  WARNING: no recognized gene-symbol column found; "
                  f"falling back to first column '{sym_col}'")
        else:
            print(f"  gene-symbol col → '{sym_col}'")

        # Detect LFC column
        for col in ("log2FoldChange", "logFC", "avg_log2FC", "Log2FC",
                    "log2FC", "LFC", "FoldChange"):
            if col in df.columns:
                df["_lfc"] = pd.to_numeric(df[col], errors="coerce")
                print(f"  LFC col → '{col}'")
                break
        else:
            df["_lfc"] = pd.to_numeric(df.iloc[:, -1], errors="coerce")
            print("  LFC col → last column (fallback)")

        # Case-insensitive index so lookups don't fail on case mismatches
        # (e.g. gene supplied as "Il6" vs. file storing "IL6" or "il6")
        df["_symbol_lower"] = df[sym_col].astype(str).str.lower()
        df = df.set_index("_symbol_lower")
        expr_dfs[label] = df

    return expr_dfs


# ═════════════════════════════════════════════════════════════════════════════
#  8. plot_gene_heatmap  (Panel C)
# ═════════════════════════════════════════════════════════════════════════════
def plot_gene_heatmap(
    ax: plt.Axes, expr_dfs: dict
) -> None:
    """
    Z-scored log2FC heatmap of three gene modules across four conditions.
    Conditions ordered: Mouse scab | Mouse wound | Spalax scab | Spalax wound.
    """
    MODULES = {
        "Inflammatory\neffectors": [g for g, c in GENE_PRIMARY.items() if c == "immune_activation"],
        "Feedback\nregulators":    [g for g, c in GENE_PRIMARY.items() if c == "feedback_regulation"],
        "Repair /\nremodeling":    [g for g, c in GENE_PRIMARY.items() if c == "tissue_repair"],
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
            gene_key = gene.lower()
            vals = []
            for cond in COL_ORDER:
                df = expr_dfs[cond]
                v = df.at[gene_key, "_lfc"] if gene_key in df.index else np.nan
                vals.append(v)

            n_found = sum(1 for v in vals if not np.isnan(v))
            if n_found >= 1:           # include if found in ≥1 condition
                rows.append(vals)
                row_labels.append(gene)
                row_modules.append(mod_name)
            else:
                print(f"  WARNING: {gene} (lookup key='{gene_key}') – not found in any "
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
        ax.get_xticklabels(), rotation=0, ha="center", fontsize=9.5
    )
    yticklabels = ax.set_yticklabels(
        ax.get_yticklabels(), fontstyle="italic", fontsize=6.4, va="center"
    )

    # ── Color-code gene names by module (matches Panel B's palette), and
    #    flag genes with a documented dual/secondary regulatory role ─────────
    for tick, gene, mod_name in zip(yticklabels, row_labels, row_modules):
        cat = MOD_CAT[mod_name]
        tick.set_color(COLORS[cat])
        tick.set_fontweight("bold")
        if gene in GENE_SECONDARY:
            tick.set_text(tick.get_text() + " \u2020")   # dagger = dual role

    # Legend explaining the gene-name color coding — placed below the
    # (now horizontal, compact) x-axis labels so it stays inside Panel C's
    # own margin and never intrudes on Panel A/B above.
    patches = [
        mpatches.Patch(color=COLORS[k], label=CAT_LABELS[k])
        for k in ("immune_activation", "feedback_regulation", "tissue_repair")
    ]
    leg = ax.legend(
        handles=patches, fontsize=7, loc="upper center",
        bbox_to_anchor=(0.5, -0.14), ncol=3, framealpha=0.85,
        handlelength=1.1, borderpad=0.5, columnspacing=1.2,
        title="Gene category  (\u2020 = documented dual/secondary role)",
        title_fontsize=7,
    )

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
    expr_dfs = load_expression_files(DATA_DIR)

    # ── Figure skeleton ───────────────────────────────────────────────────────
    fig = plt.figure(figsize=(18, 14))

    # Row 0: A (GO dot plot) | B (gene–GO network — gets extra width so the
    #        right-hand GO term labels fit, but not so much that A's own
    #        y-axis labels get squeezed illegible)
    # Row 1: C (full width, heatmap) spans both columns
    gs = gridspec.GridSpec(
        2, 2, figure=fig,
        hspace=0.25, wspace=0.22,
        left=0.13, right=0.97, top=0.96, bottom=0.13,
        width_ratios=[3, 4],
    )
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, :])

    # ── Draw panels ───────────────────────────────────────────────────────────
    print("\n── Panel A: GO enrichment ──────────────────────────────────────")
    plot_go_enrichment(axA, go_df)

    print("\n── Panel B: Gene–GO network ────────────────────────────────────")
    plot_gene_go_network(axB, go_df)

    print("\n── Panel C: Expression heatmap ─────────────────────────────────")
    plot_gene_heatmap(axC, expr_dfs)

    #print("\n── Panel D: KEGG enrichment ─────────────────────────────────────")
    #plot_kegg_enrichment(axD, kegg_df)

    # ── Bold panel labels ─────────────────────────────────────────────────────
    for ax, label in [(axA, "A"), (axB, "B"), (axC, "C"),
                       #(axD, "D")
                       ]:
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