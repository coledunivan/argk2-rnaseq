#!/usr/bin/env python3
"""
Figure 6 supplementary panels: G×E core architecture.

Two panels:
  Panel A: In-degree histogram across the 90 G×E genes — how many hub TFs
           regulate each G×E gene. Visualizes the "co-regulatory density"
           claim: the G×E set is heavily targeted, not randomly distributed.

  Panel B: 13 × 13 hub-TF Jaccard heatmap of target overlap restricted to
           G×E genes. Reveals the dense zip-2 / ceh-60 / daf-16 / pha-4 /
           pqm-1 co-regulatory module.

Reads:
  outputs/GRN_tables/nodes_full90.tsv
  outputs/GRN_tables/edges_combined_full90.tsv
  outputs/DESeq2/DESeq2_effect_classification.csv

Writes:
  outputs/supplementary/Figure6_panelB_GxE_indegree.{pdf,png}
  outputs/supplementary/Figure6_panelC_TF_coregulation.{pdf,png}
  outputs/supplementary/Figure6_GxE_core_stats.csv

Usage:
  python3 scripts/13_Supp_GxE_core_analysis.py
"""

from __future__ import annotations
import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

# ----- Style (matches Supp Fig 4 + qPCR panels) -----
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 11,
    "axes.linewidth": 0.9,
    "axes.labelsize": 11,
    "axes.titlesize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "pdf.fonttype": 42,  # editable text in PDF
    "ps.fonttype": 42,
})

NODES_TSV = "outputs/GRN_tables/nodes_full90.tsv"
EDGES_TSV = "outputs/GRN_tables/edges_combined_full90.tsv"
EFFECT_CSV = "outputs/DESeq2/DESeq2_effect_classification.csv"

OUT_DIR = Path("outputs/supplementary")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def main() -> None:
    # ---- Load ---------------------------------------------------------------
    for f in (NODES_TSV, EDGES_TSV, EFFECT_CSV):
        if not Path(f).exists():
            raise FileNotFoundError(
                f"Required input missing: {f}\n"
                "Run 01_DESeq2_and_TF_enrichment.R + 07_GRN_inputs.R first."
            )

    nodes = pd.read_csv(NODES_TSV, sep="\t")
    edges = pd.read_csv(EDGES_TSV, sep="\t")
    eff   = pd.read_csv(EFFECT_CSV)

    # Identify G×E targets via the nodes table (cosmid IDs already match
    # the edges file's casing convention). Cross-check against the effect
    # classification CSV for consistency.
    gxe_nodes = nodes[(nodes["node_class"] == "effector") &
                      (nodes["effect_category"] == "GxE")]
    gxe_ids = set(gxe_nodes["node_id"].astype(str))
    gxe_syms_lc = set(
        eff.loc[eff["effect_class"] == "GxE", "gene_symbol"].astype(str).str.lower()
    )
    print(f"G×E genes (effect_class CSV):              {len(gxe_syms_lc)}")
    print(f"G×E effector nodes in GRN nodes table:     {len(gxe_ids)}")
    print(f"Edges in GRN:                              {len(edges)}")
    print(f"Hub TFs (unique sources):                  {edges['source'].nunique()}")

    # Use the node-table set as ground truth for matching edges (case-aligned)
    gxe_genes = gxe_ids

    # ---- Per-G×E-gene in-degree (count of distinct hub TFs targeting it) ----
    edges_gxe = edges[edges["target"].isin(gxe_genes)].copy()
    in_deg = (
        edges_gxe.groupby("target")["source"].nunique()
        .reindex(sorted(gxe_genes), fill_value=0)
        .sort_values(ascending=False)
    )

    # Distribution including the genes with zero TF edges in the network
    bins = list(range(int(in_deg.max()) + 2))
    counts, edges_bin = np.histogram(in_deg.values, bins=bins)
    n_zero        = int((in_deg == 0).sum())
    n_core5       = int((in_deg >= 5).sum())
    n_core10      = int((in_deg >= 10).sum())
    median_deg    = float(np.median(in_deg.values))

    print(f"\nIn-degree distribution:")
    print(f"  G×E genes with 0 hub-TF edges:     {n_zero}")
    print(f"  G×E genes with ≥5 hub-TF edges:    {n_core5}")
    print(f"  G×E genes with ≥10 hub-TF edges:   {n_core10}  (the 'dense core')")
    print(f"  Median in-degree (incl. zeros):    {median_deg}")
    print(f"  Median in-degree (≥1 only):        "
          f"{float(np.median(in_deg[in_deg > 0])):.1f}")

    # ---- PANEL B: In-degree histogram --------------------------------------
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    colors = []
    for k in range(len(counts)):
        # Colour zero-edge bin grey (unmapped in network), then graduate up
        if k == 0:
            colors.append("#B0B0B0")
        elif k >= 10:
            colors.append("#C2185B")     # dense core
        elif k >= 5:
            colors.append("#E89441")     # convergence zone
        else:
            colors.append("#4A90C2")     # sparsely regulated
    ax.bar(range(len(counts)), counts, color=colors, edgecolor="black",
           linewidth=0.7, width=0.85)
    for k, c in enumerate(counts):
        if c > 0:
            ax.text(k, c + 0.4, str(c), ha="center", va="bottom",
                    fontsize=9, color="#333")
    ax.axvline(median_deg, color="black", linestyle="--", linewidth=1.0,
               alpha=0.6, label=f"median = {int(median_deg)}")
    ax.set_xlabel("Number of hub TFs regulating the G×E gene")
    ax.set_ylabel("Number of G×E genes")
    ax.set_title("G×E genes are densely co-regulated by hub TFs",
                 loc="left", fontweight="bold")
    ax.set_xticks(range(len(counts)))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Inline legend
    from matplotlib.patches import Patch
    legend_handles = [
        Patch(facecolor="#B0B0B0", edgecolor="black",
              label=f"No TF edge in GRN  (n={counts[0]})"),
        Patch(facecolor="#4A90C2", edgecolor="black",
              label=f"1–4 TFs  (n={int(counts[1:5].sum())})"),
        Patch(facecolor="#E89441", edgecolor="black",
              label=f"5–9 TFs (convergence)  (n={int(counts[5:10].sum())})"),
        Patch(facecolor="#C2185B", edgecolor="black",
              label=f"≥10 TFs (dense core)  (n={int(counts[10:].sum())})"),
    ]
    ax.legend(handles=legend_handles, loc="upper right", frameon=False,
              handlelength=1.2, fontsize=9)
    plt.tight_layout()

    pdf_a = OUT_DIR / "Figure6_panelB_GxE_indegree.pdf"
    png_a = OUT_DIR / "Figure6_panelB_GxE_indegree.png"
    fig.savefig(pdf_a, bbox_inches="tight")
    fig.savefig(png_a, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved: {pdf_a}")

    # ---- PANEL C: TF–TF Jaccard heatmap (G×E target overlap) ---------------
    hub_tfs = sorted(edges_gxe["source"].unique())
    print(f"\nBuilding {len(hub_tfs)}×{len(hub_tfs)} TF co-regulation heatmap...")

    # For each TF, set of G×E targets
    tf_targets = {tf: set(edges_gxe.loc[edges_gxe["source"] == tf, "target"])
                  for tf in hub_tfs}

    # Symmetric Jaccard matrix; diagonal = 1
    n = len(hub_tfs)
    jacc = np.ones((n, n), dtype=float)
    counts_mat = np.zeros((n, n), dtype=int)
    for i, ti in enumerate(hub_tfs):
        for j, tj in enumerate(hub_tfs):
            si, sj = tf_targets[ti], tf_targets[tj]
            counts_mat[i, j] = len(si & sj)
            if i == j:
                jacc[i, j] = 1.0
            else:
                union = len(si | sj)
                jacc[i, j] = (len(si & sj) / union) if union else 0.0

    # Hierarchical clustering on (1 − Jaccard); reorder rows/cols
    dist = 1.0 - jacc
    np.fill_diagonal(dist, 0.0)
    cond = squareform(dist, checks=False)
    link = linkage(cond, method="average")
    order = leaves_list(link)
    hub_tfs_ord = [hub_tfs[k] for k in order]
    jacc_ord = jacc[np.ix_(order, order)]
    counts_ord = counts_mat[np.ix_(order, order)]

    # Plot
    fig, ax = plt.subplots(figsize=(7.5, 6.8))
    im = ax.imshow(jacc_ord, cmap="magma_r", aspect="equal",
                   vmin=0, vmax=jacc_ord[~np.eye(n, dtype=bool)].max())
    # Annotate each cell with intersection count (raw shared G×E targets)
    for i in range(n):
        for j in range(n):
            v = counts_ord[i, j]
            if v == 0:
                continue
            txt_color = "white" if jacc_ord[i, j] > 0.35 else "#202020"
            ax.text(j, i, str(v), ha="center", va="center",
                    fontsize=8, color=txt_color)
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(hub_tfs_ord, rotation=45, ha="right",
                       fontstyle="italic")
    ax.set_yticklabels(hub_tfs_ord, fontstyle="italic")
    # Bold zip-2 axis labels for emphasis
    for label in ax.get_xticklabels():
        if label.get_text() == "zip-2":
            label.set_fontweight("bold")
            label.set_color("#C2185B")
    for label in ax.get_yticklabels():
        if label.get_text() == "zip-2":
            label.set_fontweight("bold")
            label.set_color("#C2185B")

    cbar = plt.colorbar(im, ax=ax, shrink=0.7, pad=0.02)
    cbar.set_label("Jaccard index (G×E target overlap)", fontsize=10)
    ax.set_title("Hub-TF co-regulation of G×E genes\n"
                 "(cells annotated with # shared G×E targets)",
                 loc="left", fontweight="bold", pad=10)
    plt.tight_layout()

    pdf_b = OUT_DIR / "Figure6_panelC_TF_coregulation.pdf"
    png_b = OUT_DIR / "Figure6_panelC_TF_coregulation.png"
    fig.savefig(pdf_b, bbox_inches="tight")
    fig.savefig(png_b, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {pdf_b}")

    # ---- Stats CSV ----------------------------------------------------------
    stats_rows = []
    for tf in hub_tfs:
        targs = tf_targets[tf]
        stats_rows.append({
            "TF": tf,
            "n_GxE_targets": len(targs),
            "GxE_targets":   ";".join(sorted(targs)),
        })
    pd.DataFrame(stats_rows).sort_values(
        "n_GxE_targets", ascending=False
    ).to_csv(OUT_DIR / "Figure6_GxE_core_stats.csv", index=False)

    # Per-gene core membership — match via lowercase gene_symbol
    eff_lc = eff.copy()
    eff_lc["gene_symbol_lc"] = eff_lc["gene_symbol"].astype(str).str.lower()
    core_df = (
        pd.DataFrame({"gene": in_deg.index, "n_hub_TFs": in_deg.values})
        .merge(
            eff_lc[["gene_symbol_lc", "gene_symbol", "lfc_treat_N2",
                    "lfc_treat_RB", "lfc_gxe", "padj_gxe"]],
            left_on="gene", right_on="gene_symbol_lc", how="left"
        )
        .drop(columns="gene_symbol_lc")
    )
    core_df["category"] = pd.cut(
        core_df["n_hub_TFs"], bins=[-0.5, 0.5, 4.5, 9.5, 100],
        labels=["unmapped", "sparse_1-4", "convergence_5-9", "dense_core_10+"]
    )
    core_df.to_csv(OUT_DIR / "Figure6_GxE_per_gene_indegree.csv", index=False)
    print(f"Saved: {OUT_DIR}/Figure6_GxE_per_gene_indegree.csv")
    print(f"Saved: {OUT_DIR}/Figure6_GxE_core_stats.csv")

    # ---- Headline summary --------------------------------------------------
    print("\n" + "=" * 60)
    print("HEADLINE NUMBERS for Methods / Results text")
    print("=" * 60)
    # Broader TF coverage — count G×E genes regulated by ANY TF (full TFLink),
    # not restricted to the 15 hub TFs rendered in Figure 6. This is the
    # number reported in the thesis (81/90).
    tf_targets_path = Path("outputs/DESeq2/DESeq2_TF_targets_DE.csv")
    if tf_targets_path.exists():
        tfd = pd.read_csv(tf_targets_path)
        gxe_with_any_tf = (
            (tfd["effect_class"] == "GxE") & (tfd["n_TFs"] >= 1)
        ).sum()
        print(f"  G×E genes with ≥1 TF (full TFLink):  {gxe_with_any_tf} / 90"
              f"   ← thesis figure")
    print(f"  Hub TFs in the rendered GRN:        {n}")
    print(f"  G×E genes with ≥1 hub-TF edge:      {(in_deg > 0).sum()} / 90"
          f"   (hub-TF subset only)")
    print(f"  Median in-degree (regulated only):  "
          f"{float(np.median(in_deg[in_deg > 0])):.1f}")
    print(f"  G×E genes with ≥5 hub-TF edges:     {n_core5}")
    print(f"  G×E genes with ≥10 hub-TF edges:    {n_core10}  (the 'dense core')")
    # Strongest co-regulator pair (off-diagonal)
    iu = np.triu_indices(n, k=1)
    pair_idx = np.argmax(counts_mat[iu])
    i_max, j_max = iu[0][pair_idx], iu[1][pair_idx]
    print(f"  Tightest co-regulator pair (raw count): "
          f"{hub_tfs[i_max]} ↔ {hub_tfs[j_max]} "
          f"(shared = {counts_mat[i_max, j_max]} G×E targets, "
          f"Jaccard = {jacc[i_max, j_max]:.2f})")


if __name__ == "__main__":
    main()
