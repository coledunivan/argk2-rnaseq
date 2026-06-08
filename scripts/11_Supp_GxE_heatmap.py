#!/usr/bin/env python3
"""
Supplemental Figure: Full G×E vertical heatmap, no TF-module grouping.

Changes from prior version:
  - Removes the primary TF/module grouping color bar.
  - Uses a vertical heatmap: genes are rows, contrasts are columns.
  - Column order:
      1. Genotype: argk-2−/− untreated vs N2 untreated
      2. Treatment: N2 +CoCl2 vs N2 untreated
      3. Stress genotype: argk-2−/− +CoCl2 vs N2 +CoCl2
      4. G×E interaction: non-additive component
  - Keeps all G×E genes, including TFLink-regulated and unregulated genes.
  - Sorts genes by G×E interaction value so sub-additive and supra-additive
    patterns are visually organized without module grouping.

INPUTS:
  DESeq2_effect_classification.csv
  data/reference/RNAseq_to_TF_Targets.csv   # only used for common gene-name lookup

OUTPUTS:
  Figure_S_full_GxE_heatmap_vertical_no_modules.pdf
  Figure_S_full_GxE_heatmap_vertical_no_modules.png
"""

import sys
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    from _contrasts import validate_contrast_columns, provenance_footer
except ImportError:
    validate_contrast_columns = None
    provenance_footer = None
    print("[warn] _contrasts.py not found. Continuing without contrast validation/provenance footer.")

# ---- Configuration ---------------------------------------------------------
def _resolve_input(filename, search_paths=None):
    """Find an input CSV at one of several canonical locations."""
    from pathlib import Path as _P
    if search_paths is None:
        search_paths = [".", "outputs/DESeq2", "../outputs/DESeq2"]
    for sp in search_paths:
        cand = _P(sp) / filename
        if cand.exists():
            return str(cand)
    return None

DESEQ_FILE  = _resolve_input("DESeq2_effect_classification.csv") \
              or "DESeq2_effect_classification.csv"
TFLINK_FILE = "data/reference/RNAseq_to_TF_Targets.csv"
import os
_OUT_DIR = "outputs/supplementary"
os.makedirs(_OUT_DIR, exist_ok=True)
OUT_PDF     = os.path.join(_OUT_DIR, "Figure_S_full_GxE_heatmap_vertical_no_modules.pdf")
OUT_PNG     = os.path.join(_OUT_DIR, "Figure_S_full_GxE_heatmap_vertical_no_modules.png")

HEATMAP_VMAX = 4.0

# These are the actual columns in the DESeq2 classification table.
# The order here controls the left-to-right order in the heatmap.
HEATMAP_COLS = [
    "lfc_geno_unt",    # argk-2 untreated vs N2 untreated
    "lfc_treat_N2",    # N2 +CoCl2 vs N2 untreated
    "lfc_geno_trt",    # argk-2 +CoCl2 vs N2 +CoCl2
    "lfc_gxe",         # non-additive interaction term
]

# Short labels to fit as column headers.
# The caption/figure title defines each one explicitly.
HEATMAP_LABELS = [
    "Genotype",
    "Treatment",
    "Stress\ngenotype",
    "G×E",
]

FULL_LABELS = {
    "Genotype": "argk-2−/− untreated vs N2 untreated",
    "Treatment": "N2 +CoCl₂ vs N2 untreated",
    "Stress genotype": "argk-2−/− +CoCl₂ vs N2 +CoCl₂",
    "G×E": "non-additive interaction term",
}

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 10,
    "axes.linewidth": 0.9,
    "pdf.fonttype": 42,
})


def load_gxe_genes():
    """Load G×E gene list from DESeq2 output."""
    if not Path(DESEQ_FILE).exists():
        sys.exit(f"[fatal] {DESEQ_FILE} not found. Place it next to this script.")

    print(f"[ok] reading G×E classification from {DESEQ_FILE}")
    genome = pd.read_csv(DESEQ_FILE)

    for col in HEATMAP_COLS + ["padj_gxe", "lfc_treat_RB"]:
        if col in genome.columns:
            genome[col] = pd.to_numeric(genome[col], errors="coerce")

    missing = [c for c in HEATMAP_COLS if c not in genome.columns]
    if missing:
        sys.exit(f"[fatal] Missing required contrast columns: {missing}")

    if validate_contrast_columns is not None:
        validate_contrast_columns(genome)

    gxe = genome[genome["effect_class"] == "GxE"].copy()
    print(f"   {len(gxe)} G×E genes from DESeq2")
    return gxe, genome


def load_tflink_for_names():
    """Load TFLink table for common-name lookup only."""
    if not Path(TFLINK_FILE).exists():
        print(f"[warn] {TFLINK_FILE} not found. Using IDs as labels.")
        return None
    return pd.read_csv(TFLINK_FILE, low_memory=False)


def get_common_name(cosmid, tf_df=None):
    """Resolve cosmid/gene ID to common gene name via TFLink Name.Target column."""
    if tf_df is None or "Gene" not in tf_df.columns or "Name.Target" not in tf_df.columns:
        return cosmid
    hits = tf_df[tf_df["Gene"] == cosmid]["Name.Target"].dropna()
    if len(hits) == 0:
        return cosmid
    name = hits.iloc[0]
    if name == "-" or pd.isna(name):
        return cosmid
    return str(name).split(";")[0].split(",")[0].strip().lower()


def main():
    gxe, genome = load_gxe_genes()
    tf_raw = load_tflink_for_names()

    gene_col = "gene_symbol" if "gene_symbol" in gxe.columns else "Gene"
    if gene_col not in gxe.columns:
        sys.exit("[fatal] Could not find a gene ID column: expected 'gene_symbol' or 'Gene'.")

    # Build plotting table.
    rows = []
    for _, row in gxe.iterrows():
        gene_id = row[gene_col]
        common = get_common_name(gene_id, tf_raw)
        label = common if common != gene_id else gene_id
        rows.append({
            "gene_id": gene_id,
            "label": label,
            **{k: row.get(k, np.nan) for k in HEATMAP_COLS},
        })

    df = pd.DataFrame(rows)

    # No TF/module grouping. Sort by G×E interaction to preserve the biological
    # sub-additive → supra-additive structure without using module labels.
    df = df.sort_values("lfc_gxe", ascending=True).reset_index(drop=True)

    mat = df[HEATMAP_COLS].fillna(0).values
    mat = np.clip(mat, -HEATMAP_VMAX, HEATMAP_VMAX)

    n_genes = len(df)
    print(f"\nHeatmap: {n_genes} G×E genes, no TF-module grouping")

    # Figure size scales with gene count; tall enough for readable gene labels.
    fig_h = max(9, n_genes * 0.18)
    fig_w = 5.8
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    im = ax.imshow(
        mat,
        cmap=plt.cm.RdBu_r,
        aspect="auto",
        vmin=-HEATMAP_VMAX,
        vmax=HEATMAP_VMAX,
        interpolation="nearest",
    )

    # Column labels at top.
    ax.set_xticks(np.arange(len(HEATMAP_COLS)))
    ax.set_xticklabels(HEATMAP_LABELS, fontsize=10)
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", length=0, pad=6)

    # Gene labels on y-axis.
    ax.set_yticks(np.arange(n_genes))
    ax.set_yticklabels(df["label"], fontsize=7)
    ax.tick_params(axis="y", length=0, pad=2)

    # Light grid to make the 4 contrast columns easy to read.
    ax.set_xticks(np.arange(-0.5, len(HEATMAP_COLS), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_genes, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.35)
    ax.tick_params(which="minor", bottom=False, left=False)

    for s in ["top", "right", "left", "bottom"]:
        ax.spines[s].set_visible(False)

    title = (
        f"Directional response of all {n_genes} G×E genes\n"
        "Genotype = argk-2−/− untreated vs N2 untreated; "
        "Treatment = N2 +CoCl₂ vs N2 untreated; "
        "Stress genotype = argk-2−/− +CoCl₂ vs N2 +CoCl₂"
    )
    ax.set_title(title, fontsize=10.5, fontweight="bold", pad=36)

    # Colorbar.
    cax = fig.add_axes([0.86, 0.25, 0.035, 0.5])
    cb = fig.colorbar(im, cax=cax, orientation="vertical", ticks=[-HEATMAP_VMAX, 0, HEATMAP_VMAX])
    cb.set_label("log$_2$FC", fontsize=9)
    cb.ax.tick_params(labelsize=8)

    # Footer/caption for the G×E column.
    fig.text(
        0.12,
        0.02,
        "G×E = genotype-by-treatment interaction term; genes sorted by G×E effect size. TF-module grouping removed.",
        ha="left",
        va="bottom",
        fontsize=7,
        color="0.35",
    )

    if provenance_footer is not None:
        provenance_footer(
            fig,
            contrast_keys=HEATMAP_COLS,
            script_name="figure_supplemental_full_gxe_vertical_no_modules.py",
        )

    plt.savefig(OUT_PDF, dpi=400, bbox_inches="tight")
    plt.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
    print(f"\n[done] {OUT_PDF}  +  {OUT_PNG}")


if __name__ == "__main__":
    main()
