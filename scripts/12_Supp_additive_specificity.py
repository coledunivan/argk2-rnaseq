#!/usr/bin/env python3
"""
Supplemental Figure — CoCl₂ response is maintained for additive genes
and collapses specifically for G×E genes.

Shows that the transcriptional collapse observed in the G×E effectors
is NOT a global consequence of CoCl₂ exposure in the argk-2⁻/⁻ background.
Genes classified as "Treatment" or "Additive" (where genotype and treatment
effects sum normally) maintain wild-type-like CoCl₂ induction, while
"GxE" genes show the characteristic sub-additive collapse or supra-additive
runaway.

Reads: DESeq2_effect_classification.csv
Produces: Figure_S_additive_vs_gxe.pdf + .png
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

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

GENOME_FILE = _resolve_input("DESeq2_effect_classification.csv") \
              or "DESeq2_effect_classification.csv"
import os
_OUT_DIR = "outputs/supplementary"
os.makedirs(_OUT_DIR, exist_ok=True)
OUT_PDF = os.path.join(_OUT_DIR, "Figure_S_additive_vs_gxe.pdf")
OUT_PNG = os.path.join(_OUT_DIR, "Figure_S_additive_vs_gxe.png")

LFC_ARGK2_TRT_ALIASES = ["lfc_treat_RB", "lfc_treat_argk2", "lfc_treat_argk_2"]

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica","Arial","Liberation Sans","DejaVu Sans"],
    "font.size": 11,
    "axes.linewidth": 0.9,
    "pdf.fonttype": 42,
})


def main():
    if not Path(GENOME_FILE).exists():
        print(f"[error] {GENOME_FILE} not found"); return

    genome = pd.read_csv(GENOME_FILE)

    # Resolve argk-2 treatment LFC column
    alias = next((a for a in LFC_ARGK2_TRT_ALIASES if a in genome.columns), None)
    if alias is None:
        print(f"[warn] no treatment-in-argk2 LFC column found — using lfc_geno_trt as proxy")
        if "lfc_geno_trt" in genome.columns:
            genome["lfc_treat_RB"] = genome["lfc_geno_trt"]
        else:
            print("[error] no suitable column"); return
    elif alias != "lfc_treat_RB":
        genome["lfc_treat_RB"] = genome[alias]

    for col in ["lfc_treat_N2","lfc_treat_RB","lfc_gxe","padj_gxe"]:
        if col in genome.columns:
            genome[col] = pd.to_numeric(genome[col], errors="coerce")

    df = genome[genome["lfc_treat_N2"].notna() & genome["lfc_treat_RB"].notna()].copy()

    fig, axes = plt.subplots(1, 3, figsize=(14, 5))

    lim = 9.0

    # --- Panel A: Scatter of Treatment/Additive genes ---
    ax = axes[0]
    cats_ok = ["Treatment","Additive"]
    ok = df[df["effect_class"].isin(cats_ok)]
    ax.plot([-lim, lim], [-lim, lim], color="#666", lw=1.0, ls="-", zorder=1)
    ax.scatter(ok["lfc_treat_N2"].clip(-lim, lim),
               ok["lfc_treat_RB"].clip(-lim, lim),
               s=8, c="#4A90C2", alpha=0.5, edgecolors="none", zorder=2)
    # Correlation
    from scipy.stats import pearsonr
    mask = ok["lfc_treat_N2"].notna() & ok["lfc_treat_RB"].notna()
    r_ok, p_ok = pearsonr(ok.loc[mask, "lfc_treat_N2"], ok.loc[mask, "lfc_treat_RB"])
    # Median deviation from diagonal
    dev_ok = (ok["lfc_treat_RB"] - ok["lfc_treat_N2"]).dropna()
    ax.text(0.05, 0.95,
            f"Treatment + Additive\nn = {len(ok)}\n"
            f"median deviation = {dev_ok.median():.2f}",
            transform=ax.transAxes, ha="left", va="top", fontsize=10,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#4A90C2",
                      alpha=0.9))
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("CoCl$_2$ response in N2 (log$_2$FC)")
    ax.set_ylabel("CoCl$_2$ response in argk-2$^{-/-}$ (log$_2$FC)")
    ax.set_title("Additive genes: CoCl$_2$ response\nsymmetric around y = x",
                 fontsize=11.5, fontweight="bold")
    for s in ["top","right"]: ax.spines[s].set_visible(False)
    ax.set_aspect("equal")

    # --- Panel B: Scatter of GxE genes ---
    ax = axes[1]
    gxe = df[df["effect_class"] == "GxE"]
    ax.plot([-lim, lim], [-lim, lim], color="#666", lw=1.0, ls="-", zorder=1)
    ax.scatter(gxe["lfc_treat_N2"].clip(-lim, lim),
               gxe["lfc_treat_RB"].clip(-lim, lim),
               s=12, c="#E07B3D", alpha=0.6, edgecolors="none", zorder=2)
    dev_gxe = (gxe["lfc_treat_RB"] - gxe["lfc_treat_N2"]).dropna()
    ax.text(0.05, 0.95,
            f"G$\\times$E genes\nn = {len(gxe)}\n"
            f"median deviation = {dev_gxe.median():.2f}",
            transform=ax.transAxes, ha="left", va="top", fontsize=10,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#E07B3D",
                      alpha=0.9))
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("CoCl$_2$ response in N2 (log$_2$FC)")
    ax.set_ylabel("CoCl$_2$ response in argk-2$^{-/-}$ (log$_2$FC)")
    ax.set_title("G$\\times$E genes: CoCl$_2$ response\nsuppressed below y = x",
                 fontsize=11.5, fontweight="bold")
    for s in ["top","right"]: ax.spines[s].set_visible(False)
    ax.set_aspect("equal")

    # --- Panel C: Box/violin comparing deviation from additivity ---
    ax = axes[2]
    # Deviation = lfc_treat_RB - lfc_treat_N2 (additive = ~0)
    for cat_set, label, color, xpos in [
        (["Treatment","Additive"], "Additive\n+ Treatment", "#4A90C2", 0),
        (["GxE"], "G$\\times$E", "#E07B3D", 1),
    ]:
        sub = df[df["effect_class"].isin(cat_set)].copy()
        dev = (sub["lfc_treat_RB"] - sub["lfc_treat_N2"]).dropna()
        # Violin
        parts = ax.violinplot([dev.values], positions=[xpos],
                              showmedians=True, showextrema=False)
        for pc in parts["bodies"]:
            pc.set_facecolor(color)
            pc.set_alpha(0.4)
        parts["cmedians"].set_color(color)
        parts["cmedians"].set_linewidth(2)
        # Median text
        med = dev.median()
        ax.text(xpos + 0.15, med, f"median = {med:.2f}",
                fontsize=9, color=color, va="center")

    ax.axhline(0, color="#666", lw=0.8, ls="--", zorder=0)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Additive\n+ Treatment", "G$\\times$E"], fontsize=11)
    ax.set_ylabel("Deviation from additive expectation\n"
                   "(argk-2 CoCl$_2$ response $-$ N2 CoCl$_2$ response, log$_2$FC)",
                   fontsize=9)
    ax.set_title("Suppression is G$\\times$E-specific,\nnot a global CoCl$_2$ effect",
                 fontsize=11.5, fontweight="bold")
    for s in ["top","right"]: ax.spines[s].set_visible(False)

    # Mann-Whitney test between the two distributions
    from scipy.stats import mannwhitneyu
    add_dev = (df[df["effect_class"].isin(["Treatment","Additive"])]["lfc_treat_RB"] -
               df[df["effect_class"].isin(["Treatment","Additive"])]["lfc_treat_N2"]).dropna()
    gxe_dev = (df[df["effect_class"]=="GxE"]["lfc_treat_RB"] -
               df[df["effect_class"]=="GxE"]["lfc_treat_N2"]).dropna()
    _, mw_p = mannwhitneyu(add_dev, gxe_dev, alternative="two-sided")
    ax.text(0.5, 0.97, f"Mann-Whitney p = {mw_p:.1e}",
            transform=ax.transAxes, ha="center", va="top",
            fontsize=9.5, fontweight="bold", color="#333")

    for i, lbl in enumerate("abc"):
        axes[i].text(-0.12, 1.06, lbl, transform=axes[i].transAxes,
                     fontsize=18, fontweight="bold")

    plt.tight_layout(w_pad=2.0)
    plt.savefig(OUT_PDF, dpi=400, bbox_inches="tight")
    plt.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
    print(f"[done] {OUT_PDF}  +  {OUT_PNG}")


if __name__ == "__main__":
    main()
