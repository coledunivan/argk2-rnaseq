#!/usr/bin/env python3
"""
qPCR reanalysis pipeline - argk-2 oxidative stress study.

Takes Bio-Rad CFX "Quantification Cq Results.csv" files from multiple runs,
normalizes to tba-1 (DeltaCq), computes three key contrasts via DeltaDeltaCq,
and builds a publishable multi-panel figure.

The three contrasts match the RNA-seq main effects in Figure 6:

    Contrast 1  "N2 treated vs N2 untreated"
                = treatment main effect at reference genotype
                = DeltaCq(N2 treated) - DeltaCq(N2 untreated)

    Contrast 2  "argk-2(-/-) untreated vs N2 untreated"
                = genotype main effect at baseline
                = DeltaCq(mutant untreated) - DeltaCq(N2 untreated)

    Contrast 3  "argk-2(-/-) treated vs N2 treated"   <- KEY: genotype under stress
                = genotype effect under treatment = genotype + GxE interaction
                = DeltaCq(mutant treated) - DeltaCq(N2 treated)

This is the correct framing for testing the RNA-seq prediction that argk-2
mutants respond differently to oxidative stress than wild type.

Outputs:
    qPCR_combined_deltadeltaCq.csv       per-gene per-contrast log2FC table
    qPCR_figure_panel.pdf                publishable multi-panel figure
    qPCR_figure_panel.png                raster preview
"""
import os
import re
import glob
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy import stats

# ============================================================================
# CONFIGURATION
# ============================================================================

QPCR_ROOT = Path("data/qpcr") if Path("data/qpcr").exists() else Path(".")
REF_GENE  = "tba-1"

# Minimum biological runs required for a gene x contrast to appear in the
# main figure. Single-run measurements go to the supplementary table only.
MIN_RUNS_FOR_FIG = 2

# Genes to exclude from the figure entirely
GENES_EXCLUDE = ["fbxa-79", "nhr-58"]

# Rename targets that were mislabeled on qPCR plates
TARGET_RENAMES  = {"t22f3.3": "t22f3.11"}

# Per-run dots >this many IQRs from the group median are hidden from the
# bar-chart overlay (NOT removed from aggregate stats).
DOT_OUTLIER_IQR = 2.5

# Map all strain labels used across plates to a canonical name
STRAIN_ALIASES = {
    "n2":           "N2",
    "fau5":         "argk-2",    # lab name for argk-2(-/-) (RB2060 genetic background)
    "rb2060":       "argk-2",
    "argk-2":       "argk-2",
    "argk2":        "argk-2",
}

# Condition parser: split "N2 treated" into ("N2", "treated")
def parse_sample(sample):
    """Return (strain, treatment) or (None, None) if not a biological sample."""
    if not isinstance(sample, str):
        return None, None
    s = sample.strip().lower()
    # Reject standard curve dilution points
    if re.fullmatch(r"d\d+", s):
        return None, None
    # Reject pilot Greek-letter replicates
    if " alpha" in s or " beta" in s or " gamma" in s:
        # These are the 'oac-7' pilot run - merge back
        m = re.match(r"^(n2|fau5?)\s+([ct])\s", s)
        if m:
            strain = m.group(1)
            treatment = "treated" if m.group(2) == "t" else "untreated"
        else:
            return None, None
    else:
        # Normal parsing: "strain" or "strain treated"
        parts = s.split()
        if not parts:
            return None, None
        strain_raw = parts[0]
        treatment = "treated" if "treated" in s else "untreated"
        strain = strain_raw

    canonical = STRAIN_ALIASES.get(strain)
    if canonical is None:
        return None, None
    return canonical, treatment

# ============================================================================
# STEP 1: LOAD ALL Cq FILES
# ============================================================================

def load_all_cq(root):
    """Return long-format DataFrame of every Cq well across every run.

    Columns: run, target, strain, treatment, replicate, cq
    """
    rows = []
    for run_dir in sorted(root.iterdir()):
        if not run_dir.is_dir():
            continue
        # Older CFX exports use "Cq Results" (with space), newer exports use
        # "Cq_Results" (underscore). Accept both.
        cq_files = list(run_dir.glob("*Cq Results*.csv")) \
                 + list(run_dir.glob("*Cq_Results*.csv"))
        if not cq_files:
            continue
        cq_path = cq_files[0]
        run_name = run_dir.name

        df = pd.read_csv(cq_path)
        if "Cq" not in df.columns:
            continue

        for _, row in df.iterrows():
            target = str(row.get("Target", "")).strip()
            sample = str(row.get("Sample", "")).strip()
            # Some Bio-Rad CFX plate layouts record the condition in
            # "Biological Set Name" instead of "Sample" (e.g., the ColeD_41726
            # run of 04/2026, where the FAU5-treated wells have an empty
            # Sample column and the label in Biological Set Name). Fall back
            # to that column when Sample is missing.
            if sample.lower() in ("", "nan"):
                bsn = str(row.get("Biological Set Name", "")).strip()
                if bsn.lower() not in ("", "nan"):
                    sample = bsn
            cq_raw = row.get("Cq")

            if not target or not sample or target.lower() in ("nan", ""):
                continue
            if sample.lower() == "nan":
                continue

            # Parse Cq value
            try:
                cq = float(cq_raw)
                if np.isnan(cq):
                    continue
            except (TypeError, ValueError):
                continue

            strain, treatment = parse_sample(sample)
            if strain is None:
                continue  # standard curve point or unparseable

            rows.append({
                "run":       run_name,
                "target":    target.lower(),
                "strain":    strain,
                "treatment": treatment,
                "sample":    sample,
                "well":      row.get("Well", ""),
                "cq":        cq,
            })
    return pd.DataFrame(rows)

# ============================================================================
# STEP 2: QC - FLAG OUTLIER WELLS WITHIN TECHNICAL REPLICATES
# ============================================================================

def deduplicate_runs(df):
    """Some folders in the qPCR dataset are email forwards of the same plate
    already in another folder (the '_No_subject*' and 'Fw__*' folders). These
    have identical Cq values - if left in, they inflate n_runs with the same
    measurement counted twice.

    Match runs by the tuple of (sorted wells -> cq values) and keep only the
    first occurrence of each unique signature. Uses run-name length as a
    tiebreaker so the folder with the descriptive name wins over '_No_subject_'.
    """
    # Build a signature per run from well-level Cq values
    sigs = {}
    for run, grp in df.groupby("run"):
        sig = tuple(sorted(
            (r["target"], r["strain"], r["treatment"], r["well"], round(r["cq"], 6))
            for _, r in grp.iterrows()
        ))
        sigs[run] = sig

    # Group runs by signature, prefer descriptive names
    from collections import defaultdict
    sig_to_runs = defaultdict(list)
    for run, sig in sigs.items():
        sig_to_runs[sig].append(run)

    keep_runs = set()
    for sig, runs in sig_to_runs.items():
        # Prefer runs whose names do NOT start with '_No_subject' or 'Fw__'
        def _priority(r):
            if r.startswith("_No_subject"): return 2
            if r.startswith("Fw__"):        return 1
            return 0
        runs_sorted = sorted(runs, key=lambda r: (_priority(r), -len(r)))
        keep_runs.add(runs_sorted[0])
        if len(runs) > 1:
            dropped = [r for r in runs if r != runs_sorted[0]]
            print(f"  DEDUP: keeping '{runs_sorted[0]}', "
                  f"dropping {dropped} (identical Cq content)")

    return df[df["run"].isin(keep_runs)].copy()

def flag_outliers(df, sd_threshold=0.5, max_drop=2):
    """Within each (run, target, strain, treatment) group, iteratively drop
    the most extreme well (furthest from mean) until either SD < threshold
    or we've dropped max_drop wells. Stops at 2 dropped wells so we don't
    annihilate small groups.

    Technical replicate outlier filtering is standard qPCR QC - a few bad
    wells from edge effects, bubbles, or pipetting errors can ruin a
    technical triplicate.
    """
    df = df.copy()
    df["keep"] = True
    for key, grp in df.groupby(["run","target","strain","treatment"]):
        if len(grp) < 3:
            continue
        working_idx = list(grp.index)
        dropped = 0
        while dropped < max_drop and len(working_idx) >= 2:
            working_cqs = df.loc[working_idx, "cq"].values
            sd = np.std(working_cqs, ddof=1)
            if sd < sd_threshold:
                break
            m = np.mean(working_cqs)
            worst_pos = int(np.argmax(np.abs(working_cqs - m)))
            df.loc[working_idx[worst_pos], "keep"] = False
            working_idx.pop(worst_pos)
            dropped += 1
    return df

# ============================================================================
# STEP 3: COMPUTE DELTA CQ (TARGET - REF) WITHIN EACH RUN/CONDITION
# ============================================================================

def compute_delta_cq(df, ref_gene=REF_GENE):
    """For each (run, strain, treatment) group:
      - Mean Cq of tba-1 (after outlier removal)
      - For each target gene, per-well DeltaCq = Cq_target - mean_Cq_ref
    This is the standard qPCR normalization step: each target well is
    referenced to the average of its housekeeping wells on the same plate
    and condition.

    Returns long-format DataFrame: run, target, strain, treatment, dcq.
    """
    df = df[df["keep"]].copy()
    out = []
    for (run, strain, treatment), grp in df.groupby(["run","strain","treatment"]):
        ref_rows = grp[grp["target"] == ref_gene]
        if ref_rows.empty:
            continue
        ref_mean = ref_rows["cq"].mean()

        for _, row in grp.iterrows():
            if row["target"] == ref_gene:
                continue
            dcq = row["cq"] - ref_mean
            out.append({
                "run":       run,
                "target":    row["target"],
                "strain":    strain,
                "treatment": treatment,
                "dcq":       dcq,
            })
    return pd.DataFrame(out)

# ============================================================================
# STEP 4: COMPUTE DELTA DELTA CQ FOR EACH OF THE 3 CONTRASTS
# ============================================================================

CONTRASTS = [
    {
        "name":    "N2 treated\nvs N2 untreated",
        "short":   "treat_N2",
        "num":     ("N2", "treated"),
        "den":     ("N2", "untreated"),
        "meaning": "Treatment effect in wild type",
    },
    {
        "name":    "argk-2$^{-/-}$ untreated\nvs N2 untreated",
        "short":   "geno_unt",
        "num":     ("argk-2", "untreated"),
        "den":     ("N2", "untreated"),
        "meaning": "Baseline genotype effect",
    },
    {
        "name":    "argk-2$^{-/-}$ treated\nvs N2 treated",
        "short":   "geno_trt",
        "num":     ("argk-2", "treated"),
        "den":     ("N2", "treated"),
        "meaning": "Genotype effect under treatment (genotype + G$\\times$E)",
    },
]

def compute_ddcq_contrasts(dcq_df):
    """For each (run, target, contrast), compute:
      - mean DeltaCq of numerator condition
      - mean DeltaCq of denominator condition
      - DeltaDeltaCq = num - den
      - log2FC = -DeltaDeltaCq   (since lower Cq = higher expression)
    Each run is treated as an independent biological experiment.

    The analytical convention: log2FC > 0 means numerator is UP relative to
    denominator. DeltaDeltaCq < 0 means the numerator had lower Cq (more
    mRNA), so log2FC = -DeltaDeltaCq.

    Returns long-format DataFrame with columns:
      run, target, contrast, ddcq, log2fc, n_num, n_den
    """
    out = []
    for (run, target), grp in dcq_df.groupby(["run","target"]):
        for c in CONTRASTS:
            num = grp[(grp["strain"] == c["num"][0]) &
                      (grp["treatment"] == c["num"][1])]["dcq"]
            den = grp[(grp["strain"] == c["den"][0]) &
                      (grp["treatment"] == c["den"][1])]["dcq"]
            if num.empty or den.empty:
                continue
            ddcq = num.mean() - den.mean()
            log2fc = -ddcq
            out.append({
                "run":      run,
                "target":   target,
                "contrast": c["short"],
                "ddcq":     ddcq,
                "log2fc":   log2fc,
                "n_num":    len(num),
                "n_den":    len(den),
                "dcq_num_vals": list(num.values),
                "dcq_den_vals": list(den.values),
            })
    return pd.DataFrame(out)

# ============================================================================
# STEP 5: AGGREGATE ACROSS RUNS (per-gene per-contrast mean and CI)
# ============================================================================

def aggregate_across_runs(ddcq_df):
    """Per (target, contrast), compute:
      - mean log2FC across runs
      - SEM across runs
      - t-test vs 0 (Welch's one-sample)
      - n_runs, n_total_wells
    Each run contributes one log2FC value, so we're testing whether the
    per-run log2FC estimates differ significantly from zero.
    """
    rows = []
    for (target, contrast), grp in ddcq_df.groupby(["target","contrast"]):
        vals = grp["log2fc"].values
        if len(vals) == 0:
            continue
        mean = np.mean(vals)
        n = len(vals)
        if n >= 2:
            sem = np.std(vals, ddof=1) / np.sqrt(n)
            t_stat, p_val = stats.ttest_1samp(vals, 0.0)
        else:
            sem = np.nan
            p_val = np.nan

        n_total_wells = int(grp["n_num"].sum() + grp["n_den"].sum())
        rows.append({
            "target":         target,
            "contrast":       contrast,
            "log2fc_mean":    mean,
            "log2fc_sem":     sem,
            "n_runs":         n,
            "n_wells_total":  n_total_wells,
            "p_value":        p_val,
            "runs":           ";".join(sorted(set(grp["run"]))),
        })
    return pd.DataFrame(rows)

def bh_adjust(pvals):
    """Benjamini-Hochberg FDR correction."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    finite_mask = np.isfinite(p)
    finite = p[finite_mask]
    if len(finite) == 0:
        return np.full(n, np.nan)
    order = np.argsort(finite)
    ranked = finite[order]
    adj = ranked * len(finite) / np.arange(1, len(finite) + 1)
    # Enforce monotonicity
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    out_finite = np.empty_like(finite)
    out_finite[order] = adj
    out = np.full(n, np.nan)
    out[finite_mask] = out_finite
    return out

def star(p):
    """Significance stars. Returns '' for non-significant or NaN."""
    if p is None or not np.isfinite(p):
        return ""
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return ""

# ============================================================================
# STEP 6: BUILD THE FIGURE
# ============================================================================

CONTRAST_COLORS = {
    "treat_N2": "#4A90C2",  # blue   - treatment effect in WT
    "geno_unt": "#5DAA5D",  # green  - baseline genotype
    "geno_trt": "#E07B3D",  # orange - genotype under stress (GxE)
}
CONTRAST_ORDER = ["treat_N2", "geno_unt", "geno_trt"]
CONTRAST_LABELS = {c["short"]: c["name"] for c in CONTRASTS}

COMPARISON_CSV = "qPCR_vs_RNAseq_comparison.csv"


def _load_comparison():
    """Load qPCR vs RNA-seq comparison table if it exists. Returns a DataFrame
    with only the rows where both platforms have a log2FC, or None if the
    file is missing/empty."""
    # Search canonical locations (10_qPCR_compare.R writes to outputs/qpcr/)
    candidates = [
        os.path.join("outputs", "qpcr", COMPARISON_CSV),
        os.path.join("data", "qpcr", COMPARISON_CSV),
        COMPARISON_CSV,
    ]
    csv_path = next((p for p in candidates if os.path.exists(p)), None)
    if csv_path is None:
        return None
    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        print(f"[warn] could not read {csv_path}: {e}")
        return None
    if df.empty:
        return None
    df = df[df["qpcr_log2fc"].notna() & df["rnaseq_log2fc"].notna()].copy()
    if df.empty:
        return None
    return df


def _draw_panel_c(ax, comp_df):
    """Render the qPCR vs RNA-seq scatter comparison panel onto ax."""
    from scipy import stats as _stats

    # Plot bounds - symmetric around zero with a little padding
    all_vals = np.concatenate([
        comp_df["qpcr_log2fc"].values,
        comp_df["rnaseq_log2fc"].values,
    ])
    amax = np.nanmax(np.abs(all_vals))
    if not np.isfinite(amax) or amax < 1:
        amax = 1
    amax = amax * 1.15

    # Diagonal reference line (perfect agreement, y=x)
    ax.plot([-amax, amax], [-amax, amax],
            color="#888", lw=0.9, ls="--", zorder=1,
            label="y = x")

    # Zero reference lines
    ax.axhline(0, color="#ccc", lw=0.6, zorder=0)
    ax.axvline(0, color="#ccc", lw=0.6, zorder=0)

    # Scatter by contrast - use the same colors as panels a/b
    for contrast in CONTRAST_ORDER:
        sub = comp_df[comp_df["contrast"] == contrast]
        if sub.empty:
            continue
        ax.errorbar(
            sub["rnaseq_log2fc"].values,
            sub["qpcr_log2fc"].values,
            yerr=sub["qpcr_sem"].fillna(0).values,
            fmt="o",
            markersize=7,
            color=CONTRAST_COLORS[contrast],
            markeredgecolor="#222",
            markeredgewidth=0.6,
            ecolor="#555",
            elinewidth=0.7,
            capsize=2,
            label=CONTRAST_LABELS[contrast].replace("\n", " "),
            zorder=3,
        )

    # Gene labels with collision avoidance (requires adjustText; falls back to fixed offset)
    label_thresh = 0.75
    texts = []
    for _, row in comp_df.iterrows():
        if (abs(row["qpcr_log2fc"]) >= label_thresh or
                abs(row["rnaseq_log2fc"]) >= label_thresh):
            t = ax.text(
                row["rnaseq_log2fc"], row["qpcr_log2fc"],
                f"$\\it{{{row['gene']}}}$",
                fontsize=7.5, color="#222", zorder=4,
            )
            texts.append(t)
    if texts:
        try:
            from adjustText import adjust_text
            adjust_text(
                texts, ax=ax,
                arrowprops=dict(arrowstyle="-", color="#aaa", lw=0.5),
                expand=(1.2, 1.4),
                force_text=(0.3, 0.5),
            )
        except ImportError:
            print("[warn] adjustText not installed — labels may overlap. "
                  "Fix with: pip install adjustText")

    # Spearman correlation across all points
    x = comp_df["rnaseq_log2fc"].values
    y = comp_df["qpcr_log2fc"].values
    if len(x) >= 3:
        r, p = _stats.spearmanr(x, y)
        n = len(x)
        txt = f"Spearman $\\rho$ = {r:.2f}\n$n$ = {n}, $p$ = {p:.1g}"
    else:
        txt = f"$n$ = {len(x)} (too few for $\\rho$)"

    ax.text(
        0.03, 0.97, txt,
        transform=ax.transAxes,
        fontsize=8.5, ha="left", va="top",
        bbox=dict(boxstyle="round,pad=0.35",
                  facecolor="white", edgecolor="#888", lw=0.5),
        zorder=5,
    )

    ax.set_xlim(-amax, amax)
    ax.set_ylim(-amax, amax)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"RNA-seq log$_2$(fold change)"
                  + "\n(DESeq2, shrunk)", fontsize=9.5)
    ax.set_ylabel(r"qPCR log$_2$(fold change)"
                  + r"   ($\Delta\Delta$C$_q$ normalized to $\it{tba\text{-}1}$)",
                  fontsize=9.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", length=3)
    ax.legend(loc="lower right", bbox_to_anchor=(0.99, 0.01),
              bbox_transform=ax.transAxes,
              fontsize=6.5, frameon=False,
              handletextpad=0.3, labelspacing=0.25)


def make_supplement_figure(comp_df, out_pdf, out_png):
    """Per-gene qPCR vs RNA-seq scatter panels for supplementary figure."""
    from scipy import stats as _stats

    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica","Arial","Liberation Sans","DejaVu Sans"],
        "font.size": 9,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    genes = [g for g in comp_df["gene"].unique() if
             comp_df[comp_df["gene"] == g].shape[0] >= 2]
    n = len(genes)
    ncols = min(n, 4)
    nrows = int(np.ceil(n / ncols))

    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(3.2 * ncols, 3.4 * nrows),
                             squeeze=False)

    panel_labels = "abcdefghij"

    for idx, gene in enumerate(genes):
        ax = axes[idx // ncols][idx % ncols]
        sub = comp_df[comp_df["gene"] == gene].copy()

        # Axis limits - symmetric, gene-specific
        all_vals = np.concatenate([sub["qpcr_log2fc"].values,
                                   sub["rnaseq_log2fc"].values])
        amax = max(np.nanmax(np.abs(all_vals)) * 1.3, 1.0)

        ax.plot([-amax, amax], [-amax, amax],
                color="#aaa", lw=0.8, ls="--", zorder=1)
        ax.axhline(0, color="#ddd", lw=0.5, zorder=0)
        ax.axvline(0, color="#ddd", lw=0.5, zorder=0)

        for contrast in CONTRAST_ORDER:
            row = sub[sub["contrast"] == contrast]
            if row.empty:
                continue
            ax.errorbar(
                row["rnaseq_log2fc"].values,
                row["qpcr_log2fc"].values,
                yerr=row["qpcr_sem"].fillna(0).values,
                fmt="o", markersize=8,
                color=CONTRAST_COLORS[contrast],
                markeredgecolor="#222", markeredgewidth=0.6,
                ecolor="#555", elinewidth=0.8, capsize=2.5,
                label=CONTRAST_LABELS[contrast].replace("\n", " "),
                zorder=3,
            )

        # Spearman on the 3 contrast points
        x = sub["rnaseq_log2fc"].values
        y = sub["qpcr_log2fc"].values
        if len(x) >= 3:
            rho, p = _stats.spearmanr(x, y)
            p_str = f"p={p:.3f}" if p >= 0.001 else "p<0.001"
            stat_txt = f"Spearman ρ = {rho:.2f}\n{p_str}  (n={len(x)})"
        else:
            stat_txt = f"n={len(x)} (insufficient for ρ)"

        ax.text(0.05, 0.97, stat_txt,
                transform=ax.transAxes, fontsize=7,
                ha="left", va="top",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                          edgecolor="#bbb", lw=0.5))

        ax.set_xlim(-amax, amax)
        ax.set_ylim(-amax, amax)
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(f"$\\it{{{gene}}}$", fontsize=10, pad=4)
        ax.set_xlabel(r"RNA-seq log$_2$FC (DESeq2)", fontsize=7.5)
        ax.set_ylabel(r"qPCR log$_2$FC ($\Delta\Delta$Cq)", fontsize=7.5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(axis="both", length=3, labelsize=7.5)

        # Panel label (a, b, c…)
        ax.text(-0.14, 1.08, panel_labels[idx], transform=ax.transAxes,
                fontsize=13, fontweight="bold", ha="left", va="top")

    # Hide any unused axes
    for idx in range(n, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    # Shared legend below panels
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels,
               loc="lower center",
               ncol=len(CONTRAST_ORDER),
               fontsize=8, frameon=False,
               bbox_to_anchor=(0.5, -0.04),
               handletextpad=0.4, columnspacing=1.2)

    fig.suptitle("Supplementary Figure — Per-gene qPCR vs RNA-seq validation",
                 fontsize=10, fontweight="bold", y=1.01)

    plt.tight_layout(rect=[0, 0.06, 1, 1])
    plt.savefig(out_pdf, dpi=400, bbox_inches="tight")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"[done] saved {out_pdf}")
    print(f"[done] saved {out_png}")
    plt.close()


def make_figure(agg_df, ddcq_df, out_pdf, out_png, subtitle=None):
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica","Arial","Liberation Sans","DejaVu Sans"],
        "font.size": 10,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    comp_df = _load_comparison()
    have_panel_c = comp_df is not None
    if have_panel_c:
        # Filter to only genes in this figure's subset, drop missing RNA-seq values
        gene_order_set = set(agg_df["target"].unique())
        comp_df = comp_df[comp_df["gene"].isin(gene_order_set)].dropna(
            subset=["rnaseq_log2fc", "qpcr_log2fc"]
        ).copy()
        have_panel_c = len(comp_df) >= 3
        if have_panel_c:
            print(f"[info] panel C: {len(comp_df)} points for genes {sorted(gene_order_set)}")
        else:
            print(f"[info] panel C: too few points after filtering — skipping")
    else:
        print(f"[info] no {COMPARISON_CSV} found - skipping panel C")

    # Order genes by mean absolute log2FC in the key G x E contrast (geno_trt),
    # descending - most mutant-specific responders at the top
    geno_trt = agg_df[agg_df["contrast"] == "geno_trt"].copy()
    if geno_trt.empty:
        gene_order = sorted(agg_df["target"].unique())
    else:
        geno_trt["abs_fc"] = np.abs(geno_trt["log2fc_mean"])
        gene_order = geno_trt.sort_values("abs_fc", ascending=False)["target"].tolist()
        missing = [g for g in sorted(agg_df["target"].unique()) if g not in gene_order]
        gene_order = gene_order + missing

    n_genes = len(gene_order)
    n_contrasts = len(CONTRAST_ORDER)

    # Figure layout:
    #   Without panel C: 2 stacked rows [bar, dual-heatmap]
    #   With panel C:    left column = [bar, dual-heatmap], right column = [scatter]
    if have_panel_c:
        bar_width_in = max(7.5, 0.75 * n_genes + 1.8)
        fig_w = bar_width_in + 4.8
        fig_h = 8.2
        fig = plt.figure(figsize=(fig_w, fig_h))
        gs = fig.add_gridspec(
            2, 2,
            width_ratios=[bar_width_in, 4.8],
            height_ratios=[4.5, 2.0],
            wspace=0.22, hspace=0.45,
            left=0.065, right=0.975, top=0.945, bottom=0.11,
        )
        ax_bar = fig.add_subplot(gs[0, 0])
        # Split bottom-left cell into qPCR heatmap | RNA-seq heatmap
        gs_hm  = gs[1, 0].subgridspec(1, 2, width_ratios=[1, 1], wspace=0.08)
        ax_hm  = fig.add_subplot(gs_hm[0, 0])   # qPCR
        ax_hm2 = fig.add_subplot(gs_hm[0, 1])   # RNA-seq
        ax_sc  = fig.add_subplot(gs[:, 1])
    else:
        fig = plt.figure(figsize=(max(9, 0.85 * n_genes + 2.5), 8.2))
        gs = fig.add_gridspec(
            2, 1,
            height_ratios=[4.5, 2.0],
            hspace=0.45,
            left=0.085, right=0.955, top=0.945, bottom=0.11,
        )
        ax_bar = fig.add_subplot(gs[0])
        gs_hm  = gs[1].subgridspec(1, 2, width_ratios=[1, 1], wspace=0.08)
        ax_hm  = fig.add_subplot(gs_hm[0, 0])
        ax_hm2 = fig.add_subplot(gs_hm[0, 1])
        ax_sc  = None

    # -------- Panel A: grouped bar chart --------
    bar_w = 0.26
    x_base = np.arange(n_genes)

    # Jitter width for per-run overlay dots
    jitter_w = bar_w * 0.3

    for i, contrast in enumerate(CONTRAST_ORDER):
        offset = (i - 1) * bar_w
        means, sems, ps, ns = [], [], [], []
        for gene in gene_order:
            row = agg_df[(agg_df["target"] == gene) &
                         (agg_df["contrast"] == contrast)]
            if row.empty:
                means.append(np.nan); sems.append(np.nan); ps.append(np.nan); ns.append(0)
            else:
                means.append(row["log2fc_mean"].iloc[0])
                sems.append(row["log2fc_sem"].iloc[0])
                ps.append(row["padj_BH"].iloc[0] if "padj_BH" in row.columns
                          else row["p_value"].iloc[0])
                ns.append(int(row["n_runs"].iloc[0]))

        xs = x_base + offset
        means_arr = np.array(means)
        sems_arr  = np.array(sems)

        bars = ax_bar.bar(
            xs, means_arr,
            width=bar_w,
            color=CONTRAST_COLORS[contrast],
            edgecolor="#222",
            lw=0.6,
            label=CONTRAST_LABELS[contrast].replace("\n"," "),
            zorder=2,
        )

        # Error bars (only where SEM is finite)
        finite_mask = np.isfinite(sems_arr) & np.isfinite(means_arr)
        if finite_mask.any():
            ax_bar.errorbar(
                xs[finite_mask], means_arr[finite_mask],
                yerr=sems_arr[finite_mask],
                fmt="none", ecolor="#222", elinewidth=0.8, capsize=2.5, zorder=3,
            )

        # Per-run dots overlaid on each bar (IQR outliers hidden)
        for xi, gene in zip(xs, gene_order):
            per_run = ddcq_df[(ddcq_df["target"] == gene) &
                              (ddcq_df["contrast"] == contrast)]
            if per_run.empty:
                continue
            vals = per_run["log2fc"].values
            if len(vals) >= 4:
                q1, q3 = np.percentile(vals, [25, 75])
                iqr = q3 - q1
                keep_mask = (vals >= q1 - DOT_OUTLIER_IQR * iqr) & \
                            (vals <= q3 + DOT_OUTLIER_IQR * iqr)
            else:
                keep_mask = np.ones(len(vals), dtype=bool)
            rng = np.random.default_rng(hash(gene + contrast) & 0xFFFFFFFF)
            jitters = rng.uniform(-jitter_w, jitter_w, size=keep_mask.sum())
            ax_bar.scatter(
                xi + jitters,
                vals[keep_mask],
                s=12, c="white", edgecolors="#222", linewidths=0.5,
                zorder=4,
            )

        # Significance stars above bars
        for xi, m, p in zip(xs, means_arr, ps):
            if not np.isfinite(m):
                continue
            s = star(p)
            if s:
                y_off = 0.15 if m >= 0 else -0.35
                ax_bar.text(
                    xi, m + y_off, s,
                    ha="center", va="bottom" if m >= 0 else "top",
                    fontsize=10, fontweight="bold", color="#222",
                )

    ax_bar.axhline(0, color="#222", lw=0.8, zorder=1)
    ax_bar.set_xticks(x_base)
    ax_bar.set_xticklabels(gene_order, rotation=40, ha="right", fontsize=10,
                           style="italic")
    ax_bar.set_ylabel(r"log$_2$(fold change)   [qPCR, $\Delta\Delta$C$_q$ normalized to "
                      + r"$\it{tba\text{-}1}$]",
                      fontsize=10)
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)
    ax_bar.tick_params(axis="y", length=3)
    ax_bar.tick_params(axis="x", length=3, pad=2)

    # Legend
    leg = ax_bar.legend(
        loc="lower left",
        bbox_to_anchor=(0.01, 0.01),
        bbox_transform=ax_bar.transAxes,
        borderaxespad=0.3,
        fontsize=9,
        frameon=True,
        framealpha=0.85,
        edgecolor="#ccc",
        title="Contrast",
        title_fontsize=9.5,
    )
    leg.get_title().set_fontweight("bold")

    # Panel A label
    ax_bar.text(-0.06, 1.04, "a", transform=ax_bar.transAxes,
                fontsize=18, fontweight="bold",
                ha="left", va="top")

    # -------- Panel B: qPCR heatmap (left) | RNA-seq heatmap (right) --------

    # Build qPCR matrix from agg_df.
    # Note: p-value matrices are intentionally omitted here — _draw_heatmap()
    # annotates cells with numeric log2FC only; significance is already shown
    # in panel A via stars above bars. If heatmap star annotations are added
    # in future, re-introduce pv_qpcr/pv_rna and pass them into _draw_heatmap.
    hm_qpcr = np.full((n_contrasts, n_genes), np.nan)
    for i, contrast in enumerate(CONTRAST_ORDER):
        for j, gene in enumerate(gene_order):
            row = agg_df[(agg_df["target"] == gene) &
                         (agg_df["contrast"] == contrast)]
            if not row.empty:
                hm_qpcr[i, j] = row["log2fc_mean"].iloc[0]

    # Build RNA-seq matrix from comp_df
    hm_rna = np.full((n_contrasts, n_genes), np.nan)
    if comp_df is not None:
        for i, contrast in enumerate(CONTRAST_ORDER):
            for j, gene in enumerate(gene_order):
                match = comp_df[(comp_df["gene"] == gene) &
                                (comp_df["contrast"] == contrast)]
                if not match.empty and not np.isnan(match["rnaseq_log2fc"].iloc[0]):
                    hm_rna[i, j] = match["rnaseq_log2fc"].iloc[0]

    # Shared color scale across both heatmaps
    all_vals = np.concatenate([hm_qpcr.ravel(), hm_rna.ravel()])
    vmax = np.nanmax(np.abs(all_vals))
    vmax = max(vmax, 1.0)
    vmax = min(vmax, 6.0)

    def _draw_heatmap(ax, matrix, title, show_ylabels):
        im = ax.imshow(
            matrix, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
            aspect="auto", interpolation="nearest",
        )
        ax.set_xticks(np.arange(n_genes))
        ax.set_xticklabels(gene_order, rotation=40, ha="right",
                           fontsize=9, style="italic")
        ax.set_yticks(np.arange(n_contrasts))
        if show_ylabels:
            ax.set_yticklabels([CONTRAST_LABELS[c] for c in CONTRAST_ORDER],
                               fontsize=8.5, linespacing=1.1)
        else:
            ax.set_yticklabels([])
        ax.tick_params(axis="both", length=0, pad=3)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_title(title, fontsize=9, fontweight="bold", pad=4)
        # Annotate cells: log2FC value only
        for i in range(n_contrasts):
            for j in range(n_genes):
                v = matrix[i, j]
                if np.isnan(v):
                    continue
                txt_color = "white" if abs(v) > vmax * 0.5 else "#222"
                ax.text(j, i, f"{v:.1f}", ha="center", va="center",
                        fontsize=7, color=txt_color)
        return im

    im_qpcr = _draw_heatmap(ax_hm,  hm_qpcr, "qPCR (\u0394\u0394Cq)",     show_ylabels=True)
    im_rna  = _draw_heatmap(ax_hm2, hm_rna,  "RNA-seq (DESeq2)", show_ylabels=False)

    # Shared colorbar below the two heatmaps
    cbar_ref = im_qpcr
    if have_panel_c:
        cbar_ax = fig.add_axes([0.065, 0.03, 0.30, 0.013])
        cb = fig.colorbar(cbar_ref, cax=cbar_ax, orientation="horizontal")
    else:
        cbar_ax = fig.add_axes([0.085, 0.03, 0.40, 0.013])
        cb = fig.colorbar(cbar_ref, cax=cbar_ax, orientation="horizontal")
    cb.set_label(r"log$_2$ FC", fontsize=8.5, labelpad=3)
    cb.ax.tick_params(labelsize=7.5, length=2)
    cb.outline.set_linewidth(0.5)

    # Panel B label
    ax_hm.text(-0.10, 1.28, "b", transform=ax_hm.transAxes,
               fontsize=18, fontweight="bold",
               ha="left", va="top")

    # -------- Panel C: qPCR vs RNA-seq scatter (optional) --------
    if have_panel_c:
        _draw_panel_c(ax_sc, comp_df)
        ax_sc.text(-0.12, 1.04, "c", transform=ax_sc.transAxes,
                   fontsize=18, fontweight="bold",
                   ha="left", va="top")

    if subtitle:
        fig.suptitle(subtitle, fontsize=11, fontweight="bold",
                     x=0.5, y=0.995, ha="center", va="top")

    plt.savefig(out_pdf, dpi=400, bbox_inches="tight")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"[done] saved {out_pdf}")
    print(f"[done] saved {out_png}")

# ============================================================================
# MAIN
# ============================================================================

def main():
    # ── Fast-path: load from pre-computed CSVs if they exist ──────────────────
    # Search canonical locations (data/qpcr/ ships with the repo; outputs/qpcr/
    # is where 09 writes when run from raw; repo root is a legacy fallback).
    def _find_qpcr_csv(filename):
        for d in ["data/qpcr", "outputs/qpcr", "."]:
            p = Path(d) / filename
            if p.exists():
                return p
        return Path(filename)  # default (will fail .exists() check below)

    agg_path  = _find_qpcr_csv("qPCR_combined_deltadeltaCq.csv")
    ddcq_path = _find_qpcr_csv("qPCR_per_run_ddcq.csv")

    if agg_path.exists() and ddcq_path.exists():
        print(f"[step 1-5] Pre-computed results found — skipping raw-data pipeline.")
        print(f"  agg_path  = {agg_path}")
        print(f"  ddcq_path = {ddcq_path}")
        agg_df  = pd.read_csv(agg_path)
        ddcq_df = pd.read_csv(ddcq_path)
        print(f"  Loaded {len(agg_df)} rows from {agg_path}")
        print(f"  Loaded {len(ddcq_df)} rows from {ddcq_path}")
        # Fix any mislabeled target names and persist the correction to disk
        # so that the R comparison script reads canonical names (e.g. "t22f3.11"
        # not "t22f3.3"). Without re-saving, R would read the old names and fail
        # to resolve those genes to cosmid IDs in the count matrix.
        if TARGET_RENAMES:
            agg_df["target"]  = agg_df["target"].replace(TARGET_RENAMES)
            ddcq_df["target"] = ddcq_df["target"].replace(TARGET_RENAMES)
            agg_df.to_csv(agg_path, index=False)
            ddcq_df.to_csv(ddcq_path, index=False)
            print(f"  Renamed targets and re-saved CSVs: {TARGET_RENAMES}")
    else:
        print("[step 1] Loading Cq files from all runs...")
        cq_df = load_all_cq(QPCR_ROOT)
        print(f"  Total wells: {len(cq_df)}")
        if len(cq_df) == 0 or "run" not in cq_df.columns:
            print()
            print("=" * 64)
            print("No raw Bio-Rad Cq files found under data/qpcr/. Skipping 13a.")
            print(f"  Searched: {QPCR_ROOT}")
            print("  Expected: per-run subdirectories with Cq*.csv files")
            print("  (the deposited pre-computed qPCR CSVs in data/qpcr/ are")
            print("   sufficient for downstream comparison via 10_qPCR_compare.R)")
            print("=" * 64)
            return
        print(f"  Runs:        {cq_df['run'].nunique()}")
        print(f"  Targets:     {sorted(cq_df['target'].unique())}")
        print(f"  Strains:     {sorted(cq_df['strain'].unique())}")

        print("\n[step 1b] Deduplicating runs (email-forwarded copies)...")
        cq_df = deduplicate_runs(cq_df)
        print(f"  Runs after dedup: {cq_df['run'].nunique()}")

        print("\n[step 2] Flagging technical outlier wells...")
        cq_df = flag_outliers(cq_df)
        n_dropped = (~cq_df["keep"]).sum()
        print(f"  Wells flagged and dropped: {n_dropped}")

        print("\n[step 3] Computing DeltaCq (target - tba-1)...")
        dcq_df = compute_delta_cq(cq_df, ref_gene=REF_GENE)
        print(f"  DeltaCq rows: {len(dcq_df)}")
        if dcq_df.empty:
            print("ERROR: no DeltaCq rows produced - check REF_GENE and parse_sample()")
            return

        print("\n[step 4] Computing DeltaDeltaCq contrasts per run...")
        ddcq_df = compute_ddcq_contrasts(dcq_df)
        print(f"  DeltaDeltaCq rows: {len(ddcq_df)}")
        ddcq_save = ddcq_df.drop(columns=["dcq_num_vals","dcq_den_vals"])

        print("\n[step 5] Aggregating across runs...")
        agg_df = aggregate_across_runs(ddcq_df)
        agg_df["padj_BH"] = bh_adjust(agg_df["p_value"].values)
        agg_df = agg_df.sort_values(["contrast","target"]).reset_index(drop=True)

        # Apply renames BEFORE saving so the on-disk CSVs use the canonical
        # names. The R comparison script reads these files directly and must
        # see "t22f3.11" (not "t22f3.3") to resolve the cosmid correctly.
        if TARGET_RENAMES:
            ddcq_save["target"] = ddcq_save["target"].replace(TARGET_RENAMES)
            agg_df["target"]    = agg_df["target"].replace(TARGET_RENAMES)
            print(f"  Renamed targets before saving: {TARGET_RENAMES}")

        ddcq_save.to_csv("qPCR_per_run_ddcq.csv", index=False)
        print(f"  Wrote qPCR_per_run_ddcq.csv")
        agg_df.to_csv("qPCR_combined_deltadeltaCq.csv", index=False)
        print(f"  Wrote qPCR_combined_deltadeltaCq.csv (all rows)")

    # Figure subset: require >= MIN_RUNS_FOR_FIG biological runs per row.
    # Also require: each gene must have >=MIN_RUNS for the KEY G x E contrast
    # (geno_trt), otherwise it has no anchor for ordering or interpretation.
    agg_fig = agg_df[agg_df["n_runs"] >= MIN_RUNS_FOR_FIG].copy()
    keyc = "geno_trt"
    genes_with_key = set(agg_fig[agg_fig["contrast"] == keyc]["target"])
    agg_fig = agg_fig[agg_fig["target"].isin(genes_with_key)]

    # Also filter per-run dots to same gene set
    ddcq_fig = ddcq_df[ddcq_df["target"].isin(genes_with_key)].copy()

    n_dropped_genes = (agg_df["target"].nunique() - len(genes_with_key))
    print(f"\n  For the figure: {len(genes_with_key)} genes with >={MIN_RUNS_FOR_FIG} "
          f"runs in the key G x E contrast")
    print(f"  ({n_dropped_genes} genes relegated to supplementary - "
          f"insufficient runs)")

    print("\nPer-gene/contrast summary (figure subset):")
    disp = agg_fig[["target","contrast","log2fc_mean","log2fc_sem",
                    "n_runs","p_value","padj_BH"]].copy()
    disp["log2fc_mean"] = disp["log2fc_mean"].round(2)
    disp["log2fc_sem"]  = disp["log2fc_sem"].round(2)
    disp["p_value"]     = disp["p_value"].apply(
        lambda x: f"{x:.2g}" if np.isfinite(x) else "")
    disp["padj_BH"]     = disp["padj_BH"].apply(
        lambda x: f"{x:.2g}" if np.isfinite(x) else "")
    print(disp.to_string(index=False))

    # Apply exclusions
    agg_fig  = agg_fig[~agg_fig["target"].isin(GENES_EXCLUDE)].copy()
    ddcq_fig = ddcq_fig[~ddcq_fig["target"].isin(GENES_EXCLUDE)].copy()
    if GENES_EXCLUDE:
        print(f"\n  Excluded genes: {GENES_EXCLUDE}")

    print("\n[step 6] Building figure...")
    import os
    _qpcr_out = "outputs/qpcr"
    os.makedirs(_qpcr_out, exist_ok=True)
    make_figure(agg_fig, ddcq_fig,
                os.path.join(_qpcr_out, "qPCR_figure_panel.pdf"),
                os.path.join(_qpcr_out, "qPCR_figure_panel.png"))

    print("\n[step 7] Building supplementary per-gene correlation figure...")
    comp_all = _load_comparison()
    if comp_all is not None:
        make_supplement_figure(comp_all,
                               os.path.join(_qpcr_out, "qPCR_supplement_pergene.pdf"),
                               os.path.join(_qpcr_out, "qPCR_supplement_pergene.png"))
    else:
        print("  [skip] no comparison CSV found")

if __name__ == "__main__":
    main()
