#!/usr/bin/env python3
"""
Figure_GRN_publishable: predicted TF-centered regulatory architecture
for the 90 G×E genes from the argk-2 × CoCl2 dataset.

Scope of this analysis
----------------------
This is a STRUCTURAL / ARCHITECTURAL GRN, not a quantitative-predictive
model. Specifically:

    * 13 TFs (and 2 upstream kinases shown as glyphs) are pre-selected
      from prior biological knowledge and partitioned by hand into four
      modules (M0–M3).
    * TF -> effector edges are taken directly from TFLink — they are
      "predicted" only in the sense that they are database-curated
      regulatory relationships. No edges are inferred here.
    * Each G×E effector is assigned to its 'primary module' by simple
      majority vote over its TFLink-supported regulating TFs (alphabetical
      tie-break). No probabilistic / Bayesian inference is performed.
    * Panel B visualises every TF -> effector edge so multi-module
      convergence is readable directly off the wiring.
    * The convergence statistic printed under Panel B quantifies how many
      G×E targets sit at multi-module convergence points.

The convergence pattern is consistent with non-additive G×E behaviour
emerging from coordinated baseline shifts across several TF programs
that share targets. That is a MECHANISTIC HYPOTHESIS, not a tested
quantitative prediction. No held-out validation, ML inference, or
null-permutation test against shuffled TF -> module labels is performed
in this script.
"""
from __future__ import annotations
import sys
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyArrowPatch
import numpy as np
import pandas as pd

def _resolve_input(filename, search_paths=None):
    """Find an input TSV/CSV at one of several canonical locations."""
    from pathlib import Path as _P
    if search_paths is None:
        search_paths = [".", "outputs/GRN_tables", "../outputs/GRN_tables"]
    for sp in search_paths:
        cand = _P(sp) / filename
        if cand.exists():
            return str(cand)
    return None

NODES_FILE = _resolve_input("nodes_full90.tsv") or "nodes_full90.tsv"
EDGES_FILE = _resolve_input("edges_combined_full90.tsv") or "edges_combined_full90.tsv"
import os
OUT_DIR = "outputs/figures/Figure6"
os.makedirs(OUT_DIR, exist_ok=True)
OUT_PDF = os.path.join(OUT_DIR, "Figure_GRN_publishable_v2.pdf")
OUT_PNG = os.path.join(OUT_DIR, "Figure_GRN_publishable_v2.png")
HEATMAP_VMAX = 4.0

# Panel A: only label the top N effectors per module (sorted by |lfc_gxe|).
# All dots are still drawn — this preserves the quantitative density pattern
# while keeping labels readable. Set to None to label everything.
LABEL_TOP_N_PER_MODULE = 8

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 10,
    "axes.linewidth": 0.8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.facecolor": "white",
    "figure.facecolor": "white",
})

MODULES = {
    # PQM-1 binds DAE elements (IIS/metabolic), CEH-60 is a PBX developmental
    # regulator, ELT-2 is the intestinal GATA master — all three set baseline
    # transcriptional tone across growth, metabolic, and stress-response genes.
    "M0": {"tag": "developmental & metabolic program", "color": "#2B7A9E",
           "members": ["pqm-1", "ceh-60", "elt-2"]},
    # FOS-1 here: OR=1.83, padj=0.020 — more enriched than any M0 TF.
    # AP-1 is an acute xenobiotic / oxidative stress responder, not a baseline
    # intestinal regulator.
    "M1": {"tag": "xenobiotic / stress response", "color": "#D96A2F",
           "members": ["zip-2", "skn-1", "cebp-1", "fos-1"]},
    # BLMP-1 is G×E-depleted (OR=0.49) — noted in panel C but NOT visually
    # distinguished in the network (its TFLink regulatory edges are real;
    # depletion reflects target-set breadth, not absence of involvement).
    "M2": {"tag": "lipid / developmental",       "color": "#4F9E5B",
           "members": ["nhr-28", "nhr-77", "blmp-1"]},
    # AKT-2 and AAK-2: Ser/Thr kinases, not DNA-binding TFs → see
    # M3_SIGNALING_KINASES below. DAF-16/PHA-4 have large target sets that
    # include 49/43 G×E genes; Fisher OR≈1 reflects target-set breadth, not
    # absence of regulatory involvement.
    "M3": {"tag": "IIS / FOXO axis",             "color": "#6E4F9E",
           "members": ["daf-16", "pha-4", "hif-1"]},
}

# AKT-2 and AAK-2: upstream Ser/Thr kinases that phosphorylate DAF-16 and
# retain it in the cytoplasm. They are NOT transcription factors — they do not
# bind DNA. Shown as diamond glyphs inside M3 with no edges to effectors.
M3_SIGNALING_KINASES = ["akt-2", "aak-2"]

MODULE_CENTERS = {
    "M1": np.array([-3.0, 0.0]),
    "M0": np.array([ 0.0, 2.0]),
    "M2": np.array([ 0.0,-2.0]),
    "M3": np.array([ 3.0, 0.0]),
}

LOCAL_COORDS = {
    "M0": {"pqm-1": (-0.82,0.08),"ceh-60": (-0.20,-0.38),"fos-1": (0.40,0.12),"elt-2": (0.12,-0.74)},
    "M1": {"zip-2": (-0.74,-0.70),"skn-1": (-0.84,0.10),"cebp-1": (0.14,-0.72)},
    "M2": {"nhr-77": (-0.48,0.72),"nhr-28": (0.38,0.72),"blmp-1": (0.0,-0.12)},
    "M3": {"hif-1": (-0.06,0.24),"daf-16": (0.40,0.40),"akt-2": (0.80,0.04),"pha-4": (0.28,-0.44),"aak-2": (0.72,-0.74)},
}

# Wider shelves for all unlabeled GxE effectors
B_EFFECTOR_SLOTS = {
    "M0": [(-1.12, 1.02), (-0.64, 1.18), (-0.14, 1.20), (0.36, 1.18), (0.86, 1.00),
           (-1.06, 0.62), (-0.54, 0.74), (-0.02, 0.80), (0.50, 0.72), (0.98, 0.58)],
    "M1": [(-2.18, 0.78), (-1.72, 0.42), (-2.10, 0.06), (-1.58, -0.26), (-2.16, -0.62),
           (-1.66, -0.92), (-2.34, 0.38), (-2.34, -0.20), (-2.36, -0.80)],
    "M2": [(-0.96, -0.82), (-0.46, -0.98), (0.00, -1.06), (0.46, -0.98), (0.96, -0.82),
           (-0.92, -1.36), (-0.40, -1.52), (0.12, -1.58), (0.62, -1.48), (1.02, -1.28)],
    "M3": [(1.06, 0.86), (1.58, 0.50), (2.06, 0.12), (1.26, -0.28), (1.78, -0.64),
           (1.36, -0.98), (2.22, 0.54), (2.28, -0.12), (2.18, -0.72)],
}

# ---- Geometry constants ------------------------------------------------
HALO_RADIUS   = 1.42   # shared halo radius for panels A and B
TF_CIRCLE_R   = 0.17   # radius of TF circle on the perimeter
EFF_DOT_R     = 0.076  # radius of effector dot (both panels)
_PHI          = (1 + 5**0.5) / 2   # golden ratio for Fibonacci packing


# ---- Interior effector packing ----------------------------------------
def _fibonacci_layout(n, max_r):
    """Pack n items inside a circle of radius max_r using Fibonacci spiral.
    Returns a list of 2-d numpy offsets from the circle centre."""
    if n == 0:  return []
    if n == 1:  return [np.array([0.0, 0.0])]
    return [
        max_r * np.sqrt(i / n) * np.array([
            np.cos(2 * np.pi * i / _PHI**2),
            np.sin(2 * np.pi * i / _PHI**2),
        ])
        for i in range(n)
    ]


def build_effector_positions_inside(chosen_effectors, target_to_tfs,
                                    halo_radius=HALO_RADIUS):
    """Compute interior positions for effectors inside their primary halo.
    Used for BOTH panel A (labeled) and panel B (unlabeled dots), so the
    effectors occupy the same spots in both panels for visual continuity.
    """
    # Group by primary module, sorted so multi-module effectors go first
    grouped = {M: [] for M in MODULES}
    for eff in chosen_effectors:
        pm = primary_module(eff, target_to_tfs) or "M0"
        n_mods = len({tf_mod(t) for t in target_to_tfs.get(eff, set())
                      if tf_mod(t)})
        n_tfs  = len(target_to_tfs.get(eff, set()))
        grouped[pm].append((eff, n_mods, n_tfs))

    eff_pos = {}
    for M, items in grouped.items():
        ctr   = MODULE_CENTERS[M]
        items_sorted = [e for e, _, _ in
                        sorted(items, key=lambda x: (-x[1], -x[2], x[0]))]
        # Pack within 72 % of halo radius so they sit clearly inside the ring
        spots = _fibonacci_layout(len(items_sorted), halo_radius * 0.72)
        for eff, spot in zip(items_sorted, spots):
            eff_pos[eff] = ctr + spot
    return eff_pos


# ---- TF perimeter placement -------------------------------------------
def build_tf_perimeter_positions(halo_radius=HALO_RADIUS):
    """Place TFs in an outward-facing arc on each module's circumference."""
    tf_pos  = {}
    ARC_HALF = np.pi * 0.44
    for M, minfo in MODULES.items():
        ctr = MODULE_CENTERS[M]
        tfs = minfo["members"]
        n   = len(tfs)
        outward = (np.arctan2(ctr[1], ctr[0])
                   if np.hypot(*ctr) > 0.01 else np.pi / 2)
        if n == 1:
            angles = [outward]
        else:
            angles = [outward - ARC_HALF + 2 * ARC_HALF * i / (n - 1)
                      for i in range(n)]
        for tf, angle in zip(tfs, angles):
            tf_pos[tf] = ctr + halo_radius * np.array([np.cos(angle),
                                                        np.sin(angle)])
    return tf_pos


def build_kinase_positions(halo_radius=HALO_RADIUS):
    """Place AKT-2 and AAK-2 as diamond nodes on the inward side of M3's halo
    (the side facing the network centre) so they are visually distinct from
    the outward-arc TFs and look upstream of the DAF-16/PHA-4 nodes."""
    ctr  = MODULE_CENTERS["M3"]
    inward = np.arctan2(ctr[1], ctr[0]) + np.pi  # opposite to outward
    positions = {}
    for i, k in enumerate(M3_SIGNALING_KINASES):
        angle = inward + (i - 0.5) * 0.55          # small angular offset each
        positions[k] = ctr + halo_radius * 0.70 * np.array(   # inside halo
            [np.cos(angle), np.sin(angle)])
    return positions


    """Place TFs in an outward-facing arc on each module's circumference.

    TFs are distributed within ±80° of the outward direction (away from
    the network centre). This prevents any TF from being placed on the
    inward-facing side of a halo, which previously put M2 TFs in the
    centre of the network where they collided with other modules' content.
    """
    tf_pos = {}
    ARC_HALF = np.pi * 0.44   # ±79°  — outward arc only
    for M, minfo in MODULES.items():
        ctr = MODULE_CENTERS[M]
        tfs = minfo["members"]
        n   = len(tfs)
        outward = (np.arctan2(ctr[1], ctr[0])
                   if np.hypot(*ctr) > 0.01 else np.pi / 2)
        if n == 1:
            angles = [outward]
        else:
            angles = [outward - ARC_HALF + 2 * ARC_HALF * i / (n - 1)
                      for i in range(n)]
        for tf, angle in zip(tfs, angles):
            tf_pos[tf] = ctr + halo_radius * np.array([np.cos(angle),
                                                        np.sin(angle)])
    return tf_pos


# ---- Label alignment helper -------------------------------------------
def _label_align(pos, ctr):
    """(ha, va) for a label that should float radially OUTWARD from ctr."""
    angle = np.arctan2(pos[1] - ctr[1], pos[0] - ctr[0])
    if   -np.pi/4  <= angle <  np.pi/4:   return "left",   "center"
    elif  np.pi/4  <= angle <  3*np.pi/4: return "center", "bottom"
    elif -3*np.pi/4 <= angle < -np.pi/4:  return "center", "top"
    else:                                   return "right",  "center"

def read_table(path: str) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        sys.exit(f"[error] missing {path}")
    try:
        df = pd.read_csv(p, sep="\t")
        if df.shape[1] == 1:
            df = pd.read_csv(p)
    except Exception:
        df = pd.read_csv(p)
    return df

def standardize_nodes(nodes: pd.DataFrame) -> pd.DataFrame:
    rename_map = {}
    id_col = next((c for c in ["node_id","name","id","node","label"] if c in nodes.columns), None)
    if id_col is None:
        sys.exit(f"[error] Could not find node id column. Columns: {list(nodes.columns)}")
    if id_col != "node_id":
        rename_map[id_col] = "node_id"
    disp_col = next((c for c in ["display_label","node_label","pretty_label","label","name"] if c in nodes.columns and c != id_col), None)
    if disp_col and disp_col != "display_label":
        rename_map[disp_col] = "display_label"
    nodes = nodes.rename(columns=rename_map)
    if "display_label" not in nodes.columns:
        nodes["display_label"] = nodes["node_id"]
    if "effect_category" not in nodes.columns:
        nodes["effect_category"] = ""
    if "lfc_gxe" not in nodes.columns:
        nodes["lfc_gxe"] = 0.0
    nodes["node_id"] = nodes["node_id"].astype(str).str.strip().str.lower()
    nodes["display_label"] = nodes["display_label"].astype(str).str.strip()
    nodes["effect_category"] = nodes["effect_category"].astype(str).fillna("")
    nodes["lfc_gxe"] = pd.to_numeric(nodes["lfc_gxe"], errors="coerce").fillna(0.0)
    if "node_class" not in nodes.columns:
        regs = {r for m in MODULES.values() for r in m["members"]}
        nodes["node_class"] = np.where(nodes["node_id"].isin(regs), "TF", "effector")
    return nodes

def standardize_edges(edges: pd.DataFrame) -> pd.DataFrame:
    rename_map = {}
    if "source" not in edges.columns:
        for c in ["from","src"]:
            if c in edges.columns:
                rename_map[c] = "source"; break
    if "target" not in edges.columns:
        for c in ["to","dst"]:
            if c in edges.columns:
                rename_map[c] = "target"; break
    edges = edges.rename(columns=rename_map)
    if "source" not in edges.columns or "target" not in edges.columns:
        sys.exit(f"[error] Could not find source/target columns. Columns: {list(edges.columns)}")
    if "edge_type" not in edges.columns:
        edges["edge_type"] = "regulatory"
    edges["source"] = edges["source"].astype(str).str.strip().str.lower()
    edges["target"] = edges["target"].astype(str).str.strip().str.lower()
    edges["edge_type"] = edges["edge_type"].astype(str).fillna("regulatory")
    return edges

def load_data():
    return standardize_nodes(read_table(NODES_FILE)), standardize_edges(read_table(EDGES_FILE))

def tf_mod(tf: str):
    for M, info in MODULES.items():
        if tf in info["members"]:
            return M
    return None

def regulatory_edges_only(edges):
    keep = edges["edge_type"].str.contains("regulatory", case=False, na=False)
    return edges.loc[keep, ["source","target","edge_type"]].drop_duplicates().copy()

def build_tf_target_maps(reg_edges):
    regs = {r for m in MODULES.values() for r in m["members"]}
    target_to_tfs = {}
    for _, r in reg_edges.iterrows():
        a, b = r["source"], r["target"]
        a_reg, b_reg = a in regs, b in regs
        if a_reg and not b_reg:
            src, tgt = a, b
        elif b_reg and not a_reg:
            src, tgt = b, a
        else:
            continue
        target_to_tfs.setdefault(tgt, set()).add(src)
    return target_to_tfs

def choose_effectors(nodes, target_to_tfs, min_inputs=2):
    info = nodes.set_index("node_id")
    rows = []
    for eff, tfs in target_to_tfs.items():
        if eff not in info.index or info.loc[eff, "node_class"] != "effector":
            continue
        if len(tfs) < min_inputs:
            continue
        effect_cat = str(info.loc[eff, "effect_category"])
        n_mods = len({tf_mod(tf) for tf in tfs if tf_mod(tf) is not None})
        lfc = float(info.loc[eff, "lfc_gxe"])
        rows.append((eff, effect_cat == "GxE", n_mods, len(tfs), abs(lfc)))
    rows = sorted(rows, key=lambda x: (-int(x[1]), -x[2], -x[3], -x[4], x[0]))
    return [r[0] for r in rows]

def primary_module(eff, target_to_tfs):
    counts = {}
    for tf in target_to_tfs.get(eff, set()):
        mod = tf_mod(tf)
        if mod is not None:
            counts[mod] = counts.get(mod, 0) + 1
    return max(counts.items(), key=lambda kv: (kv[1], kv[0]))[0] if counts else None


def compute_convergence_stats(chosen_effectors, target_to_tfs):
    """Quantify multi-module convergence at G×E targets.

    The mechanistic claim of the figure is that small baseline shifts across
    several TF programs can compound at shared targets. That claim is only
    architecturally meaningful if many G×E targets actually receive TFLink
    input from more than one module. This function produces the headline
    number rendered under Panel B and printed to stdout.

    Returns
    -------
    dict with keys:
        total           – G×E targets with at least one in-set TF input
        by_n            – {n_modules: count} histogram, n in 1..4
        multi_module    – count receiving input from >= 2 modules
        three_plus      – count receiving input from >= 3 modules
        four            – count receiving input from all 4 modules
    """
    by_n = {1: 0, 2: 0, 3: 0, 4: 0}
    total = 0
    for eff in chosen_effectors:
        mods = {tf_mod(tf) for tf in target_to_tfs.get(eff, set())
                if tf_mod(tf) is not None}
        n = len(mods)
        if n >= 1:
            total += 1
            by_n[min(n, 4)] += 1
    multi_module = by_n[2] + by_n[3] + by_n[4]
    three_plus   = by_n[3] + by_n[4]
    return {
        "total": total,
        "by_n": by_n,
        "multi_module": multi_module,
        "three_plus": three_plus,
        "four": by_n[4],
    }


def build_tf_positions():
    pos = {}
    for M in MODULES:
        ctr = MODULE_CENTERS[M]
        for tf, local in LOCAL_COORDS[M].items():
            pos[tf] = ctr + np.array(local)
    return pos

def build_effector_positions(chosen_effectors, target_to_tfs):
    grouped = {M: [] for M in MODULES}
    for eff in chosen_effectors:
        pm = primary_module(eff, target_to_tfs) or "M0"
        tfs = target_to_tfs.get(eff, set())
        grouped[pm].append((eff, len({tf_mod(tf) for tf in tfs if tf_mod(tf) is not None}), len(tfs)))
    eff_pos = {}
    for M in grouped:
        ranked = sorted(grouped[M], key=lambda x: (-x[1], -x[2], x[0]))
        slots = B_EFFECTOR_SLOTS[M]
        for i, (eff, _, _) in enumerate(ranked):
            if i < len(slots):
                eff_pos[eff] = np.array(slots[i])
            else:
                base = np.array(slots[-1])
                eff_pos[eff] = base + np.array([0.0, -0.26*(i-len(slots)+1)])
    return eff_pos

def color_for_lfc(val):
    cmap = plt.cm.RdBu_r
    val = float(np.clip(val, -HEATMAP_VMAX, HEATMAP_VMAX))
    return cmap((val + HEATMAP_VMAX) / (2 * HEATMAP_VMAX))

def draw_module_halo(ax, center, radius, color, M):
    ax.add_patch(Circle(center, radius=radius, facecolor=color,
                        edgecolor=color, alpha=0.09, lw=1.0, zorder=0))
    # Chip label points OUTWARD from the network origin so M2 (at the
    # bottom) gets a label below its halo rather than above it (which
    # would land in the middle of the network).
    ctr = np.asarray(center, dtype=float)
    if np.hypot(*ctr) > 0.01:
        outward = ctr / np.hypot(*ctr)
    else:
        outward = np.array([0.0, 1.0])
    lp = ctr + outward * (radius + 0.18)
    angle = np.arctan2(outward[1], outward[0])
    if   abs(angle) < np.pi/4:         ha, va = "left",   "center"
    elif abs(angle) > 3 * np.pi / 4:   ha, va = "right",  "center"
    elif angle > 0:                     ha, va = "center", "bottom"
    else:                               ha, va = "center", "top"
    ax.text(lp[0], lp[1], M, ha=ha, va=va,
            fontsize=10.2, fontweight="bold", color=color,
            bbox=dict(boxstyle="round,pad=0.06", fc="white", ec=color,
                      lw=0.5, alpha=0.95), zorder=10)

def draw_tf_perimeter(ax, p, module_ctr, color, label):
    """TF circle on halo edge with outward-facing label.
    All TF circles use the same solid style — the dashed treatment was
    removed because non-significant Fisher enrichment reflects target-set
    breadth, not absence of regulatory involvement."""
    ax.add_patch(Circle(p, radius=TF_CIRCLE_R, facecolor="white",
                        edgecolor=color, lw=1.4, zorder=6))
    d   = p - module_ctr
    d_n = d / (np.linalg.norm(d) + 1e-9)
    lp  = p + d_n * (TF_CIRCLE_R + 0.15)
    ha, va = _label_align(p, module_ctr)
    ax.text(lp[0], lp[1], label, ha=ha, va=va, fontsize=8.0,
            bbox=dict(boxstyle="round,pad=0.05", fc="white", ec="none",
                      alpha=0.96),
            zorder=9)


def draw_kinase_node(ax, p, module_ctr, color, label):
    """Diamond glyph for upstream signaling kinases (AKT-2, AAK-2).
    These are Ser/Thr kinases, not DNA-binding TFs — no regulatory edges
    to effectors are drawn from them. The diamond visually distinguishes
    them from circular TF nodes."""
    from matplotlib.patches import RegularPolygon
    ax.add_patch(RegularPolygon(p, numVertices=4, radius=TF_CIRCLE_R * 1.1,
                                orientation=np.pi / 4,   # ◇ orientation
                                facecolor="white", edgecolor=color,
                                lw=1.2, zorder=6))
    d   = p - module_ctr
    d_n = d / (np.linalg.norm(d) + 1e-9)
    lp  = p + d_n * (TF_CIRCLE_R + 0.18)
    ha, va = _label_align(p, module_ctr)
    ax.text(lp[0], lp[1], label, ha=ha, va=va, fontsize=7.6,
            style="italic",
            bbox=dict(boxstyle="round,pad=0.04", fc="white", ec="none",
                      alpha=0.96),
            zorder=9)


def draw_eff_dot(ax, p, fill, size=EFF_DOT_R):
    """Small coloured dot — used in both panels (labelled in A, plain in B)."""
    ax.add_patch(Circle(p, radius=size, facecolor=fill,
                        edgecolor="#555", lw=0.60, zorder=4))

def draw_eff_labeled(ax, p, fill, label, module_ctr=None):
    """Effector dot with gene-name label pointing radially outward from the
    module centre — avoids the central crowding that occurs when all labels
    sit on top of their dots."""
    draw_eff_dot(ax, p, fill)
    if module_ctr is not None:
        d   = p - module_ctr
        dn  = d / (np.linalg.norm(d) + 1e-9)
        lp  = p + dn * (EFF_DOT_R + 0.10)
        ha, va = _label_align(p, module_ctr)
    else:
        lp  = p + np.array([0.0, EFF_DOT_R + 0.08])
        ha, va = "center", "bottom"
    ax.text(lp[0], lp[1], label, ha=ha, va=va, fontsize=5.8,
            bbox=dict(boxstyle="round,pad=0.03", fc="white", ec="none",
                      alpha=0.90),
            zorder=6)

def draw_edge(ax, p1, p2, color, lw=0.60, alpha=0.18):
    ax.plot([p1[0], p2[0]], [p1[1], p2[1]],
            color=color, lw=lw, alpha=alpha, zorder=1)

def init_axis(ax, title=None, xlim=(-4.8, 4.8), ylim=(-3.85, 3.85)):
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal")
    ax.axis("off")
    if title:
        # Title at top of axes in axes-fraction coords so it never collides
        # with data-space content (the M0 chip floats at the top of the
        # data area and previously overlapped a title placed in data coords).
        ax.set_title(title, fontsize=11.5, fontweight="bold",
                     loc="left", pad=4)

# ============================================================
# PANEL A — effectors inside their module halos (no TFs)
# ============================================================
def panel_A(ax, nodes, eff_pos_inside, chosen_effectors, target_to_tfs):
    """All G×E dots drawn for density; only top LABEL_TOP_N_PER_MODULE
    per module (by |lfc_gxe|) receive a text label.  This keeps M0
    readable even when it contains 50+ effectors."""
    init_axis(ax, "A. G×E targets organised into four modules")
    info = nodes.set_index("node_id")

    # ---- decide which effectors get labels ----
    label_set = set()
    if LABEL_TOP_N_PER_MODULE is None:
        label_set = set(chosen_effectors)
    else:
        by_mod = {}
        for eff in chosen_effectors:
            pm = primary_module(eff, target_to_tfs)
            if pm:
                by_mod.setdefault(pm, []).append(eff)
        for effs in by_mod.values():
            ranked = sorted(
                effs,
                key=lambda e: abs(float(info.loc[e, "lfc_gxe"])
                                  if e in info.index else 0),
                reverse=True,
            )
            label_set.update(ranked[:LABEL_TOP_N_PER_MODULE])

    for M, m in MODULES.items():
        draw_module_halo(ax, MODULE_CENTERS[M], radius=HALO_RADIUS,
                         color=m["color"], M=M)

    texts, xs, ys = [], [], []
    for eff in chosen_effectors:
        if eff not in eff_pos_inside:
            continue
        lfc  = float(info.loc[eff, "lfc_gxe"]) if eff in info.index else 0.0
        fill = color_for_lfc(lfc)
        p    = eff_pos_inside[eff]
        # Labelled dots get a slightly larger circle so adjustText arrow
        # lands visibly on the dot.
        r    = EFF_DOT_R * 1.15 if eff in label_set else EFF_DOT_R
        draw_eff_dot(ax, p, fill, size=r)

        if eff in label_set:
            xs.append(p[0]); ys.append(p[1])
            label = info.loc[eff, "display_label"] if eff in info.index else eff
            t = ax.text(p[0], p[1] + r + 0.04, label,
                        ha="center", va="bottom", fontsize=5.8,
                        bbox=dict(boxstyle="round,pad=0.03", fc="white",
                                  ec="none", alpha=0.93),
                        zorder=6)
            texts.append(t)

    if texts:
        try:
            from adjustText import adjust_text
            adjust_text(
                texts, x=xs, y=ys, ax=ax,
                expand=(1.2, 1.4),
                force_text=(0.25, 0.35),
                force_static=(0.15, 0.20),
                only_move={"text": "xy", "static": "xy"},
                arrowprops=dict(arrowstyle="-", color="#999999",
                                lw=0.30, alpha=0.75),
            )
        except Exception:
            pass

    cax  = ax.inset_axes([0.80, 0.02, 0.14, 0.022])
    norm = mpl.colors.Normalize(vmin=-HEATMAP_VMAX, vmax=HEATMAP_VMAX)
    cb   = mpl.colorbar.ColorbarBase(cax, cmap=plt.cm.RdBu_r, norm=norm,
                                     orientation='horizontal')
    cb.ax.tick_params(labelsize=6.4, length=2)
    cb.set_label("log$_2$FC (G×E)", fontsize=6.6, labelpad=1)


# ============================================================
# PANEL B — TFs on module perimeter, effectors inside (unlabelled)
# ============================================================
def panel_B(ax, nodes, tf_perim_pos, eff_pos_inside,
            target_to_tfs, chosen_effectors, conv_stats=None):
    """Add the TF regulatory layer. TF circles sit on the circumference of
    each module halo so they can be identified without crowding the interior.
    Edges run from each TF on the perimeter to the effectors it regulates
    inside the halo. Effector dots are the same interior positions as panel A
    but without labels — the reader already knows which gene is which.

    If conv_stats is provided (from compute_convergence_stats), a small
    annotation in the lower-left of the panel reports the convergence
    headline number — i.e. the empirical count of G×E targets receiving
    TFLink input from >=2 modules. This is the structural fact behind the
    'compounding at shared targets' mechanism."""
    # Wider x-limits to give room for TF labels outside the M1 / M3 halos
    init_axis(ax, "B. TF regulatory layer", xlim=(-5.4, 5.4), ylim=(-4.0, 4.0))
    info = nodes.set_index("node_id")

    # 1. Module halos
    for M, m in MODULES.items():
        draw_module_halo(ax, MODULE_CENTERS[M], radius=HALO_RADIUS,
                         color=m["color"], M=M)

    # 2. Edges: TF perimeter → effector interior
    for eff in chosen_effectors:
        if eff not in eff_pos_inside:
            continue
        for tf in sorted(target_to_tfs.get(eff, set())):
            M = tf_mod(tf)
            if M is None or tf not in tf_perim_pos:
                continue
            draw_edge(ax, tf_perim_pos[tf], eff_pos_inside[eff],
                      MODULES[M]["color"])

    # 3. Effector dots (unlabelled) — same positions as panel A
    for eff in chosen_effectors:
        if eff not in eff_pos_inside:
            continue
        lfc  = float(info.loc[eff, "lfc_gxe"]) if eff in info.index else 0.0
        draw_eff_dot(ax, eff_pos_inside[eff], color_for_lfc(lfc))

    # 4. TF circles — uniform solid style for all TFs
    for M, m in MODULES.items():
        ctr = MODULE_CENTERS[M]
        for tf in m["members"]:
            if tf not in tf_perim_pos:
                continue
            label = info.loc[tf, "display_label"] if tf in info.index else tf
            draw_tf_perimeter(ax, tf_perim_pos[tf], ctr, m["color"], label)

    # 5. Signaling kinases (AKT-2, AAK-2) as diamond glyphs inside M3 halo
    kinase_pos = build_kinase_positions(HALO_RADIUS)
    for k in M3_SIGNALING_KINASES:
        if k not in kinase_pos:
            continue
        label = info.loc[k, "display_label"] if k in info.index else k
        draw_kinase_node(ax, kinase_pos[k], MODULE_CENTERS["M3"],
                         MODULES["M3"]["color"], label)

    # Legend — upper-left avoids collision with M1/M2 perimeter labels
    handles = [
        Line2D([0],[0], marker='o', markersize=6.5,
               markerfacecolor='white', markeredgewidth=1.4,
               markeredgecolor='#555', linestyle='None',
               label='Transcription factor'),
        Line2D([0],[0], marker='D', markersize=5.5,
               markerfacecolor='white', markeredgewidth=1.2,
               markeredgecolor='#6E4F9E', linestyle='None',
               label='Signaling kinase'),
        Line2D([0],[0], marker='o', markersize=4.6,
               markerfacecolor='#bbbbbb', markeredgewidth=0.7,
               markeredgecolor='#555', linestyle='None',
               label='G×E effector'),
        Line2D([0,1],[0,0], color='#888', lw=0.8, label='TF → effector'),
    ]
    leg = ax.legend(handles=handles, loc="upper left",
                    bbox_to_anchor=(0.01, 0.99),
                    frameon=True, fontsize=7.6,
                    borderpad=0.36, handlelength=1.2)
    leg.get_frame().set_alpha(0.97)

    cax  = ax.inset_axes([0.82, 0.04, 0.11, 0.020])
    norm = mpl.colors.Normalize(vmin=-HEATMAP_VMAX, vmax=HEATMAP_VMAX)
    cb   = mpl.colorbar.ColorbarBase(cax, cmap=plt.cm.RdBu_r, norm=norm,
                                     orientation='horizontal')
    cb.ax.tick_params(labelsize=6.4, length=2)
    cb.set_label("log$_2$FC", fontsize=6.6, labelpad=1)

    # Convergence headline — placed in lower-left axes-fraction coords so it
    # sits opposite the colorbar (lower-right) and below the legend (upper-
    # left). Numbers are computed from the actual TFLink edge file at
    # runtime — not hard-coded.
    if conv_stats is not None:
        total = conv_stats["total"]
        multi = conv_stats["multi_module"]
        three = conv_stats["three_plus"]
        msg = (f"Multi-module convergence: "
               f"{multi}/{total} G×E targets receive TFLink input "
               f"from ≥2 modules ({three} from ≥3).")
        ax.text(0.015, 0.04, msg,
                transform=ax.transAxes,
                ha="left", va="bottom",
                fontsize=6.6, color="#222",
                bbox=dict(boxstyle="round,pad=0.20", fc="white",
                          ec="#888", lw=0.5, alpha=0.96),
                zorder=20)

def simple_box(ax, x, y, text, edgecolor, fontsize=8.4, pad=0.15):
    ax.text(x, y, text, ha="center", va="center", fontsize=fontsize,
            color="#111",
            bbox=dict(boxstyle=f"round,pad={pad}", fc="white",
                      ec=edgecolor, lw=0.78, alpha=0.96), zorder=20)

def panel_C(ax):
    """ARGK-2 setpoint model — clean schematic with no overlapping elements.

    Layout (data-unit coordinates, xlim=±4.8, ylim=-3.65 to 3.85):

         M0 (0, 2.5)
              ↓
    M1 (-3.1, 0.4)  ←→  core (0, 0.4)  ←→  M3 (3.1, 0.4)
              ↑
         M2 (0, -1.6)

    ────────────── separator ──────────────

    [N2 baseline]  →  [ARGK-2 loss]  →  [CoCl₂]      y = -2.75

    [ARGK-2 sets resting engagement …]                 y = -3.35
    """
    init_axis(ax, "C. ARGK-2 setpoint model", ylim=(-3.75, 3.85))

    # ---- Module halos (light background only, no chip — this is a schematic)
    for M, m in MODULES.items():
        ax.add_patch(Circle(MODULE_CENTERS[M], radius=1.08,
                            facecolor=m["color"], edgecolor=m["color"],
                            alpha=0.07, lw=0.7, zorder=0))

    # ---- Box positions (chosen to give ≥0.5 unit gap between every pair)
    m0_xy  = (0.0,  2.50)
    m1_xy  = (-3.1, 0.40)
    core   = (0.0,  0.40)
    m3_xy  = (3.1,  0.40)
    m2_xy  = (0.0, -1.60)

    simple_box(ax, *m0_xy,
               "M0\ndevelopmental &\nmetabolic program",
               MODULES["M0"]["color"], fontsize=7.6)
    simple_box(ax, *m1_xy,
               "M1\nxenobiotic /\nstress response",
               MODULES["M1"]["color"], fontsize=7.6)
    simple_box(ax, *core,
               "shared G×E\ncore",
               "#C77B42", fontsize=8.6, pad=0.14)
    simple_box(ax, *m3_xy,
               "M3\nIIS / FOXO axis",
               MODULES["M3"]["color"], fontsize=7.2)
    # M2: BLMP-1 is G×E-depleted (OR=0.49) — note explicitly so the
    # schematic is not misleading about M2's involvement.
    simple_box(ax, *m2_xy,
               "M2\nlipid / developmental\n(NHR-28, NHR-77, BLMP-1)",
               MODULES["M2"]["color"], fontsize=7.2)

    # ---- Arrows from each module box edge to the core box edge.
    # Endpoints are estimated from fontsize / number of lines so arrows
    # start/end at box edges, not box centres.
    # M0 (3 lines, ~0.50 half-height) → core (2 lines, ~0.38 half-height)
    ax.add_patch(FancyArrowPatch((0.0, 2.00), (0.0, 0.78),
                 arrowstyle='->', mutation_scale=9, lw=0.9,
                 color=MODULES["M0"]["color"], alpha=0.85, zorder=8))
    # M1 (2 lines, ~0.38 half-height; ~1.10 half-width) → core (~0.50 half-width)
    ax.add_patch(FancyArrowPatch((-2.00, 0.40), (-0.55, 0.40),
                 arrowstyle='->', mutation_scale=9, lw=0.9,
                 color=MODULES["M1"]["color"], alpha=0.85, zorder=8))
    # M3 (~0.95 half-width) → core
    ax.add_patch(FancyArrowPatch((2.15, 0.40), (0.55, 0.40),
                 arrowstyle='->', mutation_scale=9, lw=0.9,
                 color=MODULES["M3"]["color"], alpha=0.85, zorder=8))
    # M2 (3 lines, ~0.50 half-height) → core bottom
    ax.add_patch(FancyArrowPatch((0.0, -1.10), (0.0, -0.02),
                 arrowstyle='->', mutation_scale=9, lw=0.9,
                 color=MODULES["M2"]["color"], alpha=0.85, zorder=8))

    # ---- Separator
    ax.axhline(y=-2.12, xmin=0.04, xmax=0.96,
               color="#cccccc", lw=0.7, ls="--", zorder=1)

    # ---- Bottom strip: genotype progression
    y_strip = -2.75
    simple_box(ax, -3.25, y_strip,
               "N2 baseline\nbalanced / unprimed",
               "#AAAAAA", fontsize=7.6)
    simple_box(ax,  0.00, y_strip,
               "ARGK-2 loss\nbaseline pre-engagement",
               "#C77B42", fontsize=7.6)
    simple_box(ax,  3.25, y_strip,
               "CoCl$_2$\nreorganized response",
               "#AAAAAA", fontsize=7.6)
    # Arrows between strip boxes
    ax.add_patch(FancyArrowPatch((-2.18, y_strip), (-1.15, y_strip),
                 arrowstyle='->', mutation_scale=8, lw=0.8, color="#666"))
    ax.add_patch(FancyArrowPatch(( 1.15, y_strip), ( 2.18, y_strip),
                 arrowstyle='->', mutation_scale=8, lw=0.8, color="#666"))

    # ---- Bottom summary (framed as a hypothesis, not a result)
    ax.text(0.0, -3.42,
            "Proposed: ARGK-2 may set resting engagement of four "
            "convergent TF modules",
            ha="center", va="center", fontsize=6.8, color="#111",
            bbox=dict(boxstyle="round,pad=0.14", fc="white",
                      ec="#C77B42", lw=0.78, alpha=0.96), zorder=20)

def main():
    nodes, edges = load_data()
    reg_edges = regulatory_edges_only(edges)
    target_to_tfs = build_tf_target_maps(reg_edges)
    # min_inputs=1 so every G×E gene with at least one TFLink connection
    # appears in panel A, even those regulated by only one of the 15 TFs.
    chosen_effectors = choose_effectors(nodes, target_to_tfs, min_inputs=1)
    if not chosen_effectors:
        sys.exit("[error] No effectors passed the filter.")

    # Positions shared across panels A and B for visual continuity
    eff_pos_inside = build_effector_positions_inside(
        chosen_effectors, target_to_tfs, halo_radius=HALO_RADIUS)
    tf_perim_pos = build_tf_perimeter_positions(halo_radius=HALO_RADIUS)

    # Empirical multi-module convergence — the headline structural fact
    # under the convergence claim. Computed from the TFLink edges only.
    conv_stats = compute_convergence_stats(chosen_effectors, target_to_tfs)

    print(f"[ok] {len(chosen_effectors)} effectors across "
          f"{len(MODULES)} modules")
    for M in MODULES:
        n = sum(1 for e in chosen_effectors
                if eff_pos_inside.get(e) is not None
                and primary_module(e, target_to_tfs) == M)
        print(f"     {M}: {n} effectors, "
              f"{len(MODULES[M]['members'])} TFs on perimeter")
    print(f"[convergence] of {conv_stats['total']} G×E targets with "
          f"in-set TF input:")
    print(f"     1 module : {conv_stats['by_n'][1]}")
    print(f"     2 modules: {conv_stats['by_n'][2]}")
    print(f"     3 modules: {conv_stats['by_n'][3]}")
    print(f"     4 modules: {conv_stats['by_n'][4]}")
    print(f"     >=2 mods : {conv_stats['multi_module']} "
          f"({100*conv_stats['multi_module']/max(conv_stats['total'],1):.0f}%)")
    print(f"     >=3 mods : {conv_stats['three_plus']} "
          f"({100*conv_stats['three_plus']/max(conv_stats['total'],1):.0f}%)")

    fig = plt.figure(figsize=(7.1, 8.0), constrained_layout=False)
    gs  = fig.add_gridspec(3, 1, height_ratios=[0.92, 1.12, 0.94],
                           hspace=0.08, left=0.05, right=0.985,
                           top=0.985, bottom=0.03)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[1, 0])
    axC = fig.add_subplot(gs[2, 0])

    panel_A(axA, nodes, eff_pos_inside, chosen_effectors, target_to_tfs)
    panel_B(axB, nodes, tf_perim_pos, eff_pos_inside,
            target_to_tfs, chosen_effectors, conv_stats=conv_stats)
    panel_C(axC)

    plt.savefig(OUT_PDF, dpi=700, bbox_inches="tight")
    plt.savefig(OUT_PNG, dpi=700, bbox_inches="tight")
    print(f"[ok] wrote {OUT_PDF}")
    print(f"[ok] wrote {OUT_PNG}")

if __name__ == "__main__":
    main()
