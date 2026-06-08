#!/usr/bin/env bash
# =============================================================================
# run_all.sh -- argk-2 x CoCl2 RNA-seq pipeline driver
# =============================================================================
#
# Runs every script in dependency order from the repo root. Every script writes
# directly to its canonical outputs/<subdir>/ location, so a partial run still
# leaves a usable, organized output tree.
#
#   00 QC report                -- QC plots (PCA, distance heatmap, dispersion)
#   01 DESeq2 + TF enrichment   -- canonical model fit + 5 contrasts +
#                                  effect classification (90 GxE genes) +
#                                  TF GxE target enrichment (zip-2 OR=17.83)
#   02 Figure 1                 -- PCA, KO validation, volcanos, Venn
#   03 Figure 2 (KEGG)          -- KEGG functional composition (4-contrast)
#   04 Figure 3 (GO)            -- GO BP enrichment + immune effect map
#   05 Figure 4 (TF bubbles)    -- TF target activity across 4 contrasts
#   06 Figure 5 (GxE)           -- GxE scatter + heatmap + bars (90/68/22)
#   07 GRN inputs               -- nodes_full90.tsv + edges_combined_full90.tsv
#   08 Figure 6 (GRN)           -- TF regulatory network rendering
#   09 qPCR analysis            -- DeltaDeltaCq from Bio-Rad Cq files (skips if absent)
#   10 qPCR comparison          -- qPCR vs RNA-seq Spearman + figure
#   11 Supp GxE heatmap         -- full 90-gene GxE heatmap
#   12 Supp additive specificity-- Mann-Whitney p=9.2e-15
#   13 Supp GxE core analysis   -- in-degree histogram + TF Jaccard heatmap
# =============================================================================
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$REPO_ROOT"

mkdir -p outputs/logs outputs/QC outputs/DESeq2 outputs/GRN_tables \
         outputs/figures outputs/qpcr outputs/supplementary

run() {
  local label="$1" cmd="$2"
  echo
  echo "--- $label ---"
  if eval "$cmd" 2>&1 | tee "outputs/logs/${label// /_}.log"; then
    echo "[ok] $label"
  else
    echo "[FAIL] $label  -- see outputs/logs/${label// /_}.log"
    exit 1
  fi
}

echo "=== argk-2 x CoCl2 pipeline driver ==="
echo "Repo:   $REPO_ROOT"
echo "Date:   $(date -Iseconds)"
echo "R:      $(Rscript --version 2>&1 | head -1)"
echo "Python: $(python3 --version)"

run "00 QC report"                 "Rscript scripts/00_QC_report.R || true"
run "01 DESeq2 and TF enrichment"  "Rscript scripts/01_DESeq2_and_TF_enrichment.R"
run "02 Figure 1"                  "Rscript scripts/02_Figure1.R"
run "03 Figure 2 KEGG"             "Rscript scripts/03_Figure2_KEGG.R"
run "04 Figure 3 GO"               "Rscript scripts/04_Figure3_GO.R"
run "05 Figure 4 TF bubbles"       "Rscript scripts/05_Figure4_TF.R"
run "06 Figure 5 GxE"              "Rscript scripts/06_Figure5_GxE.R"
run "07 GRN inputs"                "Rscript scripts/07_GRN_inputs.R"
run "08 Figure 6 GRN"              "python3 scripts/08_Figure6_GRN.py"
run "10 qPCR comparison"           "Rscript scripts/10_qPCR_compare.R"
run "09 qPCR analysis"             "python3 scripts/09_qPCR_analysis.py"
run "11 Supp GxE heatmap"          "python3 scripts/11_Supp_GxE_heatmap.py"
run "12 Supp additive specificity" "python3 scripts/12_Supp_additive_specificity.py"
run "13 Supp GxE core analysis"    "python3 scripts/13_Supp_GxE_core_analysis.py"

echo
echo "=== Pipeline complete. ==="
echo
# R sometimes drops an empty Rplots.pdf at the repo root when a graphics
# device is left implicitly open by a package (most often pheatmap on macOS).
# It's harmless, but we clean it up so the repo root stays tidy.
rm -f Rplots.pdf

echo "Outputs in:"
echo "  outputs/QC/             QC plots (PCA, distance heatmap, dispersion)"
echo "  outputs/DESeq2/         contrast tables + effect classification + dds.RData"
echo "  outputs/GRN_tables/     nodes/edges TSVs for Figure 6"
echo "  outputs/figures/        publication-quality PDFs (main figures)"
echo "  outputs/supplementary/  supplementary figures + tables"
echo "  outputs/qpcr/           qPCR DeltaDeltaCq tables + figures"
echo "  outputs/logs/           per-step run logs"
echo
echo "Main figure PDFs:"
find outputs/figures -name '*.pdf' | sort
