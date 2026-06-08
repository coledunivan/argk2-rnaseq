#!/usr/bin/env Rscript
# =============================================================================
# requirements.R — install everything the pipeline needs
# =============================================================================
#
# Usage:
#   Rscript requirements.R
#
# Strategy:
#   1. Bioconductor packages via BiocManager.
#   2. CRAN packages via install.packages().
#   3. ashr + apeglm both come from Bioconductor.
# =============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

bioc_pkgs <- c(
  "DESeq2",         # core differential expression — pipeline pinned to v1.42
  "apeglm",         # LFC shrinkage for single coefficients
  "ashr",           # LFC shrinkage for linear-combination contrasts
  "org.Ce.eg.db",   # C. elegans annotation DB
  "GO.db",          # GO term hierarchy + names
  "AnnotationDbi",  # annotation interface
  "clusterProfiler",# KEGG/GO enrichment (OPTIONAL — offline helpers in _utils.R
                    # cover the same functionality)
  "biomaRt",        # OPTIONAL — Ensembl queries
  "limma"           # used in Figure 1 for batch-correction of PCA only
)

cran_pkgs <- c(
  "dplyr", "tidyr", "readr", "tibble", "stringr", "purrr", "ggplot2",
  "ggrepel", "patchwork", "scales",
  "VennDiagram",    # Figure 1 Venn (replaces ggVennDiagram)
  "pheatmap",       # not currently used by the deposited scripts but a
                    # frequent companion in RNA-seq pipelines
  "RColorBrewer",
  "data.table",
  "matrixStats"
)

cat("Installing Bioconductor packages...\n")
BiocManager::install(bioc_pkgs, ask = FALSE, update = FALSE)

cat("Installing CRAN packages...\n")
install.packages(cran_pkgs, repos = "https://cloud.r-project.org")

cat("\nVerifying installation...\n")
for (p in c(bioc_pkgs, cran_pkgs)) {
  ok <- requireNamespace(p, quietly = TRUE)
  cat(sprintf("  %-22s %s\n", p, if (ok) "OK" else "MISSING"))
}
cat("\nDone. If any packages show MISSING, see the project README.\n")
