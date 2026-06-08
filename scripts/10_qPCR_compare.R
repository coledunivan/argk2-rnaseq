#!/usr/bin/env Rscript
# ==============================================================================
# qPCR -> RNA-seq Comparison Panel (Panel C of qPCR figure)
# ==============================================================================
#
# This script regenerates the DESeq2 contrasts that match the three qPCR
# contrasts exactly, then pulls out the log2 fold-changes for every gene
# assayed by qPCR. The output is a small table that qpcr_analysis.py reads
# to render a side-by-side qPCR-vs-RNA-seq comparison panel.
#
# The three matched contrasts:
#
#   treat_N2  = treatment main effect at reference genotype
#             = log2(N2 treated / N2 untreated)
#             = DESeq2 coefficient "treatment_treated_vs_untreated"
#
#   geno_unt  = genotype main effect at untreated baseline
#             = log2(argk-2 untreated / N2 untreated)
#             = DESeq2 coefficient "genotype_RB2060_vs_N2"
#
#   geno_trt  = genotype effect UNDER TREATMENT (genotype + GxE interaction)
#             = log2(argk-2 treated / N2 treated)
#             = linear combination: genotype_RB2060_vs_N2 + interaction
#
# These are the SAME three contrasts the qPCR analysis computes, so the
# comparison is apples-to-apples (same numerator and denominator conditions
# on both platforms).
#
# INPUTS (place in working directory):
#   data/reference/RNAseq_to_TF_Targets.csv   - merged counts + TFLink edges
#   data/Sample_Metadata_Table.csv  - sample metadata
#   qPCR_combined_deltadeltaCq.csv    - produced by qpcr_analysis.py
#
# OUTPUTS:
#   qPCR_vs_RNAseq_comparison.csv     - joined table for plotting
#
# USAGE:
#   Rscript qpcr_rnaseq_comparison.R
#
# NOTES:
# 1. LFC shrinkage matters here - without it, low-count genes that the qPCR
#    measured cleanly will show inflated RNA-seq LFCs and the comparison will
#    look artificially worse than it is. apeglm is used for single-coefficient
#    contrasts; ashr is used for manual linear combinations.
#
# 2. qPCR target names are WormBase common names ("ugt-31"), while the count
#    matrix rows are cosmid IDs ("Y75B8A.23"). We resolve this by looking up
#    each qPCR gene in the Name.Target column of the TFLink file (which
#    provides cosmid <-> common name pairs) before pulling from DESeq2.
#
# 3. Some qPCR genes may not be in the RNA-seq count matrix at all (if they
#    were dropped during upstream filtering). These are flagged as NA in the
#    output rather than silently skipped.
# ==============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
})

# ---- 1. Load data ------------------------------------------------------------

cat("Loading data...\n")

merged <- read.csv("data/reference/RNAseq_to_TF_Targets.csv",
                   stringsAsFactors = FALSE, check.names = FALSE)
colnames(merged) <- trimws(colnames(merged))
merged <- merged[, colnames(merged) != ""]

meta <- read.csv("data/Sample_Metadata_Table.csv", stringsAsFactors = FALSE)
colnames(meta) <- trimws(colnames(meta))
meta$SampleLabel <- trimws(meta$SampleLabel)

qpcr_csv <- {
  candidates <- c("qPCR_combined_deltadeltaCq.csv",
                  "data/qpcr/qPCR_combined_deltadeltaCq.csv",
                  "outputs/qpcr/qPCR_combined_deltadeltaCq.csv")
  found <- candidates[file.exists(candidates)]
  if (length(found) == 0) {
    stop("Cannot find qPCR_combined_deltadeltaCq.csv. Searched: ",
         paste(candidates, collapse = ", "),
         ".\nRun scripts/09_qPCR_analysis.py first (needs raw Bio-Rad Cq files in data/qpcr/<run>/),\n",
         "or provide the pre-computed CSV at one of the above paths.")
  }
  found[1]
}
message("Reading qPCR data from: ", qpcr_csv)
if (!file.exists(qpcr_csv)) {
  stop("Cannot find ", qpcr_csv, " in working directory. ",
       "Run qpcr_analysis.py first to generate it.")
}
qpcr <- read.csv(qpcr_csv, stringsAsFactors = FALSE)

# Apply the same target name corrections as qpcr_analysis.py (TARGET_RENAMES).
# The on-disk CSV may have been generated before the rename was applied, so
# old plate names like "t22f3.3" must be corrected to "t22f3.11" here to
# ensure cosmid resolution succeeds.
TARGET_RENAMES <- c("t22f3.3" = "t22f3.11")
qpcr$target <- dplyr::recode(qpcr$target, !!!TARGET_RENAMES)
if (any(names(TARGET_RENAMES) %in% unique(qpcr$target))) {
  cat("  WARNING: some targets were still using old names after recode — check TARGET_RENAMES\n")
} else {
  cat("  Applied target renames:", paste(names(TARGET_RENAMES), "->",
      TARGET_RENAMES, collapse = ", "), "\n")
}

cat("  Merged file rows:", nrow(merged), "\n")
cat("  qPCR gene x contrast rows:", nrow(qpcr), "\n")

# ---- 2. Build count matrix (dedup by Gene, lossless) ------------------------

sample_cols <- grep("^A[0-9]+$", colnames(merged), value = TRUE)
counts_dedup <- merged[!duplicated(merged$Gene), ]
count_mat_raw <- as.matrix(counts_dedup[, sample_cols])
rownames(count_mat_raw) <- counts_dedup$Gene
storage.mode(count_mat_raw) <- "double"
cat("  Unique genes in count matrix:", nrow(count_mat_raw), "\n")

# Common-name -> cosmid lookup from TFLink Name.Target column
# Some genes have multiple cosmid IDs under one common name; take the first
# (which in this dataset is the correct/canonical one for all qPCR targets).
name_to_cosmid <- merged %>%
  filter(!is.na(Name.Target), Name.Target != "-", Name.Target != "") %>%
  mutate(Name.Target = tolower(trimws(Name.Target))) %>%
  distinct(Name.Target, Gene) %>%
  group_by(Name.Target) %>%
  dplyr::slice(1) %>%
  ungroup()
name_lookup <- setNames(name_to_cosmid$Gene, name_to_cosmid$Name.Target)

# ---- 3. Metadata subset + factors -------------------------------------------

meta_sub <- meta[meta$genotype %in% c("N2", "RB2060"), ]
meta_sub <- droplevels(meta_sub)
rownames(meta_sub) <- meta_sub$SampleLabel

missing <- setdiff(meta_sub$SampleLabel, colnames(count_mat_raw))
if (length(missing) > 0) {
  cat("WARNING: dropping samples missing from counts:", 
      paste(missing, collapse = ", "), "\n")
  meta_sub <- meta_sub[!meta_sub$SampleLabel %in% missing, ]
}

count_mat <- count_mat_raw[, meta_sub$SampleLabel]
count_mat <- round(count_mat)
storage.mode(count_mat) <- "integer"

meta_sub$genotype   <- factor(meta_sub$genotype,  levels = c("N2", "RB2060"))
meta_sub$treatment  <- factor(meta_sub$treatment, levels = c("untreated", "treated"))
meta_sub$experiment <- factor(meta_sub$experiment)

cat("\nDesign:\n")
print(with(meta_sub, table(genotype, treatment, experiment)))
cat("Total samples:", nrow(meta_sub), "\n\n")

# ---- 4. DESeq2 --------------------------------------------------------------

cat("Running DESeq2...\n")
dds <- tryCatch({
  d <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData   = meta_sub,
    design    = ~ genotype * treatment + genotype:experiment + treatment:experiment
  )
  DESeq(d)
}, error = function(e) {
  cat("Full model failed:", conditionMessage(e),
      "\n  Falling back to additive batch model\n")
  d <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData   = meta_sub,
    design    = ~ experiment + genotype * treatment
  )
  DESeq(d)
})
cat("Result names:", paste(resultsNames(dds), collapse = ", "), "\n\n")

# ---- 5. Extract the three qPCR-matched contrasts with LFC shrinkage ---------

have_apeglm <- requireNamespace("apeglm", quietly = TRUE)
have_ashr   <- requireNamespace("ashr",   quietly = TRUE)
if (!have_apeglm) cat("WARNING: apeglm not installed - falling back to unshrunk LFCs\n")
if (!have_ashr)   cat("WARNING: ashr not installed - falling back to unshrunk LFCs\n")

shrink <- function(coef_name = NULL, contrast_list = NULL, type = "apeglm") {
  pkg_ok <- (type == "apeglm" && have_apeglm) || (type == "ashr" && have_ashr)
  if (pkg_ok) {
    res <- tryCatch({
      if (!is.null(coef_name)) {
        lfcShrink(dds, coef = coef_name, type = type)
      } else {
        lfcShrink(dds, contrast = contrast_list, type = type)
      }
    }, error = function(e) {
      cat("  shrink failed (", type, "):", conditionMessage(e), 
          "- using unshrunk\n"); NULL
    })
    if (!is.null(res)) return(as.data.frame(res))
  }
  if (!is.null(coef_name)) {
    return(as.data.frame(results(dds, name = coef_name)))
  } else {
    return(as.data.frame(results(dds, contrast = contrast_list)))
  }
}

cat("Extracting qPCR-matched contrasts with shrinkage...\n")

# Verify that the coefficient names we rely on are actually present.
# If the full model fell back to the additive model, the coefficient names
# will still match (DESeq2 uses the same convention), but an unexpected
# model failure could produce different names silently.
rn <- resultsNames(dds)
required_coefs <- c("treatment_treated_vs_untreated",
                    "genotype_RB2060_vs_N2",
                    "genotypeRB2060.treatmenttreated")
missing_coefs <- setdiff(required_coefs, rn)
if (length(missing_coefs) > 0) {
  cat("FATAL: the following expected DESeq2 coefficients are missing:\n")
  cat("  ", paste(missing_coefs, collapse = ", "), "\n")
  cat("Available coefficients:\n  ", paste(rn, collapse = "\n   "), "\n")
  stop("Cannot extract qPCR-matched contrasts. Check model design and factor levels.")
}

# Contrast 1: treat_N2 = treatment main effect at N2 baseline
res_treat_N2 <- shrink(coef_name = "treatment_treated_vs_untreated",
                       type = "apeglm")

# Contrast 2: geno_unt = genotype main effect at untreated baseline  
res_geno_unt <- shrink(coef_name = "genotype_RB2060_vs_N2",
                       type = "apeglm")

# Contrast 3: geno_trt = genotype effect under treatment
#             = genotype_RB2060_vs_N2 + interaction term
#             Use ashr because this is a manual linear combination
res_geno_trt <- shrink(
  contrast_list = list(c("genotype_RB2060_vs_N2",
                         "genotypeRB2060.treatmenttreated")),
  type = "ashr")

# ---- 6. Resolve qPCR common names to cosmid IDs -----------------------------

qpcr_genes <- unique(qpcr$target)
cat("\nResolving", length(qpcr_genes), "qPCR gene names to cosmid IDs...\n")

resolve_gene <- function(g) {
  g_clean <- tolower(trimws(g))
  # Direct cosmid match first (in case the qPCR name IS the cosmid)
  if (g_clean %in% tolower(rownames(count_mat))) {
    idx <- which(tolower(rownames(count_mat)) == g_clean)[1]
    return(rownames(count_mat)[idx])
  }
  # Common name lookup via TFLink
  if (g_clean %in% names(name_lookup)) {
    cosmid <- name_lookup[[g_clean]]
    if (cosmid %in% rownames(count_mat)) return(cosmid)
  }
  NA_character_
}

qpcr_resolved <- data.frame(
  qpcr_name = qpcr_genes,
  cosmid    = sapply(qpcr_genes, resolve_gene),
  stringsAsFactors = FALSE
)
unresolved <- qpcr_resolved$qpcr_name[is.na(qpcr_resolved$cosmid)]
if (length(unresolved) > 0) {
  cat("  WARNING: could not resolve these qPCR gene names to count matrix:\n")
  cat("  ", paste(unresolved, collapse = ", "), "\n")
  cat("  (check spelling, or they may be genes filtered out upstream)\n")
}
resolved <- qpcr_resolved$qpcr_name[!is.na(qpcr_resolved$cosmid)]
cat("  Resolved:", length(resolved), "/", length(qpcr_genes), "\n")

# ---- 7. Build the joined output table ---------------------------------------
#
# For each qPCR gene x contrast row, add the matching RNA-seq LFC and padj
# from DESeq2. Rows where the gene cannot be resolved get NA values.

get_rnaseq <- function(cosmid, contrast) {
  if (is.na(cosmid)) return(c(lfc = NA_real_, padj = NA_real_, baseMean = NA_real_))
  res <- switch(contrast,
                "treat_N2" = res_treat_N2,
                "geno_unt" = res_geno_unt,
                "geno_trt" = res_geno_trt,
                NULL)
  if (is.null(res) || !(cosmid %in% rownames(res))) {
    return(c(lfc = NA_real_, padj = NA_real_, baseMean = NA_real_))
  }
  row <- res[cosmid, ]
  c(lfc      = as.numeric(row$log2FoldChange),
    padj     = as.numeric(row$padj),
    baseMean = as.numeric(row$baseMean))
}

out_rows <- vector("list", nrow(qpcr))
for (i in seq_len(nrow(qpcr))) {
  q <- qpcr[i, ]
  cosmid <- qpcr_resolved$cosmid[qpcr_resolved$qpcr_name == q$target][1]
  rs <- get_rnaseq(cosmid, q$contrast)
  out_rows[[i]] <- data.frame(
    gene             = q$target,
    cosmid           = ifelse(is.na(cosmid), "", cosmid),
    contrast         = q$contrast,
    qpcr_log2fc      = q$log2fc_mean,
    qpcr_sem         = q$log2fc_sem,
    qpcr_n_runs      = q$n_runs,
    qpcr_pvalue      = q$p_value,
    qpcr_padj_BH     = q$padj_BH,
    rnaseq_log2fc    = unname(rs["lfc"]),
    rnaseq_padj      = unname(rs["padj"]),
    rnaseq_baseMean  = unname(rs["baseMean"]),
    stringsAsFactors = FALSE
  )
}
out_df <- do.call(rbind, out_rows)

# Direction agreement flag - only meaningful where both platforms have a value
out_df$direction_agrees <- with(out_df,
  ifelse(is.na(qpcr_log2fc) | is.na(rnaseq_log2fc),
         NA,
         sign(qpcr_log2fc) == sign(rnaseq_log2fc)))

{
  out_dir <- "outputs/qpcr"
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  out_path <- file.path(out_dir, "qPCR_vs_RNAseq_comparison.csv")
  write.csv(out_df, out_path, row.names = FALSE)
  cat("Wrote: ", out_path, "\n", sep = "")
}

cat("\n========== SUMMARY ==========\n")
cat("\nDone.\n")
cat("  Total rows:", nrow(out_df), "\n")
cat("  Rows with both qPCR and RNA-seq values:",
    sum(!is.na(out_df$qpcr_log2fc) & !is.na(out_df$rnaseq_log2fc)), "\n")

comparable <- out_df[!is.na(out_df$qpcr_log2fc) & !is.na(out_df$rnaseq_log2fc), ]
if (nrow(comparable) >= 3) {
  cat("\nSpearman correlation (qPCR vs RNA-seq log2FC) by contrast:\n")
  cat("  [Spearman rho used here to match the figure statistic]\n")
  for (c in c("treat_N2", "geno_unt", "geno_trt")) {
    sub <- comparable[comparable$contrast == c, ]
    if (nrow(sub) >= 3) {
      rho <- cor(sub$qpcr_log2fc, sub$rnaseq_log2fc, method = "spearman")
      cat(sprintf("  %-10s n=%2d  rho=%+.3f\n", c, nrow(sub), rho))
    } else {
      cat(sprintf("  %-10s n=%2d  (too few for correlation)\n", c, nrow(sub)))
    }
  }
  all_rho <- cor(comparable$qpcr_log2fc, comparable$rnaseq_log2fc, method = "spearman")
  cat(sprintf("  %-10s n=%2d  rho=%+.3f\n", "OVERALL", nrow(comparable), all_rho))

  agrees <- sum(comparable$direction_agrees, na.rm = TRUE)
  cat(sprintf("\nDirection agreement: %d/%d rows agree on sign (%.0f%%)\n",
              agrees, nrow(comparable), 100 * agrees / nrow(comparable)))
}

cat("\nRun qpcr_analysis.py next to render the comparison panel.\n")
