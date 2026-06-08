#!/usr/bin/env Rscript

# DESeq2 Differential Expression for ARGK-2 Integrated GRN

#
# INPUTS (place in working directory):
#   - data/reference/Cleaned_Annotated_Counts.csv            (gene x sample count matrix)
#   - data/Sample_Metadata_Table.csv        (sample metadata)
#   - data/reference/RNAseq_to_TF_Targets.csv.gz       (TFLink TF-target pairs)
#
# MODEL: ~ genotype * treatment + genotype:experiment + treatment:experiment
#   Falls back to ~ experiment + genotype * treatment if rank-deficient
#
# OUTPUTS:
#   DESeq2_treatment_in_N2.csv           contrast results
#   DESeq2_treatment_in_RB2060.csv       contrast results
#   DESeq2_genotype_untreated.csv        contrast results
#   DESeq2_genotype_treated.csv          contrast results
#   DESeq2_GxE_interaction.csv           interaction term results
#   DESeq2_effect_classification.csv     per-gene effect class + LFCs + padj
#   DESeq2_TF_targets_DE.csv             DE TF targets with connectivity
#   DESeq2_TF_summary.csv                TF node label stats for GRN
#
# USAGE: Rscript DESeq2_GRN_analysis.R
# ==============================================================================

library(DESeq2)
library(dplyr)

# Output directory for all DESeq2 tables + dds object (so the repo root stays clean)
DESEQ_OUT <- "outputs/DESeq2"
if (!dir.exists(DESEQ_OUT)) dir.create(DESEQ_OUT, recursive = TRUE)

# ---- 1. Load data ------------------------------------------------------------
#
# The counts and the TFLink edges are merged into a single file:
# `data/reference/RNAseq_to_TF_Targets.csv`. It has one row per TF->target edge, and
# the same gene appears many times (once per edge). The count columns A1..A32
# are identical across all rows for a given gene, so we deduplicate by Gene.

cat("Loading data...\n")

merged <- read.csv("data/reference/RNAseq_to_TF_Targets.csv",
                   stringsAsFactors = FALSE, check.names = FALSE)
colnames(merged) <- trimws(colnames(merged))
# Drop any nameless leading index column from Excel/pandas export
merged <- merged[, colnames(merged) != ""]

cat("  Merged file rows:", nrow(merged), "\n")
cat("  Merged file columns:", ncol(merged), "\n")

# Identify sample columns (A followed by digits)
sample_cols <- grep("^A[0-9]+$", colnames(merged), value = TRUE)
cat("  Sample columns found:", paste(sample_cols, collapse = ", "), "\n")

# ---- 1a. Build the count matrix (deduplicated by Gene) ----------------------

# Take the first row per gene for the count matrix. Counts are identical
# across duplicate rows, so this is lossless for expression data.
counts_dedup <- merged[!duplicated(merged$Gene), ]
count_mat_raw <- as.matrix(counts_dedup[, sample_cols])
rownames(count_mat_raw) <- counts_dedup$Gene
storage.mode(count_mat_raw) <- "double"

cat("  Unique genes in count matrix:", nrow(count_mat_raw), "\n")

# Save a common-name map (Name.Target is TFLink's common name for the target)
# Take the first non-"-" common name per gene
name_map_df <- merged %>%
  filter(!is.na(Name.Target), Name.Target != "-", Name.Target != "") %>%
  distinct(Gene, Name.Target) %>%
  group_by(Gene) %>%
  dplyr::slice(1) %>%
  ungroup()
gene_symbols_lookup <- setNames(name_map_df$Name.Target, name_map_df$Gene)

# For the DESeq2 output, gene_symbol is just the Gene (cosmid) identifier;
# the downstream Python figure script will pull the display name from the
# Name.Target column of the TFLink file itself.
gene_symbols <- rownames(count_mat_raw)

# ---- 1b. Build the TF edge table --------------------------------------------

# Define TFs of interest
key_tfs <- c("zip-2","daf-16","pqm-1","ceh-60","cebp-1",
             "fos-1","skn-1","elt-2","nhr-28","nhr-77",
             "blmp-1","pha-4","sma-9","unc-62","sta-1")

tf_edges <- merged %>%
  filter(Name.TF %in% key_tfs) %>%
  select(TF_name = Name.TF,
         target_cosmid = Gene,
         target_common_name = Name.Target) %>%
  distinct()

cat("  Key TF-target edges:", nrow(tf_edges), "\n")
cat("  Unique target genes:", length(unique(tf_edges$target_cosmid)), "\n")

# ---- 1c. Load metadata ------------------------------------------------------

meta <- read.csv("data/Sample_Metadata_Table.csv", stringsAsFactors = FALSE)
colnames(meta) <- trimws(colnames(meta))
meta$SampleLabel <- trimws(meta$SampleLabel)

# ---- 2. Prepare count matrix -------------------------------------------------

# Subset metadata to N2 and RB2060 only, then drop unused factor levels
# so ghost levels (MAH172, RB2598) don't interfere with the design matrix
meta_sub <- meta[meta$genotype %in% c("N2", "RB2060"), ]
meta_sub <- droplevels(meta_sub)
rownames(meta_sub) <- meta_sub$SampleLabel

# Verify count matrix columns match metadata
missing_in_counts <- setdiff(meta_sub$SampleLabel, colnames(count_mat_raw))
if (length(missing_in_counts) > 0) {
  cat("WARNING: These samples are in metadata but missing from count matrix:\n")
  cat(paste(missing_in_counts, collapse = ", "), "\n")
  cat("\nDiagnostic - sample columns actually present in counts:\n")
  cat(paste(colnames(count_mat_raw), collapse = ", "), "\n")
  meta_sub <- meta_sub[!meta_sub$SampleLabel %in% missing_in_counts, ]
}

count_mat <- count_mat_raw[, meta_sub$SampleLabel]
count_mat <- round(count_mat)
storage.mode(count_mat) <- "integer"

# Set factor levels (N2 and untreated as reference)
meta_sub$genotype   <- factor(meta_sub$genotype,  levels = c("N2", "RB2060"))
meta_sub$treatment  <- factor(meta_sub$treatment, levels = c("untreated", "treated"))
meta_sub$experiment <- factor(meta_sub$experiment)

gene_symbols <- rownames(count_mat)

cat("\nDesign:\n")
print(with(meta_sub, table(genotype, treatment, experiment)))
cat("Total samples:", nrow(meta_sub), "\n")

# ---- 3. DESeq2 ---------------------------------------------------------------

# Option A: Full model with experiment interactions
cat("\n--- Trying full model: ~ genotype * treatment + genotype:experiment + treatment:experiment ---\n")

model_used <- "full"
tryCatch({
  dds <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData   = meta_sub,
    design    = ~ genotype * treatment + genotype:experiment + treatment:experiment
  )
  dds <- DESeq(dds)
  cat("Full model converged.\n")
}, error = function(e) {
  cat("Full model failed:", conditionMessage(e), "\n")
  cat("--- Falling back: ~ experiment + genotype * treatment ---\n")
  model_used <<- "additive_batch"
  dds <<- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData   = meta_sub,
    design    = ~ experiment + genotype * treatment
  )
  dds <<- DESeq(dds)
  cat("Additive batch model converged.\n")
})

cat("Model used:", model_used, "\n")
cat("Result names:", paste(resultsNames(dds), collapse = ", "), "\n\n")

# ---- 4. Extract contrasts (with apeglm/ashr LFC shrinkage) ------------------
#
# Shrinkage is essential here: the design is small (16 samples) and we have
# many low-count genes, so unshrunk LFCs explode for genes with high variance.
# - apeglm shrinkage works on single coefficients
# - ashr shrinkage works on linear combinations (manual contrasts)
# Both produce shrunk LFCs that are conservative and stable.

cat("Extracting contrasts with LFC shrinkage...\n")

# Make sure the shrinkage packages are available
have_apeglm <- requireNamespace("apeglm", quietly = TRUE)
have_ashr   <- requireNamespace("ashr",   quietly = TRUE)
if (!have_apeglm || !have_ashr) {
  stop("Required shrinkage packages missing: ",
       if (!have_apeglm) "apeglm " else "",
       if (!have_ashr)   "ashr"    else "",
       "\n  Falling back to unshrunk LFCs INFLATES effect_class=='GxE' counts ",
       "(e.g. 126 instead of the published 90) and silently breaks every ",
       "downstream figure that reads DESeq2_effect_classification.csv.\n",
       "  Install both before running this script:\n",
       "    BiocManager::install(c('apeglm','ashr'))")
}

shrink_or_unshrunk <- function(dds, coef_name = NULL, contrast_list = NULL,
                               type = "apeglm", alpha = 0.05) {
  # Try shrinkage first; if package missing or call fails, fall back to plain results()
  pkg_ok <- (type == "apeglm" && have_apeglm) ||
    (type == "ashr"   && have_ashr)
  if (pkg_ok) {
    res <- tryCatch({
      if (!is.null(coef_name)) {
        lfcShrink(dds, coef = coef_name, type = type)
      } else {
        lfcShrink(dds, contrast = contrast_list, type = type)
      }
    }, error = function(e) {
      cat("  shrink failed (", type, "):", conditionMessage(e),
          "- using unshrunk\n")
      NULL
    })
    if (!is.null(res)) return(as.data.frame(res))
  }
  # Fallback to unshrunk
  if (!is.null(coef_name)) {
    return(as.data.frame(results(dds, name = coef_name, alpha = alpha)))
  } else {
    return(as.data.frame(results(dds, contrast = contrast_list, alpha = alpha)))
  }
}

# Treatment effect in N2 (reference genotype)
# Coefficient name: treatment_treated_vs_untreated
res_treat_N2 <- shrink_or_unshrunk(dds,
                                   coef_name = "treatment_treated_vs_untreated",
                                   type = "apeglm")
res_treat_N2$gene_symbol <- gene_symbols

# Treatment effect in RB2060 = treatment + interaction (manual contrast - use ashr)
res_treat_RB <- shrink_or_unshrunk(
  dds,
  contrast_list = list(c("treatment_treated_vs_untreated",
                         "genotypeRB2060.treatmenttreated")),
  type = "ashr")
res_treat_RB$gene_symbol <- gene_symbols

# Genotype effect in untreated (reference treatment) - apeglm
res_geno_unt <- shrink_or_unshrunk(dds,
                                   coef_name = "genotype_RB2060_vs_N2",
                                   type = "apeglm")
res_geno_unt$gene_symbol <- gene_symbols

# Genotype effect in treated = genotype + interaction (manual contrast - ashr)
res_geno_trt <- shrink_or_unshrunk(
  dds,
  contrast_list = list(c("genotype_RB2060_vs_N2",
                         "genotypeRB2060.treatmenttreated")),
  type = "ashr")
res_geno_trt$gene_symbol <- gene_symbols

# G x E interaction term directly - apeglm
int_name <- grep("RB2060.*treated", resultsNames(dds), value = TRUE)[1]
res_interaction <- shrink_or_unshrunk(dds, coef_name = int_name, type = "apeglm")
res_interaction$gene_symbol <- gene_symbols

# Quick LFC distribution diagnostic so we can see if shrinkage worked
cat("\nLFC magnitude diagnostics (post-shrinkage):\n")
for (nm in c("res_treat_N2","res_treat_RB","res_geno_unt","res_geno_trt","res_interaction")) {
  obj <- get(nm)
  vals <- abs(obj$log2FoldChange)
  vals <- vals[!is.na(vals)]
  cat(sprintf("  %-18s  median=%.2f  90%%=%.2f  max=%.2f\n",
              nm, median(vals), quantile(vals, 0.9), max(vals)))
}

# ---- 5. Save contrast results ------------------------------------------------

write.csv(res_treat_N2, file.path(DESEQ_OUT, "DESeq2_treatment_in_N2.csv"),     row.names = TRUE)
write.csv(res_treat_RB, file.path(DESEQ_OUT, "DESeq2_treatment_in_RB2060.csv"), row.names = TRUE)
write.csv(res_geno_unt, file.path(DESEQ_OUT, "DESeq2_genotype_untreated.csv"),  row.names = TRUE)
write.csv(res_geno_trt, file.path(DESEQ_OUT, "DESeq2_genotype_treated.csv"),    row.names = TRUE)
write.csv(res_interaction, file.path(DESEQ_OUT, "DESeq2_GxE_interaction.csv"),     row.names = TRUE)

# ---- 6. Effect classification ------------------------------------------------

cat("Classifying genes...\n")

LFC_THRESH <- 1.0
PADJ_THRESH <- 0.05

classify <- data.frame(
  gene_symbol   = gene_symbols,
  lfc_treat_N2  = res_treat_N2$log2FoldChange,
  padj_treat_N2 = res_treat_N2$padj,
  lfc_treat_RB  = res_treat_RB$log2FoldChange,
  padj_treat_RB = res_treat_RB$padj,
  lfc_geno_unt  = res_geno_unt$log2FoldChange,
  padj_geno_unt = res_geno_unt$padj,
  lfc_geno_trt  = res_geno_trt$log2FoldChange,
  padj_geno_trt = res_geno_trt$padj,
  lfc_gxe       = res_interaction$log2FoldChange,
  padj_gxe      = res_interaction$padj,
  stringsAsFactors = FALSE
)

# Significance flags
sig_tn2 <- !is.na(classify$padj_treat_N2) & classify$padj_treat_N2 < PADJ_THRESH &
  abs(classify$lfc_treat_N2) > LFC_THRESH
sig_trb <- !is.na(classify$padj_treat_RB) & classify$padj_treat_RB < PADJ_THRESH &
  abs(classify$lfc_treat_RB) > LFC_THRESH
sig_gut <- !is.na(classify$padj_geno_unt) & classify$padj_geno_unt < PADJ_THRESH &
  abs(classify$lfc_geno_unt) > LFC_THRESH
sig_gtr <- !is.na(classify$padj_geno_trt) & classify$padj_geno_trt < PADJ_THRESH &
  abs(classify$lfc_geno_trt) > LFC_THRESH
sig_gxe <- !is.na(classify$padj_gxe) & classify$padj_gxe < PADJ_THRESH

sig_treat <- sig_tn2 | sig_trb
sig_geno  <- sig_gut | sig_gtr

# Assign effect class
classify$effect_class <- "NS"
classify$effect_class[sig_treat & !sig_geno]            <- "Treatment"
classify$effect_class[!sig_treat & sig_geno]            <- "Genotype"
classify$effect_class[sig_treat & sig_geno & sig_gxe]   <- "GxE"
classify$effect_class[sig_treat & sig_geno & !sig_gxe]  <- "Additive"

write.csv(classify, file.path(DESEQ_OUT, "DESeq2_effect_classification.csv"), row.names = FALSE)

cat("\nEffect classification (|log2FC| >", LFC_THRESH, ", padj <", PADJ_THRESH, "):\n")
print(table(classify$effect_class))

# ---- 7. Filter to TF targets + add connectivity -----------------------------

cat("\nBuilding TF target DE list...\n")
all_tf_genes <- unique(tf_edges$target_cosmid)
tf_targets_de <- classify[classify$gene_symbol %in% all_tf_genes &
                            classify$effect_class != "NS", ]

# TF connectivity per gene
tf_per_gene <- tf_edges %>%
  group_by(target_cosmid) %>%
  summarise(n_TFs = n_distinct(TF_name),
            TF_list = paste(sort(unique(TF_name)), collapse = ";"),
            .groups = "drop")

# Common name mapping
name_map <- tf_edges %>%
  filter(target_common_name != "-") %>%
  select(target_cosmid, target_common_name) %>%
  distinct() %>%
  group_by(target_cosmid) %>%
  dplyr::slice(1) %>%
  ungroup()

tf_targets_de <- tf_targets_de %>%
  left_join(tf_per_gene, by = c("gene_symbol" = "target_cosmid")) %>%
  left_join(name_map, by = c("gene_symbol" = "target_cosmid")) %>%
  arrange(desc(n_TFs))

write.csv(tf_targets_de, file.path(DESEQ_OUT, "DESeq2_TF_targets_DE.csv"), row.names = FALSE)

cat("  DE TF targets:", nrow(tf_targets_de), "\n")
print(table(tf_targets_de$effect_class))

# ---- 8. TF-level summary for GRN node labels --------------------------------

cat("\nTF summary (DE targets and GxE proportions):\n")
de_genes  <- classify$gene_symbol[classify$effect_class != "NS"]
gxe_genes <- classify$gene_symbol[classify$effect_class == "GxE"]

tf_stats <- tf_edges %>%
  group_by(TF_name) %>%
  summarise(
    total_targets = n_distinct(target_cosmid),
    DE_targets    = sum(target_cosmid %in% de_genes),
    GxE_targets   = sum(target_cosmid %in% gxe_genes),
    GxE_pct       = ifelse(DE_targets > 0,
                           round(GxE_targets / DE_targets * 100, 1), 0),
    .groups = "drop"
  ) %>%
  arrange(desc(DE_targets))

print(as.data.frame(tf_stats))
write.csv(tf_stats, file.path(DESEQ_OUT, "DESeq2_TF_summary.csv"), row.names = FALSE)

# ---- 8b. Fisher exact test: G x E enrichment per TF -------------------------
#
# 2x2 contingency table for each TF, comparing its target set against the
# rest of the testable genome:
#                         GxE  | non-GxE DE
#   This TF's targets:    a    | b
#   All other genes:      c    | d
#
# A significant enrichment means this TF preferentially regulates G x E
# interaction genes more than expected by chance.

cat("\n---- G x E enrichment by TF (Fisher exact test) ----\n")

testable    <- classify$gene_symbol[!is.na(classify$padj_gxe)]
gxe_set     <- intersect(gxe_genes, testable)
nonNS_de    <- intersect(de_genes, testable)
non_gxe_de  <- setdiff(nonNS_de, gxe_set)

bg_total <- length(testable)

fisher_rows <- lapply(tf_stats$TF_name, function(tf) {
  tf_targets <- unique(tf_edges$target_cosmid[tf_edges$TF_name == tf])
  tf_targets <- intersect(tf_targets, testable)
  
  a <- sum(tf_targets %in% gxe_set)
  b <- sum(tf_targets %in% non_gxe_de)
  c <- length(gxe_set)    - a
  d <- length(non_gxe_de) - b
  
  if ((a + b) == 0) {
    return(data.frame(TF_name = tf, a = a, b = b, c = c, d = d,
                      odds_ratio = NA_real_, p_value = NA_real_,
                      stringsAsFactors = FALSE))
  }
  
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE),
                    alternative = "greater")
  data.frame(TF_name = tf, a = a, b = b, c = c, d = d,
             odds_ratio = unname(ft$estimate),
             p_value    = ft$p.value,
             stringsAsFactors = FALSE)
})

fisher_df <- do.call(rbind, fisher_rows)
fisher_df$padj_BH <- p.adjust(fisher_df$p_value, method = "BH")
fisher_df <- fisher_df[order(fisher_df$p_value), ]

# Pretty print
fisher_print <- fisher_df
fisher_print$odds_ratio <- round(fisher_print$odds_ratio, 2)
fisher_print$p_value    <- signif(fisher_print$p_value,  3)
fisher_print$padj_BH    <- signif(fisher_print$padj_BH,  3)
print(fisher_print)
write.csv(fisher_df, file.path(DESEQ_OUT, "DESeq2_TF_GxE_enrichment.csv"), row.names = FALSE)

# Merge BH-adjusted p-value back into tf_stats so build_figure.py can show it
tf_stats <- merge(tf_stats, fisher_df[, c("TF_name","odds_ratio","p_value","padj_BH")],
                  by = "TF_name", all.x = TRUE)
tf_stats <- tf_stats[order(-tf_stats$DE_targets), ]
write.csv(tf_stats, file.path(DESEQ_OUT, "DESeq2_TF_summary.csv"), row.names = FALSE)

# ---- 8c. STRING PPI query for the GRN TFs -----------------------------------
#
# Pulls high-confidence (combined score >= 700) protein-protein interactions
# among the 12 TFs in the figure, from the STRING REST API for C. elegans
# (taxon 6239). Writes STRING_PPI_edges.csv which build_figure.py reads.

cat("\n---- Querying STRING for PPI edges among GRN TFs ----\n")

grn_tfs <- c("zip-2","daf-16","pqm-1","ceh-60","cebp-1",
             "fos-1","skn-1","elt-2","nhr-28","nhr-77","blmp-1","pha-4")

string_url <- paste0(
  "https://string-db.org/api/tsv/network?",
  "identifiers=", paste(grn_tfs, collapse = "%0d"),
  "&species=6239",
  "&required_score=700",
  "&network_type=physical"
)

ppi_df <- tryCatch({
  read.delim(url(string_url), stringsAsFactors = FALSE)
}, error = function(e) {
  cat("  STRING query failed:", conditionMessage(e), "\n")
  cat("  (Network issue? PPI dashed lines will be skipped in the figure.)\n")
  NULL
})

if (!is.null(ppi_df) && nrow(ppi_df) > 0) {
  # Standardize column names across STRING API versions
  ppi_df$preferredName_A <- tolower(ppi_df$preferredName_A)
  ppi_df$preferredName_B <- tolower(ppi_df$preferredName_B)
  
  # Restrict to pairs where BOTH endpoints are in our GRN TF list
  ppi_df <- ppi_df[ppi_df$preferredName_A %in% grn_tfs &
                     ppi_df$preferredName_B %in% grn_tfs, ]
  
  # Keep one row per unordered pair (A-B == B-A), highest score
  pair_key <- apply(ppi_df[, c("preferredName_A","preferredName_B")], 1,
                    function(x) paste(sort(x), collapse = "|"))
  ppi_df$pair_key <- pair_key
  ppi_df <- ppi_df[order(-ppi_df$score), ]
  ppi_df <- ppi_df[!duplicated(ppi_df$pair_key), ]
  
  ppi_out <- data.frame(
    TF_A  = ppi_df$preferredName_A,
    TF_B  = ppi_df$preferredName_B,
    score = ppi_df$score,
    stringsAsFactors = FALSE
  )
  write.csv(ppi_out, file.path(DESEQ_OUT, "STRING_PPI_edges.csv"), row.names = FALSE)
  cat("  Found", nrow(ppi_out), "high-confidence PPI edges among GRN TFs\n")
  print(ppi_out)
} else {
  cat("  No PPI edges returned. No CSV written.\n")
}

# ---- 9. Summary --------------------------------------------------------------

cat("\n========== SUMMARY ==========\n")
cat("Model:", model_used, "\n")
cat("Thresholds: |log2FC| >", LFC_THRESH, ", padj <", PADJ_THRESH, "\n\n")

cat("DE genes per contrast:\n")
cat("  Treatment in N2:     ", sum(sig_tn2), "\n")
cat("  Treatment in RB2060: ", sum(sig_trb), "\n")
cat("  Genotype untreated:  ", sum(sig_gut), "\n")
cat("  Genotype treated:    ", sum(sig_gtr), "\n")
cat("  GxE interaction:     ", sum(sig_gxe), "\n\n")

cat("Effect classification:\n")
print(table(classify$effect_class))
cat("\nTF targets (DE only):\n")
print(table(tf_targets_de$effect_class))

cat("\n========== OUTPUT FILES ==========\n")
cat("  DESeq2_treatment_in_N2.csv\n")
cat("  DESeq2_treatment_in_RB2060.csv\n")
cat("  DESeq2_genotype_untreated.csv\n")
cat("  DESeq2_genotype_treated.csv\n")
cat("  DESeq2_GxE_interaction.csv\n")
cat("  DESeq2_effect_classification.csv  <- main file for figure rebuild\n")
cat("  DESeq2_TF_targets_DE.csv          <- DE targets with TF connectivity\n")
cat("  DESeq2_TF_summary.csv             <- TF node label stats for GRN\n")
cat("\nRun rebuild_GRN_figures.py next to regenerate Panel E and F.\n")
# ---- Save full dds object for downstream figure scripts (e.g. Figure 2) ----
save(dds, file = file.path(DESEQ_OUT, "DESeq2_dds_interaction_model.RData"))
cat("Wrote:", file.path(DESEQ_OUT, "DESeq2_dds_interaction_model.RData"), "
")
