# ============================================================================
# Figure 1 -- Acute cobalt (II) chloride exposure induces transcriptional
# responses in wildtype and arginine kinase loss-of-function mutants.
#
# Cole Dunivan -- ARGK-2 thesis manuscript
#
# Panels:
#   A: PCA of VST-transformed counts, with inter-experiment batch effect
#      removed via limma::removeBatchEffect (visualization only)
#   B: Knockout validation -- log2FC of target argk transcript in its
#      respective mutant vs N2 (untreated baseline)
#   C: Volcano plots -- argk-2(-/-) and argk-4(-/-) vs N2, in each treatment
#   D: Three-set Venn diagrams (untreated, blue; CoCl2-treated, red)
#
# Pipeline decisions documented in methods:
#   * argk-1 (MAH172) excluded -- incomplete knockdown in raw counts.
#   * Protein-coding filter: REMOVED -- new counts file (data/RNASEQ61125.csv)
#     does not carry a biotype/GENENAME column. If biotype annotation is
#     needed, join an external annotation table before the DESeq2 fit.
#   * PCA-only batch correction: limma::removeBatchEffect() is applied to
#     the VST matrix using experiment as the batch variable, with the
#     genotype*treatment design preserved. This is for visualization
#     only and does not propagate to DESeq2 testing -- those tests
#     handle batch correctly via the `experiment` term in the design
#     formula and are computed from the raw counts.
#   * DEGs for Panel D: BH-adjusted p < 0.05 (no LFC cutoff). LFC > 1
#     is applied only on Panel C volcano highlighting. Update figure
#     caption to match.
#
# CHANGES vs previous version:
#   [FIX 1] COUNTS_FILE gene identifier column renamed from "gene_symbol"
#           to "gene" to match the new counts file header.
#   [FIX 2] Biotype filter removed -- new counts file has no GENENAME /
#           biotype column. The stopifnot() guard and all GENENAME-dependent
#           filtering lines have been dropped.
#   [FIX 3] ANNOT_COLS no longer includes "GENENAME" (absent from new file).
#           "ENTREZID", "ENSEMBL", "WORMBASE" are also absent; any_of()
#           handles this gracefully but the list is updated for clarity.
#   [FIX 4] KO_TARGETS gene_id updated to include the "CELE_" prefix that
#           the new counts file uses (e.g. CELE_W10C8.5, CELE_F46H5.3).
#           Without this Panel B would silently return NA for both targets.
# ============================================================================


# ---- 0. Setup --------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
  library(stringr); library(purrr); library(ggplot2)
  library(DESeq2)
  library(limma)         # removeBatchEffect (for PCA only)
  library(ggrepel)
  library(VennDiagram)   # apt-available alternative to ggVennDiagram
  library(grid)
  library(patchwork)
  library(scales)
})
# VennDiagram writes a noisy log via futile.logger; silence it
suppressPackageStartupMessages(library(futile.logger))
flog.threshold(ERROR)

set.seed(2025)

# ---- Config ----
COUNTS_FILE   <- "data/RNASEQ61125.csv"
METADATA_FILE <- "data/Sample_Metadata_Table.csv"
OUT_DIR       <- "outputs/figures/Figure1"
DEG_DIR       <- file.path(OUT_DIR, "DEGs")

PADJ_CUTOFF   <- 0.05
LFC_CUTOFF    <- 1
MIN_COUNT_SUM <- 10

EXCLUDE_STRAINS  <- c("MAH172")
# [FIX 2] EXCLUDE_BIOTYPES removed -- no biotype column in new counts file.
# [FIX 3] ANNOT_COLS updated -- none of these columns exist in new file;
#          kept so select(-any_of(ANNOT_COLS)) remains a no-op gracefully.
ANNOT_COLS       <- c("ENTREZID", "ENSEMBL", "WORMBASE")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DEG_DIR, showWarnings = FALSE, recursive = TRUE)


# ---- Style -----------------------------------------------------------------

COL_UP   <- "#cb181d"
COL_DOWN <- "#2171b5"
COL_NS   <- "grey70"
COL_NA   <- "grey85"

COL_UNTR <- c(low = "#f7fbff", high = "#2171b5")
COL_TR   <- c(low = "#fff5f0", high = "#cb181d")

GROUP_COLORS <- c(
  "N2_untreated"     = "#3CB39C",
  "N2_treated"       = "#97C04F",
  "RB2060_untreated" = "#6FAEDB",
  "RB2060_treated"   = "#3FC5D8",
  "RB2598_untreated" = "#F2A6BC",
  "RB2598_treated"   = "#E66FA8"
)

GROUP_LABEL_EXPRS <- c(
  "N2_untreated"     = "'reference untreated'",
  "N2_treated"       = "'reference treated'",
  "RB2060_untreated" = "italic('argk-2')^'-/-'*' untreated'",
  "RB2060_treated"   = "italic('argk-2')^'-/-'*' treated'",
  "RB2598_untreated" = "italic('argk-4')^'-/-'*' untreated'",
  "RB2598_treated"   = "italic('argk-4')^'-/-'*' treated'"
)

GROUP_LEVELS <- c("N2_untreated",     "N2_treated",
                  "RB2060_untreated", "RB2060_treated",
                  "RB2598_untreated", "RB2598_treated")

GENOTYPE_LABELS <- c(
  "N2"     = "N2",
  "RB2060" = "argk-2(-/-)",
  "RB2598" = "argk-4(-/-)"
)

# [FIX 4] gene_id values updated with CELE_ prefix to match new counts file.
KO_TARGETS <- tibble(
  strain     = c("RB2060",       "RB2598"),
  gene_label = c("argk-2",       "argk-4"),
  gene_id    = c("CELE_W10C8.5", "CELE_F46H5.3")
)


# ---- 1. Load and align data ------------------------------------------------

# Use robust CSV loader from shared utils (handles CR-only line endings)
# Source the utils file from one of these standard locations:
for (utils_path in c("scripts/_utils.R", "_utils.R", "../scripts/_utils.R")) {
  if (file.exists(utils_path)) { source(utils_path); break }
}
if (!exists("read_csv_robust")) {
  stop("Could not source scripts/_utils.R. Run this script from the repo root.")
}

raw_counts <- read_csv_robust(COUNTS_FILE, check.names = FALSE,
                              stringsAsFactors = FALSE)

# [FIX 1] Guard updated: new file uses "gene" not "gene_symbol"; no GENENAME.
stopifnot("gene" %in% colnames(raw_counts))

# [FIX 1] Renamed gene column; [FIX 2] biotype filter block removed entirely.
raw_counts <- raw_counts %>%
  filter(!is.na(gene),
         !duplicated(gene))

# Coerce sample columns to numeric (CR line endings can make them character)
sample_cols <- setdiff(colnames(raw_counts), c("gene", ANNOT_COLS))
for (cc in sample_cols) {
  raw_counts[[cc]] <- as.numeric(as.character(raw_counts[[cc]]))
}

cat("Genes loaded: ", nrow(raw_counts), "\n", sep = "")
cat("(Note: biotype filter not applied -- no GENENAME column in new counts file.)\n")

# [FIX 1] column renamed from "gene_symbol" to "gene"
counts <- raw_counts %>%
  column_to_rownames("gene") %>%
  select(-any_of(ANNOT_COLS))

coldata <- read_csv_robust(METADATA_FILE)
id_col_fig1 <- intersect(c("SampleLabel", "library_name", "Sample", "sample"),
                         colnames(coldata))[1]
if (is.na(id_col_fig1)) id_col_fig1 <- colnames(coldata)[1]
rownames(coldata) <- trimws(coldata[[id_col_fig1]])
coldata <- coldata[!coldata$genotype %in% EXCLUDE_STRAINS, , drop = FALSE]

# Harmonize sample naming convention (rename GEO library names to A-labels if needed)
counts_mat_temp <- as.matrix(counts)
counts_mat_temp <- harmonize_sample_names(counts_mat_temp, coldata,
                                          target_col = id_col_fig1,
                                          alt_col = setdiff(c("SampleLabel", "GEO_library_name"),
                                                            id_col_fig1)[1])
counts <- as.data.frame(counts_mat_temp, check.names = FALSE)

shared_samples <- intersect(colnames(counts), rownames(coldata))
if (length(shared_samples) == 0L) {
  stop("No overlapping sample names between counts and metadata.")
}
counts  <- counts[, shared_samples]
coldata <- coldata[shared_samples, , drop = FALSE]

counts <- as.matrix(counts)
mode(counts) <- "integer"

coldata$genotype   <- relevel(factor(coldata$genotype),   ref = "N2")
coldata$treatment  <- relevel(factor(coldata$treatment),  ref = "untreated")
coldata$experiment <- factor(coldata$experiment)

stopifnot(all(colnames(counts) == rownames(coldata)))

cat("\nLoaded ", nrow(counts), " genes x ", ncol(counts), " samples.\n", sep = "")
cat("Sample breakdown:\n")
print(table(coldata$genotype, coldata$treatment))


# ---- 2. DESeq2 models ------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = coldata,
  design    = ~ genotype + treatment + experiment +
    genotype:treatment +
    genotype:experiment +
    treatment:experiment
)
dds <- dds[rowSums(counts(dds)) >= MIN_COUNT_SUM, ]
dds <- DESeq(dds)

make_within_dds <- function(tr_label) {
  cd <- coldata[coldata$treatment == tr_label, , drop = FALSE]
  cd$genotype <- relevel(droplevels(cd$genotype), ref = "N2")
  cc <- counts[, rownames(cd)]
  d  <- DESeqDataSetFromMatrix(countData = cc, colData = cd, design = ~ genotype)
  d  <- d[rowSums(counts(d)) >= MIN_COUNT_SUM, ]
  DESeq(d)
}
dds_un <- make_within_dds("untreated")
dds_tr <- make_within_dds("treated")

dds_n2 <- {
  cd <- coldata[coldata$genotype == "N2", , drop = FALSE]
  cd$treatment <- relevel(droplevels(cd$treatment), ref = "untreated")
  cc <- counts[, rownames(cd)]
  d  <- DESeqDataSetFromMatrix(countData = cc, colData = cd, design = ~ treatment)
  d  <- d[rowSums(counts(d)) >= MIN_COUNT_SUM, ]
  DESeq(d)
}


# ---- 3. Panel A -- PCA -----------------------------------------------------

vsd <- vst(dds, blind = FALSE)

# Batch correction (VISUALIZATION ONLY).
mat <- assay(vsd)
mm  <- model.matrix(~ genotype * treatment,
                    data = as.data.frame(colData(vsd)))
assay(vsd) <- limma::removeBatchEffect(mat,
                                       batch  = colData(vsd)$experiment,
                                       design = mm)

pcaData <- plotPCA(vsd, intgroup = c("genotype", "treatment", "experiment"),
                   returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

cat("\nANOVA on principal components (after batch correction):\n")
cat("  PC1:\n"); print(summary(aov(PC1 ~ genotype * treatment + experiment, pcaData))[[1]])
cat("  PC2:\n"); print(summary(aov(PC2 ~ genotype * treatment + experiment, pcaData))[[1]])

pcaData <- pcaData %>%
  mutate(
    group = factor(paste(genotype, treatment, sep = "_"),
                   levels = GROUP_LEVELS)
  )

panelA <- ggplot(pcaData, aes(x = PC1, y = PC2, colour = group)) +
  geom_point(size = 3.5, alpha = 0.95) +
  scale_colour_manual(
    name   = NULL,
    values = GROUP_COLORS,
    labels = parse(text = GROUP_LABEL_EXPRS[GROUP_LEVELS]),
    breaks = GROUP_LEVELS
  ) +
  guides(colour = guide_legend(ncol = 2, byrow = TRUE)) +
  labs(x = sprintf("PC1: %d%% variance", percentVar[1]),
       y = sprintf("PC2: %d%% variance", percentVar[2])) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor     = element_blank(),
    legend.position      = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background    = element_rect(fill = alpha("white", 0.7),
                                        colour = NA),
    legend.key.size      = unit(0.4, "cm"),
    legend.text          = element_text(size = 10)
  )

ggsave(file.path(OUT_DIR, "Figure1_PanelA_PCA.pdf"),
       panelA, width = 7.5, height = 5.5)


# ---- 4. Panel B -- Knockout validation -------------------------------------

ko_rows <- vector("list", nrow(KO_TARGETS))
for (i in seq_len(nrow(KO_TARGETS))) {
  strain  <- KO_TARGETS$strain[i]
  gene_id <- KO_TARGETS$gene_id[i]   # [FIX 4] now "CELE_W10C8.5" etc.
  label   <- KO_TARGETS$gene_label[i]
  
  res <- as.data.frame(results(dds_un,
                               contrast = c("genotype", strain, "N2"),
                               alpha    = PADJ_CUTOFF))
  if (!gene_id %in% rownames(res)) {
    warning(sprintf("Target %s not found in dds_un results for %s.",
                    gene_id, strain))
    next
  }
  ko_rows[[i]] <- tibble(
    gene_label = label,
    strain     = strain,
    log2FC     = res[gene_id, "log2FoldChange"],
    padj       = res[gene_id, "padj"]
  )
}

ko_lfc <- bind_rows(ko_rows)
if (nrow(ko_lfc) == 0L) {
  stop("Knockout validation: no target genes found. ",
       "Check KO_TARGETS$gene_id against rownames(counts).")
}
ko_lfc <- ko_lfc %>%
  mutate(
    stars = case_when(
      is.na(padj) ~ "ns",
      padj < 1e-4 ~ "***",
      padj < 1e-3 ~ "**",
      padj < 0.05 ~ "*",
      TRUE        ~ "ns"
    ),
    star_y = ifelse(log2FC < 0, log2FC - 0.4, log2FC + 0.4)
  )

cat("\nKnockout validation (log2FC, padj):\n")
print(ko_lfc, n = Inf)

panelB <- ggplot(ko_lfc, aes(x = gene_label, y = log2FC, fill = gene_label)) +
  geom_col(width = 0.55, colour = "black") +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
  geom_text(aes(label = stars, y = star_y), size = 5.5, fontface = "bold") +
  scale_fill_manual(
    values = c("argk-2" = GROUP_COLORS[["RB2060_untreated"]],
               "argk-4" = GROUP_COLORS[["RB2598_untreated"]]),
    guide  = "none"
  ) +
  labs(x = "Target transcript",
       y = expression(log[2]~"fold change (mutant vs N2)")) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank())

ggsave(file.path(OUT_DIR, "Figure1_PanelB_KO_validation.pdf"),
       panelB, width = 4, height = 4)


# ---- 5. Panel C -- Volcano plots -------------------------------------------

contrast_grid <- expand.grid(
  strain    = c("RB2060", "RB2598"),
  treatment = c("untreated", "treated"),
  stringsAsFactors = FALSE
) %>%
  mutate(
    pretty_strain = unname(GENOTYPE_LABELS[strain]),
    label    = paste0(pretty_strain, " vs N2 (", treatment, ")"),
    file_tag = paste0(strain, "_vs_N2_", treatment)
  )

make_volcano <- function(res, title_text) {
  df <- as.data.frame(res) %>%
    rownames_to_column("Gene") %>%
    mutate(
      sig = case_when(
        is.na(padj)                                       ~ "NA",
        padj < PADJ_CUTOFF & log2FoldChange >  LFC_CUTOFF ~ "Up",
        padj < PADJ_CUTOFF & log2FoldChange < -LFC_CUTOFF ~ "Down",
        TRUE                                              ~ "ns"
      ),
      sig = factor(sig, levels = c("Up", "Down", "ns", "NA"))
    )
  
  top_lab <- df %>%
    filter(sig %in% c("Up", "Down")) %>%
    arrange(padj) %>%
    slice_head(n = 5)
  
  ggplot(df, aes(x = log2FoldChange, y = -log10(padj), colour = sig)) +
    geom_point(alpha = 0.55, size = 1) +
    scale_colour_manual(
      values = c("Up" = COL_UP, "Down" = COL_DOWN,
                 "ns" = COL_NS, "NA" = COL_NA),
      labels = c("Up (padj<0.05, log2FC>1)",
                 "Down (padj<0.05, log2FC<-1)",
                 "Not significant", "NA"),
      drop = FALSE
    ) +
    geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF),
               linetype = "dashed", linewidth = 0.3) +
    geom_hline(yintercept = -log10(PADJ_CUTOFF),
               linetype = "dashed", linewidth = 0.3) +
    geom_text_repel(data = top_lab, aes(label = Gene),
                    size = 3, colour = "black",
                    box.padding = 0.5, max.overlaps = Inf,
                    segment.colour = "grey50",
                    show.legend = FALSE) +
    labs(title = title_text,
         x = expression(log[2]~"fold change"),
         y = expression(-log[10]~"(adjusted p-value)"),
         colour = "Significance") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12),
          panel.grid.minor = element_blank())
}

volcano_plots <- vector("list", nrow(contrast_grid))
for (i in seq_len(nrow(contrast_grid))) {
  row      <- contrast_grid[i, ]
  this_dds <- if (row$treatment == "untreated") dds_un else dds_tr
  res <- results(this_dds,
                 contrast = c("genotype", row$strain, "N2"),
                 alpha    = PADJ_CUTOFF)
  
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene_symbol")
  write.csv(res_df,
            file = file.path(DEG_DIR,
                             paste0(row$strain, "_vs_N2_",
                                    row$treatment, "_DEGs.csv")),
            row.names = FALSE)
  
  p <- make_volcano(res, row$label)
  volcano_plots[[i]] <- p
  ggsave(file.path(OUT_DIR,
                   paste0("Figure1_PanelC_Volcano_", row$file_tag, ".pdf")),
         p, width = 5, height = 4.5)
}

res_n2 <- results(dds_n2,
                  contrast = c("treatment", "treated", "untreated"),
                  alpha    = PADJ_CUTOFF)
n2_df <- as.data.frame(res_n2) %>%
  rownames_to_column("gene_symbol")
write.csv(n2_df,
          file = file.path(DEG_DIR, "N2_treated_vs_untreated_DEGs.csv"),
          row.names = FALSE)

panelC <- wrap_plots(volcano_plots, ncol = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")


# ---- 6. Panel D -- Venn diagrams -------------------------------------------
#
# Two Venns, each internally consistent (all circles = same comparison type):
#
#   UNTREATED (2-set, blue):
#     Each circle = genotype effect vs N2 within the untreated condition.
#     N2 is the reference and has no circle (N2 unt vs N2 unt = 0 DEGs).
#       argk-2 circle: RB2060 untreated vs N2 untreated
#       argk-4 circle: RB2598 untreated vs N2 untreated
#
#   TREATED (3-set, red):
#     Each circle = treatment response (treated vs untreated) within that
#     genotype. All three circles are the same comparison type.
#       N2 circle    : N2 treated vs N2 untreated      (wildtype stress response)
#       argk-2 circle: RB2060 treated vs RB2060 untreated (mutant stress response)
#       argk-4 circle: RB2598 treated vs RB2598 untreated (mutant stress response)
#     Overlap = conserved stress-response genes.
#     Unique to N2 = genes N2 responds to that the mutants fail to mount.

# DEG filter: padj < 0.05 AND |log2FC| > 1 (both thresholds match figure caption)
deg_filter <- function(df) {
  df %>%
    filter(!is.na(padj), !is.na(log2FoldChange),
           padj < PADJ_CUTOFF,
           abs(log2FoldChange) > LFC_CUTOFF) %>%
    pull(gene_symbol) %>%
    unique()
}
read_deg <- function(file) {
  read.csv(file.path(DEG_DIR, file)) %>%
    deg_filter()
}

# --- Write within-mutant treatment DEG files (needed for treated Venn) ------
# Same approach as dds_n2: fit a per-strain model, then extract treatment contrast.
make_within_strain_dds <- function(strain_label) {
  cd <- coldata[coldata$genotype == strain_label, , drop = FALSE]
  cd$treatment <- relevel(droplevels(cd$treatment), ref = "untreated")
  cc <- counts[, rownames(cd)]
  d  <- DESeqDataSetFromMatrix(countData = cc, colData = cd,
                               design = ~ treatment)
  d  <- d[rowSums(counts(d)) >= MIN_COUNT_SUM, ]
  DESeq(d)
}

for (strain in c("RB2060", "RB2598")) {
  dds_s <- make_within_strain_dds(strain)
  res_s <- results(dds_s,
                   contrast = c("treatment", "treated", "untreated"),
                   alpha    = PADJ_CUTOFF)
  df_s  <- as.data.frame(res_s) %>% rownames_to_column("gene_symbol")
  write.csv(df_s,
            file      = file.path(DEG_DIR,
                                  paste0(strain, "_treated_vs_untreated_DEGs.csv")),
            row.names = FALSE)
}

# --- Full-model treatment contrasts (overwrite the within-strain files) ----
# The within-strain DESeq2 fits above use only 8 samples each and drop the
# experiment batch term, producing inflated DEG counts (e.g. argk-4 came out
# with >3000 treatment DEGs in those fits). Re-derive each genotype's
# treatment response from the SAME full-factorial model used in scripts
# 01, 03, and 09 (apeglm for the reference-genotype coefficient, ashr for
# the linear-combination contrasts), so Panel D's Venn is statistically
# sound and consistent with the rest of the figure suite.
message("Re-extracting per-genotype treatment DEGs from full-factorial model with shrinkage...")
write_full_model_trt_deg <- function(out_file, coef = NULL, contrast_list = NULL) {
  res <- if (!is.null(coef)) {
    DESeq2::lfcShrink(dds, coef = coef, type = "apeglm", quiet = TRUE)
  } else {
    DESeq2::lfcShrink(dds, contrast = contrast_list, type = "ashr", quiet = TRUE)
  }
  df <- as.data.frame(res) %>% tibble::rownames_to_column("gene_symbol")
  write.csv(df, file = file.path(DEG_DIR, out_file), row.names = FALSE)
  n_sig <- sum(!is.na(df$padj) & df$padj < PADJ_CUTOFF &
               abs(df$log2FoldChange) > LFC_CUTOFF)
  message(sprintf("  %-45s %5d DEGs (padj<%.2f, |LFC|>%.0f)",
                  out_file, n_sig, PADJ_CUTOFF, LFC_CUTOFF))
}
write_full_model_trt_deg("N2_treated_vs_untreated_DEGs.csv",
                         coef = "treatment_treated_vs_untreated")
write_full_model_trt_deg("RB2060_treated_vs_untreated_DEGs.csv",
                         contrast_list = list(c("treatment_treated_vs_untreated",
                                                "genotypeRB2060.treatmenttreated")))
write_full_model_trt_deg("RB2598_treated_vs_untreated_DEGs.csv",
                         contrast_list = list(c("treatment_treated_vs_untreated",
                                                "genotypeRB2598.treatmenttreated")))

# --- Load gene lists --------------------------------------------------------
genes_argk2_un  <- read_deg("RB2060_vs_N2_untreated_DEGs.csv")      # argk-2 unt vs N2 unt
genes_argk4_un  <- read_deg("RB2598_vs_N2_untreated_DEGs.csv")      # argk-4 unt vs N2 unt
genes_n2_stress <- read_deg("N2_treated_vs_untreated_DEGs.csv")      # N2 tx vs N2 unt
genes_argk2_tx  <- read_deg("RB2060_treated_vs_untreated_DEGs.csv")  # argk-2 tx vs argk-2 unt
genes_argk4_tx  <- read_deg("RB2598_treated_vs_untreated_DEGs.csv")  # argk-4 tx vs argk-4 unt

cat("\nDEG set sizes (padj < 0.05, |log2FC| > 1):\n")
cat(sprintf("  argk-2(-/-) vs N2 untreated          : %d\n", length(genes_argk2_un)))
cat(sprintf("  argk-4(-/-) vs N2 untreated          : %d\n", length(genes_argk4_un)))
cat(sprintf("  N2 treated vs N2 untreated           : %d\n", length(genes_n2_stress)))
cat(sprintf("  argk-2(-/-) treated vs untreated     : %d\n", length(genes_argk2_tx)))
cat(sprintf("  argk-4(-/-) treated vs untreated     : %d\n", length(genes_argk4_tx)))

# Single 3-set Venn — logically coherent: each circle is the *treatment*
# response DEG set (treated vs untreated) for one genotype. All three
# circles measure the same kind of contrast, so overlaps and uniques are
# directly interpretable as conserved vs genotype-specific stress responses.
venn_tr <- list(
  "N2"          = genes_n2_stress,
  "argk-2(-/-)" = genes_argk2_tx,
  "argk-4(-/-)" = genes_argk4_tx
)

make_venn <- function(venn_list, low, high, title_text, subtitle_text = NULL) {
  # VennDiagram::venn.diagram returns a gList of grobs. wrap_elements(panel=)
  # turns it into a patchwork-compatible plot we can layout next to ggplots.
  n <- length(venn_list)
  fills <- if (n == 2) c(low, high) else c(low, high, "grey75")
  vd <- venn.diagram(
    x          = venn_list,
    filename   = NULL,
    fill       = fills,
    alpha      = 0.55,
    lwd        = 1,
    col        = "black",
    fontfamily      = "sans",
    cat.fontfamily  = "sans",
    cex             = 1.2,
    cat.cex         = 1.05,
    cat.default.pos = "outer",
    main            = title_text,
    main.cex        = 1.2,
    main.fontfamily = "sans",
    sub             = subtitle_text,
    sub.cex         = 0.85,
    sub.fontfamily  = "sans",
    margin          = 0.08,
    disable.logging = TRUE
  )
  # Clean up the per-call VennDiagram .log file that's silently created
  unlink(list.files(pattern = "^VennDiagram.*\\.log$"))
  patchwork::wrap_elements(panel = grid::gTree(children = vd))
}

p_venn_tr <- make_venn(
  venn_tr, COL_TR["low"], COL_TR["high"],
  title_text    = "Treatment response across genotypes",
  subtitle_text = "treated vs untreated within each genotype | padj < 0.05, |log\u2082FC| > 1"
)

ggsave(file.path(OUT_DIR, "Figure1_PanelD_Venn_Treatment.pdf"),
       p_venn_tr, width = 6.5, height = 5.5)

panelD <- p_venn_tr


# ---- 7. Composite figure ---------------------------------------------------

figure1 <- (panelA | panelB) /
  panelC /
  panelD +
  plot_layout(heights = c(1.0, 1.4, 1.0)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave(file.path(OUT_DIR, "Figure1_Composite.pdf"),
       figure1, width = 12, height = 14)
ggsave(file.path(OUT_DIR, "Figure1_Composite.png"),
       figure1, width = 12, height = 14, dpi = 600)

message("Done. All Figure 1 outputs written to: ", normalizePath(OUT_DIR))