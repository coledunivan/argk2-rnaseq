# =============================================================================
# Fig3_combined.R — recreate Figure 3 end-to-end in R
#
# Builds a 3-panel figure from raw counts:
#   A — GO BP enrichment dot plot: N2 treated vs untreated (top 20 terms)
#   B — GO BP enrichment dot plot: argk-2-/- vs N2 (treated) (top 9 terms)
#   C — Immune effect map: scatter of genotype effect vs treatment effect
#       for immune-annotated genes, colored by effect category
#
# DESeq2 model matches Fig4_combined.R and Fig5_combined.R so contrast values
# are consistent across all three figures.
#
# Required packages: DESeq2, clusterProfiler, org.Ce.eg.db, GO.db,
#                    AnnotationDbi, ggplot2, patchwork, dplyr, tibble,
#                    tidyr, purrr, scales, stringr
#
# To install missing Bioconductor packages:
#   if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#   BiocManager::install(c("DESeq2", "clusterProfiler", "org.Ce.eg.db",
#                          "GO.db", "AnnotationDbi"))
#
# Outputs (in ./Fig3_outputs/):
#   Fig3_recreated.pdf, Fig3_recreated.png   — combined 3-panel figure
#   Fig3_panel_A_GO.csv                      — full Panel A GO enrichment table
#   Fig3_panel_B_GO.csv                      — full Panel B GO enrichment table
#   Fig3_panel_C_data.csv                    — immune genes + categories + LFCs
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(scales)
  library(stringr)
})
# Offline GO helpers (load_go_offline, enrichGO_offline) live in _utils.R
for (utils_path in c("scripts/_utils.R","_utils.R","../scripts/_utils.R")) {
  if (file.exists(utils_path)) { source(utils_path); break }
}
stopifnot(exists("enrichGO_offline"))

# Force dplyr verbs (DESeq2 + AnnotationDbi mask select/filter)
select  <- dplyr::select
filter  <- dplyr::filter
mutate  <- dplyr::mutate
summarize <- dplyr::summarize
arrange <- dplyr::arrange
group_by <- dplyr::group_by
distinct <- dplyr::distinct
rename <- dplyr::rename

# ---- Inputs (auto-detect) ---------------------------------------------------
find_first <- function(patterns, label) {
  hits <- character(0)
  for (p in patterns) hits <- c(hits, Sys.glob(p))
  hits <- unique(hits)
  if (!length(hits)) stop(sprintf("Could not find %s file. Tried: %s",
                                  label, paste(patterns, collapse = ", ")))
  hits[1]
}
rna_file  <- find_first(c("data/reference/RNAseq_to_TF_Targets*.csv", "RNAseq_to_TF_Targets*.csv", "*TF_Targets*.csv"), "counts")
meta_file <- find_first(c("data/Sample_Metadata_Table*.csv", "Sample_Metadata_Table*.csv", "*Metadata*.csv"),  "metadata")
message("Using counts:   ", rna_file)
message("Using metadata: ", meta_file)

go_padj_cut <- 0.10   # GO enrichment input cutoff (matches GOEnrichment.R)
de_padj_cut <- 0.05   # significance cutoff for Panel C categorization
min_sum     <- 10
out_dir     <- "outputs/figures/Figure3"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

PANEL_A_TOP <- 20
PANEL_B_TOP <- 9

# Immune GO term regex (matches Immune_Effect_Map.R)
immune_regex <- paste0(
  "(immune|innate|host defen[cs]e|antimicrob|pathogen|bacteri|fungi?|virus|viral|",
  "antiviral|antibacteri|antifung|lysozyme|lys-\\d+|clec|c[- ]?type lectin|lectin|",
  "irg-\\d+|pals-\\d+|nlp-\\d+|spp-\\d+|defensin|pattern[- ]?recognition|pgrp|",
  "peptidoglycan|tlr|toll[- ]like|nlr|nod[- ]like|interferon|ifn-|complement|",
  "response to bacteri|response to fung|response to virus)"
)

# ---- Load -------------------------------------------------------------------
rna_raw <- read.csv(rna_file,  check.names = FALSE, stringsAsFactors = FALSE)
meta_df <- read.csv(meta_file, stringsAsFactors = FALSE)

count_cols <- grep("^A\\d+", colnames(rna_raw), value = TRUE)
stopifnot(length(count_cols) > 0)

counts <- rna_raw %>%
  select(Gene, all_of(count_cols)) %>%
  distinct() %>%
  column_to_rownames("Gene")

meta_df <- meta_df %>%
  filter(SampleLabel %in% colnames(counts)) %>%
  column_to_rownames("SampleLabel") %>%
  mutate(
    genotype  = factor(genotype,  levels = c("N2", "RB2060", "RB2598")),
    treatment = factor(treatment, levels = c("untreated", "treated"))
  ) %>%
  filter(!is.na(genotype), !is.na(treatment), !is.na(experiment))

shared <- intersect(colnames(counts), rownames(meta_df))
if (!length(shared)) stop("No overlapping samples between counts and metadata.")
counts  <- counts[, shared, drop = FALSE]
meta_df <- meta_df[shared, , drop = FALSE]

# ---- DESeq2 -----------------------------------------------------------------
message("Running DESeq2…")
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = meta_df,
  design    = ~ genotype * treatment + genotype:experiment + treatment:experiment
)
dds <- dds[rowSums(counts(dds)) > min_sum, ]
dds <- DESeq(dds)

dds_simple <- dds
dds_simple$genotype  <- relevel(dds_simple$genotype,  ref = "N2")
dds_simple$treatment <- relevel(dds_simple$treatment, ref = "untreated")
design(dds_simple) <- ~ genotype + treatment + genotype:treatment
dds_simple <- DESeq(dds_simple)

# ---- Contrast helpers (match Fig 4/5) ---------------------------------------
get_genotype_within_treatment <- function(dds_obj, genotype_alt = "RB2060", trt = "untreated") {
  rn <- resultsNames(dds_obj)
  main_coef <- rn[grepl("^genotype", rn) & grepl(genotype_alt, rn) & !grepl("treatment", rn)]
  inter     <- rn[grepl("genotype", rn)  & grepl(genotype_alt, rn) &  grepl("treatment.*treated", rn)]
  if (trt == "untreated") results(dds_obj, name = main_coef)
  else                    results(dds_obj, list(c(main_coef, inter)))
}
get_treatment_within_genotype <- function(dds_obj, genotype_level = "N2") {
  rn <- resultsNames(dds_obj)
  trt_pos <- rn[grepl("^treatment", rn) & grepl("treated", rn)]
  inter   <- rn[grepl("genotype", rn) & grepl(genotype_level, rn) & grepl("treatment.*treated", rn)]
  if (genotype_level == "N2") results(dds_obj, name = trt_pos)
  else                        results(dds_obj, list(c(trt_pos, inter)))
}

message("Extracting contrasts…")
res_treat_N2  <- get_treatment_within_genotype(dds_simple, "N2")
res_treat_RB  <- get_treatment_within_genotype(dds_simple, "RB2060")
res_geno_unt  <- get_genotype_within_treatment(dds_simple, "RB2060", "untreated")
res_geno_trt  <- get_genotype_within_treatment(dds_simple, "RB2060", "treated")
inter_name    <- grep("genotype.*RB2060.*treatment.*treated",
                      resultsNames(dds_simple), value = TRUE)
res_gxe       <- results(dds_simple, name = inter_name)

# ---- Gene → ENTREZ mapping (CSV-based, offline) -----------------------------
# Build cosmid → common-name (SYMBOL) map from local TFLink file,
# then SYMBOL → ENTREZID from the offline org_Ce_eg_GO_map.csv.
message("Mapping gene IDs to ENTREZ (offline via TFLink + GO map)…")
go_map_full <- load_go_offline("data/reference/org_Ce_eg_GO_map.csv")
sym_to_entrez_map <- go_map_full %>%
  dplyr::distinct(SYMBOL, ENTREZID) %>%
  dplyr::filter(!is.na(SYMBOL), !is.na(ENTREZID))
sym_to_entrez_map$SYMBOL_lc <- tolower(sym_to_entrez_map$SYMBOL)

tfl_raw <- read.csv("data/reference/RNAseq_to_TF_Targets.csv",
                    stringsAsFactors = FALSE, check.names = FALSE)
tfl_raw <- tfl_raw[, colnames(tfl_raw) != "" & !is.na(colnames(tfl_raw)), drop = FALSE]
cosmid_to_sym <- data.frame(
  gene      = trimws(as.character(tfl_raw$Gene)),
  SYMBOL_lc = tolower(trimws(as.character(tfl_raw$Name.Target))),
  stringsAsFactors = FALSE
)
cosmid_to_sym <- cosmid_to_sym[
  cosmid_to_sym$SYMBOL_lc != "" & cosmid_to_sym$SYMBOL_lc != "-" &
    !is.na(cosmid_to_sym$SYMBOL_lc), , drop = FALSE]
cosmid_to_sym <- cosmid_to_sym[!duplicated(cosmid_to_sym$gene), ]

map_to_entrez <- function(gene_names) {
  cleaned <- sub("^CELE_", "", as.character(gene_names))
  # Direct cosmid match in TFLink → SYMBOL → ENTREZID
  step1 <- data.frame(gene = cleaned, stringsAsFactors = FALSE)
  step1 <- merge(step1, cosmid_to_sym, by = "gene", all.x = TRUE,
                 sort = FALSE)
  step2 <- merge(step1, sym_to_entrez_map[, c("SYMBOL_lc","ENTREZID")],
                 by = "SYMBOL_lc", all.x = TRUE, sort = FALSE)
  # Also try matching cosmid directly as ENTREZ symbol (some gene_symbol values
  # are already common-name lowercased)
  needs_fallback <- is.na(step2$ENTREZID)
  if (any(needs_fallback)) {
    fb_keys <- tolower(step2$gene[needs_fallback])
    fb <- sym_to_entrez_map$ENTREZID[match(fb_keys, sym_to_entrez_map$SYMBOL_lc)]
    step2$ENTREZID[needs_fallback] <- fb
  }
  # Reorder back to input
  result <- step2$ENTREZID[match(cleaned, step2$gene)]
  names(result) <- gene_names
  result
}

all_genes <- rownames(counts)
gene_to_entrez <- map_to_entrez(all_genes)
n_mapped <- sum(!is.na(gene_to_entrez))
message(sprintf("ENTREZ mapping: %d / %d genes mapped (%.0f%%)",
                n_mapped, length(all_genes), 100 * n_mapped / length(all_genes)))

# Universe for enrichGO = all genes with ENTREZ in our DE table
universe_entrez <- na.omit(unique(gene_to_entrez))

# ---- GO enrichment helper ---------------------------------------------------
run_go_BP <- function(de_genes, label) {
  entrez <- na.omit(unique(gene_to_entrez[de_genes]))
  message(sprintf("[%s] %d DE genes → %d unique ENTREZ", label,
                  length(de_genes), length(entrez)))
  if (length(entrez) < 5) {
    warning("Not enough mapped DE genes for ", label)
    return(NULL)
  }
  ego_df <- enrichGO_offline(
    gene         = entrez,
    universe     = universe_entrez,
    go_map       = go_map_full,
    ontology     = "BP",
    pvalueCutoff = 0.10,
    minGSSize    = 10,
    maxGSSize    = 500
  )
  if (is.null(ego_df) || nrow(ego_df) == 0) {
    warning("No significant GO BP terms for ", label)
    return(NULL)
  }
  # Replicate clusterProfiler's `readable=TRUE`: map ENTREZ in geneID -> SYMBOL
  ez2sym <- setNames(sym_to_entrez_map$SYMBOL, as.character(sym_to_entrez_map$ENTREZID))
  ego_df$geneID <- vapply(strsplit(ego_df$geneID, "/"), function(ids) {
    syms <- ez2sym[ids]; syms[is.na(syms)] <- ids[is.na(syms)]
    paste(syms, collapse = "/")
  }, character(1))
  df <- ego_df %>%
    mutate(
      GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x)
        as.numeric(x[1]) / as.numeric(x[2]))
    )
  df
}

# ---- Panel A: GO for N2 treated vs untreated --------------------------------
de_treat_N2 <- as.data.frame(res_treat_N2) %>%
  rownames_to_column("Gene") %>%
  filter(!is.na(padj), padj < go_padj_cut) %>%
  pull(Gene)
go_A <- run_go_BP(de_treat_N2, "Panel A: N2 Treated vs Untreated")
if (!is.null(go_A)) write.csv(go_A, file.path(out_dir, "Fig3_panel_A_GO.csv"),
                              row.names = FALSE)

# ---- Panel B: GO for argk-2-/- vs N2 (treated) ------------------------------
de_geno_trt <- as.data.frame(res_geno_trt) %>%
  rownames_to_column("Gene") %>%
  filter(!is.na(padj), padj < go_padj_cut) %>%
  pull(Gene)
go_B <- run_go_BP(de_geno_trt, "Panel B: argk-2-/- vs N2 (treated)")
if (!is.null(go_B)) write.csv(go_B, file.path(out_dir, "Fig3_panel_B_GO.csv"),
                              row.names = FALSE)

# ---- GO dotplot builder -----------------------------------------------------
make_go_dotplot <- function(go_df, title, top_n) {
  d <- go_df %>%
    arrange(desc(GeneRatio_num)) %>%
    slice_head(n = top_n) %>%
    mutate(Description = factor(Description, levels = rev(Description)))

  ggplot(d, aes(x = GeneRatio_num, y = Description,
                size = Count, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "#C0392B", high = "#2C7FB8",
                         name = "p.adjust",
                         labels = function(x) format(x, scientific = TRUE, digits = 2)) +
    scale_size_continuous(range = c(2.5, 7), name = "Count") +
    labs(x = "GeneRatio", y = NULL, title = title) +
    theme_bw(base_size = 10) +
    theme(
      plot.title       = element_text(face = "italic", size = 10, hjust = 0.5),
      axis.text.y      = element_text(size = 8),
      axis.text.x      = element_text(size = 8),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray92", linewidth = 0.3),
      legend.title     = element_text(size = 7.5),
      legend.text      = element_text(size = 7),
      legend.key.size  = unit(0.35, "cm"),
      legend.spacing.y = unit(0.05, "cm")
    )
}

panel_A <- if (!is.null(go_A)) {
  make_go_dotplot(go_A, "N2 Treated vs. Untreated", PANEL_A_TOP) +
    labs(tag = "A") +
    theme(plot.tag = element_text(size = 16, face = "bold"),
          plot.tag.position = c(0.01, 0.99))
} else {
  ggplot() + theme_void() + labs(title = "Panel A: no GO terms", tag = "A")
}

panel_B <- if (!is.null(go_B)) {
  make_go_dotplot(go_B, expression(italic("argk")*"-2"^"-/-"*" vs. N2 (Treated)"),
                  PANEL_B_TOP) +
    labs(tag = "B") +
    theme(plot.tag = element_text(size = 16, face = "bold"),
          plot.tag.position = c(0.01, 0.99))
} else {
  ggplot() + theme_void() + labs(title = "Panel B: no GO terms", tag = "B")
}

# =============================================================================
# Panel C — Immune effect map
# =============================================================================
message("Building immune effect map…")

# Master per-gene table: 4 main LFCs + 4 padjs + interaction
master <- data.frame(
  Gene          = rownames(res_treat_N2),
  lfc_treat_N2  = res_treat_N2$log2FoldChange,
  lfc_treat_RB  = res_treat_RB$log2FoldChange,
  lfc_geno_unt  = res_geno_unt$log2FoldChange,
  lfc_geno_trt  = res_geno_trt$log2FoldChange,
  lfc_gxe       = res_gxe$log2FoldChange,
  padj_treat_N2 = res_treat_N2$padj,
  padj_treat_RB = res_treat_RB$padj,
  padj_geno_unt = res_geno_unt$padj,
  padj_geno_trt = res_geno_trt$padj,
  padj_gxe      = res_gxe$padj,
  stringsAsFactors = FALSE
)

# Classify each gene
classified <- master %>%
  mutate(
    sig_treat = (!is.na(padj_treat_N2) & padj_treat_N2 < de_padj_cut) |
                (!is.na(padj_treat_RB) & padj_treat_RB < de_padj_cut),
    sig_geno  = (!is.na(padj_geno_unt) & padj_geno_unt < de_padj_cut) |
                (!is.na(padj_geno_trt) & padj_geno_trt < de_padj_cut),
    sig_inter = !is.na(padj_gxe) & padj_gxe < de_padj_cut,
    category  = case_when(
      sig_inter | (sig_geno & sig_treat) ~ "Interaction",
      sig_geno  & !sig_treat             ~ "Genotype only",
      sig_treat & !sig_geno              ~ "Treatment only",
      TRUE                                ~ NA_character_
    ),
    mean_geno_lfc  = rowMeans(cbind(lfc_geno_unt, lfc_geno_trt), na.rm = TRUE),
    mean_treat_lfc = rowMeans(cbind(lfc_treat_N2, lfc_treat_RB), na.rm = TRUE)
  ) %>%
  filter(!is.na(category))

# Look up GO BP annotations for ALL classified genes, filter to immune terms
# (CSV-backed equivalent of the original AnnotationDbi::select on org.Ce.eg.db + GO.db)
classified$ENTREZ <- gene_to_entrez[classified$Gene]
entrez_classified <- na.omit(unique(classified$ENTREZ))

go_tbl <- go_map_full %>%
  filter(ENTREZID %in% entrez_classified, ONTOLOGY == "BP", !is.na(GO))

# GO term names from GO.db (apt-available), fall back to GO IDs as names
if (requireNamespace("GO.db", quietly = TRUE)) {
  term_tbl <- suppressMessages(AnnotationDbi::select(
    GO.db::GO.db,
    keys     = unique(go_tbl$GO),
    columns  = "TERM",
    keytype  = "GOID"
  ))
} else {
  term_tbl <- data.frame(GOID = unique(go_tbl$GO),
                         TERM = unique(go_tbl$GO),
                         stringsAsFactors = FALSE)
}

go_annot <- go_tbl %>%
  left_join(term_tbl, by = c("GO" = "GOID")) %>%
  filter(!is.na(TERM))

# ENTREZ IDs that have at least one immune-related GO BP term
immune_entrez <- go_annot %>%
  filter(str_detect(TERM, regex(immune_regex, ignore_case = TRUE))) %>%
  pull(ENTREZID) %>%
  unique()

immune_df <- classified %>%
  filter(ENTREZ %in% immune_entrez) %>%
  distinct(Gene, .keep_all = TRUE)

message(sprintf("Immune-annotated genes for Panel C: %d", nrow(immune_df)))
message(sprintf("  by category: %s",
                paste(names(table(immune_df$category)),
                      "=", table(immune_df$category), collapse = ", ")))

write.csv(immune_df, file.path(out_dir, "Fig3_panel_C_data.csv"), row.names = FALSE)

# Plot
pal_effect <- c(
  "Genotype only"  = "#009E73",
  "Treatment only" = "#0072B2",
  "Interaction"    = "#D55E00"
)

panel_C <- ggplot(immune_df,
                  aes(x = mean_geno_lfc, y = mean_treat_lfc, color = category)) +
  geom_hline(yintercept = 0, linetype = "22", linewidth = 0.3, color = "gray60") +
  geom_vline(xintercept = 0, linetype = "22", linewidth = 0.3, color = "gray60") +
  geom_point(alpha = 0.9, size = 2.2, stroke = 0) +
  scale_color_manual(values = pal_effect, name = NULL,
                     breaks = c("Genotype only", "Treatment only", "Interaction")) +
  labs(
    x = expression("Genotype effect (log"[2]*"FC)"),
    y = expression("Treat. effect (log"[2]*"FC)"),
    tag = "C"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.tag = element_text(size = 16, face = "bold"),
    plot.tag.position = c(0.01, 0.99),
    legend.position = c(0.99, 0.99),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = alpha("white", 0.95),
                                     color = NA),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray92", linewidth = 0.3),
    axis.title = element_text(size = 9)
  )

# =============================================================================
# Combine: A | (B / C)
# =============================================================================
right_col <- panel_B / panel_C + plot_layout(heights = c(1, 1))
final     <- panel_A | right_col
final     <- final + plot_layout(widths = c(1.15, 1))

ggsave(file.path(out_dir, "Fig3_recreated.pdf"), final,
       width = 13, height = 8, device = cairo_pdf)
ggsave(file.path(out_dir, "Fig3_recreated.png"), final,
       width = 13, height = 8, dpi = 300)

message("Done. Outputs in: ", out_dir)
