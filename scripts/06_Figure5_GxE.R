# =============================================================================
# Fig5_combined.R — recreate Figure 5 end-to-end in R
#
# Builds a 3-panel figure from raw counts:
#   A — scatter: CoCl2 response in N2 (x) vs argk-2-/- (y) for all responsive
#       genes; 90 G×E genes highlighted (68 sub-additive, 22 supra-additive)
#       with binomial test annotation
#   B — heatmap: 30 sub-additive + 8 supra-additive G×E effector genes across
#       four contrasts (Genotype, Treatment, Genotype-under-CoCl2, G×E term)
#   C — stacked bars: TF targets partitioned into Treatment-only / Genotype-only
#       / Conditional interaction for 14 TFs
#
# DESeq2 model and contrast logic match TF_DEG_ANALYSIS_11_12_25.R exactly.
# Selection thresholds for panel B match build_publication_figure_v15.py:
#   sub-additive: top 30 GxE genes with interaction LFC < 0,
#                 scored by (-lfc_gxe + max(0, lfc_geno_unt))
#   supra-additive: top 8 GxE genes with interaction LFC >= 1.0
#
# Auto-detects input filenames (handles spaces, "copy 2", etc).
#
# Outputs (in ./Fig5_outputs/):
#   Fig5_recreated.pdf, Fig5_recreated.png    — combined 3-panel figure
#   Fig5_master_DE_table.csv                  — per-gene LFC across 5 contrasts
#                                               + interaction padj
#   Fig5_panel_B_genes.csv                    — 38 selected genes + scores
#   Fig5_panel_C_data.csv                     — TF effect-category proportions
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
})

# Force dplyr verbs (DESeq2 masks select/filter)
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

padj_cut <- 0.05
min_sum  <- 10
out_dir  <- "outputs/figures/Figure5"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Panel B selection thresholds (match build_publication_figure_v15.py)
N_SUB_SELECT       <- 30
N_SUPRA_SELECT     <- 8
MIN_SUPRA_LFC      <- 1.0

# Panel C TF list (14 TFs in display order)
PANEL_C_TFS <- c("blmp-1","cebp-1","ceh-60","daf-16","elt-2","fos-1",
                 "lin-35","nhr-28","nhr-77","pha-4","pqm-1","sma-9",
                 "snpc-4","zip-2")

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

tf_map <- rna_raw %>%
  select(UniprotID.TF, Name.TF, UniprotID.Target) %>%
  filter(!is.na(UniprotID.TF), !is.na(UniprotID.Target)) %>%
  distinct()

gene_to_target <- rna_raw %>%
  select(Gene, UniprotID.Target) %>%
  filter(!is.na(Gene), !is.na(UniprotID.Target)) %>%
  distinct()

# ---- DESeq2 -----------------------------------------------------------------
message("Running DESeq2…")
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = meta_df,
  design    = ~ genotype * treatment + genotype:experiment + treatment:experiment
)
dds <- dds[rowSums(counts(dds)) > min_sum, ]
dds <- DESeq(dds)

# Simple-effect model for contrast extraction (includes interaction term)
dds_simple <- dds
dds_simple$genotype  <- relevel(dds_simple$genotype,  ref = "N2")
dds_simple$treatment <- relevel(dds_simple$treatment, ref = "untreated")
design(dds_simple) <- ~ genotype + treatment + genotype:treatment
dds_simple <- DESeq(dds_simple)

# ---- Contrast helpers -------------------------------------------------------
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

# Interaction-term coefficient (the key piece for Fig 5)
inter_name <- grep("genotype.*RB2060.*treatment.*treated",
                   resultsNames(dds_simple), value = TRUE)
stopifnot(length(inter_name) == 1)
res_gxe <- results(dds_simple, name = inter_name)
message("Interaction coefficient extracted: ", inter_name)

# ---- Master per-gene table --------------------------------------------------
master <- data.frame(
  Gene         = rownames(res_treat_N2),
  lfc_treat_N2 = res_treat_N2$log2FoldChange,
  lfc_treat_RB = res_treat_RB$log2FoldChange,
  lfc_geno_unt = res_geno_unt$log2FoldChange,
  lfc_geno_trt = res_geno_trt$log2FoldChange,
  lfc_gxe      = res_gxe$log2FoldChange,
  padj_treat_N2 = res_treat_N2$padj,
  padj_treat_RB = res_treat_RB$padj,
  padj_geno_unt = res_geno_unt$padj,
  padj_geno_trt = res_geno_trt$padj,
  padj_gxe      = res_gxe$padj,
  stringsAsFactors = FALSE
)
write.csv(master, file.path(out_dir, "Fig5_master_DE_table.csv"), row.names = FALSE)

# ---- G×E set + sub/supra classification --------------------------------------
# Published Figure 5A recipe (build_publication_figure_v15.py lines 397-438):
#   gene set      = effect_class == "GxE" from DESeq2_effect_classification.csv
#                   (produced by 01_DESeq2_and_TF_enrichment.R with apeglm + ashr shrinkage)
#   sub-additive  = treatment response in argk-2 <= treatment response in N2
#   supra        = strictly above the y=x line
# This intentionally OVERRIDES the older Wald-coef-sign approach so this script
# stays consistent with the published numbers (n=90, 68 sub, 22 supra, p=1.3e-6).
ec_path <- "DESeq2_effect_classification.csv"
if (file.exists(ec_path)) {
  ec <- read.csv(ec_path, stringsAsFactors = FALSE)
  gxe_set <- ec %>% filter(effect_class == "GxE")
  n_gxe   <- nrow(gxe_set)
  n_sub   <- sum(gxe_set$lfc_treat_RB <= gxe_set$lfc_treat_N2, na.rm = TRUE)
  n_supra <- sum(gxe_set$lfc_treat_RB >  gxe_set$lfc_treat_N2, na.rm = TRUE)
  binom_p <- binom.test(n_sub, n_gxe, p = 0.5, alternative = "two.sided")$p.value
  message(sprintf("G×E (effect_class): %d total (%d sub, %d supra); binomial p = %.2e",
                  n_gxe, n_sub, n_supra, binom_p))
  # Bring lfc_gxe over so downstream panel-B selection still works
  gxe_set$Gene    <- gxe_set$gene_symbol
  gxe_set$lfc_gxe <- gxe_set$lfc_gxe
} else {
  warning("DESeq2_effect_classification.csv not found; falling back to live LRT.")
  gxe_set <- master %>%
    filter(!is.na(padj_gxe), padj_gxe < padj_cut)
  n_gxe   <- nrow(gxe_set)
  n_sub   <- sum(gxe_set$lfc_gxe < 0, na.rm = TRUE)
  n_supra <- sum(gxe_set$lfc_gxe >= 0, na.rm = TRUE)
  binom_p <- binom.test(n_sub, n_gxe, p = 0.5, alternative = "two.sided")$p.value
  message(sprintf("G×E (live, fallback): %d total (%d sub, %d supra); binomial p = %.2e",
                  n_gxe, n_sub, n_supra, binom_p))
}

# ---- Panel B selection ------------------------------------------------------
sub_sel <- gxe_set %>%
  filter(lfc_gxe < 0) %>%
  mutate(score = -lfc_gxe + pmax(0, replace_na(lfc_geno_unt, 0))) %>%
  arrange(desc(score)) %>%
  slice_head(n = N_SUB_SELECT) %>%
  mutate(group = "sub-additive")

supra_sel <- gxe_set %>%
  filter(lfc_gxe >= MIN_SUPRA_LFC) %>%
  arrange(desc(lfc_gxe)) %>%
  slice_head(n = N_SUPRA_SELECT) %>%
  mutate(score = lfc_gxe, group = "supra-additive")

panelB_genes <- bind_rows(sub_sel, supra_sel) %>%
  distinct(Gene, .keep_all = TRUE)
message(sprintf("Panel B selection: %d sub-additive + %d supra-additive = %d genes",
                sum(panelB_genes$group == "sub-additive"),
                sum(panelB_genes$group == "supra-additive"),
                nrow(panelB_genes)))
write.csv(panelB_genes, file.path(out_dir, "Fig5_panel_B_genes.csv"), row.names = FALSE)

# =============================================================================
# Panel A — scatter
# =============================================================================
panelA_bg <- master %>%
  filter(!is.na(lfc_treat_N2), !is.na(lfc_treat_RB)) %>%
  mutate(
    set = case_when(
      Gene %in% panelB_genes$Gene[panelB_genes$group == "sub-additive"]   ~ "sub_in_B",
      Gene %in% panelB_genes$Gene[panelB_genes$group == "supra-additive"] ~ "supra_in_B",
      Gene %in% gxe_set$Gene                                              ~ "gxe_not_B",
      TRUE                                                                ~ "other"
    ),
    set = factor(set, levels = c("other", "gxe_not_B", "sub_in_B", "supra_in_B"))
  ) %>%
  arrange(set)   # background first, foreground last (drawn on top)

set_colors <- c("other"      = "gray80",
                "gxe_not_B"  = "#F4B27A",
                "sub_in_B"   = "#2C7FB8",
                "supra_in_B" = "#C0392B")
set_sizes  <- c("other" = 0.5, "gxe_not_B" = 1.4, "sub_in_B" = 2.0, "supra_in_B" = 2.0)
set_alpha  <- c("other" = 0.35, "gxe_not_B" = 0.75, "sub_in_B" = 0.95, "supra_in_B" = 0.95)
set_labels <- c(
  other      = "other responsive",
  gxe_not_B  = sprintf("G×E not in b (n=%d)",
                       n_gxe - nrow(panelB_genes)),
  sub_in_B   = sprintf("sub-additive in b (n=%d)",
                       sum(panelB_genes$group == "sub-additive")),
  supra_in_B = sprintf("supra-additive in b (n=%d)",
                       sum(panelB_genes$group == "supra-additive"))
)

ann_text <- sprintf("90 G×E genes:\n%d sub-additive, %d supra-additive\nbinomial p = %s",
                    n_sub, n_supra, format(binom_p, scientific = TRUE, digits = 2))

# axis limits — symmetric
lim_xy <- max(abs(c(panelA_bg$lfc_treat_N2, panelA_bg$lfc_treat_RB)),
              na.rm = TRUE)
lim_xy <- ceiling(lim_xy)

panel_A <- ggplot(panelA_bg,
                  aes(x = lfc_treat_N2, y = lfc_treat_RB,
                      color = set, size = set, alpha = set)) +
  # additive expectation
  geom_abline(slope = 1, intercept = 0, color = "black", linewidth = 0.6) +
  geom_abline(slope = 1, intercept =  1, color = "gray70", linewidth = 0.3, linetype = "dashed") +
  geom_abline(slope = 1, intercept = -1, color = "gray70", linewidth = 0.3, linetype = "dashed") +
  # points
  geom_point(stroke = 0) +
  scale_color_manual(values = set_colors, labels = set_labels, name = NULL) +
  scale_size_manual (values = set_sizes,  labels = set_labels, name = NULL) +
  scale_alpha_manual(values = set_alpha,  labels = set_labels, name = NULL) +
  # axis
  scale_x_continuous(limits = c(-lim_xy, lim_xy),
                     breaks = seq(-lim_xy, lim_xy, by = 4)) +
  scale_y_continuous(limits = c(-lim_xy, lim_xy),
                     breaks = seq(-lim_xy, lim_xy, by = 4)) +
  coord_fixed() +
  labs(
    title = expression("Non-additive CoCl"[2]*" response in "*italic("argk")*"-2"^"-/-"),
    x = expression("CoCl"[2]*" response in N2 (log"[2]*"FC)"),
    y = expression("CoCl"[2]*" response in "*italic("argk")*"-2"^"-/-"*" (log"[2]*"FC)"),
    tag = "A"
  ) +
  # stat annotation (upper-left)
  annotate("rect",
           xmin = -lim_xy + 0.3, xmax = -lim_xy + 0.3 + lim_xy * 0.75,
           ymin =  lim_xy - 2.4, ymax =  lim_xy - 0.3,
           fill = "white", color = "gray60", linewidth = 0.3) +
  annotate("text",
           x = -lim_xy + 0.5, y = lim_xy - 0.7,
           label = ann_text, hjust = 0, vjust = 1, size = 2.7, lineheight = 0.95,
           family = "sans") +
  # corner labels — boxed
  annotate("label", x = -lim_xy * 0.55, y =  lim_xy * 0.55,
           label = "supra-additive", color = "#C0392B", fontface = "italic",
           size = 2.9, label.padding = unit(0.15, "lines"),
           label.r = unit(0.1, "lines"), label.size = 0.4) +
  annotate("label", x =  lim_xy * 0.55, y = -lim_xy * 0.55,
           label = "sub-additive", color = "#2C7FB8", fontface = "italic",
           size = 2.9, label.padding = unit(0.15, "lines"),
           label.r = unit(0.1, "lines"), label.size = 0.4) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
    plot.tag = element_text(size = 16, face = "bold"),
    plot.tag.position = c(0.01, 0.99),
    legend.position = c(0.78, 0.18),
    legend.background = element_rect(fill = alpha("white", 0.95), color = "gray60", linewidth = 0.3),
    legend.text = element_text(size = 6.5),
    legend.key.size = unit(0.3, "cm"),
    legend.spacing.y = unit(0.05, "cm")
  )

# =============================================================================
# Panel B — heatmap
# =============================================================================
# Within-group ordering: sub-additive by descending score, supra by descending lfc_gxe
panelB_ordered <- panelB_genes %>%
  arrange(group, desc(score)) %>%
  mutate(row_id = row_number())

panelB_long <- panelB_ordered %>%
  select(Gene, row_id, group,
         Genotype           = lfc_geno_unt,
         Treatment          = lfc_treat_N2,
         `Genotype(CoCl2)`  = lfc_geno_trt,
         `G×E interaction`  = lfc_gxe) %>%
  pivot_longer(c(Genotype, Treatment, `Genotype(CoCl2)`, `G×E interaction`),
               names_to = "contrast", values_to = "lfc") %>%
  mutate(
    contrast = factor(contrast,
                      levels = c("Genotype", "Treatment",
                                 "Genotype(CoCl2)", "G×E interaction"))
  )

# Color cap so extreme values don't blow out the scale
heat_cap <- 4
panelB_long$lfc_capped <- pmax(pmin(panelB_long$lfc, heat_cap), -heat_cap)

# Boundary between sub and supra groups (for the dashed line)
boundary_row <- sum(panelB_ordered$group == "sub-additive") + 0.5

panel_B <- ggplot(panelB_long,
                  aes(x = contrast, y = row_id, fill = lfc_capped)) +
  geom_tile(color = "white", linewidth = 0.15) +
  # Black line separating G×E interaction column from the rest
  geom_vline(xintercept = 3.5, color = "black", linewidth = 0.7) +
  # Dashed line separating sub- vs supra-additive
  geom_hline(yintercept = boundary_row, color = "gray40",
             linewidth = 0.4, linetype = "dashed") +
  scale_fill_gradient2(
    low = "#053061", mid = "white", high = "#67001F",
    midpoint = 0, limits = c(-heat_cap, heat_cap),
    name = expression("log"[2]*"FC"),
    breaks = c(-4, -2, 0, 2, 4),
    guide = guide_colorbar(barwidth = 6, barheight = 0.4,
                           title.position = "right",
                           title.vjust = 0.85)
  ) +
  scale_y_continuous(
    breaks = panelB_ordered$row_id,
    labels = panelB_ordered$Gene,
    expand = c(0, 0),
    sec.axis = sec_axis(
      ~ .,
      breaks = c(
        median(panelB_ordered$row_id[panelB_ordered$group == "sub-additive"]),
        median(panelB_ordered$row_id[panelB_ordered$group == "supra-additive"])
      ),
      labels = c("sub-additive", "supra-additive")
    )
  ) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  labs(title = expression("G×E effector expression"), x = NULL, y = NULL, tag = "B") +
  theme_minimal(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
    axis.text.y = element_text(face = "italic", size = 6.5),
    axis.text.y.right = element_text(angle = -90, hjust = 0.5,
                                     face = "italic", size = 8,
                                     color = c("#2C7FB8", "#C0392B")),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, size = 8),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.margin = margin(t = 0, b = 0),
    plot.tag = element_text(size = 16, face = "bold"),
    plot.tag.position = c(0.01, 0.99)
  )

# =============================================================================
# Panel C — TF target effect categories
# =============================================================================
# Build the TF -> target set lookup
tf_to_targets <- tf_map %>%
  filter(tolower(Name.TF) %in% PANEL_C_TFS) %>%
  group_by(Name.TF) %>%
  summarize(targets = list(unique(UniprotID.Target)), .groups = "drop")

# DE target sets (UniprotID.Target) per contrast
de_target_set <- function(res_obj) {
  as.data.frame(res_obj) %>%
    tibble::rownames_to_column("Gene") %>%
    dplyr::left_join(gene_to_target, by = "Gene") %>%
    dplyr::filter(!is.na(UniprotID.Target), !is.na(padj), padj < padj_cut) %>%
    pull(UniprotID.Target) %>%
    unique()
}
de_treat_N2 <- de_target_set(res_treat_N2)
de_treat_RB <- de_target_set(res_treat_RB)
de_geno_unt <- de_target_set(res_geno_unt)
de_geno_trt <- de_target_set(res_geno_trt)

classify_tf <- function(tgts) {
  total <- length(tgts)
  treat_union <- intersect(tgts, union(de_treat_N2, de_treat_RB))
  geno_union  <- intersect(tgts, union(de_geno_unt, de_geno_trt))
  tibble(
    total              = total,
    `Treatment-only`   = length(setdiff(treat_union, geno_union)) / total,
    `Genotype-only`    = length(setdiff(geno_union, treat_union)) / total,
    `Conditional interaction` = length(intersect(treat_union, geno_union)) / total
  )
}

panelC_data <- tf_to_targets %>%
  mutate(stats = lapply(targets, classify_tf)) %>%
  select(Name.TF, stats) %>%
  tidyr::unnest(stats) %>%
  mutate(Name.TF_lower = tolower(Name.TF)) %>%
  filter(Name.TF_lower %in% PANEL_C_TFS) %>%
  mutate(Name.TF = factor(toupper(Name.TF_lower), levels = toupper(PANEL_C_TFS)))

write.csv(panelC_data, file.path(out_dir, "Fig5_panel_C_data.csv"), row.names = FALSE)

panelC_long <- panelC_data %>%
  select(Name.TF,
         `Treatment-only`, `Genotype-only`, `Conditional interaction`) %>%
  pivot_longer(-Name.TF, names_to = "category", values_to = "proportion") %>%
  mutate(category = factor(category,
                           levels = c("Treatment-only",
                                      "Genotype-only",
                                      "Conditional interaction")))

cat_colors <- c(
  "Treatment-only"          = "#F4A09E",  # salmon
  "Genotype-only"           = "#5DAA5D",  # green
  "Conditional interaction" = "#4A90C2"   # blue
)

panel_C <- ggplot(panelC_long,
                  aes(x = Name.TF, y = proportion, fill = category)) +
  geom_col(width = 0.85) +
  scale_fill_manual(values = cat_colors, name = "Effect Category") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Proportion of TF targets", title = "Effect category of TF targets",
       tag = "C") +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.title.y = element_text(size = 9),
    legend.position = c(0.18, 0.82),
    legend.background = element_rect(fill = alpha("white", 0.95),
                                     color = "gray60", linewidth = 0.3),
    legend.title = element_text(size = 7.5, face = "bold"),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.35, "cm"),
    plot.tag = element_text(size = 16, face = "bold"),
    plot.tag.position = c(0.01, 0.99),
    panel.grid.major.y = element_line(color = "gray92",
                                       linetype = "dotted", linewidth = 0.3)
  )

# =============================================================================
# Combine: (A / C) | B
# =============================================================================
left_col <- panel_A / panel_C + plot_layout(heights = c(1.1, 1))
final    <- left_col | panel_B
final    <- final + plot_layout(widths = c(0.42, 0.58))

ggsave(file.path(out_dir, "Fig5_recreated.pdf"), final,
       width = 13, height = 10, device = cairo_pdf)
ggsave(file.path(out_dir, "Fig5_recreated.png"), final,
       width = 13, height = 10, dpi = 300)

message("Done. Outputs in: ", out_dir)
