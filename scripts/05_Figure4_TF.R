# =============================================================================
# Fig4_combined.R — recreate Figure 4 end-to-end in R
#
# Builds a single 3-panel figure from raw counts:
#   A — grouped bar: % DE TF targets in N2 (TR vs UN) vs argk-2 KO (TR vs UN)
#   B — slope graph: Genotype / Treatment / Genotype + Treatment across top TFs
#   C — 2x2 bubble grid: top-10 TFs per contrast, size = # DE targets,
#                        color = proportion DE, x = INDEX (count × prop)
#
# DESeq2 model and contrast logic match TF_DEG_ANALYSIS_11_12_25.R exactly.
# Auto-detects input filenames (handles spaces, "copy 2", etc).
#
# Outputs (in ./Fig4_outputs/):
#   Fig4_recreated.pdf, Fig4_recreated.png         — combined 3-panel figure
#   <contrast>_by_proportion.csv                   — TF stats per contrast
#   <contrast>_by_count.csv
#   <contrast>_by_index.csv
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
out_dir  <- "outputs/figures/Figure4"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

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

# Simple-effect model for contrast extraction
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

# Combined "Genotype + Treatment": RB2060 treated vs N2 untreated
get_combined_RB_TR_vs_N2_UN <- function(dds_obj, genotype_alt = "RB2060") {
  rn <- resultsNames(dds_obj)
  main_geno <- rn[grepl("^genotype", rn) & grepl(genotype_alt, rn) & !grepl("treatment", rn)]
  main_trt  <- rn[grepl("^treatment", rn) & grepl("treated", rn)]
  inter     <- rn[grepl("genotype", rn) & grepl(genotype_alt, rn) & grepl("treatment.*treated", rn)]
  # E[RB2060_treated] - E[N2_untreated] = main_geno + main_trt + interaction
  results(dds_obj, list(c(main_geno, main_trt, inter)))
}

contrast_fns <- list(
  RB2060_UN_vs_N2     = function() get_genotype_within_treatment(dds_simple, "RB2060", "untreated"),
  RB2060_TR_vs_N2     = function() get_genotype_within_treatment(dds_simple, "RB2060", "treated"),
  N2_TR_vs_UN         = function() get_treatment_within_genotype(dds_simple, "N2"),
  RB2060_TR_vs_UN     = function() get_treatment_within_genotype(dds_simple, "RB2060"),
  RB2060_TR_vs_N2UN   = function() get_combined_RB_TR_vs_N2_UN(dds_simple, "RB2060")
)

# ---- Per-contrast TF summaries ---------------------------------------------
summarize_TFs <- function(res_obj, label) {
  res_df <- as.data.frame(res_obj) %>%
    rownames_to_column("Gene") %>%
    left_join(gene_to_target, by = "Gene") %>%
    filter(!is.na(padj))
  
  sig_targets <- res_df %>% filter(padj < padj_cut) %>% pull(UniprotID.Target)
  
  tf_summary <- tf_map %>%
    filter(UniprotID.Target %in% res_df$UniprotID.Target) %>%
    group_by(UniprotID.TF, Name.TF) %>%
    summarize(
      total_targets = n(),
      deg_targets   = sum(UniprotID.Target %in% sig_targets),
      prop_deg      = ifelse(total_targets > 0, deg_targets / total_targets, 0),
      .groups = "drop"
    ) %>%
    mutate(contrast = label, index = deg_targets * prop_deg)
  
  write.csv(arrange(tf_summary, desc(prop_deg)),
            file.path(out_dir, paste0(label, "_by_proportion.csv")), row.names = FALSE)
  write.csv(arrange(tf_summary, desc(deg_targets)),
            file.path(out_dir, paste0(label, "_by_count.csv")), row.names = FALSE)
  write.csv(arrange(tf_summary, desc(index)),
            file.path(out_dir, paste0(label, "_by_index.csv")), row.names = FALSE)
  tf_summary
}

message("Summarising TFs per contrast…")
summaries <- imap(contrast_fns, ~ summarize_TFs(.x(), .y))

# =============================================================================
# Panel A — grouped bar: %DE in reference vs argk-2 KO (TR vs UN)
# =============================================================================
fig_tfs <- c("blmp-1","cebp-1","ceh-60","daf-16","elt-2","fos-1","lin-35",
             "nhr-28","nhr-77","pha-4","pqm-1","sma-9","snpc-4","zip-2")

panelA_data <- bind_rows(
  summaries$N2_TR_vs_UN     %>% mutate(grp = "% DE in reference (TR vs UN)"),
  summaries$RB2060_TR_vs_UN %>% mutate(grp = "% DE in argk-2-/- (TR vs UN)")
) %>%
  mutate(Name.TF = tolower(Name.TF)) %>%
  filter(Name.TF %in% fig_tfs) %>%
  mutate(
    Name.TF = factor(toupper(Name.TF), levels = toupper(fig_tfs)),
    grp     = factor(grp, levels = c("% DE in reference (TR vs UN)",
                                     "% DE in argk-2-/- (TR vs UN)"))
  )

panel_A <- ggplot(panelA_data, aes(x = Name.TF, y = prop_deg, fill = grp)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  scale_fill_manual(values = c("#C97187", "#5B9BD5"),
                    labels = c("% DE in reference (TR vs UN)",
                               expression("% DE in " * italic("argk")*"-2"^"-/-" * " (TR vs UN)"))) +
  scale_y_continuous(limits = c(0, 0.55), breaks = c(0, 0.25, 0.5), expand = c(0, 0)) +
  labs(x = NULL, y = "Proportion of\nTF Targets", fill = NULL, tag = "A") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.title.y = element_text(size = 11),
    legend.position = c(0.22, 0.92),
    legend.background = element_rect(fill = alpha("white", 0.95), color = "lightgray"),
    legend.text = element_text(size = 8.5),
    legend.key.size = unit(0.4, "cm"),
    plot.tag = element_text(size = 18, face = "bold"),
    plot.tag.position = c(0.01, 0.98),
    panel.grid.major.y = element_line(color = "gray90", linetype = "dotted", linewidth = 0.3)
  )

# =============================================================================
# Panel B — slope: Genotype → Treatment → Genotype + Treatment
# =============================================================================
slope_tfs <- c("zip-2","cebp-1","ceh-60","elt-2","pqm-1","nhr-28",
               "fos-1","daf-16","nhr-77","pha-4","blmp-1")

panelB_data <- bind_rows(
  summaries$RB2060_UN_vs_N2   %>% mutate(stage = "Genotype",             stage_x = 1),
  summaries$N2_TR_vs_UN       %>% mutate(stage = "Treatment",            stage_x = 2),
  summaries$RB2060_TR_vs_N2UN %>% mutate(stage = "Genotype + Treatment", stage_x = 3)
) %>%
  mutate(Name.TF = tolower(Name.TF)) %>%
  filter(Name.TF %in% slope_tfs) %>%
  mutate(Name.TF_upper = toupper(Name.TF))

# Right-side label positions: spread vertically to avoid overlap
end_pts <- panelB_data %>%
  filter(stage_x == 3) %>%
  arrange(desc(prop_deg)) %>%
  select(Name.TF, Name.TF_upper, prop_deg) %>%
  mutate(y_raw = prop_deg)

n_lab <- nrow(end_pts)
y_lo <- min(end_pts$y_raw) - 0.02
y_hi <- max(end_pts$y_raw) + 0.04
end_pts$y_lab <- seq(y_hi, y_lo, length.out = n_lab)
end_pts$x_lab <- 3.15
end_pts$x_connect <- 3.05

# Stable color order tied to slope_tfs vector
panelB_data$Name.TF <- factor(panelB_data$Name.TF, levels = slope_tfs)
slope_palette <- c("#1f77b4","#d62728","#2ca02c","#9467bd","#e377c2",
                   "#bcbd22","#17becf","#ff7f0e","#8c564b","#7f7f7f","#aec7e8")
names(slope_palette) <- slope_tfs

panel_B <- ggplot(panelB_data, aes(x = stage_x, y = prop_deg,
                                   color = Name.TF, group = Name.TF)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.4) +
  # connector lines from final data point to label
  geom_segment(data = end_pts,
               aes(x = 3, xend = x_connect, y = y_raw, yend = y_lab),
               inherit.aes = FALSE,
               color = "gray70", linewidth = 0.4) +
  geom_text(data = end_pts,
            aes(x = x_lab, y = y_lab, label = Name.TF_upper),
            inherit.aes = FALSE,
            hjust = 0, size = 3.1, color = "black") +
  scale_color_manual(values = slope_palette, guide = "none") +
  scale_x_continuous(breaks = 1:3,
                     labels = c("Genotype", "Treatment", "Genotype + Treatment"),
                     limits = c(0.8, 3.85), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.04, max(0.6, y_hi + 0.02)), expand = c(0, 0)) +
  labs(x = NULL, y = "Proportion of DE TF targets", tag = "B") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x  = element_text(size = 10),
    axis.title.y = element_text(size = 11),
    plot.tag = element_text(size = 18, face = "bold"),
    plot.tag.position = c(0.01, 0.98),
    panel.grid.major.y = element_line(color = "gray90", linetype = "dotted", linewidth = 0.3)
  )

# =============================================================================
# Panel C — 2x2 bubble grid: top-10 TFs per contrast by INDEX
# Shared size & color scales across all four contrasts so visual magnitudes
# can be compared directly. Limits are computed from the pooled top-10 of
# each contrast (i.e. the points actually plotted).
# =============================================================================
panel_c_contrasts <- list(
  summaries$RB2060_UN_vs_N2,
  summaries$N2_TR_vs_UN,
  summaries$RB2060_TR_vs_N2,
  summaries$RB2060_TR_vs_UN
)
panel_c_top10 <- bind_rows(lapply(panel_c_contrasts, function(d) {
  d %>% arrange(desc(index)) %>% slice_head(n = 10)
}))
size_limits  <- range(panel_c_top10$deg_targets, na.rm = TRUE)
color_limits <- range(panel_c_top10$prop_deg,    na.rm = TRUE)
message(sprintf("Panel C shared size limits:  %d – %d  # DE targets",
                size_limits[1],  size_limits[2]))
message(sprintf("Panel C shared color limits: %.3f – %.3f Proportion DE",
                color_limits[1], color_limits[2]))

make_bubble <- function(df, title, size_limits, color_limits, show_tag = FALSE) {
  topI <- df %>% arrange(desc(index)) %>% slice_head(n = 10) %>%
    mutate(Name.TF_upper = toupper(Name.TF))
  
  p <- ggplot(topI, aes(x = pmax(index, 0.8), y = reorder(Name.TF_upper, index))) +
    geom_point(aes(size = deg_targets, color = prop_deg),
               alpha = 0.95, stroke = 0.3) +
    scale_size_continuous(range = c(1.8, 10),
                          limits = size_limits,
                          name = "# DE targets") +
    scale_color_gradient(low = "#deebf7", high = "#08306b",
                         limits = color_limits,
                         name = "Proportion DE") +
    scale_x_log10(limits = c(0.7, 4000),
                  breaks = c(1, 10, 100, 1000),
                  labels = c("1E00", "1E01", "1E02", "1E03")) +
    labs(title = title, x = "INDEX", y = NULL) +
    theme_classic(base_size = 9) +
    theme(
      plot.title   = element_text(face = "italic", size = 9.5, hjust = 0.5),
      axis.text.y  = element_text(size = 7.5),
      axis.text.x  = element_text(size = 7),
      axis.title.x = element_text(size = 8),
      legend.title = element_text(size = 6.5),
      legend.text  = element_text(size = 6),
      legend.key.size = unit(0.3, "cm"),
      legend.spacing.y = unit(0.05, "cm"),
      legend.box.spacing = unit(0.1, "cm"),
      panel.grid.major.x = element_line(color = "gray92", linetype = "dotted", linewidth = 0.3)
    )
  if (show_tag) p <- p + labs(tag = "C") +
    theme(plot.tag = element_text(size = 18, face = "bold"),
          plot.tag.position = c(-0.05, 1.05))
  p
}

c_topleft  <- make_bubble(summaries$RB2060_UN_vs_N2,
                          expression(italic("argk")*"-2"^"-/-"*" vs reference (untreated)"),
                          size_limits, color_limits, show_tag = TRUE)
c_topright <- make_bubble(summaries$N2_TR_vs_UN,
                          "treated vs untreated (reference)",
                          size_limits, color_limits)
c_botleft  <- make_bubble(summaries$RB2060_TR_vs_N2,
                          expression(italic("argk")*"-2"^"-/-"*" vs reference (treated)"),
                          size_limits, color_limits)
c_botright <- make_bubble(summaries$RB2060_TR_vs_UN,
                          expression("treated vs untreated ("*italic("argk")*"-2"^"-/-"*")"),
                          size_limits, color_limits)

# Combine into a 2x2 and collect the (now identical) guides into one legend
panel_C <- wrap_plots(c_topleft, c_topright, c_botleft, c_botright, ncol = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right",
        legend.box = "vertical",
        legend.spacing.y = unit(0.1, "cm"))

# =============================================================================
# Combine A + B + C
# =============================================================================
left_col <- panel_A / panel_B
final <- left_col | panel_C
final <- final + plot_layout(widths = c(1.05, 1.0))

ggsave(file.path(out_dir, "Fig4_recreated.pdf"), final,
       width = 15, height = 9, device = cairo_pdf)
ggsave(file.path(out_dir, "Fig4_recreated.png"), final,
       width = 15, height = 9, dpi = 300)

message("Done. Outputs in: ", out_dir)
