#!/usr/bin/env Rscript
# =============================================================================
# 03_Figure2_KEGG.R
# =============================================================================
#
# PURPOSE
#   Render Figure 2 -- KEGG functional composition of the CoCl2 stress response
#   across the four canonical contrasts.
#
# INPUTS (all produced by 01_DESeq2_and_TF_enrichment.R, in outputs/DESeq2/)
#   DESeq2_treatment_in_N2.csv         -- N2: treated vs untreated
#   DESeq2_treatment_in_RB2060.csv     -- argk-2: treated vs untreated
#   DESeq2_genotype_untreated.csv      -- argk-2 vs N2 (untreated baseline)
#   DESeq2_genotype_treated.csv        -- argk-2 vs N2 (treated)
#   DESeq2_GxE_interaction.csv         -- genotype x treatment interaction term
#
# Plus the offline KEGG cel mappings (in data/reference/):
#   kegg_cel_pathway_to_gene.csv
#   kegg_cel_pathway_names.csv
#
# DEG-calling filter throughout: padj < 0.05 AND |log2FoldChange| > 1
# (matches the Figure 1 Venn filter and is consistent with the full pipeline)
#
# OUTPUTS
#   outputs/figures/Figure2/Figure2_KEGG_recreated.pdf   -- 3-panel figure
#   outputs/figures/Figure2/Figure2_KEGG_recreated.png   -- raster preview
#   outputs/DESeq2/KEGG_Functional_Composition_All_Contrasts.csv
#                                                        -- long-format DEG x pathway
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
  library(stringr); library(purrr); library(ggplot2)
  library(patchwork)
})
for (utils_path in c("scripts/_utils.R", "_utils.R", "../scripts/_utils.R")) {
  if (file.exists(utils_path)) { source(utils_path); break }
}
stopifnot(exists("load_kegg_offline"))

DESEQ_DIR <- "outputs/DESeq2"
FIG_DIR   <- "outputs/figures/Figure2"
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

# ---- 1. Load 03's contrast tables and label them with friendly contrast names ----
read_contrast <- function(filename, contrast_label) {
  path <- file.path(DESEQ_DIR, filename)
  if (!file.exists(path)) {
    stop("Missing ", path, "\n  Run 01_DESeq2_and_TF_enrichment.R first.")
  }
  df <- read.csv(path, stringsAsFactors = FALSE)
  # First column is the unnamed rownames col -> gene_symbol fallback
  if (colnames(df)[1] == "X" || colnames(df)[1] == "") {
    colnames(df)[1] <- "gene_symbol_rn"
  }
  if (!"gene_symbol" %in% colnames(df)) df$gene_symbol <- df$gene_symbol_rn
  df$contrast <- contrast_label
  df
}

deg_n2_trt       <- read_contrast("DESeq2_treatment_in_N2.csv",     "N2: Treated vs Untreated")
deg_rb2060_trt   <- read_contrast("DESeq2_treatment_in_RB2060.csv", "argk-2: Treated vs Untreated")
deg_geno_unt     <- read_contrast("DESeq2_genotype_untreated.csv",  "argk-2 vs N2 (Untreated)")
deg_geno_trt     <- read_contrast("DESeq2_genotype_treated.csv",    "argk-2 vs N2 (Treated)")
deg_interaction  <- read_contrast("DESeq2_GxE_interaction.csv",     "Interaction (G x E)")

message(sprintf("Loaded contrasts (rows): N2=%d, argk-2=%d, geno_unt=%d, geno_trt=%d, int=%d",
                nrow(deg_n2_trt), nrow(deg_rb2060_trt),
                nrow(deg_geno_unt), nrow(deg_geno_trt), nrow(deg_interaction)))

# ---- 2. Build KEGG gene -> pathway map ----
kegg_paths <- load_kegg_offline()
pathway_names <- as_tibble(kegg_paths$KEGGPATHID2NAME) |>
  setNames(c("pathway_id", "pathway_name"))
gene2path <- as_tibble(kegg_paths$KEGGPATHID2EXTID) |>
  setNames(c("pathway_id", "kegg_gene")) |>
  inner_join(pathway_names, by = "pathway_id") |>
  mutate(kegg_gene = str_remove(kegg_gene, "^CELE_"))
message(sprintf("KEGG: %d unique pathways, %d gene-pathway edges",
                nrow(pathway_names), nrow(gene2path)))

# ---- 3. Compose DEG x KEGG cross-table ----
# Filter: padj < 0.05 AND |log2FoldChange| > 1
make_deg_kegg <- function(df) {
  df |>
    transmute(Gene = gene_symbol,
              log2FC = log2FoldChange,
              padj   = padj,
              contrast = contrast) |>
    filter(!is.na(log2FC), !is.na(padj),
           padj < 0.05, abs(log2FC) > 1) |>
    inner_join(gene2path, by = c("Gene" = "kegg_gene"),
               relationship = "many-to-many")
}

deg_kegg_3 <- bind_rows(
  make_deg_kegg(deg_n2_trt),
  make_deg_kegg(deg_rb2060_trt),
  make_deg_kegg(deg_geno_unt)
)
deg_kegg_4 <- bind_rows(
  deg_kegg_3,
  make_deg_kegg(deg_geno_trt)
)
message(sprintf("DEGs matched to KEGG: 3-contrast=%d, 4-contrast=%d",
                nrow(deg_kegg_3), nrow(deg_kegg_4)))

# Save canonical KEGG composition CSV (3-contrast version is the headline)
write.csv(deg_kegg_3,
          file.path(DESEQ_DIR, "KEGG_Functional_Composition_All_Contrasts.csv"),
          row.names = FALSE)
message("Wrote: ", file.path(DESEQ_DIR, "KEGG_Functional_Composition_All_Contrasts.csv"))

# Reusable shortcut for Panel A code below
deg_kegg <- deg_kegg_3

# ---- 4. Panel A: KEGG violin (3 contrasts) ----
message("Building Panel A...")

pathway_order_A <- deg_kegg |>
  group_by(contrast, pathway_name) |>
  filter(n() >= 5) |>
  group_by(pathway_name) |>
  summarise(mean_abs = mean(abs(log2FC)), .groups = "drop") |>
  arrange(desc(mean_abs)) |>
  slice_head(n = 12) |>
  pull(pathway_name)

plot_A_data <- deg_kegg |>
  group_by(contrast, pathway_name) |>
  filter(n() >= 3) |>
  ungroup() |>
  filter(pathway_name %in% pathway_order_A) |>
  mutate(
    pathway_name = factor(pathway_name, levels = rev(pathway_order_A)),
    contrast = case_when(
      contrast == "argk-2 vs N2 (Untreated)"     ~ "argk-2-/- vs N2\n(Untreated)",
      contrast == "N2: Treated vs Untreated"     ~ "N2: Treated vs\nUntreated",
      contrast == "argk-2: Treated vs Untreated" ~ "argk-2-/-: Treated\nvs Untreated",
      TRUE ~ contrast
    )
  )

colors_A <- c(
  "argk-2-/- vs N2\n(Untreated)"      = "#4E79A7",
  "N2: Treated vs\nUntreated"          = "#E15759",
  "argk-2-/-: Treated\nvs Untreated"   = "#59A14F"
)
plot_A_data$contrast <- factor(plot_A_data$contrast, levels = names(colors_A))

panel_A <- ggplot(plot_A_data, aes(x = log2FC, y = pathway_name, fill = contrast)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.85, linewidth = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.3) +
  scale_fill_manual(values = colors_A) +
  facet_wrap(~ contrast, ncol = 1) +
  theme_classic(base_size = 8) +
  theme(
    axis.text.y      = element_text(size = 7, face = "bold", color = "black"),
    axis.text.x      = element_text(size = 7, color = "black"),
    axis.title.x     = element_text(size = 8, face = "bold", color = "black"),
    axis.title.y     = element_blank(),
    strip.text       = element_text(face = "bold", size = 7, color = "black"),
    strip.background = element_rect(fill = "grey95", color = NA),
    legend.position  = "none",
    plot.margin      = margin(2, 4, 2, 2),
    panel.spacing    = unit(0.15, "cm"),
    axis.line        = element_line(linewidth = 0.3),
    axis.ticks       = element_line(linewidth = 0.2)
  ) +
  labs(x = expression(log[2]~Fold~Change))

# ---- 5. Panel B: Coverage heatmap (top 15 pathways x 4 contrasts) ----
message("Building Panel B...")

panelB_data <- deg_kegg_4 |>
  group_by(pathway_name, contrast) |>
  summarise(n_deg = n_distinct(Gene), .groups = "drop")

top15_paths <- panelB_data |>
  group_by(pathway_name) |>
  summarise(total = sum(n_deg), .groups = "drop") |>
  arrange(desc(total)) |>
  slice_head(n = 15) |>
  pull(pathway_name)

col_labels <- c("N2 Trt", "argk-2 Trt", "argk-2/N2 Un", "argk-2/N2 Tr")
heatmap_df <- panelB_data |>
  filter(pathway_name %in% top15_paths) |>
  mutate(contrast_short = case_when(
    contrast == "N2: Treated vs Untreated"     ~ "N2 Trt",
    contrast == "argk-2: Treated vs Untreated" ~ "argk-2 Trt",
    contrast == "argk-2 vs N2 (Untreated)"     ~ "argk-2/N2 Un",
    contrast == "argk-2 vs N2 (Treated)"       ~ "argk-2/N2 Tr"
  )) |>
  complete(pathway_name = top15_paths, contrast_short = col_labels,
           fill = list(n_deg = 0)) |>
  mutate(
    pathway_name   = factor(pathway_name,  levels = rev(top15_paths)),
    contrast_short = factor(contrast_short, levels = col_labels)
  )

panel_B <- ggplot(heatmap_df, aes(x = contrast_short, y = pathway_name, fill = n_deg)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = n_deg), size = 2.6, color = "black") +
  scale_fill_gradient(low = "#F0F4FA", high = "#1F4E79", name = "DEGs") +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1,
                                    size = 7, color = "black"),
    axis.text.y      = element_text(size = 7, face = "bold", color = "black"),
    axis.title       = element_blank(),
    panel.grid       = element_blank(),
    plot.margin      = margin(2, 4, 2, 2),
    legend.position  = "right",
    legend.key.width = unit(0.35, "cm"),
    legend.key.height = unit(0.6, "cm"),
    legend.title     = element_text(size = 7),
    legend.text      = element_text(size = 6.5)
  )

# ---- 6. Panel C: Interaction-term DEGs by pathway ----
message("Building Panel C...")

interaction_sig <- deg_interaction |>
  filter(!is.na(log2FoldChange), !is.na(padj), padj < 0.05) |>
  transmute(Gene = gene_symbol, log2FC = log2FoldChange, padj)

interaction_kegg <- interaction_sig |>
  inner_join(gene2path, by = c("Gene" = "kegg_gene"),
             relationship = "many-to-many")

interaction_kegg_filt <- interaction_kegg |>
  group_by(pathway_name) |>
  filter(n() >= 3) |>
  ungroup()

int_pathway_order <- interaction_kegg_filt |>
  group_by(pathway_name) |>
  summarise(n = n(), .groups = "drop") |>
  arrange(desc(n)) |>
  slice_head(n = 15) |>
  pull(pathway_name)

bar_data_C <- interaction_kegg_filt |>
  filter(pathway_name %in% int_pathway_order) |>
  mutate(
    direction = ifelse(log2FC > 0, "Up in argk-2-/-", "Down in argk-2-/-"),
    pathway_name = factor(pathway_name, levels = rev(int_pathway_order))
  ) |>
  group_by(pathway_name, direction) |>
  summarise(n = n_distinct(Gene), .groups = "drop")

panel_C <- ggplot(bar_data_C, aes(x = n, y = pathway_name, fill = direction)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = c(
    "Up in argk-2-/-"   = "#B2182B",
    "Down in argk-2-/-" = "#2166AC"
  )) +
  theme_classic(base_size = 8) +
  theme(
    axis.text.y      = element_text(size = 7, face = "bold", color = "black"),
    axis.text.x      = element_text(size = 7, color = "black"),
    axis.title.x     = element_text(size = 8, face = "bold", color = "black"),
    axis.title.y     = element_blank(),
    legend.position  = "bottom",
    legend.title     = element_blank(),
    legend.text      = element_text(size = 7),
    legend.key.size  = unit(0.25, "cm"),
    plot.margin      = margin(2, 4, 2, 2),
    axis.line        = element_line(linewidth = 0.3),
    axis.ticks       = element_line(linewidth = 0.2)
  ) +
  labs(x = "Interaction DEGs")

# ---- 7. Assemble + save ----
message("Assembling Figure 2...")
right_col <- panel_B / panel_C + plot_layout(heights = c(1.1, 0.9))
fig2 <- (panel_A | right_col) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 10))

out_pdf <- file.path(FIG_DIR, "Figure2_KEGG_recreated.pdf")
out_png <- file.path(FIG_DIR, "Figure2_KEGG_recreated.png")
ggsave(out_pdf, fig2, width = 7.5, height = 9, bg = "white")
ggsave(out_png, fig2, width = 7.5, height = 9, dpi = 400, bg = "white")
message("\u2714 Saved: ", out_pdf)
message("\u2714 Saved: ", out_png)
message("\n\u2714 Figure 2 complete.")
