#!/usr/bin/env Rscript
#
# rebuild_grn_inputs.R
#
# PURPOSE
#   Rebuild nodes.tsv and edges_combined.tsv so the network includes ALL
#   90 G x E genes from DESeq2_effect_classification.csv, using only
#   TFLink regulatory edges (no STRING PPI required).
#
# WHY NO STRING
#   The Python figure script (grn2_publishable_finalish.py) calls
#   regulatory_edges_only() before building its network, so PPI edges
#   have no effect on the final figure. Removing the STRING dependency
#   eliminates the need to download multi-GB protein-interaction files.
#
# INPUTS
#   EFFECT_CLASS_FILE   DESeq2_effect_classification.csv
#   TFLINK_FILE         Your TFLink export — CSV or TSV, auto-detected
#
# OUTPUTS
#   nodes.tsv                  TFs + signaling + all 90 G×E effectors
#   edges_combined.tsv         Regulatory edges only (edge_type=regulatory)
#   TFLink_coverage_report.tsv Which G×E genes were/weren't covered
#
# DEPENDENCIES
#   install.packages(c("dplyr", "readr", "stringr"))

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})


# =========================================================================
# CONFIG
# =========================================================================
EFFECT_CLASS_FILE <- {
  candidates <- c("DESeq2_effect_classification.csv",
                  "outputs/DESeq2/DESeq2_effect_classification.csv")
  found <- candidates[file.exists(candidates)]
  if (length(found) == 0) {
    stop("Cannot find DESeq2_effect_classification.csv. ",
         "Run 01_DESeq2_and_TF_enrichment.R first.")
  }
  found[1]
}
TFLINK_FILE       <- "data/reference/RNAseq_to_TF_Targets.csv"   # CSV or TSV

# The 15 TFs that are members of MODULES in grn2_publishable_finalish.py.
# Must match the script's MODULES dict exactly (case-insensitive).
CHARACTERISTIC_TFS <- c(
  "pqm-1", "ceh-60", "fos-1", "elt-2",    # M0  developmental / metabolic
  "zip-2",  "skn-1",  "cebp-1",            # M1  immune surveillance
  "nhr-28", "nhr-77", "blmp-1",            # M2  secondary axis
  "daf-16", "pha-4",  "akt-2",             # M3  IIS / stress
  "aak-2",  "hif-1"                        # M3  (signaling anchors)
)

# Extra signaling proteins kept as nodes for completeness
SIGNALING_PROTEINS <- c(
  "aak-2", "pmk-1", "atf-7", "akt-1", "akt-2", "hif-1", "par-4"
)

{
  grn_out <- "outputs/GRN_tables"
  if (!dir.exists(grn_out)) dir.create(grn_out, recursive = TRUE)
  OUT_NODES    <- file.path(grn_out, "nodes_full90.tsv")
  OUT_EDGES    <- file.path(grn_out, "edges_combined_full90.tsv")
  OUT_COVERAGE <- file.path(grn_out, "TFLink_coverage_report.tsv")
}


# =========================================================================
# HELPERS
# =========================================================================

# Read CSV or TSV by file extension
read_auto <- function(path, ...) {
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    read_csv(path, show_col_types = FALSE, ...)
  } else {
    read_tsv(path, show_col_types = FALSE, ...)
  }
}

# Detect which columns are the TF and target gene names.
# Handles several TFLink release naming conventions, plus common
# custom export formats.
detect_tf_target_cols <- function(df) {
  tf_cands <- c(
    "Name.TF", "NameTF", "TF.Name", "TF_Name", "TF", "Name_TF",
    "TFName", "tf", "Transcription.Factor", "TF_symbol", "TF_name"
  )
  tgt_cands <- c(
    "Name.Target", "NameTarget", "Target.Name", "Target_Name",
    "Target", "Name_Target", "TargetName", "target",
    "Gene.Name", "gene_name", "Gene", "gene", "target_gene",
    "Target_symbol", "Target_name"
  )
  tf_col  <- intersect(tf_cands,  names(df))[1]
  tgt_col <- intersect(tgt_cands, names(df))[1]

  if (is.na(tf_col) || is.na(tgt_col)) {
    cat("\nAll column names in your TFLink file:\n")
    cat(paste(" ", names(df), collapse = "\n"), "\n\n")
    stop(
      "Cannot auto-detect TF / Target columns.\n",
      "Add the correct column names to tf_cands / tgt_cands in ",
      "detect_tf_target_cols() and re-run."
    )
  }
  cat(sprintf("          TF column  = '%s'\n", tf_col))
  cat(sprintf("          Target col = '%s'\n", tgt_col))
  list(tf = tf_col, target = tgt_col)
}

normalize_gene <- function(x) {
  x |> as.character() |> str_trim() |> tolower()
}


# =========================================================================
# STEP 1  Load DESeq2 classification
# =========================================================================
cat("[step 1] loading", EFFECT_CLASS_FILE, "\n")

cls <- read_auto(EFFECT_CLASS_FILE) |>
  mutate(gene_symbol_lc = normalize_gene(gene_symbol))

stopifnot("effect_class" %in% names(cls))

gxe_genes <- cls |> filter(effect_class == "GxE") |> pull(gene_symbol_lc)
cat(sprintf("         %d G×E genes found\n", length(gxe_genes)))


# =========================================================================
# STEP 1b  Convert sequence names → common names for TFLink matching
# =========================================================================
# DESeq2 row IDs are WormBase SEQUENCE NAMES (e.g. k02f2.6, y58a7a.3).
# TFLink stores targets by their PUBLIC / COMMON NAME (e.g. irg-2, gst-19).
# Without this conversion, every filter(target %in% gxe_genes) returns 0.
#
# We download the WormBase gene-ID table (~2 MB) which maps:
#   sequence_name  →  common_name
# and build a lookup vector. Genes with no common name keep their sequence
# name (they'll only appear in TFLink if TFLink also uses the sequence name).
#
# If the download fails (no internet, or URL changed), the script falls back
# to sequence names and prints a warning — you'll get 0 edges again and
# will need to run this step manually.

cat("\n[step 1b] converting WormBase sequence names → common names\n")

convert_sequence_to_common <- function(sequence_names) {
  # Build the map from the LOCAL TFLink file (RNAseq_to_TF_Targets.csv),
  # which already contains both the sequence name (Gene column) and the
  # common name (Name.Target column) for every target. This avoids the
  # unreachable WormBase URL the original version tried to hit, which
  # silently fell back to using sequence names and produced 0 edges.
  cat("  Building seq→common map from local TFLink file (offline)\n")
  tfraw <- tryCatch(read_auto(TFLINK_FILE),
                    error = function(e) NULL)
  if (is.null(tfraw) || !all(c("Gene","Name.Target") %in% names(tfraw))) {
    cat("  [warn] Local TFLink file missing Gene/Name.Target columns;",
        "falling back to sequence names (will produce 0 edges).\n")
    return(tibble(seq_name = sequence_names, tflink_name = sequence_names))
  }
  lookup <- tfraw |>
    transmute(seq_name    = tolower(trimws(as.character(Gene))),
              tflink_name = tolower(trimws(as.character(Name.Target)))) |>
    filter(!is.na(seq_name),    seq_name    != "",
           !is.na(tflink_name), tflink_name != "", tflink_name != "-") |>
    distinct(seq_name, .keep_all = TRUE)

  result <- tibble(seq_name = sequence_names) |>
    left_join(lookup, by = "seq_name") |>
    mutate(tflink_name = coalesce(tflink_name, seq_name))

  n_converted <- sum(result$tflink_name != result$seq_name, na.rm = TRUE)
  cat(sprintf("         %d / %d sequence names mapped to common names\n",
              n_converted, length(sequence_names)))
  result
}

name_map      <- convert_sequence_to_common(gxe_genes)
# gxe_names_for_tflink: the names to use when filtering TFLink targets
gxe_tflink    <- name_map$tflink_name
# Keep the original sequence names for the nodes table (DESeq2 reference)
# but note the common-name equivalent for each gene
cat(sprintf("         %d unique TFLink-format names for the 90 G×E genes\n",
            length(unique(gxe_tflink))))


# =========================================================================
# STEP 2  Load TFLink and build regulatory edges
# =========================================================================
cat("\n[step 2] loading TFLink:", TFLINK_FILE, "\n")

tflink_raw <- read_auto(TFLINK_FILE)
cols       <- detect_tf_target_cols(tflink_raw)

tflink <- tflink_raw |>
  transmute(
    tf     = normalize_gene(.data[[cols$tf]]),
    target = normalize_gene(.data[[cols$target]])
  ) |>
  filter(!is.na(tf), tf != "", !is.na(target), target != "") |>
  distinct()

cat(sprintf("         %d unique TF->target pairs loaded\n", nrow(tflink)))

# 2a  Characteristic TFs → G×E genes
#     Filter using the TFLink-format names (common names), then join back
#     the original sequence name so node IDs stay consistent with DESeq2.
gxe_reg_edges_raw <- tflink |>
  filter(tf %in% CHARACTERISTIC_TFS, target %in% gxe_tflink)
cat(sprintf("[2a]     %d edges: characteristic TFs → G×E targets\n",
            nrow(gxe_reg_edges_raw)))

# Map common names back to sequence names for node-table consistency
common_to_seq <- setNames(name_map$seq_name, name_map$tflink_name)
gxe_reg_edges <- gxe_reg_edges_raw |>
  mutate(target = coalesce(common_to_seq[target], target))

# 2b  Characteristic TFs → other TFs / signaling proteins (for TF-TF connections)
tf_tf_edges <- tflink |>
  filter(
    tf     %in% CHARACTERISTIC_TFS,
    target %in% c(CHARACTERISTIC_TFS, SIGNALING_PROTEINS),
    tf     != target
  )
cat(sprintf("[2b]     %d edges: TF → TF / signaling\n", nrow(tf_tf_edges)))

all_reg_edges <- bind_rows(gxe_reg_edges, tf_tf_edges) |> distinct()


# =========================================================================
# STEP 3  Coverage diagnostic
# =========================================================================
cat("\n[step 3] coverage report\n")

covered <- unique(gxe_reg_edges$target)   # these are now sequence names
orphans <- setdiff(gxe_genes, covered)

cat(sprintf("         %d / %d G×E genes have ≥1 edge to the 15 TFs\n",
            length(covered), length(gxe_genes)))
cat(sprintf("         %d orphans (no edge to any of the 15 TFs)\n",
            length(orphans)))

# For orphans: check which additional TFs in TFLink cover the most
# (search by TFLink common-name format)
if (length(orphans) > 0) {
  orphan_tflink_names <- name_map |>
    filter(seq_name %in% orphans) |>
    pull(tflink_name)

  extra_tfs <- tflink |>
    filter(target %in% orphan_tflink_names, !(tf %in% CHARACTERISTIC_TFS)) |>
    count(tf, name = "n_orphans_covered") |>
    arrange(desc(n_orphans_covered)) |>
    mutate(pct_orphans = paste0(round(100 * n_orphans_covered / length(orphans)), "%"))

  if (nrow(extra_tfs) > 0) {
    cat("\n  Top TFs (outside your 15) that regulate the most orphans:\n")
    print(head(extra_tfs, 12), n = 12)
    cat("  → Add the best ones to CHARACTERISTIC_TFS and re-run.\n")
  }
}

coverage_report <- tibble(
    gene_symbol     = gxe_genes,
    covered_by_15   = gxe_genes %in% covered,
    is_orphan       = gxe_genes %in% orphans
  )

write_tsv(coverage_report, OUT_COVERAGE)


# =========================================================================
# STEP 4  Build nodes table
# =========================================================================
cat("\n[step 4] assembling nodes\n")

cls_lookup <- cls |>
  select(gene_symbol_lc,
         lfc_gxe, padj_gxe,
         lfc_treat_N2, padj_treat_N2,
         lfc_geno_unt, padj_geno_unt,
         lfc_geno_trt, padj_geno_trt,
         effect_class)

make_node_row <- function(gene, cls_class) {
  d <- filter(cls_lookup, gene_symbol_lc == gene)
  # display_label: use common name if available, else sequence name
  dlabel <- coalesce(
    name_map$tflink_name[name_map$seq_name == gene][1],
    gene
  )
  if (nrow(d) == 0) {
    tibble(
      node_id         = gene,        node_class      = cls_class,
      display_label   = dlabel,
      effect_category = NA_character_,
      lfc_gxe = NA_real_,  padj_gxe      = NA_real_,
      lfc_treat_N2 = NA_real_,  padj_treat_N2 = NA_real_,
      lfc_geno_unt = NA_real_,  padj_geno_unt = NA_real_,
      lfc_geno_trt = NA_real_,  padj_geno_trt = NA_real_,
      orphan_gxe = FALSE
    )
  } else {
    tibble(
      node_id         = gene,        node_class      = cls_class,
      display_label   = dlabel,
      effect_category = d$effect_class[1],
      lfc_gxe = d$lfc_gxe[1],      padj_gxe      = d$padj_gxe[1],
      lfc_treat_N2 = d$lfc_treat_N2[1], padj_treat_N2 = d$padj_treat_N2[1],
      lfc_geno_unt = d$lfc_geno_unt[1], padj_geno_unt = d$padj_geno_unt[1],
      lfc_geno_trt = d$lfc_geno_trt[1], padj_geno_trt = d$padj_geno_trt[1],
      orphan_gxe = (cls_class == "effector") && (gene %in% orphans)
    )
  }
}

nodes <- bind_rows(
  lapply(CHARACTERISTIC_TFS, make_node_row, cls_class = "TF"),
  lapply(SIGNALING_PROTEINS, make_node_row, cls_class = "signaling"),
  lapply(gxe_genes,          make_node_row, cls_class = "effector")
) |>
  # aak-2 and hif-1 appear in both lists; keep the TF row
  arrange(match(node_class, c("TF", "signaling", "effector"))) |>
  distinct(node_id, .keep_all = TRUE)

# TF sizing: number of G×E genes each TF regulates
tf_degree <- gxe_reg_edges |>
  count(tf, name = "degree_GxE_targets")

nodes <- nodes |>
  left_join(tf_degree, by = c("node_id" = "tf")) |>
  mutate(degree_GxE_targets = coalesce(degree_GxE_targets, 0L))

cat(sprintf("         %d nodes  (%d TF · %d signaling · %d effectors · %d orphan effectors)\n",
            nrow(nodes),
            sum(nodes$node_class == "TF"),
            sum(nodes$node_class == "signaling"),
            sum(nodes$node_class == "effector"),
            sum(nodes$orphan_gxe, na.rm = TRUE)))


# =========================================================================
# STEP 5  Build edges table
# =========================================================================
cat("[step 5] assembling edges\n")

edges <- all_reg_edges |>
  transmute(
    source                = tf,
    target                = target,
    edge_type             = "regulatory",
    string_combined_score = NA_real_   # column kept for schema compatibility
  ) |>
  distinct()

cat(sprintf("         %d regulatory edges\n", nrow(edges)))


# =========================================================================
# STEP 6  Write
# =========================================================================
write_tsv(nodes, OUT_NODES)
write_tsv(edges, OUT_EDGES)

cat(sprintf("\n[done] %s  — %d rows\n", OUT_NODES,    nrow(nodes)))
cat(sprintf("[done] %s  — %d rows\n", OUT_EDGES,    nrow(edges)))
cat(sprintf("[done] %s — %d rows\n\n", OUT_COVERAGE, nrow(coverage_report)))

cat("Next steps:\n")
cat("  1. If orphan count is high, add the suggested TFs to CHARACTERISTIC_TFS\n")
cat("     (they must also be added to the MODULES dict in the Python script).\n")
cat("  2. Run:  python3 grn2_publishable_finalish.py\n")
cat("     The script reads nodes.tsv + edges_combined.tsv from the same directory.\n")
