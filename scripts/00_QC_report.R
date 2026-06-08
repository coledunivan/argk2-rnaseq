
# ---- 0) Libraries ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(pheatmap)
  library(RColorBrewer)

# Output directory for QC artifacts (so the repo root stays clean)
QC_OUT <- "outputs/QC"
if (!dir.exists(QC_OUT)) dir.create(QC_OUT, recursive = TRUE)
})

# ---- 1) FILE PATHS (EDIT IF NEEDED) ----
counts_file    <- "data/RNASEQ61125.csv"          # raw counts matrix
metadata_file  <- "data/Sample_Metadata_Table.csv"  # sample metadata
mapping_file   <- "mapping_stats.txt"               # mapping stats (tab-delimited)

# ---- 2) LOAD RAW COUNTS ----
message("Loading raw count data from: ", counts_file)

# Source shared utilities for robust CSV loading (handles CR-only line endings)
for (utils_path in c("scripts/_utils.R", "_utils.R", "../scripts/_utils.R")) {
  if (file.exists(utils_path)) { source(utils_path); break }
}
if (!exists("read_csv_robust")) {
  stop("Could not source scripts/_utils.R. Run this script from the repo root.")
}

counts_raw <- read_csv_robust(
  counts_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# If first column is unnamed or looks like an index, rename to Gene
if (colnames(counts_raw)[1] == "" || grepl("^Unnamed", colnames(counts_raw)[1])) {
  colnames(counts_raw)[1] <- "Gene"
}

# If the gene column is not called "Gene", force it:
if (!"Gene" %in% colnames(counts_raw)) {
  colnames(counts_raw)[1] <- "Gene"
}

# Remove duplicated gene rows (keep first occurrence)
counts_raw <- counts_raw %>%
  filter(!duplicated(Gene))

# Coerce sample columns to numeric (defensive against CR-mangled types)
sample_cols <- setdiff(colnames(counts_raw), "Gene")
for (cc in sample_cols) {
  counts_raw[[cc]] <- as.numeric(as.character(counts_raw[[cc]]))
}

# Build count matrix with genes as rownames
counts <- counts_raw %>%
  column_to_rownames("Gene")

message("Counts loaded: ", nrow(counts), " genes x ", ncol(counts), " samples.")

# ---- 3) LOAD METADATA ----
message("Loading metadata from: ", metadata_file)

coldata_full <- read_csv_robust(
  metadata_file,
  stringsAsFactors = FALSE
)
# Detect the sample-ID column (first one matching standard names)
id_col_qc <- intersect(c("SampleLabel", "library_name", "Sample", "sample"),
                       colnames(coldata_full))[1]
if (is.na(id_col_qc)) id_col_qc <- colnames(coldata_full)[1]
rownames(coldata_full) <- trimws(coldata_full[[id_col_qc]])
coldata <- coldata_full

# ---- 4) ALIGN COUNTS & METADATA ----
# Harmonize naming if counts use the GEO library names
counts <- harmonize_sample_names(counts, coldata,
                                 target_col = id_col_qc,
                                 alt_col    = setdiff(c("SampleLabel", "GEO_library_name"),
                                                      id_col_qc)[1])

shared_samples <- intersect(colnames(counts), rownames(coldata))

if (length(shared_samples) == 0) {
  stop("No overlapping sample IDs between counts and metadata. Check column/row names.")
}

counts  <- counts[, shared_samples]
coldata <- coldata[shared_samples, ]

# Confirm perfect alignment
stopifnot(all(colnames(counts) == rownames(coldata)))

message("Aligned samples: ", length(shared_samples), " (", paste(shared_samples, collapse = ", "), ")")

# ---- 5) CREATE DESEq2 OBJECT ----
# Uses your preferred design:
#   ~ genotype * treatment + genotype:experiment + treatment:experiment
message("Creating DESeqDataSet with design = ~ genotype * treatment + genotype:experiment + treatment:experiment")

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = coldata,
  design    = ~ genotype * treatment + genotype:experiment + treatment:experiment
)

# Optional: filter low-count genes
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]
message("Filtered to ", nrow(dds), " genes with rowSum >= 10.")

# ---- 6) RUN DESeq ----
message("Running DESeq model...")
dds <- DESeq(dds)
message("DESeq finished.")

# ---- 7) VARIANCE STABILIZING TRANSFORM (for PCA & heatmap) ----
message("Computing variance-stabilizing transform (VST)...")
vsd <- vst(dds, blind = TRUE)

# ---- 8) PCA PLOT ----
message("Generating PCA plot...")
pca_plot <- plotPCA(vsd, intgroup = c("genotype", "treatment"))

pdf(file.path(QC_OUT, "PCA_genotype_treatment.pdf"), width = 6, height = 5)
print(pca_plot + ggtitle("PCA of samples (VST)"))
dev.off()
message("Saved PCA plot to: outputs/QC/PCA_genotype_treatment.pdf")

# ---- 9) SAMPLE-TO-SAMPLE DISTANCE HEATMAP (ROBUST VERSION) ----
message("Computing sample distance heatmap...")

# 9.1 Compute distances on VST assay
vsd_mat <- assay(vsd)

# sanity check
if (anyNA(vsd_mat)) {
  warning("VST matrix contains NA values; heatmap clustering may fail.")
}

sample_dists <- dist(t(vsd_mat))
sample_dist_matrix <- as.matrix(sample_dists)

# 9.2 Use simple, consistent sample labels
sample_labels <- colnames(vsd_mat)
rownames(sample_dist_matrix) <- sample_labels
colnames(sample_dist_matrix) <- sample_labels

# 9.3 Build annotation dataframe IF genotype/treatment exist
cd_df <- as.data.frame(colData(dds))

anno_cols <- intersect(c("genotype", "treatment"), colnames(cd_df))

if (length(anno_cols) > 0) {
  annotation_df <- cd_df[, anno_cols, drop = FALSE]
  # make sure annotation rows are in the same order as columns of the matrix
  annotation_df <- annotation_df[sample_labels, , drop = FALSE]
} else {
  annotation_df <- NULL
  warning("No 'genotype' or 'treatment' columns found in colData; heatmap will be unannotated.")
}

# 9.4 Plot heatmap
pdf(file.path(QC_OUT, "Sample_distance_heatmap.pdf"), width = 7, height = 6)
pheatmap(
  sample_dist_matrix,
  annotation_col = annotation_df,
  main = "Sample-to-sample distances (VST)",
  fontsize = 9
)
dev.off()
message("Saved sample distance heatmap to: outputs/QC/Sample_distance_heatmap.pdf")

# ---- 10) DISPERSION PLOT ----
message("Generating dispersion plot...")

pdf(file.path(QC_OUT, "Dispersion_estimates.pdf"), width = 5, height = 5)
plotDispEsts(dds, main = "Dispersion estimates")
dev.off()
message("Saved dispersion plot to: outputs/QC/Dispersion_estimates.pdf")

# ---- 11) MAPPING RATE TABLE ----
# Requires a tab-delimited file "mapping_stats.txt" with:
#   sample_id  total_reads  mapped_reads
# Example:
#   A1 22000000 21000000
#   A2 23000000 22000000
# etc.
# This block will silently skip if the file isn't present.
# -----------------------------------------------
if (file.exists(mapping_file)) {
  message("Loading mapping stats from: ", mapping_file)
  
  mapping_stats <- read.table(
    mapping_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
  
  if (!all(c("sample_id", "total_reads", "mapped_reads") %in% colnames(mapping_stats))) {
    warning("mapping_stats.txt does not have the required columns: sample_id, total_reads, mapped_reads. Skipping mapping table.")
  } else {
    # Calculate mapping rate
    mapping_stats <- mapping_stats %>%
      mutate(
        mapping_rate_percent = round(mapped_reads / total_reads * 100, 2)
      )
    
    # Join to metadata for annotation
    meta_df <- as.data.frame(colData(dds))
    meta_df$sample_id <- rownames(meta_df)
    
    mapping_stats_annotated <- mapping_stats %>%
      left_join(meta_df, by = "sample_id") %>%
      relocate(sample_id, genotype, treatment)
    
    write.csv(
      mapping_stats_annotated,
      file.path(QC_OUT, "Table_SX_mapping_rates.csv"),
      row.names = FALSE
    )
    
    message("Saved mapping rate table to: outputs/QC/Table_SX_mapping_rates.csv")
  }
} else {
  warning("Mapping file '", mapping_file, "' not found. Skipping mapping rate table.")
}

message("All QC outputs generated successfully.")
