# =============================================================================
# _utils.R — shared helpers used by all R scripts in this pipeline
# =============================================================================
#
# Sourced by other scripts via:  source("scripts/_utils.R")
# Or, if running scripts from inside scripts/:  source("_utils.R")
#
# WHY THIS FILE EXISTS
#   The deposited counts file uses Classic-Mac (CR-only) line endings, which
#   leak '\r' characters into column names when read with vanilla read.csv.
#   The metadata file uses standard LF endings. The deposited counts CSV also
#   contains the 6 MAH172/argk-1 columns that were excluded from analysis,
#   while the metadata table contains only the 24 retained samples. These
#   helpers handle all of those edge cases in one place.
# =============================================================================

# ---- read_csv_robust ---------------------------------------------------------
# Read a CSV that may have CR, LF, or CRLF line endings (and/or a UTF-8 BOM
# on the first line) without leaking '\r' or '\ufeff' characters into the
# data. Wraps read.csv with line-ending and BOM normalization.
read_csv_robust <- function(path, ...) {
  txt <- readLines(path, warn = FALSE)
  if (length(txt) == 1 && grepl("\r", txt)) {
    # Classic-Mac CR-only line endings — split manually
    txt <- strsplit(txt, "\r")[[1]]
  }
  txt <- sub("\r$", "", txt)
  # Strip UTF-8 BOM (raw 3 bytes EF BB BF) from the first line if present.
  # Using raw bytes avoids locale issues with the \uFEFF codepoint in POSIX locales.
  if (length(txt) > 0) {
    raw_first <- charToRaw(txt[1])
    if (length(raw_first) >= 3 &&
        identical(as.integer(raw_first[1:3]), c(0xEF, 0xBB, 0xBF))) {
      txt[1] <- rawToChar(raw_first[-(1:3)])
    }
  }
  df <- read.csv(textConnection(txt), ...)
  # Also strip BOM/whitespace from column names defensively
  colnames(df) <- vapply(trimws(colnames(df)), function(nm) {
    r <- charToRaw(nm)
    if (length(r) >= 3 &&
        identical(as.integer(r[1:3]), c(0xEF, 0xBB, 0xBF))) {
      return(rawToChar(r[-(1:3)]))
    }
    nm
  }, character(1))
  df
}

# ---- load_counts_matrix -----------------------------------------------------
# Load a counts CSV and return an integer matrix with gene IDs as rownames.
# Aborts loudly if the file looks empty/wrong.
load_counts_matrix <- function(path, gene_col_candidates = c("gene", "gene_symbol", "Gene")) {
  counts_raw <- read_csv_robust(path, check.names = FALSE,
                                stringsAsFactors = FALSE)
  gene_col <- intersect(gene_col_candidates, colnames(counts_raw))[1]
  if (is.na(gene_col)) {
    stop("Could not find a gene-ID column in ", path,
         ". Expected one of: ", paste(gene_col_candidates, collapse = ", "),
         "\n  First few columns: ",
         paste(head(colnames(counts_raw), 5), collapse = ", "))
  }
  counts_dedup <- counts_raw[!duplicated(counts_raw[[gene_col]]), ]
  gene_ids <- counts_dedup[[gene_col]]
  count_cols <- setdiff(colnames(counts_dedup), gene_col)

  # Coerce to numeric defensively
  for (cc in count_cols) {
    counts_dedup[[cc]] <- as.numeric(as.character(counts_dedup[[cc]]))
  }

  m <- as.matrix(counts_dedup[, count_cols, drop = FALSE])
  rownames(m) <- gene_ids
  storage.mode(m) <- "integer"

  total <- sum(m, na.rm = TRUE)
  if (total == 0) {
    stop("Count matrix is all zeros after loading. Check delimiter, ",
         "encoding, and column types in ", path)
  }
  m
}

# ---- load_metadata ----------------------------------------------------------
# Load sample metadata and apply factor levels. Returns a data.frame with
# sample IDs as rownames.
load_metadata <- function(path,
                          id_candidates = c("SampleLabel", "library_name",
                                            "Sample", "sample"),
                          exclude_genotypes = c("MAH172")) {
  meta <- read_csv_robust(path, stringsAsFactors = FALSE)
  id_col <- intersect(id_candidates, colnames(meta))[1]
  if (is.na(id_col)) {
    stop("Could not find a sample-ID column in ", path,
         ". Expected one of: ", paste(id_candidates, collapse = ", "))
  }
  meta[[id_col]] <- trimws(meta[[id_col]])
  meta <- meta[!is.na(meta$genotype) &
               !is.na(meta$treatment) &
               !is.na(meta$experiment), ]
  if (length(exclude_genotypes) > 0) {
    meta <- meta[!(meta$genotype %in% exclude_genotypes), ]
  }
  meta$genotype   <- factor(meta$genotype,
                            levels = c("N2", "RB2060", "RB2598"))
  meta$treatment  <- factor(meta$treatment,
                            levels = c("untreated", "treated"))
  meta$experiment <- factor(meta$experiment)
  rownames(meta) <- meta[[id_col]]
  attr(meta, "id_col") <- id_col
  meta
}

# ---- align_counts_and_metadata ----------------------------------------------
# Drop count-matrix columns that have no metadata row, and reorder metadata
# to match the column order of the counts matrix.
align_counts_and_metadata <- function(counts, meta, verbose = TRUE) {
  shared <- intersect(colnames(counts), rownames(meta))
  if (length(shared) == 0) {
    stop("No samples shared between counts and metadata.\n",
         "  Counts columns: ", paste(head(colnames(counts), 8), collapse = ", "),
         "\n  Metadata IDs:   ", paste(head(rownames(meta), 8), collapse = ", "))
  }
  if (verbose) {
    n_drop <- ncol(counts) - length(shared)
    if (n_drop > 0) {
      dropped <- setdiff(colnames(counts), rownames(meta))
      cat("  Dropping", n_drop, "count columns not in metadata: ",
          paste(dropped, collapse = ", "), "\n")
    }
  }
  list(counts = counts[, shared, drop = FALSE],
       meta   = meta[shared, , drop = FALSE])
}

# ---- harmonize_sample_names -------------------------------------------------
# The author's pipeline uses two sample-naming conventions:
#   • A-labels (A1, A2, ...) - used in original analysis scripts, qPCR, TFLink
#   • GEO library names (RB2060_UNT_EXP1_R1, ...) - used in the GEO submission
#
# Most existing scripts assume A-labels internally because they were written
# before the GEO submission. The metadata file contains both columns
# (SampleLabel and GEO_library_name), which lets us rename whichever
# convention the counts file uses into the canonical A-label form.
#
# Returns the counts matrix with column names converted to A-labels.
harmonize_sample_names <- function(counts, meta,
                                   target_col = "SampleLabel",
                                   alt_col    = "GEO_library_name",
                                   verbose    = TRUE) {
  if (!(target_col %in% colnames(meta))) {
    if (verbose) cat("  No '", target_col,
                     "' column in metadata; skipping rename.\n", sep = "")
    return(counts)
  }
  current_cols <- colnames(counts)

  # Case 1: columns already match target_col values (e.g. A1, A2, ...)
  if (any(current_cols %in% meta[[target_col]])) {
    if (verbose) cat("  Counts columns already use ", target_col, " convention.\n",
                     sep = "")
    return(counts)
  }

  # Case 2: columns match the alt_col values — rename to target_col
  if (!is.null(alt_col) && alt_col %in% colnames(meta) &&
      any(current_cols %in% meta[[alt_col]])) {
    rename_map <- setNames(meta[[target_col]], meta[[alt_col]])
    new_cols <- rename_map[current_cols]
    n_renamed <- sum(!is.na(new_cols))
    if (verbose) cat("  Renaming ", n_renamed, " count columns from ",
                     alt_col, " -> ", target_col, " convention\n", sep = "")
    # For unmapped columns, keep the original name (will be dropped at alignment)
    new_cols[is.na(new_cols)] <- current_cols[is.na(new_cols)]
    colnames(counts) <- new_cols
    return(counts)
  }

  if (verbose) {
    cat("  Counts columns match neither ", target_col, " nor ", alt_col,
        " in metadata; passing through unchanged.\n", sep = "")
    cat("    First count columns: ", paste(head(current_cols, 3), collapse=", "), "\n")
    cat("    First ", target_col, " values: ", paste(head(meta[[target_col]], 3),
                                                     collapse = ", "), "\n", sep="")
  }
  counts
}

# ---- load_kegg_offline ------------------------------------------------------
# Drop-in replacement for clusterProfiler::download_KEGG("cel") that reads
# two pre-downloaded CSVs (produced once locally by an R session with
# Bioconductor access; see README). Returns the same shape:
#   list(KEGGPATHID2EXTID = data.frame(from=pathway, to=gene),
#        KEGGPATHID2NAME  = data.frame(from=pathway, to=name))
load_kegg_offline <- function(
  ext_file  = "data/reference/kegg_cel_pathway_to_gene.csv",
  name_file = "data/reference/kegg_cel_pathway_names.csv"
) {
  if (!file.exists(ext_file))  stop("Missing KEGG pathway-to-gene file: ", ext_file)
  if (!file.exists(name_file)) stop("Missing KEGG pathway-name file: ",   name_file)
  ext  <- read.csv(ext_file,  stringsAsFactors = FALSE)
  nm   <- read.csv(name_file, stringsAsFactors = FALSE)
  # Both ship with columns "from","to" — keep as-is to match download_KEGG output
  list(
    KEGGPATHID2EXTID = data.frame(from = ext$from, to = ext$to,
                                  stringsAsFactors = FALSE),
    KEGGPATHID2NAME  = data.frame(from = nm$from,  to = nm$to,
                                  stringsAsFactors = FALSE)
  )
}

# ---- load_go_offline --------------------------------------------------------
# Load the pre-exported gene→GO mapping for C. elegans. Returns a tibble:
#   ENTREZID, SYMBOL, GO, EVIDENCE, ONTOLOGY
load_go_offline <- function(path = "data/reference/org_Ce_eg_GO_map.csv") {
  if (!file.exists(path)) stop("Missing GO map file: ", path)
  read.csv(path, stringsAsFactors = FALSE)
}

# ---- enrichGO_offline -------------------------------------------------------
# Minimal hypergeometric GO enrichment that mimics clusterProfiler::enrichGO's
# output data frame. Uses base R's phyper(), no Bioconductor dependency.
# Returns a data frame with: ID, Description, GeneRatio, BgRatio, pvalue,
# p.adjust, qvalue, geneID, Count — matching enrichGO's @result slot columns.
#
# Args:
#   gene          : character vector of ENTREZ IDs (the "study" / DE gene set)
#   universe      : character vector of ENTREZ IDs (the background / tested set)
#   go_map        : data.frame with ENTREZID, GO, ONTOLOGY columns (from
#                   load_go_offline())
#   go_term_names : OPTIONAL data.frame with GO id → human-readable term name.
#                   If NULL and GO.db is installed, uses GO.db; otherwise GO IDs
#                   are returned as the Description.
#   ontology      : "BP" (default), "MF", or "CC"
#   pvalueCutoff  : adjusted-p cutoff for the returned set (default 0.05)
#   minGSSize     : minimum genes annotated to a term (default 10)
#   maxGSSize     : maximum genes annotated to a term (default 500)
enrichGO_offline <- function(gene, universe, go_map,
                             go_term_names = NULL,
                             ontology      = "BP",
                             pvalueCutoff  = 0.05,
                             minGSSize     = 10,
                             maxGSSize     = 500) {
  gene     <- as.character(gene);     gene     <- gene[!is.na(gene) & gene != ""]
  universe <- as.character(universe); universe <- universe[!is.na(universe) & universe != ""]
  gm <- go_map[!is.na(go_map$GO) & go_map$GO != "" &
               go_map$ONTOLOGY == ontology, , drop = FALSE]
  gm <- gm[gm$ENTREZID %in% universe, , drop = FALSE]
  gm <- gm[!duplicated(paste0(gm$ENTREZID, "|", gm$GO)), , drop = FALSE]
  # Per-GO universe and study counts
  gs_universe <- split(as.character(gm$ENTREZID), gm$GO)
  gs_universe <- gs_universe[lengths(gs_universe) >= minGSSize &
                             lengths(gs_universe) <= maxGSSize]
  if (!length(gs_universe)) return(data.frame())
  gene_in_universe <- intersect(gene, universe)
  k <- length(gene_in_universe)
  N <- length(universe)
  out <- lapply(names(gs_universe), function(go_id) {
    M     <- length(gs_universe[[go_id]])
    hits  <- intersect(gs_universe[[go_id]], gene_in_universe)
    x     <- length(hits)
    # P(X >= x) hypergeometric with x-1 (upper tail)
    p <- phyper(x - 1, M, N - M, k, lower.tail = FALSE)
    data.frame(
      ID         = go_id,
      Count      = x,
      GeneRatio  = sprintf("%d/%d", x, k),
      BgRatio    = sprintf("%d/%d", M, N),
      pvalue     = p,
      geneID     = paste(sort(hits), collapse = "/"),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  out <- out[!is.na(out$pvalue) & out$Count > 0, , drop = FALSE]
  if (!nrow(out)) return(out)
  out$p.adjust <- p.adjust(out$pvalue, method = "BH")
  out$qvalue   <- out$p.adjust   # clusterProfiler reports qvalue; BH used as proxy
  # Term name lookup
  if (is.null(go_term_names) && requireNamespace("GO.db", quietly = TRUE)) {
    term_lookup <- AnnotationDbi::select(GO.db::GO.db,
                                         keys = out$ID,
                                         columns = c("GOID", "TERM"),
                                         keytype = "GOID")
    out$Description <- term_lookup$TERM[match(out$ID, term_lookup$GOID)]
  } else if (!is.null(go_term_names)) {
    out$Description <- go_term_names[[2]][match(out$ID, go_term_names[[1]])]
  } else {
    out$Description <- out$ID
  }
  out <- out[!is.na(out$Description), , drop = FALSE]
  out <- out[order(out$pvalue), , drop = FALSE]
  # Reorder to match clusterProfiler column order
  out <- out[, c("ID","Description","GeneRatio","BgRatio","pvalue",
                 "p.adjust","qvalue","geneID","Count")]
  rownames(out) <- NULL
  out
}
