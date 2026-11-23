## =============================================================================
## PID430 - heatmaps for curated gene lists (manuscript SI Fig. 7 style)
##
## Companion to PID430_proteomics_pipeline.R (v5). Runs AFTER it, reads only its
## outputs, and computes no new statistics. One heatmap per curated list: rows =
## the listed proteins, columns = the 12 samples grouped by condition, colour =
## z-score of normalized log2 abundance per protein (the same matrix and the same
## row scaling as heatmap_significant_union.png).
##
## Inputs
##   run folder written by the main script:
##     02_processed_data/normalized_log2_abundance.csv
##     01_QC/sample_sheet.csv
##     03_differential_abundance/per_contrast_tables/*_results.csv  (symbols, adj p)
##   one or more CSV files, each a curated list. Symbols are read from a column
##   named Gene / Gene_Symbol / Symbol (any case); otherwise the first column.
##   The file name (without .csv) becomes the heatmap title.
## Outputs (in <run folder>/08_curated_heatmaps/)
##   <list>_heatmap.png          FIG_DPI
##   <list>_heatmap_values.csv   z-scores plotted + BH adj p and log2FC per contrast
##   curated_heatmap_index.csv   per list: symbols given / found / not found (named)
##   sessionInfo.txt
##
## Usage:  Rscript PID430_curated_heatmaps.R <run_folder> list1.csv [list2.csv ...]
## =============================================================================

## ---- 1. USER INPUT (used when no command-line arguments are given) ---------
RUN_DIR    <- "PID430_proteomics_results_20260923_182447"
LIST_FILES <- c("actin_actomyosin.csv", "adhesion.csv", "adhesion_cytoskeleton_combined.csv", "cell_stress_death.csv", "chromatin_modifying.csv", "ecm.csv", "hypoxia.csv", "pluripotency.csv", "tubulin_microtubule.csv")

## ---- 2. SETTINGS -----------------------------------------------------------
SHOW_MEANS   <- FALSE  # FALSE: 12 replicate columns; TRUE: 4 group-mean columns
CLUSTER_ROWS <- FALSE  # FALSE keeps the order of the curated list
Z_LIMIT      <- 2.5    # colour scale fixed at +/- this z for all lists, so panels compare
FIG_DPI      <- 300    # main script FIG_DPI
GRP_DISPLAY  <- c("0-3 A", "0-3 N", "1-5 A", "1-5 N")              # main script
GRP_COLS     <- c("#4E79A7", "#F28E2B", "#59A14F", "#E15759")      # main script

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) { RUN_DIR <- args[1]; LIST_FILES <- args[-1] }

for (p in c("pheatmap")) if (!requireNamespace(p, quietly = TRUE))
  stop("Package '", p, "' is not installed.", call. = FALSE)
library(pheatmap)

f_norm  <- file.path(RUN_DIR, "02_processed_data", "normalized_log2_abundance.csv")
f_samp  <- file.path(RUN_DIR, "01_QC", "sample_sheet.csv")
tab_dir <- file.path(RUN_DIR, "03_differential_abundance", "per_contrast_tables")
for (f in c(f_norm, f_samp, tab_dir)) if (!file.exists(f))
  stop("Not found (is this a pipeline run folder?): ", f, call. = FALSE)
for (f in LIST_FILES) if (!file.exists(f)) stop("List file not found: ", f, call. = FALSE)
out_dir <- file.path(RUN_DIR, "08_curated_heatmaps")
dir.create(out_dir, showWarnings = FALSE)

## ---- 3. Load the pipeline outputs ------------------------------------------
nm  <- read.csv(f_norm, check.names = FALSE)
mat <- as.matrix(nm[, -1]); rownames(mat) <- nm$Accession
## Only proteins quantified in all 12 samples, nothing imputed. The pipeline
## (REQUIRE_COMPLETE = TRUE) already guarantees this; stop rather than plot if
## the matrix ever contains a missing value or is not 12 samples wide.
if (ncol(mat) != 12 || anyNA(mat))
  stop("normalized_log2_abundance.csv must be complete (12 samples, no NA); found ",
       ncol(mat), " columns and ", sum(is.na(mat)), " missing values.", call. = FALSE)
ss  <- read.csv(f_samp, check.names = FALSE)
stopifnot(identical(colnames(mat), ss$sample_internal))     # same column order as the run
ss$group_display <- factor(ss$group_display, levels = GRP_DISPLAY)
ord <- order(ss$group_display); mat <- mat[, ord]; ss <- ss[ord, ]

## Per-contrast tables: symbol lookup for all tested proteins, plus adj p / log2FC
tabs <- lapply(sort(list.files(tab_dir, "_results\\.csv$", full.names = TRUE)),
               read.csv, check.names = FALSE)
names(tabs) <- sub("_results\\.csv$", "", sort(list.files(tab_dir, "_results\\.csv$")))
sym <- setNames(tabs[[1]]$Gene_Symbol, tabs[[1]]$Accession)
stopifnot(all(rownames(mat) %in% names(sym)))   # tables and matrix come from one run

if (SHOW_MEANS) {
  mat <- sapply(GRP_DISPLAY, function(g)
    rowMeans(mat[, ss$group_display == g, drop = FALSE], na.rm = TRUE))
  ann <- data.frame(Group = factor(GRP_DISPLAY, levels = GRP_DISPLAY), row.names = GRP_DISPLAY)
  gaps <- NULL
} else {
  ann <- data.frame(Group = ss$group_display, row.names = colnames(mat))
  gaps <- cumsum(table(ss$group_display))[-length(GRP_DISPLAY)]
}
ann_col <- list(Group = setNames(GRP_COLS, GRP_DISPLAY))
brk <- seq(-Z_LIMIT, Z_LIMIT, length.out = 101)
pal <- colorRampPalette(c("#4E79A7", "white", "#E15759"))(100)   # main heatmap colours

## Lists may have no header (one symbol per line) or a header row; both work.
## A first line is treated as a header only if it is a column name such as
## Gene / Gene_Symbol / Symbol, so no gene is ever lost as a "header".
read_list <- function(f) {
  x <- read.csv(f, header = FALSE, check.names = FALSE, colClasses = "character",
                blank.lines.skip = TRUE, strip.white = TRUE, fileEncoding = "UTF-8-BOM")
  hdr <- "^(gene|genes|gene_symbol|gene symbol|gene_name|symbol|symbols)$"
  col <- 1
  if (nrow(x) && any(grepl(hdr, unlist(x[1, ]), ignore.case = TRUE))) {
    col <- grep(hdr, unlist(x[1, ]), ignore.case = TRUE)[1]; x <- x[-1, , drop = FALSE]
  }
  v <- trimws(x[[col]])
  unique(v[!is.na(v) & v != ""])   # blank lines and repeated symbols dropped, order kept
}

## ---- 4. One heatmap per list -----------------------------------------------
index <- list()
for (f in LIST_FILES) {
  lname <- tools::file_path_sans_ext(basename(f))
  given <- read_list(f)
  ## Case-insensitive match; keep the curated order. A symbol measured by more
  ## than one accession gives one row per accession, labelled with the accession.
  hits <- lapply(given, function(g) rownames(mat)[toupper(sym[rownames(mat)]) %in% toupper(g)])
  found <- given[lengths(hits) > 0]; missing <- given[lengths(hits) == 0]
  acc <- unlist(hits)
  cat(sprintf("%-40s %3d symbols: %3d found (%d rows), %3d not among tested proteins\n",
              lname, length(given), length(found), length(acc), length(missing)))
  index[[lname]] <- data.frame(list = lname, file = f, n_given = length(given),
                               n_found = length(found), n_rows = length(acc),
                               n_not_found = length(missing),
                               not_found = paste(missing, collapse = ";"))
  if (length(acc) < 2) { cat("   skipped: fewer than 2 proteins found\n"); next }

  m  <- mat[acc, , drop = FALSE]
  z  <- t(scale(t(m)))                         # per-protein z-score across the 12 samples
  lab <- sym[acc]; dup <- lab %in% lab[duplicated(lab)]
  lab[dup] <- paste0(lab[dup], " (", acc[dup], ")")
  rownames(z) <- lab
  clus <- CLUSTER_ROWS && !anyNA(z)
  if (CLUSTER_ROWS && !clus) cat("   rows not clustered: missing values present\n")

  png(file.path(out_dir, paste0(lname, "_heatmap.png")),
      width = if (SHOW_MEANS) 5 else 6.5, height = max(3, 1.6 + 0.16 * nrow(z)),
      units = "in", res = FIG_DPI)
  pheatmap(z, breaks = brk, color = pal, na_col = "grey85",
           cluster_rows = clus, cluster_cols = FALSE, gaps_col = gaps,
           annotation_col = ann, annotation_colors = ann_col,
           show_colnames = !SHOW_MEANS, border_color = "grey95",
           fontsize = 8, fontsize_row = 7, cellwidth = if (SHOW_MEANS) 22 else 14,
           main = sprintf("%s (n = %d)\nrow z-score, normalized log2 abundance",
                          gsub("_", " ", lname), nrow(z)))
  dev.off()

  ## Values behind the figure, with the pipeline's own statistics per contrast.
  stats <- do.call(cbind, lapply(names(tabs), function(cn) {
    t <- tabs[[cn]][match(acc, tabs[[cn]]$Accession), ]
    setNames(data.frame(t$log2FC, t$adj_p_BH, t$significant),
             paste0(cn, c("__log2FC", "__adj_p_BH", "__significant")))
  }))
  write.csv(data.frame(Gene_Symbol = sym[acc], Accession = acc,
                       round(z, 3), stats, check.names = FALSE),
            file.path(out_dir, paste0(lname, "_heatmap_values.csv")), row.names = FALSE)
}

write.csv(do.call(rbind, index), file.path(out_dir, "curated_heatmap_index.csv"), row.names = FALSE)
writeLines(c(paste("run folder:", normalizePath(RUN_DIR)),
             paste("lists:", paste(normalizePath(LIST_FILES), collapse = "; ")),
             sprintf("SHOW_MEANS=%s CLUSTER_ROWS=%s Z_LIMIT=%s", SHOW_MEANS, CLUSTER_ROWS, Z_LIMIT),
             "", capture.output(sessionInfo())),
           file.path(out_dir, "sessionInfo.txt"))
cat("\nWritten to", out_dir, "\n")
