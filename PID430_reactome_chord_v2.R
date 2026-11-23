## =============================================================================
## PID430 - Reactome chord diagrams (manuscript Fig. 2 j-o)
##
## Companion to PID430_proteomics_pipeline.R (v5). Runs AFTER it, reads only its
## outputs, and computes no statistics. It draws, for each contrast and each
## direction subset, the Reactome pathways the pipeline already found enriched
## (BH adjusted p < 0.05) linked to the significant proteins that drive them.
##
## What each figure shows
##   - left half : the top TOP_N pathways by adjusted p (P1 = most significant)
##   - right half: every significant protein belonging to those pathways,
##                 coloured by log2FC (positive = higher in the contrast's
##                 numerator, as in the per-contrast tables)
##   - ribbons   : pathway membership
##   - legend    : full pathway name, member count and adjusted p for each P
## Significance is printed in the legend because a chord diagram cannot encode
## it; the dot plots and Reactome.csv from the main script remain the
## quantitative record.
##
## Inputs (all inside the run folder written by the main script)
##   03_differential_abundance/per_contrast_tables/<contrast>_results.csv
##   05_ontology/<contrast>/<subset>/Reactome.csv
## Outputs
##   05_ontology/<contrast>/<subset>/Reactome_chord.png   (FIG_DPI)
##   05_ontology/Reactome_chord_index.csv                 (what each figure shows)
##
## Usage:  Rscript PID430_reactome_chord.R  [run_folder]
## =============================================================================

## ---- 1. USER INPUT ----------------------------------------------------------
RUN_DIR <- "PID430_proteomics_results_YYYYMMDD_HHMMSS"   # the main script's output folder

## ---- 2. FIXED SETTINGS (keep in step with the main script) ------------------
ENRICH_PADJ <- 0.05    # main script ENRICH_PADJ; re-applied here as a guard
TOP_N       <- 10      # pathways per chord; more makes gene labels unreadable
FIG_DPI     <- 300     # main script FIG_DPI
INCLUDE_ALL <- FALSE   # TRUE also draws the combined 'all_significant' subset

args <- commandArgs(trailingOnly = TRUE)
if (length(args)) RUN_DIR <- args[1]

suppressPackageStartupMessages({
  library(dplyr); library(circlize); library(org.Mm.eg.db)
})
onto_dir <- file.path(RUN_DIR, "05_ontology")
tab_dir  <- file.path(RUN_DIR, "03_differential_abundance", "per_contrast_tables")
if (!dir.exists(onto_dir) || !dir.exists(tab_dir))
  stop("Not a pipeline run folder: ", normalizePath(RUN_DIR, mustWork = FALSE), call. = FALSE)

## Symbol -> Entrez exactly as the main script's map_entrez(): symbols ranked by
## adjusted p, first (most significant) occurrence kept, one Entrez per symbol.
## Each Entrez id in Reactome.csv therefore maps back to the same protein row,
## and hence the same log2FC, that the enrichment used.
map_entrez <- function(symbols, rank_by) {
  keep <- !is.na(symbols) & symbols != "" & symbols != "NA"
  symbols <- symbols[keep][order(rank_by[keep], na.last = TRUE)]
  suppressMessages(AnnotationDbi::select(org.Mm.eg.db, keys = unique(symbols),
                                         keytype = "SYMBOL", columns = "ENTREZID")) %>%
    filter(!is.na(ENTREZID)) %>% distinct(SYMBOL, .keep_all = TRUE)
}

## RNA orientation check, done before anything is drawn. The RNA pairwise
## contrasts put the alphabetically earlier group first (0-3 A < 0-3 N < 1-5 A <
## 1-5 N), so every numerator must sort before its denominator. A run made before
## the diagonal was re-oriented (1-5_A_vs_0-3_N) stops here instead of drawing a
## chord whose log2FC colours and colour bar are the reverse of the RNA panels.
for (f in list.files(tab_dir, "_results\\.csv$", full.names = TRUE)) {
  h <- read.csv(f, check.names = FALSE, nrows = 1)
  if (!(h$numerator[1] < h$denominator[1]))
    stop(basename(f), ": numerator '", h$numerator[1], "' / denominator '", h$denominator[1],
         "' is the reverse of the RNA orientation. This run predates the ",
         "re-oriented pipeline; re-run PID430_proteomics_pipeline.R first.", call. = FALSE)
}

fc_col <- colorRamp2(c(-3, 0, 3), c("#4E79A7", "white", "#E15759"))  # heatmap colours
index  <- list()

for (cn in sort(list.dirs(onto_dir, full.names = FALSE, recursive = FALSE))) {
  tab_file <- file.path(tab_dir, paste0(cn, "_results.csv"))
  if (!file.exists(tab_file)) next                      # e.g. 00_gene_lists
  tab <- read.csv(tab_file, check.names = FALSE)
  sig <- tab[tab$significant %in% TRUE, ]

  for (sn in list.dirs(file.path(onto_dir, cn), full.names = FALSE, recursive = FALSE)) {
    if (sn == "all_significant" && !INCLUDE_ALL) next
    rfile <- file.path(onto_dir, cn, sn, "Reactome.csv")
    if (!file.exists(rfile)) {
      cat(sprintf("   %-16s %-28s no Reactome.csv (no enriched pathways or subset not tested)\n", cn, sn))
      next
    }
    ## Same protein subset as the enrichment call: all / higher in numerator / denominator
    num <- gsub(" ", "", tab$numerator[1]); den <- gsub(" ", "", tab$denominator[1])
    s <- if (sn == paste0("higher_in_", num)) sig[sig$log2FC > 0, ] else
         if (sn == paste0("higher_in_", den)) sig[sig$log2FC < 0, ] else sig
    em <- map_entrez(s$Gene_Symbol, s$adj_p_BH)
    fc <- setNames(s$log2FC[match(em$SYMBOL, s$Gene_Symbol)], em$ENTREZID)   # first row = min adj p
    sy <- setNames(em$SYMBOL, em$ENTREZID)

    re <- read.csv(rfile, check.names = FALSE)
    re <- re[re$p.adjust < ENRICH_PADJ, ]
    re <- re[order(re$p.adjust), ]
    n_all <- nrow(re); re <- head(re, TOP_N); k <- nrow(re)
    re$P <- paste0("P", seq_len(k))

    links <- do.call(rbind, lapply(seq_len(k), function(i)
      data.frame(from = re$P[i], entrez = strsplit(re$geneID[i], "/")[[1]])))
    miss <- setdiff(links$entrez, names(sy))
    if (length(miss))
      stop(sprintf("%s/%s: %d Reactome gene id(s) not among the significant proteins (%s). ",
                   cn, sn, length(miss), paste(head(miss, 5), collapse = ",")),
           "Reactome.csv and the results table are not from the same run.", call. = FALSE)
    links$to     <- sy[links$entrez]
    links$log2FC <- fc[links$entrez]
    links$value  <- 1

    ## Genes grouped by their first (most significant) pathway, then by log2FC,
    ## which keeps ribbons from one pathway together and limits crossings.
    g <- links %>% group_by(to) %>%
      summarise(first = min(as.integer(sub("P", "", from))), fc = log2FC[1], .groups = "drop") %>%
      arrange(first, desc(fc))
    genes <- g$to
    pcols <- setNames(hcl.colors(max(k, 2), "Dark 3")[seq_len(k)], re$P)
    gcols <- setNames(fc_col(g$fc), genes)

    out <- file.path(onto_dir, cn, sn, "Reactome_chord.png")
    png(out, width = 11, height = 7, units = "in", res = FIG_DPI)
    layout(matrix(1:2, 1), widths = c(6.2, 4.8))
    par(mar = c(0.5, 0.5, 3.2, 0.5), xpd = NA)
    circos.clear()
    circos.par(start.degree = 90, clock.wise = TRUE,
               gap.after = c(rep(2, k - 1), 12, rep(0.6, length(genes) - 1), 12))
    chordDiagram(links[, c("from", "to", "value")], order = c(re$P, genes),
                 grid.col = c(pcols, gcols), transparency = 0.35,
                 annotationTrack = "grid", annotationTrackHeight = 0.05,
                 preAllocateTracks = list(track.height = 0.12))
    circos.track(track.index = 1, bg.border = NA, panel.fun = function(x, y) {
      nm <- get.cell.meta.data("sector.index")
      circos.text(CELL_META$xcenter, CELL_META$ylim[1], nm, facing = "clockwise",
                  niceFacing = TRUE, adj = c(0, 0.5),
                  cex = if (nm %in% re$P) 0.9 else 0.55,
                  font = if (nm %in% re$P) 2 else 1)
    })
    circos.clear()
    title(main = sprintf("Reactome: %s\n%s", gsub("_", " ", cn), gsub("_", " ", sn)),
          cex.main = 1, line = 1)

    ## Legend: pathway names with member count and adjusted p, then the log2FC scale.
    par(mar = c(0.5, 0, 3.2, 0.5)); plot.new()
    lab <- sprintf("%s  %s\n      n = %d, adj p = %s", re$P,
                   vapply(re$Description, function(z) paste(strwrap(z, 50), collapse = "\n      "), ""),
                   lengths(strsplit(re$geneID, "/")), formatC(re$p.adjust, format = "e", digits = 1))
    ## Drawn by hand: legend() spaces multi-line entries unevenly and they overlap.
    text(0, 0.98, sprintf("Top %d of %d pathways, BH adj p < %.2g", k, n_all, ENRICH_PADJ),
         adj = c(0, 1), font = 2, cex = 0.7)
    lh <- 0.026; y <- 0.93                                  # line height in plot units
    for (i in seq_len(k)) {
      nl <- lengths(regmatches(lab[i], gregexpr("\n", lab[i]))) + 1
      rect(0, y - 0.018, 0.03, y, col = pcols[i], border = NA)
      text(0.045, y, lab[i], adj = c(0, 1), cex = 0.62)
      y <- y - nl * lh - 0.012
    }
    xs <- seq(0.05, 0.55, length.out = 61); v <- seq(-3, 3, length.out = 61)
    rect(xs[-61], 0.05, xs[-1], 0.08, col = fc_col(v[-61]), border = NA)
    text(c(0.05, 0.30, 0.55), 0.03, c("-3", "0", "3"), cex = 0.6)
    text(0.30, 0.115, sprintf("log2FC (%s / %s), capped at +/-3",
                              tab$numerator[1], tab$denominator[1]), cex = 0.62)
    dev.off()

    index[[length(index) + 1]] <- data.frame(
      contrast = cn, gene_set = sn, n_pathways_padj = n_all, n_pathways_shown = k,
      n_proteins_shown = length(genes), n_links = nrow(links),
      pathways_shown = paste(re$Description, collapse = " | "))
    cat(sprintf("   %-16s %-28s %2d of %2d pathways, %3d proteins\n",
                cn, sn, k, n_all, length(genes)))
  }
}

if (length(index)) {
  write.csv(bind_rows(index), file.path(onto_dir, "Reactome_chord_index.csv"), row.names = FALSE)
  cat(sprintf("\n%d chord diagrams written; index: %s\n", length(index),
              file.path(onto_dir, "Reactome_chord_index.csv")))
} else cat("\nNo Reactome results found - nothing drawn.\n")
