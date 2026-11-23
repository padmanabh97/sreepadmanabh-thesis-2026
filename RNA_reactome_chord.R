## =============================================================================
## RNA-seq - Reactome ORA + chord diagrams (companion script)
##
## Runs AFTER any of these, reads only their saved results tables, and never
## re-fits DESeq2:
##   2D_only_DESeq2_060726.R            -> "<2D output folder>/CSV/<g1>_vs_<g2>_results.csv"
##   3D_only_DESeq2_060726.R            -> "<3D output folder>/CSV/<g1>_vs_<g2>_results.csv"
##   factorial_2x2x2_DESeq2_cellmeans.R -> "<run folder>/01_tables/CSV/<effect>_results.csv"
## The folder type is detected from its layout. Each folder is processed on its
## own: universes, lists and multiple-testing never mix 2D, 3D and factorial.
##
## Statistics (the only new statistics; they mirror the RNA GO ORA settings)
##   - lists      : regulation == "Up" / "Down" exactly as written by the main
##                  script (padj < 0.05 and |raw log2FC| > 1); nothing re-called
##   - universe   : genes with non-NA padj in the same table (= GO ORA universe)
##   - test       : ReactomePA::enrichPathway, organism "mouse", BH,
##                  pvalueCutoff 0.05, qvalueCutoff 0.2, gene-set size 10-500
##   - minimum    : >= 10 Entrez-mapped genes in the list (GO rule: >= 10 genes)
##   - ID mapping : Ensembl -> Entrez with org.Mm.eg.db. Ensembl ids with more
##                  than one Entrez id are dropped; when several Ensembl ids share
##                  one Entrez id, the most significant one (lowest padj) is kept.
##                  The same map is applied to list and universe, so the list is
##                  always a subset of the universe. Genes lost at each step are
##                  counted in Reactome_index.csv.
##
## Figures (per list with >= 1 pathway at BH adj p < 0.05)
##   - chord: top TOP_N pathways (P1 = smallest adj p) linked to their member
##     genes, genes coloured by log2FC; layout and colours as PID430_reactome_chord.R
##   - when the member genes of those pathways exceed MAX_GENES, each pathway
##     keeps its members with the largest |log2FC|, taken in turn across
##     pathways until MAX_GENES genes are placed. The legend prints full member
##     counts and the number shown; the full membership is in *_chord_links.csv.
##     MAX_GENES <- Inf reproduces the uncapped proteomics behaviour.
##   - dot plot: top 15 pathways (clusterProfiler::dotplot), as for GO
##
## Outputs
##   pairwise : "<folder>/Reactome/<contrast>_<Up|Down>_Reactome.csv | _chord.png | _chord_links.csv | _dotplot.png"
##   factorial: "<run folder>/07_Reactome_ORA/<effect>_<Up|Down>_..." (same files)
##   both     : Reactome_index.csv (every list, including not-run and 0-pathway
##              lists) and Reactome_methods.txt (settings + package versions)
##
## Usage:  Rscript RNA_reactome_chord.R  "<folder 1>" ["<folder 2>" ...]
##         or set RUN_DIRS below and source() the script.
## =============================================================================

## ---- 1. USER INPUT ----------------------------------------------------------
RUN_DIRS <- c("2D samples only",
              "3D samples only",
              "DESeq2_Factorial_cellmeans")

## ---- 2. FIXED SETTINGS (keep in step with the main scripts) -----------------
ORA_P        <- 0.05   # BH-adjusted p cut-off (main scripts: pvalueCutoff 0.05)
ORA_Q        <- 0.2    # q-value cut-off       (main scripts: qvalueCutoff 0.2)
MIN_GENES    <- 10     # list size needed to run the test (main scripts: 10)
MIN_GS       <- 10     # Reactome gene-set size limits (enrichPathway defaults)
MAX_GS       <- 500
TOP_N        <- 10     # pathways per chord (as PID430_reactome_chord.R)
MAX_GENES    <- 60     # gene sectors per chord before capping; Inf = no cap
FC_CAP       <- 3      # colour scale saturates at +/- FC_CAP (as proteomics)
FIG_DPI      <- 300

args <- commandArgs(trailingOnly = TRUE)
if (length(args)) RUN_DIRS <- args

suppressPackageStartupMessages({
  library(dplyr); library(circlize); library(ReactomePA)
  library(clusterProfiler); library(enrichplot); library(ggplot2)
  library(org.Mm.eg.db); library(AnnotationDbi)
})
options(warn = 1)

## Orientation text for the factorial effects (copied from EFFECT_INFO in
## factorial_2x2x2_DESeq2_cellmeans.R); pairwise contrasts are read from the name.
FACTORIAL_POSITIVE <- c(
  Oxygen_A_vs_N        = "higher in anoxia",
  Stiffness_0.3_vs_1.5 = "higher at 0.3 mg/ml",
  Dim_3D_vs_2D         = "higher in 3D",
  DimxOxygen           = "anoxia response (A - N) more positive in 3D than in 2D",
  DimxStiffness        = "stiffness response (0.3 - 1.5) more positive in 3D than in 2D",
  StiffnessxOxygen     = "anoxia response (A - N) more positive at 0.3 than at 1.5",
  DimxStiffxOxygen     = "(A-N at 0.3 minus A-N at 1.5) more positive in 3D than in 2D")
FACTORIAL_NEGATIVE <- c(
  Oxygen_A_vs_N        = "higher in normoxia",
  Stiffness_0.3_vs_1.5 = "higher at 1.5 mg/ml",
  Dim_3D_vs_2D         = "higher in 2D",
  DimxOxygen           = "anoxia response (A - N) more negative in 3D than in 2D",
  DimxStiffness        = "stiffness response (0.3 - 1.5) more negative in 3D than in 2D",
  StiffnessxOxygen     = "anoxia response (A - N) more negative at 0.3 than at 1.5",
  DimxStiffxOxygen     = "(A-N at 0.3 minus A-N at 1.5) more negative in 3D than in 2D")

orientation <- function(name, type) {
  if (type == "factorial") {
    pos <- FACTORIAL_POSITIVE[name]; neg <- FACTORIAL_NEGATIVE[name]
    if (is.na(pos)) { pos <- "log2FC > 0"; neg <- "log2FC < 0" }
    return(list(up = unname(pos), down = unname(neg),
                scale = sprintf("log2FC (%s), positive = %s", name, pos)))
  }
  g <- strsplit(name, "_vs_", fixed = TRUE)[[1]]
  if (length(g) != 2) stop("Cannot read numerator/denominator from contrast name: ", name, call. = FALSE)
  list(up = paste("higher in", g[1]), down = paste("higher in", g[2]),
       scale = sprintf("log2FC (%s / %s)", g[1], g[2]))
}

## ---- 3. Ensembl -> Entrez, one row per Entrez id ---------------------------
map_table <- function(df) {
  df <- df[!is.na(df$padj), ]
  df <- df[order(df$padj, df$pvalue, na.last = TRUE), ]
  mp <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db, keys = unique(df$gene_id),
                                               keytype = "ENSEMBL", columns = "ENTREZID"))
  mp <- mp[!is.na(mp$ENTREZID), ]
  n_ent <- table(mp$ENSEMBL)
  multi <- names(n_ent)[n_ent > 1]
  mp <- mp[!mp$ENSEMBL %in% multi, ]
  df$entrez <- mp$ENTREZID[match(df$gene_id, mp$ENSEMBL)]
  no_entrez <- sum(is.na(df$entrez)) - sum(df$gene_id %in% multi)
  m <- df[!is.na(df$entrez), ]
  dup <- duplicated(m$entrez)                    # table is sorted: first = lowest padj
  list(tab = m[!dup, ], n_tested = nrow(df), n_no_entrez = no_entrez,
       n_multi_entrez = length(multi), n_shared_entrez_dropped = sum(dup))
}

## ---- 4. chord ------------------------------------------------------------------
fc_col <- colorRamp2(c(-FC_CAP, 0, FC_CAP), c("#4E79A7", "white", "#E15759"))

draw_chord <- function(re, tab, title_txt, scale_txt, out_png, out_links) {
  k <- nrow(re); re$P <- paste0("P", seq_len(k))
  links <- do.call(rbind, lapply(seq_len(k), function(i)
    data.frame(from = re$P[i], pathway = re$Description[i],
               entrez = strsplit(re$geneID[i], "/")[[1]], stringsAsFactors = FALSE)))
  miss <- setdiff(links$entrez, tab$entrez)
  if (length(miss)) stop(title_txt, ": ", length(miss), " Reactome gene id(s) not in the list's table",
                         call. = FALSE)
  links$gene   <- tab$gene_name[match(links$entrez, tab$entrez)]
  links$log2FC <- tab$log2FoldChange[match(links$entrez, tab$entrez)]

  ## gene cap: round-robin across pathways, largest |log2FC| first within each
  genes_all <- unique(links$entrez)
  keep <- genes_all
  if (length(genes_all) > MAX_GENES) {
    per <- lapply(re$P, function(p) {
      x <- links[links$from == p, ]; x$entrez[order(-abs(x$log2FC))] })
    keep <- character(0); i <- 1L
    while (length(keep) < MAX_GENES && any(lengths(per) >= i)) {
      for (v in per) if (length(v) >= i && !(v[i] %in% keep) && length(keep) < MAX_GENES)
        keep <- c(keep, v[i])
      i <- i + 1L
    }
  }
  links$shown <- links$entrez %in% keep
  write.csv(links[, c("from", "pathway", "entrez", "gene", "log2FC", "shown")], out_links, row.names = FALSE)

  L <- links[links$shown, ]
  ## one sector label per gene: symbol, or Entrez id if none; duplicates made unique
  u <- L[!duplicated(L$entrez), ]
  lab_map <- setNames(make.unique(ifelse(is.na(u$gene) | u$gene == "", u$entrez, u$gene)), u$entrez)
  L$to <- unname(lab_map[L$entrez]); L$value <- 1
  g <- L %>% group_by(to) %>%
    summarise(first = min(as.integer(sub("P", "", from))), fc = log2FC[1], .groups = "drop") %>%
    arrange(first, desc(fc))
  genes <- g$to
  pcols <- setNames(hcl.colors(max(k, 2), "Dark 3")[seq_len(k)], re$P)
  gcols <- setNames(fc_col(pmax(pmin(g$fc, FC_CAP), -FC_CAP)), genes)
  shown_n <- vapply(re$P, function(p) sum(L$from == p), 1L)

  png(out_png, width = 11, height = 7, units = "in", res = FIG_DPI)
  layout(matrix(1:2, 1), widths = c(6.2, 4.8))
  par(mar = c(0.5, 0.5, 3.2, 0.5), xpd = NA)
  circos.clear()
  circos.par(start.degree = 90, clock.wise = TRUE,
             gap.after = c(rep(2, k - 1), 12, rep(0.6, length(genes) - 1), 12))
  chordDiagram(L[, c("from", "to", "value")], order = c(re$P, genes),
               grid.col = c(pcols, gcols), transparency = 0.35,
               annotationTrack = "grid", annotationTrackHeight = 0.05,
               preAllocateTracks = list(track.height = 0.12))
  circos.track(track.index = 1, bg.border = NA, panel.fun = function(x, y) {
    nm <- get.cell.meta.data("sector.index")
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], nm, facing = "clockwise",
                niceFacing = TRUE, adj = c(0, 0.5),
                cex = if (nm %in% re$P) 0.9 else 0.55, font = if (nm %in% re$P) 2 else 1)
  })
  circos.clear()
  title(main = title_txt, cex.main = 1, line = 1)

  par(mar = c(0.5, 0, 3.2, 0.5)); plot.new()
  lab <- sprintf("%s  %s\n      n = %d%s, adj p = %s", re$P,
                 vapply(re$Description, function(z) paste(strwrap(z, 50), collapse = "\n      "), ""),
                 re$Count, ifelse(shown_n < re$Count, sprintf(" (%d shown)", shown_n), ""),
                 formatC(re$p.adjust, format = "e", digits = 1))
  head_txt <- sprintf("Top %d of %d pathways, BH adj p < %.2g", k, attr(re, "n_all"), ORA_P)
  if (length(keep) < length(genes_all))
    head_txt <- sprintf("%s\n%d of %d genes shown (largest |log2FC| per pathway)", head_txt,
                        length(keep), length(genes_all))
  text(0, 0.98, head_txt, adj = c(0, 1), font = 2, cex = 0.7)
  lh <- 0.026; y <- if (length(keep) < length(genes_all)) 0.90 else 0.93
  for (i in seq_len(k)) {
    nl <- lengths(regmatches(lab[i], gregexpr("\n", lab[i]))) + 1
    rect(0, y - 0.018, 0.03, y, col = pcols[i], border = NA)
    text(0.045, y, lab[i], adj = c(0, 1), cex = 0.62)
    y <- y - nl * lh - 0.012
  }
  xs <- seq(0.05, 0.55, length.out = 61); v <- seq(-FC_CAP, FC_CAP, length.out = 61)
  rect(xs[-61], 0.05, xs[-1], 0.08, col = fc_col(v[-61]), border = NA)
  text(c(0.05, 0.30, 0.55), 0.03, c(-FC_CAP, 0, FC_CAP), cex = 0.6)
  text(0.30, 0.115, paste(strwrap(sprintf("%s, capped at +/-%g", scale_txt, FC_CAP), 60),
                          collapse = "\n"), cex = 0.62)
  dev.off()
  c(n_genes_members = length(genes_all), n_genes_shown = length(keep))
}

## ---- 5. per folder ------------------------------------------------------------
for (rd in RUN_DIRS) {
  if (dir.exists(file.path(rd, "01_tables", "CSV"))) {
    type <- "factorial"; csv_dir <- file.path(rd, "01_tables", "CSV"); out_dir <- file.path(rd, "07_Reactome_ORA")
  } else if (dir.exists(file.path(rd, "CSV"))) {
    type <- "pairwise";  csv_dir <- file.path(rd, "CSV");              out_dir <- file.path(rd, "Reactome")
  } else stop("Not a results folder of the pairwise or factorial scripts: ",
              normalizePath(rd, mustWork = FALSE), call. = FALSE)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  files <- sort(list.files(csv_dir, pattern = "_results\\.csv$"))
  if (!length(files)) stop("No *_results.csv in ", csv_dir, call. = FALSE)
  cat(sprintf("\n=== %s  [%s]  %d tables ===\n", rd, type, length(files)))
  index <- list()

  for (f in files) {
    name <- sub("_results\\.csv$", "", f)
    df <- read.csv(file.path(csv_dir, f), stringsAsFactors = FALSE, check.names = FALSE)
    need <- c("gene_id", "gene_name", "log2FoldChange", "pvalue", "padj", "regulation")
    if (!all(need %in% names(df))) stop(f, ": missing column(s) ", paste(setdiff(need, names(df)), collapse = ", "),
                                        call. = FALSE)
    ori <- orientation(name, type)
    mt  <- map_table(df)
    uni <- mt$tab$entrez

    for (dir_ in c("Up", "Down")) {
      tag <- paste0(name, "_", dir_)
      n_ens <- sum(df$regulation == dir_, na.rm = TRUE)
      g   <- mt$tab$entrez[mt$tab$regulation == dir_]
      status <- "run"; n_all <- 0L; shown <- c(n_genes_members = 0, n_genes_shown = 0); top_txt <- ""
      if (length(g) < MIN_GENES) {
        status <- sprintf("not run (%d Entrez genes < %d)", length(g), MIN_GENES)
      } else {
        er <- tryCatch(enrichPathway(gene = g, organism = "mouse", universe = uni,
                                     pvalueCutoff = ORA_P, qvalueCutoff = ORA_Q, pAdjustMethod = "BH",
                                     minGSSize = MIN_GS, maxGSSize = MAX_GS, readable = FALSE),
                       error = function(e) { status <<- paste("ERROR:", conditionMessage(e)); NULL })
        if (!is.null(er)) {
          re <- as.data.frame(er)
          re <- re[re$p.adjust < ORA_P, , drop = FALSE]
          re <- re[order(re$p.adjust, re$pvalue), , drop = FALSE]
          ## readable symbols for the CSV only; geneID stays Entrez for the chord
          re$geneSymbol <- vapply(strsplit(re$geneID, "/"), function(z)
            paste(mt$tab$gene_name[match(z, mt$tab$entrez)], collapse = "/"), "")
          write.csv(re, file.path(out_dir, paste0(tag, "_Reactome.csv")), row.names = FALSE)
          n_all <- nrow(re)
          if (n_all > 0) {
            top_txt <- paste(head(re$Description, 3), collapse = " | ")
            p_d <- dotplot(er, showCategory = 15) +
              labs(title = paste(tag, "- Reactome ORA"),
                   subtitle = if (dir_ == "Up") ori$up else ori$down) + theme_minimal()
            ggsave(file.path(out_dir, paste0(tag, "_Reactome_dotplot.png")), plot = p_d,
                   width = 10, height = 8, dpi = FIG_DPI)
            top <- head(re, TOP_N); attr(top, "n_all") <- n_all
            shown <- draw_chord(top, mt$tab,
                                sprintf("Reactome: %s\n%s (%s)", gsub("_", " ", name), dir_,
                                        if (dir_ == "Up") ori$up else ori$down),
                                ori$scale,
                                file.path(out_dir, paste0(tag, "_Reactome_chord.png")),
                                file.path(out_dir, paste0(tag, "_Reactome_chord_links.csv")))
          }
        }
      }
      cat(sprintf("  %-44s list %4d Ensembl -> %4d Entrez | %s | %3d pathways | chord genes %s\n",
                  tag, n_ens, length(g), status, n_all,
                  if (n_all) sprintf("%d/%d", shown[2], shown[1]) else "-"))
      index[[tag]] <- data.frame(
        folder = rd, type = type, contrast_or_effect = name, direction = dir_,
        direction_means = if (dir_ == "Up") ori$up else ori$down,
        list_ensembl = n_ens, list_entrez = length(g),
        universe_tested = mt$n_tested, universe_entrez = length(uni),
        universe_lost_no_entrez = mt$n_no_entrez, universe_lost_multi_entrez = mt$n_multi_entrez,
        universe_lost_shared_entrez = mt$n_shared_entrez_dropped,
        status = status, n_pathways_padj = n_all,
        chord_genes_members = shown[[1]], chord_genes_shown = shown[[2]],
        top3_pathways = top_txt, stringsAsFactors = FALSE)
    }
  }
  write.csv(bind_rows(index), file.path(out_dir, "Reactome_index.csv"), row.names = FALSE)

  meth <- c(
    sprintf("Reactome over-representation analysis (%s results: %s).", type, rd),
    "Gene lists: genes called Up or Down by the main script (padj < 0.05 and |raw log2FC| > 1), unchanged.",
    "Universe: all genes with a non-NA adjusted p-value in the same contrast/effect table.",
    paste("Ensembl ids were mapped to Entrez with org.Mm.eg.db; ids with more than one Entrez id were",
          "excluded, and where several Ensembl ids shared one Entrez id the most significant was kept.",
          "The same map was applied to lists and universe."),
    sprintf(paste("Test: ReactomePA::enrichPathway (organism mouse), hypergeometric, BH correction,",
                  "adjusted p < %.2g and q < %.2g, gene sets of %d-%d genes, lists with >= %d genes."),
            ORA_P, ORA_Q, MIN_GS, MAX_GS, MIN_GENES),
    sprintf(paste("Chord diagrams show the %d pathways with the smallest adjusted p and their member genes",
                  "coloured by log2FC (capped at +/-%g); when members exceeded %s genes, the genes with the",
                  "largest |log2FC| in each pathway were shown, taken in turn across pathways."),
            TOP_N, FC_CAP, format(MAX_GENES)),
    "", "Package versions:",
    vapply(c("ReactomePA", "reactome.db", "clusterProfiler", "org.Mm.eg.db", "circlize", "enrichplot"),
           function(p) sprintf("  %s %s", p, tryCatch(as.character(packageVersion(p)), error = function(e) "n/a")), ""),
    sprintf("  %s", R.version.string), sprintf("Run: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  writeLines(meth, file.path(out_dir, "Reactome_methods.txt"))
  cat(sprintf("  -> %s\n", out_dir))
}
cat("\nDone.\n")
