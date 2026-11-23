## =============================================================================
## PID430 - volcano re-plot in the reference figure format
##
## Standalone companion to PID430_proteomics_pipeline.R (v5). Runs AFTER it and
## reads only its saved tables; no statistics are recomputed and no protein
## changes category. Each point, its position and its colour come from the run's
## per-contrast table:
##     x      = log2FC
##     y      = -log10(adj_p_BH)       (BH adjusted p, as in the pipeline volcanoes)
##     colour = the pipeline's own `significant` call, split by the sign of log2FC
##              Up        = significant, log2FC > 0 (higher in the numerator)
##              Down      = significant, log2FC < 0 (higher in the denominator)
##              Unchanged = everything else (incl. adj p < 0.05 below 1.5-fold)
## Dashed lines are the pipeline's thresholds (BH adj p 0.05, 1.5-fold), and the
## script stops if any row's `significant` flag disagrees with those thresholds,
## so the lines can never be drawn at values other than the ones that made the
## calls.
##
## Inputs  (inside the run folder):
##   03_differential_abundance/per_contrast_tables/<contrast>_results.csv  (6)
##   03_differential_abundance/interaction_oxygen_by_collagen.csv  (if INCLUDE_INTERACTION)
## Outputs (new, timestamped folder; nothing existing is touched):
##   04_figures/volcano_replot_<YYYYMMDD_HHMMSS>/volcano_<contrast>.png / .tiff
##   ... /volcano_replot_index.csv, sessionInfo.txt
##
## Usage:  Rscript PID430_volcano_replot.R <run_folder>
##     or  set RUN_DIR below and source() the file.
## =============================================================================

## ---- 1. Settings ------------------------------------------------------------
args    <- commandArgs(trailingOnly = TRUE)
RUN_DIR <- if (length(args) >= 1) args[1] else
  "PID430_proteomics_results_20261015_163542"   # <- or edit this

ADJP_CUTOFF <- 0.05          # must equal the pipeline's ADJP_CUTOFF (checked below)
FC_PRIMARY  <- 1.5           # must equal the pipeline's FC_PRIMARY (checked below)
LFC_PRIMARY <- log2(FC_PRIMARY)

N_LABELS            <- 25    # label the top N significant proteins by adj p
SHOW_TITLE          <- TRUE  # "<numerator> vs <denominator>" above each panel
INCLUDE_INTERACTION <- TRUE  # also re-plot the supplementary interaction volcano
FIG_W_IN <- 7.87; FIG_H_IN <- 8.66; FIG_DPI <- 300   # reference: 2362 x 2598 px at 300 dpi

COL_UP <- "red"; COL_DOWN <- "blue"; COL_NS <- "grey50"

## ---- 2. Packages and inputs -------------------------------------------------
for (p in c("ggplot2", "ggrepel")) if (!requireNamespace(p, quietly = TRUE))
  stop("Package '", p, "' is not installed.", call. = FALSE)
suppressPackageStartupMessages({ library(ggplot2); library(ggrepel) })

tab_dir <- file.path(RUN_DIR, "03_differential_abundance", "per_contrast_tables")
f_int   <- file.path(RUN_DIR, "03_differential_abundance", "interaction_oxygen_by_collagen.csv")
if (!dir.exists(tab_dir))
  stop("Not found (is this a pipeline run folder?): ", tab_dir, call. = FALSE)
f_tabs <- sort(list.files(tab_dir, "_results\\.csv$", full.names = TRUE))
if (length(f_tabs) != 6)
  stop("Expected 6 per-contrast tables in ", tab_dir, ", found ", length(f_tabs), call. = FALSE)
if (INCLUDE_INTERACTION && !file.exists(f_int))
  stop("Interaction table not found: ", f_int, call. = FALSE)

## RNA orientation check, done before anything is written. The RNA pairwise
## contrasts put the alphabetically earlier group first (0-3 A < 0-3 N < 1-5 A <
## 1-5 N), so every numerator must sort before its denominator. A run made before
## the diagonal was re-oriented (1-5_A_vs_0-3_N) stops here instead of producing a
## mirrored panel.
for (f in f_tabs) {
  h <- read.csv(f, check.names = FALSE, nrows = 1)
  if (!(h$numerator[1] < h$denominator[1]))
    stop(basename(f), ": numerator '", h$numerator[1], "' / denominator '", h$denominator[1],
         "' is the reverse of the RNA orientation. This run predates the ",
         "re-oriented pipeline; re-run PID430_proteomics_pipeline.R first.", call. = FALSE)
}

out_dir <- file.path(RUN_DIR, "04_figures",
                     paste0("volcano_replot_", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(out_dir, recursive = TRUE)

## ---- 3. One volcano per table -----------------------------------------------
plot_volcano <- function(d, x, title, name, xlab = expression(Log[2] ~ "fold change")) {
  need <- c(x, "adj_p_BH", "significant", "Gene_Symbol", "Accession")
  miss <- setdiff(need, names(d))
  if (length(miss)) stop(name, ": missing column(s) ", paste(miss, collapse = ", "), call. = FALSE)
  d$x   <- d[[x]]
  d$sig <- as.logical(d$significant)

  ## The calls must be exactly the pipeline's rule; otherwise the dashed lines
  ## would mark thresholds that did not produce the colours.
  rule <- d$adj_p_BH < ADJP_CUTOFF & abs(d$x) >= LFC_PRIMARY
  if (anyNA(d$x) || anyNA(d$adj_p_BH) || anyNA(d$sig) || any(rule != d$sig))
    stop(name, ": 'significant' does not match adj p < ", ADJP_CUTOFF, " and >= ",
         FC_PRIMARY, "-fold on ", sum(rule != d$sig, na.rm = TRUE),
         " rows - check ADJP_CUTOFF / FC_PRIMARY against the pipeline.", call. = FALSE)

  d$y   <- -log10(d$adj_p_BH)
  d$cls <- factor(ifelse(!d$sig, "Unchanged", ifelse(d$x > 0, "Up", "Down")),
                  levels = c("Down", "Unchanged", "Up"))
  d     <- d[order(d$cls == "Unchanged", decreasing = TRUE), ]  # coloured points on top
  n_up  <- sum(d$cls == "Up"); n_dn <- sum(d$cls == "Down")

  lab <- d[d$sig, ]; lab <- head(lab[order(lab$adj_p_BH, -abs(lab$x)), ], N_LABELS)
  lab$label <- ifelse(is.na(lab$Gene_Symbol) | lab$Gene_Symbol %in% c("", "NA"),
                      lab$Accession, lab$Gene_Symbol)

  xl <- c(min(d$x) - 1.5, max(d$x) + 1.5)
  yl <- c(0, max(d$y) * 1.18 + 1)   # headroom so the n= counts sit above the data

  p <- ggplot(d, aes(x, y)) +
    geom_point(aes(colour = cls), size = 1.5, alpha = 0.75) +
    geom_vline(xintercept = c(-LFC_PRIMARY, LFC_PRIMARY), linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = -log10(ADJP_CUTOFF), linetype = "dashed", linewidth = 0.4) +
    ## Labels kept below the n= counts; a thin leader appears only where a label
    ## had to move away from its point to avoid overlapping another label.
    geom_text_repel(data = lab, aes(label = label), size = 5, colour = "black",
                    ylim = c(NA, yl[2] * 0.9), box.padding = 0.35,
                    min.segment.length = 0.4, segment.colour = "grey40",
                    segment.size = 0.3, max.overlaps = Inf, max.iter = 20000, seed = 1) +
    annotate("text", x = xl[1] + 0.25 * diff(xl), y = yl[2] * 0.97, label = paste0("n=", n_dn),
             colour = "mediumblue", fontface = "bold", size = 6) +
    annotate("text", x = xl[1] + 0.75 * diff(xl), y = yl[2] * 0.97, label = paste0("n=", n_up),
             colour = COL_UP, fontface = "bold", size = 6) +
    scale_colour_manual(values = c(Down = COL_DOWN, Unchanged = COL_NS, Up = COL_UP),
                        drop = FALSE) +
    scale_x_continuous(limits = xl) +
    scale_y_continuous(limits = yl) +
    labs(x = xlab,
         y = expression(-Log[10] ~ "adjusted" ~ italic(P)),
         colour = NULL, title = if (SHOW_TITLE) title else NULL,
         caption = paste0("total = ", nrow(d), " proteins")) +
    guides(colour = guide_legend(override.aes = list(size = 5))) +
    theme_classic(base_size = 14) +
    theme(axis.text       = element_text(size = 18),
          axis.title      = element_text(size = 12),
          axis.line       = element_line(linewidth = 0.6),
          plot.title      = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.caption    = element_text(size = 14),
          legend.position = "bottom",
          legend.text     = element_text(size = 11))

  f_out <- file.path(out_dir, paste0("volcano_", name))
  ggsave(paste0(f_out, ".png"),  p, width = FIG_W_IN, height = FIG_H_IN, dpi = FIG_DPI)
  ggsave(paste0(f_out, ".tiff"), p, width = FIG_W_IN, height = FIG_H_IN, dpi = FIG_DPI,
         compression = "lzw")
  data.frame(plot = name, n_total = nrow(d), n_up = n_up, n_down = n_dn,
             labelled = paste(lab$label, collapse = ";"))
}

index <- list()
for (f in f_tabs) {
  d  <- read.csv(f, check.names = FALSE)
  cn <- sub("_results\\.csv$", "", basename(f))
  index[[cn]] <- plot_volcano(d, "log2FC",
    sprintf("%s vs %s", d$numerator[1], d$denominator[1]), cn)
}
if (INCLUDE_INTERACTION) {
  ## x = difference of differences; Up = anoxia response larger at 0.3 mg/ml,
  ## Down = larger at 1.5 mg/ml (the pipeline's interpretation column).
  di <- read.csv(f_int, check.names = FALSE)
  index[["INTERACTION"]] <- plot_volcano(di, "interaction_log2",
    "Oxygen x collagen interaction\n(Up = anoxia response larger at 0.3 mg/ml)",
    "INTERACTION_oxygen_by_collagen",
    xlab = expression("Interaction estimate (" * log[2] * ", difference of differences)"))
}

idx <- do.call(rbind, index)
write.csv(idx, file.path(out_dir, "volcano_replot_index.csv"), row.names = FALSE)
writeLines(c(paste("run folder:", normalizePath(RUN_DIR)),
             paste("tables:", paste(basename(c(f_tabs, if (INCLUDE_INTERACTION) f_int)), collapse = "; ")),
             sprintf("thresholds: BH adj p < %s, fold change >= %s (|log2FC| >= %.7f)",
                     ADJP_CUTOFF, FC_PRIMARY, LFC_PRIMARY),
             paste("written:", Sys.time()), "", capture.output(sessionInfo())),
           file.path(out_dir, "sessionInfo.txt"))
print(idx[, c("plot", "n_total", "n_up", "n_down")], row.names = FALSE)
cat("\nWritten to", out_dir, "\n")
