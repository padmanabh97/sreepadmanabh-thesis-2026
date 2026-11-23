## =============================================================================
## Standalone re-plot of the factorial PCA (PCA_all16) - RNA-seq only
##
## Recomputes exactly the PCA of factorial_2x2x2_DESeq2_cellmeans.R (section 5):
##   same counts.txt / sample_sheet.tsv checks, same prefilter
##   (>= 10 reads in >= 2 of 16 libraries), vst(blind = TRUE),
##   plotPCA(ntop = 500, returnData = TRUE).
## Only the drawing changes:
##   - equal units on both axes (coord_fixed): 1 PC unit is the same length on
##     PC1 and PC2, so the PC2 spread is no longer visually inflated
##   - device height computed from the data range -> little empty space
##   - one point size; 0.3 mg/ml = open symbol with thick border, 1.5 mg/ml = solid
##   - larger points and labels; one label per condition (both replicates),
##     e.g. 2D_0.3_A, 3D_1.5_N, placed at the replicate midpoint;
##     the two replicates are joined by a thin line (CONNECT_REPS)
## Nothing is written to the DESeq2_Factorial_cellmeans folder except the two
## new figure files (new names; PCA_all16.* is not overwritten).
## =============================================================================

## ---- USER INPUT -------------------------------------------------------------
# setwd("D:/ESoftBio Lab/Omics Analysis/Oxo-Mech 2D vs 3D transcriptome/fragment counts 010726/Subset Analysis 060726")
counts_file       <- "counts.txt"
sample_sheet_file <- "sample_sheet.tsv"
out_dir           <- file.path("DESeq2_Factorial_cellmeans", "04_figures")
out_stem          <- "PCA_all16_compact"
PANEL_WIDTH_IN    <- 6.5      # width of the plotting panel; height follows the data
SHOW_TITLE        <- FALSE    # TRUE adds "PCA, blind VST, top 500 genes, 16 libraries"
CONNECT_REPS      <- TRUE     # thin line joining the two replicates of each condition

## ---- fixed (identical to the factorial script) ------------------------------
MIN_COUNT <- 10; MIN_SAMPLES <- 2; NTOP <- 500
EXPECTED <- list(dim = c("2D", "3D"), stiffness = c("0.3", "1.5"),
                 oxygen = c("anoxia", "normoxia"), rep = c("Rep1", "Rep2"))

suppressPackageStartupMessages({
  library(DESeq2); library(ggplot2); library(ggrepel)
})

## ---- 1. inputs (same parsing and checks as the factorial script) ------------
countData <- read.delim(counts_file, row.names = 1, skip = 1, check.names = FALSE)
stopifnot(identical(colnames(countData)[1:5], c("Chr", "Start", "End", "Strand", "Length")))
countData <- countData[, -(1:5), drop = FALSE]
colnames(countData) <- sub("\\.bam$", "", basename(colnames(countData)))

si <- read.delim(sample_sheet_file, colClasses = "character")
si <- si[, c("sample", "dim", "stiffness", "oxygen", "rep")]
for (f in names(EXPECTED)) stopifnot(all(si[[f]] %in% EXPECTED[[f]]))
stopifnot(identical(paste(si$dim, si$stiffness, si$oxygen, si$rep, sep = "_"), si$sample),
          nrow(si) == 16, !anyDuplicated(si$sample), all(si$sample %in% colnames(countData)))
rownames(si) <- si$sample
countData <- as.matrix(countData[, si$sample, drop = FALSE])
storage.mode(countData) <- "integer"
rownames(countData) <- sub("\\..*$", "", rownames(countData))
for (f in names(EXPECTED)) si[[f]] <- factor(si[[f]], levels = EXPECTED[[f]])
si$group <- factor(paste(si$dim, si$stiffness, si$oxygen, sep = "_"))

## ---- 2. PCA (identical computation) -----------------------------------------
dds <- DESeqDataSetFromMatrix(countData, si, design = ~ 0 + group)
dds <- dds[rowSums(counts(dds) >= MIN_COUNT) >= MIN_SAMPLES, ]
message("genes after prefilter: ", nrow(dds), "  (factorial run: 15884)")
vsd <- vst(dds, blind = TRUE)
pca <- plotPCA(vsd, intgroup = c("dim", "stiffness", "oxygen"), ntop = NTOP, returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"))

## ---- 3. one label per condition, at the midpoint of its two replicates ------
pca$cond <- paste(pca$dim, pca$stiffness, toupper(substr(pca$oxygen, 1, 1)), sep = "_")  # 2D_0.3_A
## symbols: shape = oxygen, fill = stiffness (0.3 open with thick border, 1.5 solid)
## all four use fillable shapes 21/24 so open and solid symbols are exactly the same size
DIM_COL <- c("2D" = "#E66101", "3D" = "#5E3C99")                 # orange / purple (PuOr, colour-blind safe)
SYM  <- c("anoxia, 0.3 mg/ml" = 21, "normoxia, 0.3 mg/ml" = 24,  # open (white fill)
          "anoxia, 1.5 mg/ml" = 21, "normoxia, 1.5 mg/ml" = 24)  # solid (fill = colour)
pca$sym  <- factor(paste0(pca$oxygen, ", ", pca$stiffness, " mg/ml"), levels = names(SYM))
pca$fill <- ifelse(pca$stiffness == "0.3", "white", DIM_COL[as.character(pca$dim)])
## sizes: ggplot draws diameter = size + stroke, so the solid symbols get a larger
## size to make all symbols the same overall diameter
SIZE_OPEN <- 4.2; STROKE_OPEN <- 1.6; EDGE <- 0.4                # EDGE = black outline width
st_open  <- STROKE_OPEN + 2 * EDGE; st_solid <- 0.6
pca$st <- ifelse(pca$stiffness == "0.3", st_open, st_solid)
pca$sz <- ifelse(pca$stiffness == "0.3", SIZE_OPEN,
                 SIZE_OPEN + (st_open - st_solid) * ggplot2::.stroke / 2 / ggplot2::.pt)
lab <- aggregate(cbind(PC1, PC2) ~ cond + dim, data = pca, FUN = mean)
## ggrepel: label rows plus the 16 points as unlabeled rows, so labels avoid points
## normoxia labels are pushed above their pair and anoxia labels below
## (LABEL_FLIP reverses this for named labels; edit if a label sits on a point)
LABEL_FLIP <- c("2D_0.3_A")       # 2D_0.3_A sits above 2D_1.5_A, so its label goes up
lab$ny <- ifelse(grepl("_N$", lab$cond), 1, -1) * ifelse(lab$cond %in% LABEL_FLIP, -1, 1) *
          0.12 * diff(range(pca$PC2))
rep_df <- rbind(data.frame(PC1 = lab$PC1, PC2 = lab$PC2, dim = lab$dim, txt = lab$cond, ny = lab$ny),
                data.frame(PC1 = pca$PC1, PC2 = pca$PC2, dim = pca$dim, txt = "", ny = 0))

## ---- 4. plot ----------------------------------------------------------------
xr <- range(pca$PC1); yr <- range(pca$PC2)
pad_x <- 0.06 * diff(xr); pad_y <- 0.16 * diff(yr)        # room for labels at the edges
xlim <- xr + c(-1, 1) * pad_x; ylim <- yr + c(-1, 1) * pad_y

p <- ggplot(pca, aes(PC1, PC2, colour = dim)) +
  { if (CONNECT_REPS) geom_line(aes(group = cond), linewidth = 0.5, alpha = 0.5, show.legend = FALSE) } +
  ## layer 1: every symbol with a black edge (solid: thin black outline;
  ##          open: a wider black ring that layer 2 covers except at its edges)
  geom_point(aes(shape = sym, fill = fill, size = sz, stroke = st), colour = "black") +
  ## layer 2: coloured thick border of the open (0.3 mg/ml) symbols
  geom_point(data = subset(pca, stiffness == "0.3"), aes(shape = sym),
             fill = "white", size = SIZE_OPEN, stroke = STROKE_OPEN) +
  geom_text_repel(data = rep_df, aes(label = txt), size = 4.2, fontface = "bold",
                  bg.color = "white", bg.r = 0.15, point.size = 8,         # white halo
                  nudge_y = rep_df$ny, box.padding = 0.5, point.padding = 0.7, force = 3, max.time = 3, max.iter = 1e5,
                  min.segment.length = 0.4,
                  segment.size = 0.3, max.overlaps = Inf, seed = 1, show.legend = FALSE) +
  scale_colour_manual(values = DIM_COL, name = NULL) +
  scale_fill_identity() + scale_size_identity() +
  scale_shape_manual(values = SYM, name = NULL) +
  guides(colour = guide_legend(order = 1, override.aes = list(shape = 15, size = 4.5)),
         shape  = guide_legend(order = 2, nrow = 2, byrow = FALSE,
                               override.aes = list(colour = c("grey25", "grey25", "black", "black"),
                                                   fill   = c("white", "white", "grey45", "grey45"),
                                                   size   = c(3.8, 3.8, 4.6, 4.6),
                                                   stroke = c(1.4, 1.4, 0.5, 0.5)))) +
  coord_fixed(ratio = 1, xlim = xlim, ylim = ylim, expand = FALSE) +   # equal units
  labs(x = sprintf("PC1 (%d%% variance)", pv[1]), y = sprintf("PC2 (%d%%)", pv[2]),
       title = if (SHOW_TITLE) "PCA, blind VST, top 500 genes, 16 libraries" else NULL) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-6, 0, 0, 0),
        legend.text = element_text(size = 12),
        panel.grid.minor = element_blank(),
        plot.margin = margin(4, 6, 2, 4))

## device size follows the data aspect, so there are no empty bands
panel_h <- PANEL_WIDTH_IN * diff(ylim) / diff(xlim)
fig_w <- PANEL_WIDTH_IN + 0.75
fig_h <- panel_h + 1.25 + if (SHOW_TITLE) 0.35 else 0
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(out_dir, paste0(out_stem, ".png")), p, width = fig_w, height = fig_h, dpi = 300)
ggsave(file.path(out_dir, paste0(out_stem, ".pdf")), p, width = fig_w, height = fig_h)
message(sprintf("saved %s.png/.pdf  (%.1f x %.1f in)", out_stem, fig_w, fig_h))
