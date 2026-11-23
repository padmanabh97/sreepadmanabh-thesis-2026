# =============================================================================
# SUBSET GSEA -- 2D SAMPLES ONLY (Oxo-Mech)
# -----------------------------------------------------------------------------
# Companion to 2D_only_DESeq2_060726.R -- same subsetting philosophy, GSEA output.
# NO GLOBAL NORMALIZATION: the sample set is filtered to dim == "2D" BEFORE the
# DESeqDataSet is built, so size factors, dispersions and Wald statistics are all
# estimated from the 8 2D libraries alone (the 3D libraries never enter this fit).
#
# Design:  ~ group   where group = paste(stiffness, oxygen)  -> 4 levels
#          combn(4, 2) = 6 pairwise contrasts.
# For each contrast it runs gseGO (GO:BP) on the FULL ranked gene list -- no
# padj/LFC cutoff -- and writes one GSEA table + one dotplot per contrast.
#
# Ranking metric (identical to the global pairwise GSEA, so results stay comparable):
#     rank = sign(log2FoldChange) * -log10(pvalue)
#
# Output: <output_dir>/GSEA/<contrast>_GSEA_GO_BP.csv  (+ _dotplot.png)
# =============================================================================

# ============================================================================
# >>>>>>>>>>>>>>>>>  SETTINGS TO ALTER TO RUN IN YOUR FOLDER  <<<<<<<<<<<<<<<<<
# ============================================================================
# 1. WORKING DIRECTORY. Point R at the folder that holds counts.txt and
#    sample_sheet.tsv. On Windows, use EITHER forward slashes OR doubled
#    backslashes -- a single backslash is an escape character in R and will
#    error. Uncomment ONE of these:
#
setwd("D:/ESoftBio Lab/Omics Analysis/Oxo-Mech 2D vs 3D transcriptome/fragment counts 010726/Subset Analysis 060726")
# setwd("D:\\ESoftBio Lab\\Omics Analysis\\Oxo-Mech 2D vs 3D transcriptome\\fragment counts 010726\\Subset Analysis 060726")
#
# 2. INPUT FILENAMES. If your counts / sample-sheet files are named differently,
#    change the defaults on the next lines (or leave as-is if they match).
counts_file       <- "counts.txt"
sample_sheet_file <- "sample_sheet.tsv"
output_dir        <- "2D samples only"          # GSEA/ subfolder is created below
#
# 3. ANNOTATION LIBRARY (org.Mm.eg.db, GO.db). On your RStudio machine these are
#    in your normal user library, so nothing to set. (On the sandbox they live in
#    a workspace-local Rlib and R_LIBS_USER must be exported before Rscript --
#    NOT needed on your machine.)
#
# 4. SIGNIFICANCE CUTOFF for reported gene sets (BH-adjusted). Default 0.05.
pval_cutoff <- 0.05
# ============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
})

set.seed(1234)   # gseGO permutations are stochastic; fix seed for reproducibility
ts <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")

dir.create(file.path(output_dir, "GSEA"), recursive = TRUE, showWarnings = FALSE)

# ====================== 1. LOAD DATA ==========================================
ts("Loading counts and sample sheet")
countData <- read.delim(counts_file, row.names = 1, skip = 1, check.names = FALSE)
colnames(countData) <- basename(colnames(countData))
colnames(countData) <- sub("\\.bam$", "", colnames(countData))

sample_info <- read.delim(sample_sheet_file, stringsAsFactors = FALSE)

# ---- SUBSET TO 2D SAMPLES ONLY (before any DESeq step) ----
sample_info_sub <- sample_info[sample_info$dim == "2D", ]
valid_samples   <- sample_info_sub$sample
count_columns   <- intersect(valid_samples, colnames(countData))

countData <- as.matrix(countData[, count_columns, drop = FALSE])
mode(countData) <- "numeric"

sampleTable <- sample_info_sub[sample_info_sub$sample %in% count_columns, ]
sampleTable$group <- factor(paste(sampleTable$stiffness, sampleTable$oxygen, sep = "_"))
rownames(sampleTable) <- sampleTable$sample

countData <- countData[, rownames(sampleTable), drop = FALSE]
stopifnot(identical(colnames(countData), rownames(sampleTable)))
ts("Subset to", ncol(countData), "2D libraries |", nlevels(sampleTable$group), "groups:",
   paste(levels(sampleTable$group), collapse = ", "))

# ====================== 2. DESeq2  (subset fit -- no global normalization) ====
dds <- DESeqDataSetFromMatrix(countData = countData, colData = sampleTable, design = ~ group)

# Prefilter (matches 2D_only_DESeq2_060726.R): keep genes with >=10 counts in >=2 samples
keep <- rowSums(counts(dds) >= 10) >= 2
dds  <- dds[keep, ]
ts("Genes retained after prefilter:", nrow(dds))

dds <- DESeq(dds)

groups      <- levels(dds$group)
comparisons <- combn(groups, 2, simplify = FALSE)
ts("Running GSEA for", length(comparisons), "pairwise contrasts")

# ====================== 3. GSEA PER CONTRAST ==================================
gsea_made <- 0L
gsea_summary <- list()

for (comp in comparisons) {
  g1 <- comp[1]; g2 <- comp[2]
  contrast_name <- paste0(g1, "_vs_", g2)

  res    <- results(dds, contrast = c("group", g1, g2), alpha = 0.05)
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)

  # ---- Build ranked gene list (ALL tested genes, no cutoff) ----
  # Metric = sign(log2FC) * -log10(pvalue). See note at bottom on raw vs
  # shrunken LFC: because only the SIGN of log2FC enters here, the choice of
  # raw vs shrunken LFC is (almost) irrelevant for THIS metric.
  rank_df <- res_df[!is.na(res_df$pvalue) & !is.na(res_df$log2FoldChange), ]
  ranks   <- sign(rank_df$log2FoldChange) * -log10(rank_df$pvalue)
  names(ranks) <- rank_df$gene_id

  # collapse duplicate gene IDs (keep max abs) and sort decreasing
  ranks <- ranks[!duplicated(names(ranks))]
  ranks <- sort(ranks, decreasing = TRUE)

  if (length(ranks) < 10) { ts("  skip", contrast_name, "- too few ranked genes"); next }

  ego <- tryCatch(
    gseGO(geneList     = ranks,
          OrgDb        = org.Mm.eg.db,
          keyType      = "ENSEMBL",
          ont          = "BP",
          minGSSize    = 10,
          maxGSSize    = 500,
          pvalueCutoff = pval_cutoff,
          verbose      = FALSE,
          seed         = TRUE),
    error = function(e) { ts("  gseGO error for", contrast_name, ":", conditionMessage(e)); NULL })

  n_sets <- if (!is.null(ego)) nrow(as.data.frame(ego)) else 0L

  if (!is.null(ego) && n_sets > 0) {
    ego_r <- setReadable(ego, OrgDb = org.Mm.eg.db, keyType = "ENSEMBL")
    write.csv(as.data.frame(ego_r),
              file.path(output_dir, "GSEA", paste0(contrast_name, "_GSEA_GO_BP.csv")),
              row.names = FALSE)
    
    p <- dotplot(ego, showCategory = 10, split = ".sign") +
      facet_grid(. ~ .sign) + ggtitle(paste(contrast_name, "- GSEA GO BP")) +
      scale_x_continuous(expand = expansion(mult = c(0.10, 0.12))) +
      coord_cartesian(clip = "off") +
      theme(axis.text.y = element_text(size = 8),
            strip.text  = element_text(size = 20),
            panel.spacing.x = unit(0.3, "lines"),
            plot.margin     = margin(4, 6, 4, 4)
            )
    
    #p <- dotplot(ego, showCategory = 15, split = ".sign") +
         #facet_grid(. ~ .sign) + ggtitle(paste(contrast_name, "- GSEA GO BP"))
   
     ggsave(file.path(output_dir, "GSEA", paste0(contrast_name, "_GSEA_dotplot.png")),
           plot = p, width = 6, height = 8, dpi = 300)
    gsea_made <- gsea_made + 1L
  }

  gsea_summary[[contrast_name]] <- data.frame(
    Contrast      = contrast_name,
    N_ranked      = length(ranks),
    N_sig_sets    = n_sets)
  ts("  ", contrast_name, "->", n_sets, "enriched GO:BP sets")
}

# ====================== 4. SUMMARY ============================================
gsea_summary_df <- bind_rows(gsea_summary) %>% arrange(desc(N_sig_sets))
write.csv(gsea_summary_df,
          file.path(output_dir, "GSEA", "GSEA_Summary_2D_Only.csv"),
          row.names = FALSE)

ts("DONE (2D only). Contrasts:", length(comparisons),
   "| dotplots written:", gsea_made)
cat("\nTop contrasts by # enriched GO:BP sets:\n")
print(utils::head(gsea_summary_df, 10))

# =============================================================================
# NOTE -- HOW GSEA CHANGES WITH RAW vs SHRUNKEN log2FC
# -----------------------------------------------------------------------------
# It depends ENTIRELY on the ranking metric you choose:
#
#  (a) THIS SCRIPT'S METRIC: sign(log2FC) * -log10(pvalue)
#      Only the SIGN of the fold change enters -- the magnitude comes from the
#      p-value, which DESeq2 computes from the raw (unshrunk) Wald statistic and
#      which shrinkage does NOT change. So switching raw <-> shrunken LFC leaves
#      this ranking essentially untouched. The only genes that could move are
#      those whose LFC sits so close to 0 that shrinkage flips its sign -- and
#      those have p ~ 1 (rank ~ 0), so they sit in the middle of the list and
#      don't affect the leading edge. => raw vs shrunken is a non-issue here.
#
#  (b) IF YOU RANK BY log2FC ITSELF (a common alternative, e.g. rank = log2FC):
#      Now it matters a lot. RAW LFCs at low counts are wildly inflated -- a gene
#      with 2 vs 0 reads can show |LFC| > 5 with no real evidence. Ranking by raw
#      LFC pushes these noisy, low-information genes to the extremes of the list,
#      where GSEA weights them most heavily -> unstable, irreproducible leading
#      edges. SHRUNKEN LFCs (apeglm/ashr) pull those low-information genes toward
#      0 in proportion to their uncertainty, so the ranking reflects genes with
#      both large AND well-estimated effects. At n=2/group this is a large effect.
#      => if ranking by fold change, ALWAYS use shrunken LFC (lfcShrink()).
#
#  RECOMMENDATION: the signed -log10(p) metric used here is the more robust
#  choice at low replication precisely because it sidesteps the raw-LFC inflation
#  problem. If you prefer a fold-change-based ranking, add:
#     library(apeglm)  # or ashr
#     res_shrunk <- lfcShrink(dds, contrast = c("group", g1, g2), type = "ashr")
#     ranks <- sort(setNames(res_shrunk$log2FoldChange, rownames(res_shrunk)),
#                   decreasing = TRUE)
#  (apeglm needs `coef=`; ashr accepts `contrast=`, so ashr is simpler for the
#   arbitrary group-vs-group contrasts used here.)
# =============================================================================
