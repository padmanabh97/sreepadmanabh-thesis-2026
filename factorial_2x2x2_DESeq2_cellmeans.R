# =============================================================================
# FACTORIAL DIFFERENTIAL EXPRESSION - Oxo-Mech RNA-seq, 2 x 2 x 2
# dimension (2D / 3D) x stiffness (0.3 / 1.5 mg/ml) x oxygen (anoxia / normoxia)
# 16 libraries, n = 2 per cell.  One DESeq2 fit, cell-means parameterisation.
# -----------------------------------------------------------------------------
# Replaces factorial_2x2x2_DESeq2_rawLFC.R and its _shrunkenLFC twin.
# Decisions locked with the user (Task 2 plan, 2026-09-23):
#   D1  DE call = padj < 0.05 AND |raw log2FC| > 1 (post-hoc, same as subset scripts)
#   D2  no rep term: 2D and 3D were each harvested + library-prepped in ONE batch,
#       so Rep1/Rep2 are not batches. Model = ~ 0 + group (8 cells, 8 residual df)
#   D3  orientation: anoxia vs normoxia, 0.3 vs 1.5, 3D vs 2D
#       (positive log2FC = higher in the first-named level)
#   D4  within-dimension "component" effects are reported as extra COLUMNS next to
#       each interaction (and main effect) so it is visible which condition carries
#       the effect. They get NO significance call, NO DE list and NO GO:
#       the 2D-only / 3D-only subset scripts remain the official within-dimension results.
#   D5  GSEA rank = sign(log2FC) * -log10(p), identical to the subset GSEA scripts,
#       with p == 0 floored to .Machine$double.xmin so the strongest genes are kept.
#
# EFFECTS TESTED (Wald, numeric contrast vectors over the 8 cell means, log2 scale)
#   Main   Oxygen_A_vs_N          mean of (A - N) over the 4 dim x stiffness cells
#   Main   Stiffness_0.3_vs_1.5   mean of (0.3 - 1.5) over the 4 dim x oxygen cells
#   Main   Dim_3D_vs_2D           mean of (3D - 2D) over the 4 stiffness x oxygen cells
#                                 !! aliased with harvest / library prep / sequencing run
#   2-way  DimxOxygen             (A - N) in 3D  minus (A - N) in 2D       (avg over stiffness)
#   2-way  DimxStiffness          (0.3 - 1.5) in 3D minus (0.3 - 1.5) in 2D (avg over oxygen)
#   2-way  StiffnessxOxygen       (A - N) at 0.3 minus (A - N) at 1.5       (avg over dim)
#   3-way  DimxStiffxOxygen       StiffnessxOxygen in 3D minus the same in 2D
#   Omnibus LRT  ~0+group vs ~dim+stiffness+oxygen (4 df): does the gene deviate from
#                additivity at all?  -> interaction_flag on every row.
#
# The script does not run any analysis you did not ask for; it only writes files
# under OUTPUT_DIR. Every step that can fail reports the failure on the console
# AND in 05_QC/run_log.txt (nothing is hidden by sink()).
# =============================================================================

# ============================================================================
# >>>>>>>>>>>>>>>>>  SETTINGS TO ALTER TO RUN IN YOUR FOLDER  <<<<<<<<<<<<<<<<<
# ============================================================================
# 1. WORKING DIRECTORY: the folder that holds counts.txt and sample_sheet.tsv.
#    Use forward slashes (or doubled backslashes) on Windows. Uncomment ONE:
# setwd("D:/ESoftBio Lab/Omics Analysis/Oxo-Mech 2D vs 3D transcriptome/fragment counts 010726/Subset Analysis 060726")
# setwd("D:\\ESoftBio Lab\\Omics Analysis\\Oxo-Mech 2D vs 3D transcriptome\\fragment counts 010726\\Subset Analysis 060726")
#
# 2. INPUT FILENAMES
counts_file       <- Sys.getenv("COUNTS_FILE",       "counts.txt")
sample_sheet_file <- Sys.getenv("SAMPLE_SHEET_FILE", "sample_sheet.tsv")
#
# 3. OUTPUT FOLDER (new name, so the old DESeq2_Factorial outputs are not overwritten)
output_dir        <- Sys.getenv("OUTPUT_DIR",        "DESeq2_Factorial_cellmeans")
#
# 4. ashr-shrunken LFC column (reported only; never used for calls). Requires the
#    CRAN package 'ashr'. Leave TRUE for the publication run.
USE_ASHR <- as.logical(Sys.getenv("USE_ASHR", "TRUE"))
# ============================================================================

# ---- Analysis constants (match the 2D/3D subset scripts; do not edit casually) ----
ALPHA        <- 0.05      # padj threshold (Wald and LRT); also results(alpha=) for IF
LFC_CUT      <- 1         # |log2FC| > 1  (2-fold), post-hoc
MIN_COUNT    <- 10        # prefilter: >= 10 reads ...
MIN_SAMPLES  <- 2         # ... in >= 2 libraries (= smallest cell size)
GO_MIN_GENES <- 10        # ORA only when a list has >= 10 genes (subset rule)
ORA_P        <- 0.05      # enrichGO pvalueCutoff (BH-adjusted)
ORA_Q        <- 0.2       # enrichGO qvalueCutoff
GSEA_P       <- 0.05      # gseGO pvalueCutoff (BH)
GSEA_MIN_GS  <- 10
GSEA_MAX_GS  <- 500
SEED         <- 1234
FIG_DPI      <- 300

EXPECTED <- list(dim = c("2D", "3D"), stiffness = c("0.3", "1.5"),
                 oxygen = c("anoxia", "normoxia"), rep = c("Rep1", "Rep2"))

# ====================== 0. PREFLIGHT ======================
pkgs <- c("DESeq2", "clusterProfiler", "enrichplot", "org.Mm.eg.db", "AnnotationDbi",
          "ggplot2", "ggrepel", "dplyr", "pheatmap", "RColorBrewer", "openxlsx")
if (USE_ASHR) pkgs <- c(pkgs, "ashr")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs))
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "),
       "\n  CRAN: install.packages(c(...));  Bioconductor: BiocManager::install(c(...))")

suppressPackageStartupMessages({
  library(DESeq2); library(clusterProfiler); library(enrichplot)
  library(org.Mm.eg.db); library(AnnotationDbi)
  library(ggplot2); library(ggrepel); library(dplyr)
  library(pheatmap); library(RColorBrewer); library(openxlsx)
})
options(warn = 1)          # print warnings as they happen, not batched at the end
set.seed(SEED)

DIR <- list(root = output_dir,
            tab  = file.path(output_dir, "01_tables"),
            csv  = file.path(output_dir, "01_tables", "CSV"),
            ora  = file.path(output_dir, "02_GO_ORA"),
            gsea = file.path(output_dir, "03_GSEA"),
            fig  = file.path(output_dir, "04_figures"),
            vol  = file.path(output_dir, "04_figures", "Volcano"),
            heat = file.path(output_dir, "04_figures", "Heatmaps"),
            intp = file.path(output_dir, "04_figures", "Interaction_plots"),
            qc   = file.path(output_dir, "05_QC"),
            meth = file.path(output_dir, "06_methods"))
invisible(lapply(DIR, dir.create, recursive = TRUE, showWarnings = FALSE))

LOG_FILE <- file.path(DIR$qc, "run_log.txt")
cat("", file = LOG_FILE)
# msg() writes to the console AND the log (named msg, not log, so base::log is intact)
msg <- function(...) {
  line <- sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...))
  cat(line); cat(line, file = LOG_FILE, append = TRUE)
}
FACTS <- list()
fact <- function(key, value, note = "") {
  FACTS[[length(FACTS) + 1]] <<- data.frame(key = key, value = as.character(value),
                                            note = note, stringsAsFactors = FALSE)
  msg("  fact  ", key, " = ", value, if (nzchar(note)) paste0("   (", note, ")") else "")
}
msg("Output folder: ", normalizePath(output_dir))

# ====================== 1. LOAD + VALIDATE INPUTS ======================
stopifnot(file.exists(counts_file), file.exists(sample_sheet_file))
first_line <- readLines(counts_file, n = 1)
if (!grepl("featureCounts", first_line))
  stop("counts file line 1 is not a featureCounts header: ", substr(first_line, 1, 80))
fact("counts_file", normalizePath(counts_file))
fact("counts_md5", unname(tools::md5sum(counts_file)))
fact("sample_sheet_md5", unname(tools::md5sum(sample_sheet_file)))

countData <- read.delim(counts_file, row.names = 1, skip = 1, check.names = FALSE)
annot_cols <- c("Chr", "Start", "End", "Strand", "Length")
if (!identical(colnames(countData)[1:5], annot_cols))
  stop("Expected featureCounts annotation columns 2-6 to be ", paste(annot_cols, collapse = ","))
countData <- countData[, -(1:5), drop = FALSE]
colnames(countData) <- sub("\\.bam$", "", basename(colnames(countData)))

sample_info <- read.delim(sample_sheet_file, colClasses = "character")
need_cols <- c("sample", "dim", "stiffness", "oxygen", "rep")
if (!all(need_cols %in% colnames(sample_info)))
  stop("sample sheet lacks columns: ", paste(setdiff(need_cols, colnames(sample_info)), collapse = ","))
sample_info <- sample_info[, need_cols]

# --- exact level sets (a typo would otherwise create a 3rd level silently) ---
for (f in names(EXPECTED)) {
  bad <- setdiff(unique(sample_info[[f]]), EXPECTED[[f]])
  if (length(bad)) stop("Unexpected ", f, " level(s): ", paste(bad, collapse = ","))
}
# --- the sample NAME must agree with its metadata (guards against mixed-up rows) ---
expected_name <- paste(sample_info$dim, sample_info$stiffness, sample_info$oxygen,
                       sample_info$rep, sep = "_")
if (!identical(expected_name, sample_info$sample))
  stop("Sample name does not match its dim/stiffness/oxygen/rep columns for: ",
       paste(sample_info$sample[expected_name != sample_info$sample], collapse = ", "))
if (anyDuplicated(sample_info$sample)) stop("Duplicated sample names in sample sheet")
if (nrow(sample_info) != 16) stop("Expected 16 samples, sample sheet has ", nrow(sample_info))

# --- every sheet sample must be in counts; report (and ignore) extra count columns ---
miss <- setdiff(sample_info$sample, colnames(countData))
if (length(miss)) stop("Samples missing from counts file: ", paste(miss, collapse = ", "))
extra <- setdiff(colnames(countData), sample_info$sample)
fact("count_columns_not_in_sheet", length(extra),
     if (length(extra)) paste(extra, collapse = ";") else "none")

rownames(sample_info) <- sample_info$sample
countData <- as.matrix(countData[, sample_info$sample, drop = FALSE])
if (any(is.na(countData)) || any(countData < 0) || any(countData != round(countData)))
  stop("Counts must be non-negative integers (fractional featureCounts output detected?)")
storage.mode(countData) <- "integer"

if (any(grepl("\\.", rownames(countData)))) {
  msg("Stripping Ensembl version suffixes from gene IDs")
  rownames(countData) <- sub("\\..*$", "", rownames(countData))
}
if (anyDuplicated(rownames(countData))) stop("Gene IDs are not unique")
stopifnot(identical(colnames(countData), rownames(sample_info)))

# ====================== 2. DESIGN ======================
sample_info$dim       <- factor(sample_info$dim,       levels = EXPECTED$dim)
sample_info$stiffness <- factor(sample_info$stiffness, levels = EXPECTED$stiffness)
sample_info$oxygen    <- factor(sample_info$oxygen,    levels = EXPECTED$oxygen)
sample_info$rep       <- factor(sample_info$rep,       levels = EXPECTED$rep)

# 8 cells in a fixed order: dim (slowest) x stiffness x oxygen (fastest)
CELLS <- expand.grid(oxygen = EXPECTED$oxygen, stiffness = EXPECTED$stiffness,
                     dim = EXPECTED$dim, stringsAsFactors = FALSE)[, c("dim", "stiffness", "oxygen")]
CELLS$cell <- paste(CELLS$dim, CELLS$stiffness, CELLS$oxygen, sep = "_")
sample_info$group <- factor(paste(sample_info$dim, sample_info$stiffness,
                                  sample_info$oxygen, sep = "_"), levels = CELLS$cell)

bal <- table(sample_info$dim, sample_info$stiffness, sample_info$oxygen)
if (!all(bal == 2)) { print(bal); stop("Design is not balanced at n = 2 per cell") }
msg("Design balanced: 8 cells x 2 libraries")
write.csv(as.data.frame(bal, responseName = "n_libraries"),
          file.path(DIR$qc, "design_balance.csv"), row.names = FALSE)

dds <- DESeqDataSetFromMatrix(countData = countData, colData = sample_info,
                              design = ~ 0 + group)
fact("genes_in_counts", nrow(dds))
keep <- rowSums(counts(dds) >= MIN_COUNT) >= MIN_SAMPLES
dds  <- dds[keep, ]
fact("genes_after_prefilter", nrow(dds), sprintf(">= %d counts in >= %d of 16 libraries",
                                                 MIN_COUNT, MIN_SAMPLES))

# Per-dimension detection counts (labelling only; NOT used for filtering or calls):
# how many of the 8 2D and 8 3D libraries have >= MIN_COUNT raw reads for each gene.
dim_chr <- as.character(colData(dds)$dim)
n_libs_2D <- rowSums(counts(dds)[, dim_chr == "2D", drop = FALSE] >= MIN_COUNT)
n_libs_3D <- rowSums(counts(dds)[, dim_chr == "3D", drop = FALSE] >= MIN_COUNT)
stopifnot(sum(dim_chr == "2D") == 8, sum(dim_chr == "3D") == 8)
fact("genes_detected_2D_only", sum(n_libs_3D == 0),
     sprintf("tested genes with 0 of 8 3D libraries at >= %d reads", MIN_COUNT))
fact("genes_detected_3D_only", sum(n_libs_2D == 0),
     sprintf("tested genes with 0 of 8 2D libraries at >= %d reads", MIN_COUNT))

dds <- DESeq(dds)
mm  <- model.matrix(design(dds), colData(dds))
resid_df <- ncol(dds) - ncol(mm)
if (resid_df != 8) stop("Residual df should be 8, got ", resid_df)
fact("residual_df", resid_df)

coef_names <- resultsNames(dds)
expected_coef <- make.names(paste0("group", CELLS$cell))
if (!identical(coef_names, expected_coef))
  stop("Unexpected coefficient names: ", paste(coef_names, collapse = ", "))

# ====================== 3. CONTRAST VECTORS ======================
# Coding over the 8 cells: d = +1 for 3D / -1 for 2D, s = +1 for 0.3 / -1 for 1.5,
# o = +1 for anoxia / -1 for normoxia. Each weight vector sums to zero.
d <- ifelse(CELLS$dim == "3D", 1, -1)
s <- ifelse(CELLS$stiffness == "0.3", 1, -1)
o <- ifelse(CELLS$oxygen == "anoxia", 1, -1)
in_ <- function(col, lev) as.numeric(CELLS[[col]] == lev)

EFFECTS <- list(
  Oxygen_A_vs_N        = o / 4,
  Stiffness_0.3_vs_1.5 = s / 4,
  Dim_3D_vs_2D         = d / 4,
  DimxOxygen           = d * o / 2,
  DimxStiffness        = d * s / 2,
  StiffnessxOxygen     = s * o / 2,
  DimxStiffxOxygen     = d * s * o
)
# Component (simple) effects: columns only, no calls (decision D4)
COMPONENTS <- list(
  Oxygen_in_2D       = o * in_("dim", "2D") / 2,
  Oxygen_in_3D       = o * in_("dim", "3D") / 2,
  Stiffness_in_2D    = s * in_("dim", "2D") / 2,
  Stiffness_in_3D    = s * in_("dim", "3D") / 2,
  Oxygen_at_0.3      = o * in_("stiffness", "0.3") / 2,
  Oxygen_at_1.5      = o * in_("stiffness", "1.5") / 2,
  Dim_under_anoxia   = d * in_("oxygen", "anoxia") / 2,
  Dim_under_normoxia = d * in_("oxygen", "normoxia") / 2,
  Oxygen_2D_0.3      = o * in_("dim", "2D") * in_("stiffness", "0.3"),
  Oxygen_2D_1.5      = o * in_("dim", "2D") * in_("stiffness", "1.5"),
  Oxygen_3D_0.3      = o * in_("dim", "3D") * in_("stiffness", "0.3"),
  Oxygen_3D_1.5      = o * in_("dim", "3D") * in_("stiffness", "1.5")
)
# which components sit next to which effect in its per-effect table
EFFECT_COMPONENTS <- list(
  Oxygen_A_vs_N        = c("Oxygen_in_2D", "Oxygen_in_3D"),
  Stiffness_0.3_vs_1.5 = c("Stiffness_in_2D", "Stiffness_in_3D"),
  Dim_3D_vs_2D         = c("Dim_under_anoxia", "Dim_under_normoxia"),
  DimxOxygen           = c("Oxygen_in_2D", "Oxygen_in_3D"),
  DimxStiffness        = c("Stiffness_in_2D", "Stiffness_in_3D"),
  StiffnessxOxygen     = c("Oxygen_at_0.3", "Oxygen_at_1.5"),
  DimxStiffxOxygen     = c("Oxygen_2D_0.3", "Oxygen_2D_1.5", "Oxygen_3D_0.3", "Oxygen_3D_1.5")
)
EFFECT_INFO <- data.frame(
  effect = names(EFFECTS),
  tier   = c("main", "main", "main", "2-way", "2-way", "2-way", "3-way"),
  question = c(
    "Anoxia vs normoxia, averaged over dimension and stiffness",
    "0.3 vs 1.5 mg/ml, averaged over dimension and oxygen",
    "3D vs 2D, averaged over stiffness and oxygen (ALIASED with harvest/library prep/sequencing run)",
    "Is the anoxia response different in 3D than in 2D?",
    "Is the stiffness (0.3 vs 1.5) response different in 3D than in 2D?",
    "Is the anoxia response different at 0.3 than at 1.5 mg/ml?",
    "Does the stiffness-dependence of the anoxia response differ between 3D and 2D?"),
  positive_log2FC_means = c(
    "higher in anoxia", "higher at 0.3 mg/ml", "higher in 3D",
    "anoxia response (A - N) is more positive in 3D than in 2D",
    "stiffness response (0.3 - 1.5) is more positive in 3D than in 2D",
    "anoxia response (A - N) is more positive at 0.3 than at 1.5",
    "(A-N at 0.3 minus A-N at 1.5) is more positive in 3D than in 2D"),
  components_shown = vapply(EFFECT_COMPONENTS, paste, character(1), collapse = "; "),
  stringsAsFactors = FALSE)

# --- algebraic checks on the weights (stop on failure) ---
for (nm in c(names(EFFECTS), names(COMPONENTS))) {
  w <- c(EFFECTS, COMPONENTS)[[nm]]
  if (abs(sum(w)) > 1e-12) stop("Contrast ", nm, " does not sum to 0")
}
chk <- function(a, b, what) if (max(abs(a - b)) > 1e-12) stop("Contrast algebra failed: ", what)
C <- COMPONENTS; E <- EFFECTS
chk(E$Oxygen_A_vs_N, (C$Oxygen_2D_0.3 + C$Oxygen_2D_1.5 + C$Oxygen_3D_0.3 + C$Oxygen_3D_1.5) / 4,
    "oxygen main = mean of 4 cell-level oxygen effects")
chk(E$Oxygen_A_vs_N, (C$Oxygen_in_2D + C$Oxygen_in_3D) / 2, "oxygen main = mean(2D, 3D)")
chk(E$DimxOxygen, C$Oxygen_in_3D - C$Oxygen_in_2D, "DimxOxygen = Oxygen_in_3D - Oxygen_in_2D")
chk(E$DimxStiffness, C$Stiffness_in_3D - C$Stiffness_in_2D, "DimxStiffness")
chk(E$StiffnessxOxygen, C$Oxygen_at_0.3 - C$Oxygen_at_1.5, "StiffnessxOxygen")
chk(E$DimxStiffxOxygen, (C$Oxygen_3D_0.3 - C$Oxygen_3D_1.5) - (C$Oxygen_2D_0.3 - C$Oxygen_2D_1.5),
    "3-way")
msg("Contrast algebra checks passed")

weights_tab <- cbind(CELLS[, c("cell", "dim", "stiffness", "oxygen")],
                     as.data.frame(EFFECTS), as.data.frame(COMPONENTS))
write.csv(weights_tab, file.path(DIR$tab, "contrast_weights.csv"), row.names = FALSE)

# ====================== 4. ANNOTATION ======================
gene_ids  <- rownames(dds)
gene_name <- mapIds(org.Mm.eg.db, keys = gene_ids, column = "SYMBOL",
                    keytype = "ENSEMBL", multiVals = "first")
gene_name <- ifelse(is.na(gene_name), gene_ids, gene_name)
names(gene_name) <- gene_ids
fact("genes_without_symbol", sum(gene_name == gene_ids))

# ====================== 5. QC ======================
vsd <- vst(dds, blind = TRUE)
write.csv(data.frame(sample = colnames(dds), sizeFactor = sizeFactors(dds),
                     library_size = colSums(counts(dds))),
          file.path(DIR$qc, "size_factors.csv"), row.names = FALSE)

png(file.path(DIR$qc, "dispersion_estimates.png"), width = 7, height = 6, units = "in", res = FIG_DPI)
plotDispEsts(dds, main = "Dispersion estimates (~ 0 + group, 16 libraries)")
invisible(dev.off())

pca <- plotPCA(vsd, intgroup = c("dim", "stiffness", "oxygen"), ntop = 500, returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"))
pca$stiffness_lab <- paste0(pca$stiffness, " mg/ml")
p_pca <- ggplot(pca, aes(PC1, PC2, colour = dim, shape = oxygen)) +
  geom_point(aes(size = stiffness_lab), alpha = 0.9) +
  scale_size_manual(values = c("0.3 mg/ml" = 3, "1.5 mg/ml" = 5), name = "stiffness") +
  geom_text_repel(aes(label = name), size = 2.8, max.overlaps = 50, show.legend = FALSE) +
  scale_colour_brewer(palette = "Set1") +
  labs(x = sprintf("PC1 (%d%% variance)", pv[1]), y = sprintf("PC2 (%d%% variance)", pv[2]),
       title = "PCA, blind VST, top 500 genes, 16 libraries") +
  theme_bw(base_size = 12)
ggsave(file.path(DIR$fig, "PCA_all16.pdf"), plot = p_pca, width = 9, height = 7)
ggsave(file.path(DIR$fig, "PCA_all16.png"), plot = p_pca, width = 9, height = 7, dpi = FIG_DPI)

sd_mat <- as.matrix(dist(t(assay(vsd))))
ann_col <- as.data.frame(colData(dds)[, c("dim", "stiffness", "oxygen")])
png(file.path(DIR$qc, "sample_distance_heatmap.png"), width = 9, height = 8, units = "in", res = FIG_DPI)
pheatmap(sd_mat, annotation_col = ann_col, annotation_row = ann_col,
         color = colorRampPalette(rev(brewer.pal(9, "Blues")))(100),
         main = "Sample-to-sample Euclidean distance (blind VST)")
invisible(dev.off())

# ====================== 6. OMNIBUS INTERACTION LRT ======================
# DESeq2 cannot nest formulas across a ~0+group full model, so both models are given
# as design matrices (full = the 8 cell means, reduced = additive, 4 columns).
mm_red  <- model.matrix(~ dim + stiffness + oxygen, colData(dds))
if (qr(mm_red)$rank != 4 || qr(cbind(mm, mm_red))$rank != qr(mm)$rank)
  stop("Reduced LRT model is not nested in the full model")
dds_lrt <- nbinomLRT(dds, full = mm, reduced = mm_red)
res_lrt <- results(dds_lrt, alpha = ALPHA)
lrt_df  <- data.frame(gene_id = rownames(res_lrt), LRT_stat = res_lrt$stat,
                      LRT_pvalue = res_lrt$pvalue, LRT_padj = res_lrt$padj)
lrt_df$interaction_flag <- !is.na(lrt_df$LRT_padj) & lrt_df$LRT_padj < ALPHA
fact("LRT_df", ncol(mm) - ncol(mm_red))
fact("LRT_genes_tested", sum(!is.na(lrt_df$LRT_padj)))
fact("LRT_genes_interacting_padj<0.05", sum(lrt_df$interaction_flag),
     "any departure from additivity; n = 2/cell, so a non-flag is NOT evidence of additivity")
png(file.path(DIR$qc, "pvalue_hist_LRT.png"), width = 6, height = 4.5, units = "in", res = FIG_DPI)
hist(lrt_df$LRT_pvalue, breaks = 50, main = "Omnibus interaction LRT p-values", xlab = "p")
invisible(dev.off())

# ====================== 7. WALD RESULTS PER EFFECT + COMPONENTS ======================
get_res <- function(w) results(dds, contrast = w, alpha = ALPHA)

comp_lfc <- lapply(COMPONENTS, function(w) {
  r <- get_res(w); data.frame(LFC = r$log2FoldChange, lfcSE = r$lfcSE)
})

res_list <- list(); summary_rows <- list()
for (eff in names(EFFECTS)) {
  msg("=== Effect: ", eff, " ===")
  res <- get_res(EFFECTS[[eff]])
  df  <- data.frame(gene_id = rownames(res), gene_name = gene_name[rownames(res)],
                    baseMean = res$baseMean, log2FoldChange = res$log2FoldChange,
                    lfcSE = res$lfcSE, stat = res$stat, pvalue = res$pvalue, padj = res$padj,
                    stringsAsFactors = FALSE)
  if (USE_ASHR) {
    sh <- lfcShrink(dds, contrast = EFFECTS[[eff]], res = res, type = "ashr", quiet = TRUE)
    df$log2FoldChange_ashr <- sh$log2FoldChange
  }
  sig <- !is.na(df$padj) & df$padj < ALPHA
  df$regulation <- "NS"
  df$regulation[sig & df$log2FoldChange >  LFC_CUT] <- "Up"
  df$regulation[sig & df$log2FoldChange < -LFC_CUT] <- "Down"
  for (cn in EFFECT_COMPONENTS[[eff]])
    df[[paste0(cn, "_log2FC")]] <- comp_lfc[[cn]]$LFC
  df$n_libs_ge10_2D <- unname(n_libs_2D[df$gene_id])
  df$n_libs_ge10_3D <- unname(n_libs_3D[df$gene_id])
  df <- merge(df, lrt_df[, c("gene_id", "LRT_padj", "interaction_flag")], by = "gene_id", sort = FALSE)
  df <- df[order(df$padj, df$pvalue, na.last = TRUE), ]
  write.csv(df, file.path(DIR$csv, paste0(eff, "_results.csv")), row.names = FALSE)
  res_list[[eff]] <- df

  png(file.path(DIR$qc, paste0("pvalue_hist_", eff, ".png")), width = 6, height = 4.5,
      units = "in", res = FIG_DPI)
  hist(df$pvalue, breaks = 50, main = paste("Wald p-values:", eff), xlab = "p")
  invisible(dev.off())

  n_up <- sum(df$regulation == "Up"); n_dn <- sum(df$regulation == "Down")
  summary_rows[[eff]] <- data.frame(
    effect = eff, tier = EFFECT_INFO$tier[EFFECT_INFO$effect == eff],
    n_tested = sum(!is.na(df$padj)), n_padj_NA = sum(is.na(df$padj)),
    independent_filter_threshold = signif(metadata(res)$filterThreshold, 4),
    n_padj_only = sum(sig), Up = n_up, Down = n_dn, Total_DE = n_up + n_dn,
    DE_with_interaction_flag = sum(df$regulation != "NS" & df$interaction_flag),
    stringsAsFactors = FALSE)
  msg(sprintf("  tested %d | padj<0.05: %d | Up %d  Down %d", sum(!is.na(df$padj)),
              sum(sig), n_up, n_dn))
}

# --- check on the fitted object: main effect = mean of its cell-level components ---
# DESeq2 results() sets a contrast's log2FC to 0 (stat 0, p 1) when every library in the cells carrying
# non-zero weight has zero reads ("contrastAllZero"). A 2-cell component can hit this for a gene expressed
# in only one or two cells (the prefilter allows a single cell); the 8-cell main effect never does.
# The identity therefore holds exactly on the coefficients for all genes, and on the reported LFCs for
# every gene not overridden in any component involved. Both are checked; overridden genes are counted.
beta   <- coef(dds)                                   # log2 scale, columns = cells in model-matrix order
stopifnot(identical(colnames(beta), colnames(mm)))
cnt    <- counts(dds)
zeroed <- function(w) {                               # genes DESeq2 overrides for contrast w
  libs <- as.vector(mm %*% as.numeric(w != 0)) != 0
  rowSums(cnt[, libs, drop = FALSE]) == 0
}
chk_fit <- function(eff_w, comp_w, comp_sum, eff_res, what) {
  d_coef <- max(abs(beta %*% eff_w - comp_sum(lapply(comp_w, function(w) beta %*% w))), na.rm = TRUE)
  if (d_coef > 1e-6) stop("Fitted-contrast check failed on coefficients (", what, "): max |diff| = ", d_coef)
  over <- Reduce(`|`, lapply(comp_w, zeroed))
  a <- eff_res$log2FoldChange[!over]
  b <- comp_sum(lapply(names(comp_w), function(cn) comp_lfc[[cn]]$LFC))[!over]
  d_rep <- max(abs(a - b), na.rm = TRUE)
  if (d_rep > 1e-6) stop("Fitted-contrast check failed on reported LFCs (", what, "): max |diff| = ", d_rep)
  fact(paste0("genes_component_LFC_set_to_0_by_DESeq2 [", what, "]"), sum(over),
       "zero reads in every library of both cells of >= 1 component; excluded from the reported-LFC check only")
}
ox  <- res_list$Oxygen_A_vs_N[match(rownames(dds), res_list$Oxygen_A_vs_N$gene_id), ]
dxo <- res_list$DimxOxygen[match(rownames(dds), res_list$DimxOxygen$gene_id), ]
chk_fit(EFFECTS$Oxygen_A_vs_N, COMPONENTS[c("Oxygen_2D_0.3", "Oxygen_2D_1.5", "Oxygen_3D_0.3", "Oxygen_3D_1.5")],
        function(x) Reduce(`+`, x) / 4, ox, "Oxygen main = mean of 4 cell-level effects")
chk_fit(EFFECTS$DimxOxygen, COMPONENTS[c("Oxygen_in_3D", "Oxygen_in_2D")],
        function(x) x[[1]] - x[[2]], dxo, "DimxOxygen = Oxygen_in_3D - Oxygen_in_2D")
msg("Fitted-contrast checks passed")
fact("cooks_filtering_applied", "no",
     "DESeq2 flags Cook's outliers only in cells with >= 3 replicates; here every cell has 2")

de_summary <- bind_rows(summary_rows)

# ====================== 8. FIGURES: VOLCANO, HEATMAP, INTERACTION PLOTS ======================
norm_counts <- counts(dds, normalized = TRUE)
vst_mat     <- assay(vsd)

for (eff in names(EFFECTS)) {
  df <- res_list[[eff]]; df <- df[!is.na(df$padj), ]
  df$neg_log10_padj <- -log10(pmax(df$padj, .Machine$double.xmin))
  top_lab <- head(df[df$regulation != "NS", ], 15)
  p_v <- ggplot(df, aes(log2FoldChange, neg_log10_padj)) +
    geom_point(aes(colour = regulation), size = 1, alpha = 0.7) +
    scale_colour_manual(values = c(Up = "red", Down = "blue", NS = "grey70")) +
    geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = -log10(ALPHA), linetype = "dashed", linewidth = 0.4) +
    geom_text_repel(data = top_lab, aes(label = gene_name), size = 2.8, max.overlaps = 30) +
    labs(title = paste("Volcano -", eff),
         subtitle = EFFECT_INFO$positive_log2FC_means[EFFECT_INFO$effect == eff],
         x = "log2 fold change (raw)", y = "-log10(padj)") +
    theme_minimal(base_size = 12)
  ggsave(file.path(DIR$vol, paste0(eff, "_volcano.png")), plot = p_v, width = 9, height = 7, dpi = FIG_DPI)

  de_ids <- head(df$gene_id[df$regulation != "NS"], 50)
  if (length(de_ids) >= 2) {
    z <- t(scale(t(vst_mat[de_ids, , drop = FALSE])))
    z <- z[, order(colData(dds)$group)]
    rownames(z) <- gene_name[de_ids]
    png(file.path(DIR$heat, paste0(eff, "_top", length(de_ids), "_heatmap.png")),
        width = 9, height = max(5, 0.16 * length(de_ids) + 3), units = "in", res = FIG_DPI)
    pheatmap(z, annotation_col = ann_col, cluster_cols = FALSE, show_colnames = TRUE,
             fontsize_row = 7, main = paste(eff, "- top", length(de_ids), "DE genes (row z-score, VST)"))
    invisible(dev.off())
  } else msg("  heatmap skipped for ", eff, ": fewer than 2 DE genes")
}

# interaction plots: top 6 genes by padj for each 2-way / 3-way effect
cell_df <- as.data.frame(colData(dds)[, c("dim", "stiffness", "oxygen")])
IPLOT <- list(DimxOxygen = c(x = "oxygen", col = "dim", facet = "stiffness"),
              DimxStiffness = c(x = "stiffness", col = "dim", facet = "oxygen"),
              StiffnessxOxygen = c(x = "oxygen", col = "stiffness", facet = "dim"),
              DimxStiffxOxygen = c(x = "oxygen", col = "stiffness", facet = "dim"))
for (eff in names(IPLOT)) {
  df <- res_list[[eff]]; ids <- head(df$gene_id[!is.na(df$padj)], 6)
  long <- do.call(rbind, lapply(ids, function(g)
    data.frame(cell_df, gene = gene_name[g], log2_norm = log2(norm_counts[g, ] + 1))))
  long$gene <- factor(long$gene, levels = unique(gene_name[ids]))
  a <- IPLOT[[eff]]
  p_i <- ggplot(long, aes(.data[[a["x"]]], log2_norm, colour = .data[[a["col"]]],
                          group = .data[[a["col"]]])) +
    geom_point(position = position_dodge(width = 0.2), size = 2) +
    stat_summary(fun = mean, geom = "line", position = position_dodge(width = 0.2)) +
    facet_grid(as.formula(paste("gene ~", a["facet"])), scales = "free_y") +
    scale_colour_brewer(palette = "Set1") +
    labs(title = paste(eff, "- 6 genes with smallest padj"),
         subtitle = "points = libraries; lines join cell means; non-parallel lines = interaction",
         y = "log2(normalised count + 1)") +
    theme_bw(base_size = 10)
  ggsave(file.path(DIR$intp, paste0(eff, "_top6_interaction.png")), plot = p_i,
         width = 7, height = 12, dpi = FIG_DPI)
}

# ====================== 9. GO ORA (Up / Down per effect) ======================
ora_rows <- list()
for (eff in names(EFFECTS)) {
  df <- res_list[[eff]]
  universe <- df$gene_id[!is.na(df$padj)]
  for (dir_ in c("Up", "Down")) {
    genes <- df$gene_id[df$regulation == dir_]
    tag <- paste0(eff, "_", dir_)
    status <- "run"; n_terms <- 0L
    if (length(genes) < GO_MIN_GENES) {
      status <- sprintf("not run (%d genes < %d)", length(genes), GO_MIN_GENES)
    } else {
      ego <- tryCatch(enrichGO(gene = genes, universe = universe, OrgDb = org.Mm.eg.db,
                               keyType = "ENSEMBL", ont = "BP", pAdjustMethod = "BH",
                               pvalueCutoff = ORA_P, qvalueCutoff = ORA_Q, readable = TRUE),
                      error = function(e) { status <<- paste("ERROR:", conditionMessage(e)); NULL })
      if (!is.null(ego)) {
        tab <- as.data.frame(ego); n_terms <- nrow(tab)
        write.csv(tab, file.path(DIR$ora, paste0(tag, "_GO_BP.csv")), row.names = FALSE)
        if (n_terms > 0) {
          p_o <- dotplot(ego, showCategory = 15) + labs(title = paste(tag, "- GO BP ORA")) +
            theme_minimal()
          ggsave(file.path(DIR$ora, paste0(tag, "_GO_BP.png")), plot = p_o,
                 width = 10, height = 8, dpi = FIG_DPI)
        }
      }
    }
    if (grepl("^ERROR", status)) msg("  !! ORA ", tag, ": ", status) else
      msg(sprintf("  ORA %s: %d genes | %s | %d terms", tag, length(genes), status, n_terms))
    ora_rows[[tag]] <- data.frame(effect = eff, direction = dir_, n_genes = length(genes),
                                  universe = length(universe), status = status,
                                  n_terms = n_terms, stringsAsFactors = FALSE)
  }
}
ora_summary <- bind_rows(ora_rows)

# ====================== 10. GSEA per effect ======================
gsea_rows <- list()
for (eff in names(EFFECTS)) {
  df <- res_list[[eff]]
  rk <- df[!is.na(df$pvalue) & !is.na(df$log2FoldChange), ]
  n_p0 <- sum(rk$pvalue == 0)
  ranks <- sign(rk$log2FoldChange) * -log10(pmax(rk$pvalue, .Machine$double.xmin))
  names(ranks) <- rk$gene_id
  stopifnot(!anyDuplicated(names(ranks)), all(is.finite(ranks)))
  ranks <- sort(ranks, decreasing = TRUE)
  status <- "run"; n_sets <- 0L
  set.seed(SEED)
  gs <- tryCatch(gseGO(geneList = ranks, OrgDb = org.Mm.eg.db, keyType = "ENSEMBL",
                       ont = "BP", minGSSize = GSEA_MIN_GS, maxGSSize = GSEA_MAX_GS,
                       pvalueCutoff = GSEA_P, pAdjustMethod = "BH", verbose = FALSE, seed = TRUE),
                 error = function(e) { status <<- paste("ERROR:", conditionMessage(e)); NULL })
  if (!is.null(gs)) {
    gs_r <- setReadable(gs, OrgDb = org.Mm.eg.db, keyType = "ENSEMBL")
    tab <- as.data.frame(gs_r); n_sets <- nrow(tab)
    write.csv(tab, file.path(DIR$gsea, paste0(eff, "_GSEA_GO_BP.csv")), row.names = FALSE)
    if (n_sets > 0) {
      p_g <- dotplot(gs, showCategory = 10, split = ".sign") + facet_grid(. ~ .sign) +
        ggtitle(paste(eff, "- GSEA GO BP")) +
        theme(axis.text.y = element_text(size = 8), strip.text = element_text(size = 20))
      ggsave(file.path(DIR$gsea, paste0(eff, "_GSEA_dotplot.png")), plot = p_g,
             width = 12, height = 8, dpi = FIG_DPI)
    }
    msg(sprintf("  GSEA %s: %d ranked (%d p=0 floored) | %d sets", eff, length(ranks), n_p0, n_sets))
  } else msg("  !! GSEA ", eff, ": ", status)
  gsea_rows[[eff]] <- data.frame(effect = eff, n_ranked = length(ranks), n_pvalue_zero_floored = n_p0,
                                 status = status, n_sets = n_sets, stringsAsFactors = FALSE)
}
gsea_summary <- bind_rows(gsea_rows)

# ====================== 11. TABLES: summary, wide sheet, workbook ======================
de_summary <- de_summary %>%
  left_join(ora_summary %>% filter(direction == "Up") %>% select(effect, ORA_Up_terms = n_terms), by = "effect") %>%
  left_join(ora_summary %>% filter(direction == "Down") %>% select(effect, ORA_Down_terms = n_terms), by = "effect") %>%
  left_join(gsea_summary %>% select(effect, GSEA_sets = n_sets), by = "effect")
write.csv(de_summary, file.path(DIR$root, "DE_Summary_Factorial.csv"), row.names = FALSE)

wide <- data.frame(gene_id = rownames(dds), gene_name = gene_name[rownames(dds)],
                   baseMean = mcols(dds)$baseMean,
                   n_libs_ge10_2D = unname(n_libs_2D), n_libs_ge10_3D = unname(n_libs_3D),
                   stringsAsFactors = FALSE)
wide <- merge(wide, lrt_df, by = "gene_id", sort = FALSE)
for (eff in names(EFFECTS)) {
  df <- res_list[[eff]]
  keepc <- intersect(c("log2FoldChange", "log2FoldChange_ashr", "lfcSE", "pvalue", "padj", "regulation"),
                     colnames(df))
  sub <- df[, c("gene_id", keepc)]; colnames(sub)[-1] <- paste0(eff, "__", keepc)
  wide <- merge(wide, sub, by = "gene_id", sort = FALSE)
}
for (cn in names(COMPONENTS)) {
  wide[[paste0("component__", cn, "__log2FC")]] <- comp_lfc[[cn]]$LFC[match(wide$gene_id, rownames(dds))]
  wide[[paste0("component__", cn, "__lfcSE")]]  <- comp_lfc[[cn]]$lfcSE[match(wide$gene_id, rownames(dds))]
}
wide <- wide[order(wide$LRT_padj, na.last = TRUE), ]

cell_means <- sapply(CELLS$cell, function(cl) rowMeans(norm_counts[, colData(dds)$group == cl, drop = FALSE]))
cell_means <- data.frame(gene_id = rownames(norm_counts), gene_name = gene_name[rownames(norm_counts)],
                         round(cell_means, 2), check.names = FALSE)

fact("ashr_shrinkage_column", if (USE_ASHR) "yes (reported only, not used for calls)" else "NOT computed (USE_ASHR = FALSE)")
facts_df <- bind_rows(FACTS)
write.csv(facts_df, file.path(DIR$root, "00_RUN_SUMMARY.csv"), row.names = FALSE)

wb <- createWorkbook()
addWorksheet(wb, "README")
writeData(wb, "README", data.frame(line = c(
  "Factorial 2x2x2 DESeq2 analysis, one fit on 16 libraries, design ~ 0 + group (8 cells, n = 2).",
  "DE call: padj < 0.05 and |raw log2FC| > 1. log2FoldChange_ashr is shown for reference only.",
  "Orientation: anoxia vs normoxia, 0.3 vs 1.5 mg/ml, 3D vs 2D; positive = higher in the first-named level.",
  "Interaction effects are differences of log2 fold changes; see effect_definitions for what 'Up' means.",
  "component__ columns are within-context effects from the same fit. They carry NO significance call;",
  "  the 2D-only and 3D-only subset analyses are the official within-dimension results.",
  "n_libs_ge10_2D / n_libs_ge10_3D = number of the 8 2D / 8 3D libraries with >= 10 raw reads (labelling only).",
  "  A 0 in either column marks expression restricted to one dimension; interpret dim-interactions for such genes with care.",
  "interaction_flag = omnibus LRT (~0+group vs ~dim+stiffness+oxygen) padj < 0.05.",
  "Dim_3D_vs_2D is aliased with harvest, library prep and sequencing run and cannot be separated from them.")))
sheets <- list(effect_definitions = EFFECT_INFO, contrast_weights = weights_tab,
               DE_summary = de_summary, ORA_summary = ora_summary, GSEA_summary = gsea_summary,
               run_facts = facts_df, all_effects_wide = wide, cell_means_normcounts = cell_means)
for (sn in names(sheets)) { addWorksheet(wb, sn); writeData(wb, sn, sheets[[sn]]) }
for (eff in names(EFFECTS)) { addWorksheet(wb, eff); writeData(wb, eff, res_list[[eff]]) }
saveWorkbook(wb, file.path(DIR$tab, "Factorial_results.xlsx"), overwrite = TRUE)

# ====================== 12. METHODS TEXT (from computed facts) ======================
F <- setNames(facts_df$value, facts_df$key)
meth <- c(
  "Factorial differential expression (RNA-seq).",
  sprintf(paste("Gene-level counts (featureCounts) from 16 libraries (2 dimensions x 2 collagen",
                "stiffnesses x 2 oxygen conditions, n = 2 per condition) were analysed in a single",
                "DESeq2 (v%s) model with one coefficient per condition (design ~ 0 + group; 8",
                "residual degrees of freedom). Genes with at least %d reads in at least %d of the 16",
                "libraries (the size of the smallest experimental group) were retained (%s of %s); the",
                "filter was applied once, before model fitting and irrespective of condition, so that it",
                "is independent of the test statistic under the null hypothesis (Bourgon et al., 2010).",
                "For each gene, the number of 2D and of 3D libraries with at least %d reads is reported",
                "to identify expression restricted to one culture dimension (%s genes detected only in",
                "2D, %s only in 3D); these counts were not used for filtering or calls."),
          as.character(packageVersion("DESeq2")), MIN_COUNT, MIN_SAMPLES,
          F["genes_after_prefilter"], F["genes_in_counts"], MIN_COUNT,
          F["genes_detected_2D_only"], F["genes_detected_3D_only"]),
  paste("Replicate number was not included as a covariate because all 2D libraries, and",
        "separately all 3D libraries, were harvested and prepared in a single batch."),
  paste("Main effects were estimated as the average log2 fold change across the four",
        "combinations of the other two factors (anoxia vs normoxia; 0.3 vs 1.5 mg/ml; 3D vs 2D).",
        "Two-way interactions were estimated as the difference between the simple effects of",
        "one factor at the two levels of another, averaged over the third, and the three-way",
        "interaction as the corresponding difference of differences. All effects were tested",
        "with Wald tests on linear contrasts of the condition coefficients."),
  sprintf(paste("Departure from additivity was tested per gene with a likelihood-ratio test",
                "against the additive model ~ dim + stiffness + oxygen (4 df); %s of %s genes had",
                "padj < %.2f. With two replicates per condition this test has limited power, and a",
                "non-significant result does not establish additivity."),
          F["LRT_genes_interacting_padj<0.05"], F["LRT_genes_tested"], ALPHA),
  sprintf(paste("DESeq2 independent filtering on the mean of normalized counts was applied to",
                "each effect (alpha = %.2f), p-values were adjusted by Benjamini-Hochberg, and genes with padj < %.2f and |log2 fold change| > %g were called",
                "differentially expressed. Cook's-distance outlier filtering does not apply at",
                "two replicates per condition. %s"), ALPHA, ALPHA, LFC_CUT,
          if (USE_ASHR) "Adaptive-shrinkage (ashr) log2 fold changes are reported alongside, but were not used for calls." else ""),
  sprintf(paste("GO Biological Process over-representation (clusterProfiler enrichGO, Ensembl",
                "keys, universe = genes tested for that effect, BH-adjusted p < %.2f, q < %.2f) was",
                "run separately for up- and down-regulated lists of >= %d genes. GSEA (gseGO, BP,",
                "gene sets of %d-%d genes, BH-adjusted p < %.2f, fixed seed) used all tested genes",
                "ranked by sign(log2 fold change) x -log10(p)."),
          ORA_P, ORA_Q, GO_MIN_GENES, GSEA_MIN_GS, GSEA_MAX_GS, GSEA_P),
  paste("Because 2D and 3D samples were processed in separate batches, the 3D vs 2D main effect",
        "includes any harvest, library-preparation and sequencing-run differences. Interaction",
        "effects compare responses within each batch and are not affected by batch shifts",
        "common to all conditions in a batch."))
writeLines(meth, file.path(DIR$meth, "Methods_factorial_RNAseq.txt"))

writeLines(capture.output(sessionInfo()), file.path(DIR$qc, "sessionInfo.txt"))
print(de_summary)
msg("=== FACTORIAL ANALYSIS COMPLETE === outputs in ", normalizePath(output_dir))
