## =============================================================================
## PID430 PROTEOMICS - CONSOLIDATED END-TO-END PIPELINE
## =============================================================================
##
## PURPOSE
##   Complete re-analysis of the 3D oxo-mechanical label-free proteomics dataset,
##   replacing the service provider's statistical analysis. Written in response to
##   Reviewer comment 1 (quantitative rigour / statistical thresholds).
##
## WHAT THE USER MUST SUPPLY
##   Exactly one thing: the path to the vendor's ANALYSIS workbook, in section 1
##   below. Everything else - group assignment, contrasts, thresholds, output
##   folder layout - is fixed in this script by design, so that a run is
##   reproducible and cannot be silently reconfigured.
##
## METHOD SUMMARY (full prose in PROTEOMICS_METHODS.txt)
##   read xlsx -> unique peptide counts from PEPTIDES sheet -> complete-case and
##   >=2 unique peptide filters -> log2 -> equal-median normalization -> ONE limma
##   fit on all 12 samples (8 residual df) -> DEqMS count-dependent variance
##   moderation -> 6 contrasts -> BH within contrast -> adj p < 0.05 and a fold
##   change of at least 1.5x in either direction (|log2FC| >= log2(1.5), computed
##   from FC_PRIMARY, never a rounded literal) -> GO/KEGG/Reactome
##   over-representation on that SAME significant list, against the TESTED
##   proteome as background.
##
## ORIENTATION - READ THIS BEFORE INTERPRETING ANY OUTPUT
##   There is NO global reference level in this contrast set, and any sentence
##   that names one is wrong. The rule is per contrast and has no exceptions:
##
##     a POSITIVE log2FC means HIGHER IN THAT CONTRAST'S NUMERATOR, which is the
##     group named FIRST in the contrast name.
##
##   So "0-3_A_vs_1-5_A" is positive when higher in 0.3 mg/ml anoxia, and
##   "0-3_N_vs_1-5_A" is positive when higher in 0.3 mg/ml normoxia. All six
##   contrasts are oriented exactly as in the RNA-seq pairwise analysis (and the
##   vendor t-test workbook): the diagonal was previously 1-5_A_vs_0-3_N, the
##   reverse of RNA, and was re-oriented to 0-3_N_vs_1-5_A so that the two layers
##   share one direction per contrast. The interaction contrast is a
##   difference of differences and has its own sign convention, stated where it
##   is defined - it is not a member of this family.
##
##   Every results table carries an explicit `higher_in_if_positive` column
##   naming the numerator group, so the direction can be read off the table
##   without consulting this header. If you ever merge these results with
##   factorial-model output from the RNA-seq or ATAC analyses, check the sign
##   convention on each side rather than assuming a global flip.
##
## AUTHORSHIP OF DECISIONS
##   Group naming, orientation, enrichment engine and the missing-package policy
##   were chosen by the project owner. Points where this script departs from the
##   two earlier R scripts are marked "DEPARTURE" with the reason.
##
## =============================================================================


## =============================================================================
## 1. USER INPUT - THE ONLY LINE YOU NEED TO EDIT
## =============================================================================

INPUT_XLSX <- "DR_TAPOMAY_PID430_PROTEOMICS_ANALYSIS.xlsx"

## Optional: where results are written. Defaults to a timestamped folder created
## next to the input file, so repeated runs never silently overwrite each other.
OUTPUT_DIR <- NULL   # NULL = auto (<input folder>/PID430_proteomics_results_<date>)


## =============================================================================
## 2. FIXED CONFIGURATION - DO NOT EDIT WITHOUT UPDATING THE METHODS TEXT
## =============================================================================
## Every value here appears in PROTEOMICS_METHODS.txt. Changing one without
## changing the other makes the manuscript inconsistent with the analysis.

## --- 2a. Design -------------------------------------------------------------
## THREE parallel naming schemes are needed and they are NOT interchangeable:
##   GRP_INTERNAL : valid R identifiers. Required because the display names
##                  ("0-3 A") contain a hyphen and a space and cannot be parsed
##                  by makeContrasts or used as data-frame column names.
##   GRP_DISPLAY  : what the project owner wants shown in tables and figures.
##   GRP_SLUG     : filesystem-safe, used in file and folder names.
## The mapping is written to the output as a sheet so it is auditable.
GRP_INTERNAL <- c("G0_3_A", "G0_3_N", "G1_5_A", "G1_5_N")
GRP_DISPLAY  <- c("0-3 A",  "0-3 N",  "1-5 A",  "1-5 N")
GRP_SLUG     <- c("0-3_A",  "0-3_N",  "1-5_A",  "1-5_N")

## How each group maps onto the vendor's abundance column headers.
## CRITICAL AND NON-OBVIOUS: only two of the vendor's four labels carry a
## stiffness prefix. The two UNPREFIXED labels (ANOXIA, NORMOXIA) are the
## 0.3 mg/ml SOFT condition. This is an undocumented convention in the vendor
## file; it is asserted explicitly here and disclosed in the Methods.
GRP_VENDOR   <- c("ANOXIA", "NORMOXIA", "STIFF-ANOXIA", "STIFF-NORMOXIA")

## Biological meaning, carried into the Methods sheet so the output is
## self-documenting even if separated from this script.
GRP_COLLAGEN <- c("0.3 mg/ml", "0.3 mg/ml", "1.5 mg/ml", "1.5 mg/ml")
GRP_OXYGEN   <- c("anoxia",    "normoxia",  "anoxia",    "normoxia")

N_REP <- 3   # biological replicates per group

## Consistent group colours across every figure.
GRP_COLS <- c("#4E79A7", "#F28E2B", "#59A14F", "#E15759")
names(GRP_COLS) <- GRP_INTERNAL

## --- 2b. Contrasts ----------------------------------------------------------
## The same six pairwise comparisons as the existing project scripts, each in the
## same orientation as the RNA-seq pairwise contrasts (the diagonal 0-3_N_vs_1-5_A
## was re-oriented from the earlier 1-5_A_vs_0-3_N to match RNA). The numerator is
## whichever group is named first. See the ORIENTATION block in the header.
## Defined as strings rather than via
## non-standard evaluation so that the internal names are used literally and a
## typo becomes an error rather than a silent misparse.
CONTRASTS <- c(
  "0-3_A_vs_1-5_A" = "G0_3_A - G1_5_A",   # collagen effect under anoxia
  "0-3_N_vs_1-5_N" = "G0_3_N - G1_5_N",   # collagen effect under normoxia
  "0-3_A_vs_0-3_N" = "G0_3_A - G0_3_N",   # oxygen effect at 0.3 mg/ml
  "1-5_A_vs_1-5_N" = "G1_5_A - G1_5_N",   # oxygen effect at 1.5 mg/ml
  "0-3_N_vs_1-5_A" = "G0_3_N - G1_5_A",   # diagonal (RNA orientation)
  "0-3_A_vs_1-5_N" = "G0_3_A - G1_5_N"    # diagonal
)

## --- 2b-ii. INTERACTION CONTRAST (supplementary) -----------------------------
## The design is fully crossed (2 collagen densities x 2 oxygen states, n = 3),
## so the interaction is estimable, but none of the six pairwise contrasts above
## tests it. Six simple effects cannot answer "does the oxygen response depend on
## matrix stiffness?" - the question a crossed design is run to answer. Comparing
## two separate significance lists ("significant under anoxia but not normoxia")
## is the classic error of interpreting a difference between significant and
## non-significant as a significant difference, and it is exactly what a reader
## given only pairwise tables will do.
##
## Reported as SUPPLEMENTARY, not promoted into the six. It is a distinct
## question, it is the least powerful contrast in the set (all four group means
## enter, so its standard error is the largest), and the manuscript's claims rest
## on the pairwise comparisons. It exists so that any interaction statement in
## the text can be supported by a test of the interaction, and so that its
## absence from the results is a finding rather than an omission.
##
##   (anoxia - normoxia | 0.3 mg/ml) - (anoxia - normoxia | 1.5 mg/ml)
## Positive = the anoxia response is LARGER at 0.3 mg/ml than at 1.5 mg/ml.
CONTRAST_INTERACTION <- c(
  "INTERACTION_oxygen_by_collagen" = "(G0_3_A - G0_3_N) - (G1_5_A - G1_5_N)"
)

## --- 2c. Filtering and statistics -------------------------------------------
MIN_UNIQUE_PEPTIDES <- 2      # >=2 unique peptides required (removes 467 proteins)
REQUIRE_COMPLETE    <- TRUE   # complete observation in all 12 samples; NO imputation
PSEUDOCOUNT         <- 1      # added to peptide counts before the DEqMS variance fit
ADJP_CUTOFF         <- 0.05   # BH-adjusted p, applied WITHIN each contrast
## Effect-size criterion, held as the EXACT log of the intended fold change.
## 0.58 is a rounded log2(1.5) = 0.5849625, and a rounded threshold is a
## threshold nobody can reproduce from the methods sentence: "at least 1.5-fold"
## and "|log2FC| >= 0.58" are different rules, and proteins in the sliver
## [0.58, 0.5849625) satisfy the printed number while failing the stated biology.
## Deriving it from the fold change keeps text and code in agreement by
## construction. The sensitivity grid is derived the same way.
FC_PRIMARY          <- 1.5    # primary effect-size criterion, as a fold change
LFC_PRIMARY         <- log2(FC_PRIMARY)          # 0.5849625
LFC_GRID            <- c(0, log2(1.5), log2(2))  # reported sensitivity grid: 1x, 1.5x, 2x

## Peptide rows count toward a protein's evidence only if BOTH hold:
##   # Proteins == 1  (unique to one master protein group), AND
##   Quan Info blank  (excludes "Shared", "No Quan Values", "Not Reliable").
## A peptide is "detected" in a sample if its status there is one of:
QUANT_STATUS <- c("High", "Peak Found")

## --- 2d. Enrichment ---------------------------------------------------------
GO_ONTOLOGIES  <- c("BP", "MF", "CC")
RUN_KEGG       <- TRUE        # needs internet at runtime; failure is caught, not fatal
KEGG_ORGANISM  <- "mmu"
REACTOME_ORGANISM <- "mouse"  # ReactomePA naming; local reactome.db, no internet needed
ENRICH_PADJ    <- 0.05
ENRICH_MIN_GS  <- 10          # minimum gene-set size
ENRICH_MAX_GS  <- 500         # maximum gene-set size (excludes uninformative terms)

## --- 2e. Reproducibility ----------------------------------------------------
set.seed(1234)   # no stochastic step exists, but fixed for strictness

## --- 2f. Figures ------------------------------------------------------------
## One resolution for every figure. Base-graphics png() devices are sized in
## INCHES (units = "in"), not pixels, so raising the resolution adds pixels
## without shrinking the figure or enlarging its text: a pixel-sized device at a
## higher res is a physically smaller figure with relatively larger text.
FIG_DPI <- 300


## =============================================================================
## 3. PREFLIGHT - PACKAGES, INPUT FILE, OUTPUT FOLDERS
## =============================================================================
## POLICY (chosen by the project owner): check and STOP with copy-pasteable
## install instructions. This script never installs or upgrades anything, because
## silently modifying a shared R library can break other projects.

cat("\n=============================================================\n")
cat("PID430 PROTEOMICS PIPELINE\n")
cat("=============================================================\n\n")
cat("[preflight] checking packages ...\n")

## RColorBrewer was previously required here and never called; requiring a
## package the script does not use can fail a preflight for no reason.
.pkgs_cran <- c("readxl", "openxlsx", "dplyr", "tidyr", "tibble", "stringr",
                "ggplot2", "ggrepel", "pheatmap")
.pkgs_bioc <- c("limma", "DEqMS", "clusterProfiler", "org.Mm.eg.db", "enrichplot",
                "AnnotationDbi", "ReactomePA", "reactome.db")

.missing_cran <- .pkgs_cran[!vapply(.pkgs_cran, requireNamespace, logical(1), quietly = TRUE)]
.missing_bioc <- .pkgs_bioc[!vapply(.pkgs_bioc, requireNamespace, logical(1), quietly = TRUE)]

if (length(.missing_cran) || length(.missing_bioc)) {
  ## INSTALL IN A FRESH R SESSION. This matters and is not boilerplate.
  ## R cannot replace a package whose namespace is currently loaded, and it
  ## cannot unload a namespace that another loaded package imports. ggplot2 is
  ## the usual casualty: DEqMS, enrichplot, DOSE, fgsea, ggrepel and a dozen
  ## others import it, so once any of them is loaded, an installer that wants to
  ## touch ggplot2 fails with
  ##     namespace 'ggplot2' is imported by ... so cannot be unloaded
  ## That error is about session state, not about anything being broken.
  ## update = FALSE / ask = FALSE stops BiocManager from volunteering to update
  ## dependencies you did not ask it to touch, which is what usually triggers it.
  msg <- c("\nRequired R packages are missing. Nothing has been installed.",
           "",
           "Start a FRESH R session with no packages loaded (in RStudio:",
           "Session > Restart R), then paste the following, then re-run this",
           "script in that new session:\n")
  if (length(.missing_cran))
    msg <- c(msg, sprintf('install.packages(c(%s))',
                          paste0('"', .missing_cran, '"', collapse = ", ")))
  if (length(.missing_bioc))
    msg <- c(msg,
             'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")',
             sprintf('BiocManager::install(c(%s), update = FALSE, ask = FALSE)',
                     paste0('"', .missing_bioc, '"', collapse = ", ")))
  msg <- c(msg, "",
           "If you see \"namespace 'ggplot2' ... cannot be unloaded\", the session",
           "was not fresh. Restart R and repeat - do not force the unload.")
  stop(paste(msg, collapse = "\n"), call. = FALSE)
}

suppressPackageStartupMessages({
  library(readxl);  library(openxlsx); library(dplyr);   library(tidyr)
  library(tibble);  library(stringr);  library(ggplot2); library(ggrepel)
  library(pheatmap); library(limma);   library(DEqMS)
  library(clusterProfiler); library(org.Mm.eg.db); library(enrichplot)
})
cat("[preflight] all packages present\n")

## --- Input file -------------------------------------------------------------
if (!file.exists(INPUT_XLSX))
  stop("Input workbook not found:\n  ", normalizePath(INPUT_XLSX, mustWork = FALSE),
       "\nEdit INPUT_XLSX in section 1.", call. = FALSE)

.sheets <- excel_sheets(INPUT_XLSX)
if (!all(c("PROTEINS", "PEPTIDES") %in% .sheets))
  stop("Workbook must contain sheets 'PROTEINS' and 'PEPTIDES'. Found: ",
       paste(.sheets, collapse = ", "), call. = FALSE)
cat("[preflight] input:", normalizePath(INPUT_XLSX), "\n")

## --- Output folder tree -----------------------------------------------------
## Every artefact lands under one root, sorted by kind. Ontology results get
## their own tree with one sub-folder per contrast, because six contrasts x
## (3 GO ontologies + KEGG + Reactome) x (all/up/down) is otherwise unnavigable.
if (is.null(OUTPUT_DIR)) {
  ## DATETIME, not date. A date-only stamp means a second run on the same day
  ## writes into the first run's folder, mixing outputs from two runs under one
  ## RUN_SUMMARY.csv with no error and no visible sign.
  OUTPUT_DIR <- file.path(dirname(normalizePath(INPUT_XLSX)),
                          paste0("PID430_proteomics_results_",
                                 format(Sys.time(), "%Y%m%d_%H%M%S")))
}
if (dir.exists(OUTPUT_DIR) && length(list.files(OUTPUT_DIR)))
  stop("Output folder already exists and is not empty: ", OUTPUT_DIR,
       "\nMove it aside or set OUTPUT_DIR explicitly. Refusing to overwrite.",
       call. = FALSE)
DIR <- list(
  root        = OUTPUT_DIR,
  qc          = file.path(OUTPUT_DIR, "01_QC"),
  processed   = file.path(OUTPUT_DIR, "02_processed_data"),
  results     = file.path(OUTPUT_DIR, "03_differential_abundance"),
  tables      = file.path(OUTPUT_DIR, "03_differential_abundance", "per_contrast_tables"),
  figures     = file.path(OUTPUT_DIR, "04_figures"),
  volcano     = file.path(OUTPUT_DIR, "04_figures", "volcano"),
  onto        = file.path(OUTPUT_DIR, "05_ontology"),
  onto_gl     = file.path(OUTPUT_DIR, "05_ontology", "00_gene_lists"),
  qual        = file.path(OUTPUT_DIR, "06_presence_absence"),
  logs        = file.path(OUTPUT_DIR, "07_logs")
)
invisible(lapply(DIR, dir.create, recursive = TRUE, showWarnings = FALSE))
cat("[preflight] output:", DIR$root, "\n\n")

## Tee all console output to a log file so a run is auditable after the fact.
.logfile <- file.path(DIR$logs, "run_log.txt")
.logcon  <- file(.logfile, open = "wt")
sink(.logcon, split = TRUE); sink(.logcon, type = "message")
on.exit({ sink(type = "message"); sink(); close(.logcon) }, add = TRUE)

## Collects one-line facts for the summary workbook. Every notable number the
## pipeline computes is registered here, so the summary is generated rather than
## hand-typed and therefore cannot drift from the actual run.
.FACTS <- list()
fact <- function(key, value, note = "") {
  .FACTS[[length(.FACTS) + 1]] <<- data.frame(item = key, value = as.character(value),
                                              note = note, stringsAsFactors = FALSE)
  invisible(NULL)
}
fact("input_file", basename(INPUT_XLSX))
## Thresholds recorded at full precision so the methods sentence quotes the value
## the code used, not a rounded transcription of it.
fact("threshold_adj_p", sprintf("%.2f", ADJP_CUTOFF), "Benjamini-Hochberg, within contrast")
fact("threshold_fold_change", sprintf("%.1f", FC_PRIMARY), "primary effect-size criterion")
fact("threshold_log2FC_exact", sprintf("%.7f", LFC_PRIMARY),
     "= log2(1.5), computed not rounded; do not quote as 0.58")
fact("threshold_min_unique_peptides", MIN_UNIQUE_PEPTIDES,
     "applied to PD protein-level '# Unique Peptides', not the per-sample minimum")
fact("run_date", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
## ORIENTATION: there is NO single global slogan, and one must not be reinstated.
## Positive log2FC = higher in that contrast's OWN numerator. That is the only
## statement true of all six. The shorthand "higher in the 0.3 mg/ml or anoxic
## member" fails for the diagonal 0-3_N_vs_1-5_A, whose numerator (G0_3_N) is the
## 0.3 mg/ml NORMOXIC member. Each results table carries
## its own numerator / denominator / higher_in_if_positive columns; each volcano
## subtitle states it per panel. Quote those, never a global rule.
fact("orientation", "per-contrast: positive log2FC = higher in that contrast's numerator",
     "No global stiffness/oxygen slogan applies; see numerator column of each table")


## =============================================================================
## 4. READ THE PROTEIN TABLE
## =============================================================================
## DEPARTURE from deqms_log2_median-norm_noimpute_300826.R, which read a manually
## exported CSV from a local Desktop path. Two problems that fixes: the CSV was an
## untracked manual step in the provenance chain, and read.csv's default
## check.names=TRUE mangles Proteome Discoverer headers ("# PSMs" -> "X..PSMs"),
## so downstream code silently depended on the mangling. Reading the xlsx with
## .name_repair="minimal" keeps the headers verbatim.

cat("== 4. Reading PROTEINS sheet ==\n")
prot <- read_excel(INPUT_XLSX, sheet = "PROTEINS", .name_repair = "minimal")
prot <- as.data.frame(prot, check.names = FALSE)
cat(sprintf("   %d protein rows x %d columns\n", nrow(prot), ncol(prot)))
## Section 2 states the sheet shape and the master/species composition. Register
## rather than assert: if a re-export changes either, the methods text is wrong.
fact("proteins_sheet_cols", ncol(prot))
if ("Master" %in% colnames(prot))
  fact("proteins_all_master", sum(grepl("Master Protein", prot$Master)),
       sprintf("of %d rows in the PROTEINS sheet", nrow(prot)))
.spcol <- intersect(c("Species", "Organism"), colnames(prot))
if (length(.spcol))
  fact("proteins_species_observed",
       paste(sort(unique(na.omit(prot[[.spcol[1]]]))), collapse = "; "),
       sprintf("from the '%s' column of the export", .spcol[1]))

## Locate the 12 abundance columns and order them by group, not by sheet order.
## Matching on the vendor label rather than on position means a re-exported file
## with reordered columns still maps correctly - and if a label is missing we get
## a hard error instead of a silent mis-assignment.
ab_all <- grep("^Abundance: F", colnames(prot), value = TRUE)
if (length(ab_all) != length(GRP_VENDOR) * N_REP)
  stop("Expected ", length(GRP_VENDOR) * N_REP, " 'Abundance: F' columns, found ",
       length(ab_all), call. = FALSE)

ab_cols <- character(0); samp_internal <- character(0); samp_group <- character(0)
for (g in seq_along(GRP_VENDOR)) {
  ## Anchor the vendor label between ", " and "," so that "ANOXIA" cannot also
  ## match "STIFF-ANOXIA" - a substring collision that would silently merge two
  ## groups. This is the single most dangerous line in the parsing step.
  hit <- grep(paste0(", ", GRP_VENDOR[g], ","), ab_all, fixed = TRUE, value = TRUE)
  if (length(hit) != N_REP)
    stop("Group '", GRP_VENDOR[g], "': expected ", N_REP,
         " abundance columns, found ", length(hit), call. = FALSE)
  ## Order replicates by their trailing index so F10..F12 sort numerically.
  hit <- hit[order(as.integer(sub(".*,\\s*(\\d+)$", "\\1", hit)))]
  ab_cols       <- c(ab_cols, hit)
  samp_internal <- c(samp_internal, paste0(GRP_INTERNAL[g], "_", seq_len(N_REP)))
  samp_group    <- c(samp_group, rep(GRP_INTERNAL[g], N_REP))
}

## Explicit sample sheet. The design exists ONLY in the abundance headers; writing
## it out makes it checkable by eye once instead of being re-parsed on trust.
sample_sheet <- data.frame(
  sample_internal = samp_internal,
  group_internal  = samp_group,
  group_display   = GRP_DISPLAY[match(samp_group, GRP_INTERNAL)],
  collagen        = GRP_COLLAGEN[match(samp_group, GRP_INTERNAL)],
  oxygen          = GRP_OXYGEN[match(samp_group, GRP_INTERNAL)],
  vendor_column   = ab_cols,
  stringsAsFactors = FALSE
)
write.csv(sample_sheet, file.path(DIR$qc, "sample_sheet.csv"), row.names = FALSE)
cat("   sample sheet written; group assignment:\n")
print(table(sample_sheet$group_display))

## Abundance matrix. Vendor values are LINEAR (not log) and NOT normalized.
mat_raw <- as.matrix(prot[, ab_cols, drop = FALSE])
mode(mat_raw) <- "numeric"
rownames(mat_raw) <- prot$Accession
colnames(mat_raw) <- samp_internal

## Annotation kept separate from the numeric matrix and keyed on Accession.
anno <- data.frame(
  Accession   = as.character(prot$Accession),
  Gene_Symbol = as.character(prot$`Gene Symbol`),
  Description = as.character(prot$Description),
  Entrez      = as.character(prot$`Entrez Gene ID`),
  Ensembl     = as.character(prot$`Ensembl Gene ID`),
  n_peptides  = suppressWarnings(as.numeric(prot$`# Peptides`)),
  n_unique    = suppressWarnings(as.numeric(prot$`# Unique Peptides`)),
  n_psms      = suppressWarnings(as.numeric(prot$`# PSMs`)),
  fdr_conf    = as.character(prot$`Protein FDR Confidence: Combined`),
  exp_qvalue  = suppressWarnings(as.numeric(prot$`Exp. q-value: Combined`)),
  stringsAsFactors = FALSE
)
stopifnot(!any(duplicated(anno$Accession)))   # Accession must be a unique key
fact("proteins_in_export", nrow(prot))


## =============================================================================
## 5. QC AND DESIGN AUDIT
## =============================================================================
## Nothing equivalent exists in the vendor deliverable. Two of the reviewer's
## complaints are really about undisclosed data properties, so these numbers are
## computed and reported rather than assumed.

cat("\n== 5. QC ==\n")

## --- 5a. Identification confidence -----------------------------------------
## Documents the 1% SEARCH-level FDR that is already baked into the export.
## This is NOT the statistical FDR of section 10 - see the Methods for why the
## distinction matters. Cross-tabulating the two confidence columns shows the
## threshold was applied and lets us prove later that the peptide filter
## subsumes it entirely.
conf_tab <- table(
  `Protein FDR Confidence` = anno$fdr_conf,
  `Exp. q-value`           = cut(anno$exp_qvalue, c(-Inf, 0.01, 0.05, Inf),
                                 labels = c("<0.01", "0.01-0.05", ">0.05")),
  useNA = "ifany"
)
write.csv(as.data.frame(conf_tab),
          file.path(DIR$qc, "identification_confidence_crosstab.csv"), row.names = FALSE)
print(conf_tab)

## --- 5b. Detection status and match-between-runs ----------------------------
## Three states, not two: "High" = MS/MS identified in that run; "Peak Found" =
## quantified from an MS1 peak with NO MS/MS in that run, i.e. Minora
## match-between-runs; "Not Found" = absent. The Peak Found fraction must be
## disclosed - a reviewer inspecting the export can compute it independently.
fis_cols <- grep("^Found in Sample: ", colnames(prot), value = TRUE)
stopifnot(length(fis_cols) == 12)
fis <- as.matrix(prot[, fis_cols, drop = FALSE])
det_counts <- table(as.vector(fis))
print(det_counts)
for (nm in names(det_counts)) fact(paste0("detection_cells_", nm), det_counts[[nm]])

## Per-protein reliance on match-between-runs, carried into the results tables
## as a flag column rather than reported only in aggregate.
mbr_frac <- rowMeans(fis == "Peak Found", na.rm = TRUE)
names(mbr_frac) <- prot$Accession

## --- 5c. Loading / injection variation --------------------------------------
## Justifies normalization quantitatively instead of asserting it is needed.
lib_tot <- colSums(mat_raw, na.rm = TRUE)
lib_med <- apply(log2(mat_raw), 2, median, na.rm = TRUE)
cat(sprintf("   per-sample total abundance spans %.2f-fold\n", max(lib_tot) / min(lib_tot)))
cat(sprintf("   per-sample raw log2 median spans %.3f log2 (%.2f-fold)\n",
            diff(range(lib_med)), 2^diff(range(lib_med))))
## STAGE LABELS ARE LOAD-BEARING. These are computed on the RAW, UNFILTERED
## matrix (all proteins in the export, before the completeness and peptide
## filters). Section 8 reports a similar-looking median spread on the FILTERED
## matrix immediately before normalization. The two are close on this dataset
## (1.278 vs 1.273 on the frozen run) and are easily mistaken for one number, but
## quoting one and attributing it to the other stage is a methods error a
## reviewer can catch from the supplementary tables. Every loading fact below
## carries its stage in the note field; the methods must quote one stage only.
fact("loading_span_total_fold", sprintf("%.2f", max(lib_tot) / min(lib_tot)),
     "RAW unfiltered matrix, all proteins in the export")
fact("loading_span_median_log2", sprintf("%.3f", diff(range(lib_med))),
     "RAW unfiltered matrix, all proteins in the export - NOT the filtered matrix")

## --- 5d. Missingness --------------------------------------------------------
na_cell    <- sum(is.na(mat_raw))
na_persamp <- colSums(is.na(mat_raw))
na_perprot <- rowSums(is.na(mat_raw))
cat(sprintf("   missing cells: %d of %d (%.1f%%)\n",
            na_cell, length(mat_raw), 100 * na_cell / length(mat_raw)))
fact("missing_cells_pct", sprintf("%.1f", 100 * na_cell / length(mat_raw)))

write.csv(data.frame(sample = samp_internal,
                     group = GRP_DISPLAY[match(samp_group, GRP_INTERNAL)],
                     total_abundance = lib_tot, median_log2 = lib_med,
                     n_missing = na_persamp),
          file.path(DIR$qc, "per_sample_qc.csv"), row.names = FALSE)

## --- 5e. QC figure ----------------------------------------------------------
## Base graphics deliberately: it writes a fixed multi-panel layout to one file
## without needing a grid-assembly package.
png(file.path(DIR$qc, "QC_overview.png"), width = 12, height = 9, units = "in", res = FIG_DPI)
op <- par(mfrow = c(2, 2), mar = c(9, 5, 4, 2))
boxplot(log2(mat_raw), outline = FALSE, las = 2, col = GRP_COLS[samp_group],
        main = "Raw log2 abundance", ylab = "log2 abundance", names = samp_internal)
bp <- barplot(na_persamp, las = 2, col = GRP_COLS[samp_group],
              main = "Missing values per sample", ylab = "n missing",
              ylim = c(0, max(na_persamp) * 1.25))
text(bp, na_persamp, labels = na_persamp, pos = 3, cex = 0.8, xpd = NA)
tb <- table(na_perprot)
bp2 <- barplot(tb, col = "#4E79A7", main = "Missing values per protein",
               xlab = "n missing (of 12)", ylab = "n proteins",
               ylim = c(0, max(tb) * 1.15))
text(bp2, tb, labels = tb, pos = 3, cex = 0.7, xpd = NA)
## Colour by a NAMED status map, not by position. table() returns its levels in
## alphabetical order, so indexing a colour vector by seq_along() silently
## reassigns colours whenever a status is absent from an export. Unlisted
## statuses fall back to grey rather than borrowing another status's colour.
DET_COLS <- c("High" = "#59A14F", "Peak Found" = "#F28E2B", "Not Found" = "#E15759")
.det_col <- unname(DET_COLS[names(det_counts)])
.det_col[is.na(.det_col)] <- "#BAB0AC"
barplot(det_counts, col = .det_col,
        main = "Detection status (protein x sample cells)", ylab = "n cells", las = 1)
par(op); dev.off()
cat("   QC_overview.png written\n")


## =============================================================================
## 6. UNIQUE PEPTIDE COUNTS FROM THE PEPTIDES SHEET
## =============================================================================
## Logic reproduced from unique_peptide_count_table_v2_300826.R, which is correct.
## DEPARTURE: columns are resolved BY NAME from the real header row instead of by
## hardcoded index (2/3/11/20/14/27:38). The indices are right for this file, but
## a re-export with one extra column would silently shift them and the failure
## would be invisible - wrong counts, no error.
##
## SHEET STRUCTURE (verified): 62,553 rows in nested blocks.
##   row 1                       top-level protein header
##   then, repeating per protein: 1 master-protein row  (col "Master" == "Master Protein")
##                                1 peptide sub-header  (starts "Abundances (Grouped) CV")
##                                N peptide data rows
## Resolves exactly: 1 + 5,239 + 5,239 + 52,074 = 62,553, nothing unclassified.

cat("\n== 6. Unique peptide counts ==\n")
pep_raw <- read_excel(INPUT_XLSX, sheet = "PEPTIDES", col_names = FALSE,
                      col_types = "text", .name_repair = "minimal")
pep_raw <- as.data.frame(pep_raw, check.names = FALSE)
cat(sprintf("   PEPTIDES sheet: %d rows x %d cols\n", nrow(pep_raw), ncol(pep_raw)))
## Section 2 of the methods quotes both sheet shapes.
fact("peptides_sheet_rows", nrow(pep_raw),
     "raw nested block layout: master-protein rows interleaved with peptide rows")
fact("peptides_sheet_cols", ncol(pep_raw))

## Find the peptide header row (the one naming "Annotated Sequence") and use it
## to resolve column positions by name.
pep_hdr_row <- which(apply(pep_raw[1:6, ], 1,
                           function(r) any(r == "Annotated Sequence", na.rm = TRUE)))[1]
if (is.na(pep_hdr_row))
  stop("Could not locate the peptide header row in the PEPTIDES sheet.", call. = FALSE)
pep_hdr <- as.character(unlist(pep_raw[pep_hdr_row, ]))

col_of <- function(nm, hdr = pep_hdr) {
  i <- which(hdr == nm)
  if (length(i) != 1) stop("PEPTIDES sheet: expected exactly one '", nm,
                           "' column, found ", length(i), call. = FALSE)
  i
}
## POSITIONAL, and unavoidably so: the master-protein rows sit ABOVE the peptide
## header in each block, so these two columns have no entry in pep_hdr and cannot
## be resolved by name. Positional access in a sheet whose layout could change on
## a re-export is a real risk, so the positions are asserted against the content
## actually found there rather than trusted.
C_MASTER   <- 2                              # "Master" flag column (no name in pep header)
C_ACC      <- 3                              # Accession, populated on master rows only

.mflag <- as.character(pep_raw[[C_MASTER]])
.nmast <- sum(!is.na(.mflag) & .mflag == "Master Protein")
if (.nmast == 0) {
  .cand <- which(vapply(seq_len(min(8L, ncol(pep_raw))),
                        function(j) any(as.character(pep_raw[[j]]) == "Master Protein",
                                        na.rm = TRUE), logical(1)))
  stop(sprintf(paste0("PEPTIDES column %d holds no 'Master Protein' rows, so the block ",
                      "structure cannot be walked and every peptide count would be 0.\n",
                      "Columns that DO contain 'Master Protein': %s.\n",
                      "The sheet layout has changed - update C_MASTER and C_ACC."),
               C_MASTER,
               if (length(.cand)) paste(.cand, collapse = ", ") else "none in first 8"),
       call. = FALSE)
}
cat(sprintf("   ASSERTED: PEPTIDES col %d is the Master flag (%d master rows)\n",
            C_MASTER, .nmast))
C_QUAN     <- col_of("Quan Info")
C_NPROT    <- col_of("# Proteins")           # NOT "# Protein Groups" - different meaning
C_SEQ      <- col_of("Annotated Sequence")   # does not encode modification state
## Detection columns, ordered to MATCH ab_cols rather than taken in file order.
##
## These twelve columns later receive the sample names in samp_internal
## (colnames(found01) <- samp_internal). samp_internal was built by walking the
## groups in GRP_VENDOR order and sorting replicates by trailing index, so it
## follows the REGROUPED abundance order. Taking the peptide columns in whatever
## order they appear in the PEPTIDES sheet and then stamping those names on top
## is only correct if the two orders happen to coincide. On the delivered export
## they do, but nothing enforced it, and a re-export that grouped the sheets
## differently would mislabel every per-sample peptide count while producing a
## file that looks entirely normal. The twelve-sample minimum is invariant to
## column permutation, so DEqMS was never at risk; the per-sample columns of
## unique_peptide_counts.csv were.
##
## Both header families carry the same "F<n>:" sample index, so the mapping is
## recoverable rather than assumed. Note the anchored pattern with the trailing
## colon: "Found in Sample Group: ANOXIA" columns sit alongside these and carry
## no F-index, so an unanchored grep would pull in four per-group summaries.
.found_hdr <- grep("^Found in Sample: ", pep_hdr, value = TRUE)
.f_pep     <- sub("^.*\\]\\s*F(\\d+):.*$", "\\1", .found_hdr)
.f_ab      <- sub("^Abundance:\\s*F(\\d+):.*$", "\\1", ab_cols)
if (length(.found_hdr) != 12 || any(.f_pep == .found_hdr))
  stop(sprintf(paste0("Expected 12 'Found in Sample: [S<n>] F<n>: ...' columns on the ",
                      "PEPTIDES sheet with a parseable F-index; found %d, of which %d ",
                      "did not match the expected header format. Detection columns ",
                      "cannot be aligned to the abundance columns - inspect the sheet ",
                      "header before trusting any per-sample peptide count."),
               length(.found_hdr), sum(.f_pep == .found_hdr)), call. = FALSE)
if (!setequal(.f_pep, .f_ab))
  stop(sprintf(paste0("The PEPTIDES detection columns and the PROTEINS abundance ",
                      "columns describe different sample sets. Peptide F-indices: %s. ",
                      "Abundance F-indices: %s."),
               paste(sort(as.integer(.f_pep)), collapse = ","),
               paste(sort(as.integer(.f_ab)), collapse = ",")), call. = FALSE)

## Reorder into ab_cols order, so colnames(found01) <- samp_internal is correct
## by construction rather than by coincidence.
C_FOUND <- vapply(.found_hdr[match(.f_ab, .f_pep)], col_of,
                  integer(1), USE.NAMES = FALSE)
stopifnot(length(C_FOUND) == 12, !anyNA(C_FOUND), !anyDuplicated(C_FOUND))
cat(sprintf("   ASSERTED: %d detection columns aligned to abundance order (F%s)\n",
            length(C_FOUND), paste(.f_ab, collapse = ", F")))

## Walk the blocks, tagging each peptide row with the protein it belongs to.
## A loop is used rather than a vectorised fill because the block structure is
## sequential and this runs once in a few seconds.
n <- nrow(pep_raw)
acc_vec <- character(n); is_pep <- logical(n)
cur <- NA_character_; all_acc <- character(0)
mcol <- as.character(pep_raw[[C_MASTER]]); acol <- as.character(pep_raw[[C_ACC]])
for (i in seq_len(n)) {
  mv <- mcol[i]
  if (!is.na(mv) && mv == "Master Protein") {           # start of a new protein block
    cur <- acol[i]; all_acc <- c(all_acc, cur); next
  }
  if (!is.na(mv) && startsWith(mv, "Abundances (Grouped) CV")) next   # sub-header
  if (i == 1) next                                                    # top-level header
  if (!is.na(cur)) { acc_vec[i] <- cur; is_pep[i] <- TRUE }
}
cat(sprintf("   master proteins %d | peptide rows %d\n", length(all_acc), sum(is_pep)))
## STOP, not warn. The two sheets are two views of one search result; if the
## master-protein sets do not agree, the block walk has mis-parsed and the
## peptide counts are attached to the wrong proteins. A warning here scrolls past
## in a long run and the resulting counts look entirely plausible.
if (length(unique(all_acc)) != nrow(prot)) {
  .miss <- setdiff(prot$Accession, unique(all_acc))
  .extra <- setdiff(unique(all_acc), prot$Accession)
  stop(sprintf(paste0("Master protein count in PEPTIDES (%d) differs from PROTEINS (%d).\n",
                      "  in PROTEINS but not PEPTIDES: %d%s\n",
                      "  in PEPTIDES but not PROTEINS: %d%s\n",
                      "The block walk has mis-parsed, or the two sheets are from ",
                      "different searches. Peptide counts cannot be trusted."),
               length(unique(all_acc)), nrow(prot),
               length(.miss),  if (length(.miss))  paste0(" (e.g. ", .miss[1], ")")  else "",
               length(.extra), if (length(.extra)) paste0(" (e.g. ", .extra[1], ")") else ""),
       call. = FALSE)
}
cat(sprintf("   ASSERTED: master protein sets agree between sheets (%d)\n", nrow(prot)))

pepd <- pep_raw[is_pep, , drop = FALSE]
acc_p <- acc_vec[is_pep]

## --- 6a. EXPORT-CONVENTION DIAGNOSTICS --------------------------------------
## Two conventions are assumed by the filter below, and if either is wrong the
## ENTIRE DEqMS precision covariate is wrong - silently, because a wrong count
## vector still produces a plausible variance trend and plausible p-values.
## Neither can be verified from the vendor's documentation, so they are dumped
## to disk on every run and must be read once against the frozen workbook.
##
##   (a) Quan Info blank = the usable, quantified state. Proteome Discoverer's
##       vocabulary for this column is Unique / Not Unique / No Quan Values /
##       Redundant / Filtered Out / Excluded by Method, with the ordinary
##       quantified case left EMPTY. If this export instead writes "Unique" for
##       usable peptides, the filter below keeps nothing and every count is 0.
##   (b) Annotated Sequence carries no modification markup, so grouping on it
##       collapses modification variants of one backbone. If it instead encodes
##       modifications (e.g. "sequence [Oxidation]"), variants stay separate and
##       every count is inflated.
quan_tab <- as.data.frame(table(`Quan Info` = ifelse(is.na(pepd[[C_QUAN]]), "<NA>",
                                                     ifelse(pepd[[C_QUAN]] == "", "<blank>",
                                                            pepd[[C_QUAN]]))),
                          stringsAsFactors = FALSE)
names(quan_tab) <- c("quan_info_value", "n_peptide_rows")
quan_tab <- quan_tab[order(-quan_tab$n_peptide_rows), ]
write.csv(quan_tab, file.path(DIR$qc, "peptide_quan_info_states.csv"), row.names = FALSE)
cat("   Quan Info states in PEPTIDES (see 01_QC/peptide_quan_info_states.csv):\n")
print(utils::head(quan_tab, 10), row.names = FALSE)
fact("peptide_quan_info_states_observed",
     paste(quan_tab$quan_info_value, collapse = "; "),
     paste0("blank/NA is the quantified state; the named states are excluded. ",
            "Confirms the 'Quan Info blank = usable' convention for this export"))

## Annotated Sequence in this export uses Proteome Discoverer's flanking-residue
## format: "[R].AAEEAEEAR.[E]", where the bracketed single residues are the
## preceding and following amino acids in the parent protein, NOT modifications.
## Those flanks must be stripped before testing for modification markup,
## otherwise the test flags 100% of rows and says nothing.
seq_chr  <- as.character(pepd[[C_SEQ]])
## The flank can hold MORE THAN ONE residue when the preceding/following residue
## is ambiguous across the protein group: "[KR].EGVMVAK.[ED]". A single-residue
## pattern leaves those brackets in place and mis-reports them as modifications
## (1,104 rows on this export, all of them benign). Quantify one-or-more.
seq_core <- sub("^\\[[A-Z-]+\\]\\.", "", seq_chr)     # drop leading  [R]. or [KR].
seq_core <- sub("\\.\\[[A-Z-]+\\]$", "", seq_core)    # drop trailing .[E] or .[ED]
## Some rows carry a flank at only one end; strip a bare terminal dot too.
seq_core <- sub("^\\.", "", sub("\\.$", "", seq_core))
## Modification markup would survive as brackets, parentheses or lowercase
## residues inside the core sequence.
## Any residual bracket, parenthesis, dot or lowercase letter in the core is real
## markup (or an unrecognized flank form) and must be reported.
seq_mod  <- grepl("[][().]", seq_core) | grepl("[a-z]", seq_core)
write.csv(data.frame(row = seq_len(min(30L, length(seq_chr))),
                     annotated_sequence = utils::head(seq_chr, 30),
                     core_sequence      = utils::head(seq_core, 30),
                     has_modification_markup = utils::head(seq_mod, 30)),
          file.path(DIR$qc, "peptide_sequence_examples.csv"), row.names = FALSE)
cat(sprintf(paste0("   Annotated Sequence: flanking-residue format; %d of %d rows ",
                   "(%.2f%%) carry modification markup in the core\n"),
            sum(seq_mod, na.rm = TRUE), length(seq_mod),
            100 * mean(seq_mod, na.rm = TRUE)))
fact("peptide_rows_with_sequence_markup", sum(seq_mod, na.rm = TRUE),
     paste0("modification markup inside the core sequence, flanks stripped; ",
            "0 confirms modification variants collapse correctly"))
if (any(seq_mod, na.rm = TRUE))
  warning("Annotated Sequence carries modification markup on ", sum(seq_mod, na.rm = TRUE),
          " rows after stripping flanking residues. The modification-variant ",
          "collapse assumes it does not. Unique peptide counts may be inflated - ",
          "inspect 01_QC/peptide_sequence_examples.csv before quoting any count.",
          call. = FALSE)

## The flanking residues are part of the grouping key below. One peptide sequence
## occurring at two positions in the same protein would carry different flanks
## and be counted twice. Quantified here rather than assumed away.
.core_key <- paste(acc_p, seq_core)
.full_key <- paste(acc_p, seq_chr)
.n_flank_split <- length(unique(.full_key)) - length(unique(.core_key))
fact("peptide_variants_split_by_flanking_residue", .n_flank_split,
     "same core peptide counted twice due to differing flanks; 0 means no effect")
if (.n_flank_split > 0)
  cat(sprintf(paste0("   NOTE: %d (protein, peptide) pairs split by flanking residues; ",
                     "counts use the annotated form\n"), .n_flank_split))

## Row-level filter: unique to one protein AND reliably quantified.
np   <- suppressWarnings(as.numeric(pepd[[C_NPROT]]))
qi   <- pepd[[C_QUAN]]
pass <- (np == 1) & (is.na(qi) | qi == "")
pass[is.na(pass)] <- FALSE          # conservative: anything unparseable is excluded
cat(sprintf("   peptide rows passing unique+reliable filter: %d of %d\n",
            sum(pass), nrow(pepd)))
fact("peptide_rows_passing_filter", sum(pass),
     sprintf("of %d peptide rows; unique to one protein and Quan Info blank", nrow(pepd)))
## Section 4(c) of the methods quotes the uniqueness count on its own, separate
## from the combined filter above, to define what "unique" means in this export.
fact("peptide_rows_unique_to_one_protein", sum(np == 1, na.rm = TRUE),
     sprintf("of %d peptide rows; PD '# Proteins' == 1, before the Quan Info condition",
             nrow(pepd)))
fact("peptide_rows_total", nrow(pepd))

## If convention (a) is wrong, this fires rather than producing an all-zero
## count vector that would look like a real result.
if (sum(pass) < 0.1 * nrow(pepd))
  stop(sprintf(paste0("Only %d of %d peptide rows (%.1f%%) passed the unique+quantified ",
                      "filter. This strongly suggests the 'Quan Info blank = usable' ",
                      "convention does not hold for this export. Read ",
                      "01_QC/peptide_quan_info_states.csv and correct the filter ",
                      "before trusting any peptide count."),
               sum(pass), nrow(pepd), 100 * sum(pass) / nrow(pepd)), call. = FALSE)

## Per-row detection matrix across the 12 samples.
found <- as.matrix(pepd[pass, C_FOUND, drop = FALSE])
found01 <- matrix(as.integer(found %in% QUANT_STATUS),
                  nrow = nrow(found), ncol = ncol(found))
colnames(found01) <- samp_internal

## Collapse modification variants: one row per (protein, base sequence), detected
## in a sample if ANY modified form of it was. Without this, the oxidised and
## unmodified forms of one backbone would count as two independent peptides.
row_df <- data.frame(Accession = acc_p[pass],
                     Sequence  = as.character(pepd[[C_SEQ]])[pass],
                     found01, check.names = FALSE, stringsAsFactors = FALSE)
pep_lvl <- row_df %>%
  group_by(Accession, Sequence) %>%
  summarise(across(all_of(samp_internal), ~ as.integer(any(.x == 1))), .groups = "drop")

## Sum peptide-level flags per protein per sample.
prot_cnt <- pep_lvl %>%
  group_by(Accession) %>%
  summarise(across(all_of(samp_internal), sum), .groups = "drop") %>%
  as.data.frame()

## Every protein gets a row: one with no surviving peptide evidence becomes 0,
## rather than vanishing from the join (which a plain aggregate would allow).
prot_cnt <- data.frame(Accession = anno$Accession, stringsAsFactors = FALSE) %>%
  left_join(prot_cnt, by = "Accession")
prot_cnt[is.na(prot_cnt)] <- 0

## DESIGN DECISION - the minimum is taken across ALL TWELVE samples, not per
## contrast. Rationale: spectraCounteBayes fits the variance-count trend ONCE and
## reads a SINGLE count per protein, so per-contrast counts are structurally
## incompatible with the pooled 12-sample fit that provides the 8 residual df.
## The two choices genuinely diverge - do not quote a number here, read the
## registered fact `per_contrast_count_would_differ_max` from RUN_SUMMARY.csv of
## the run in hand, which counts the tested proteins whose per-contrast minimum
## would differ from the twelve-sample minimum in the contrast where the
## divergence is largest. So this is a real choice, not a cosmetic one; it is
## disclosed in the Methods. Using the twelve-sample minimum also matches the
## scope of the model actually being fitted.
prot_cnt$MIN_across_12 <- apply(prot_cnt[, samp_internal], 1, min)
prot_cnt$count_for_DEqMS <- prot_cnt$MIN_across_12 + PSEUDOCOUNT
write.csv(prot_cnt, file.path(DIR$processed, "unique_peptide_counts.csv"), row.names = FALSE)
cat(sprintf("   proteins with MIN == 0 across 12 samples: %d\n",
            sum(prot_cnt$MIN_across_12 == 0)))


## =============================================================================
## 7. FILTERING - NO IMPUTATION AT ANY POINT
## =============================================================================

cat("\n== 7. Filtering ==\n")

## TWO DIFFERENT UNIQUE-PEPTIDE QUANTITIES ARE IN PLAY. They are not
## interchangeable and the methods text must name both:
##
##   (i)  anno$n_unique - Proteome Discoverer's protein-level "# Unique Peptides",
##        a single number per protein computed across the whole experiment.
##        THIS is the inclusion filter (>= MIN_UNIQUE_PEPTIDES).
##   (ii) prot_cnt$MIN_across_12 - the minimum, across the 12 samples, of the
##        per-sample unique peptide count reconstructed from the PEPTIDES sheet.
##        THIS (plus the pseudocount) is the DEqMS precision covariate.
##
## (ii) <= (i) always, and usually strictly: a protein seen with 5 unique
## peptides overall may have only 1 in the worst sample. A reader given the
## single sentence "proteins supported by fewer than two unique peptides were
## removed" will assume one quantity governs both roles. It does not, and the
## consequence is concrete: proteins pass the filter on (i) while entering the
## variance model at (ii) = 0 or 1. The cross-tabulation below is written out so
## the disclosure is quantified from the run rather than asserted in prose.
.mn <- prot_cnt$MIN_across_12[match(anno$Accession, prot_cnt$Accession)]
pep_defs <- data.frame(
  quantity = c("PD protein-level # Unique Peptides (inclusion filter)",
               "per-sample minimum across 12 samples (DEqMS covariate)"),
  median   = c(median(anno$n_unique, na.rm = TRUE), median(.mn, na.rm = TRUE)),
  n_below_2 = c(sum(anno$n_unique < 2, na.rm = TRUE), sum(.mn < 2, na.rm = TRUE)),
  n_zero    = c(sum(anno$n_unique == 0, na.rm = TRUE), sum(.mn == 0, na.rm = TRUE)),
  stringsAsFactors = FALSE)
print(pep_defs, row.names = FALSE)
write.csv(pep_defs, file.path(DIR$qc, "unique_peptide_definitions.csv"), row.names = FALSE)

keep_complete <- if (REQUIRE_COMPLETE) complete.cases(mat_raw) else rep(TRUE, nrow(mat_raw))
keep_pep      <- !is.na(anno$n_unique) & anno$n_unique >= MIN_UNIQUE_PEPTIDES
keep          <- keep_complete & keep_pep

cascade <- data.frame(
  step = c("proteins in export", "quantified in all 12 samples",
           paste0(">= ", MIN_UNIQUE_PEPTIDES, " unique peptides")),
  n_retained = c(nrow(mat_raw), sum(keep_complete), sum(keep)),
  n_removed  = c(NA, nrow(mat_raw) - sum(keep_complete), sum(keep_complete) - sum(keep)),
  stringsAsFactors = FALSE
)
print(cascade)
write.csv(cascade, file.path(DIR$qc, "filter_cascade.csv"), row.names = FALSE)
for (i in seq_len(nrow(cascade))) fact(paste0("filter_", i, "_", cascade$step[i]),
                                       cascade$n_retained[i])

## --- 7a. The relaxed-completeness counterfactual -----------------------------
## The methods section defends complete-case analysis over imputation with a
## quantitative claim: the proteins imputation would recover are structurally
## different, not marginal members of the same population. That claim was
## previously computed by hand and pasted into the text. Compute it here so it
## is regenerated with the run and can be quoted from RUN_SUMMARY.csv.
stopifnot(identical(colnames(mat_raw), sample_sheet$sample_internal))
.grp_idx  <- split(seq_len(ncol(mat_raw)), sample_sheet$group_internal)
.ok_2of3  <- Reduce(`&`, lapply(.grp_idx, function(ix)
                    rowSums(!is.na(mat_raw[, ix, drop = FALSE])) >= 2))
.added    <- .ok_2of3 & !keep_complete          # gained by relaxing 12/12 to 2-of-3
.add_pep  <- anno$n_unique[.added]
.add_abun <- log2(apply(mat_raw[.added, , drop = FALSE], 1, median, na.rm = TRUE))
.add_pass <- sum(.add_pep >= MIN_UNIQUE_PEPTIDES, na.rm = TRUE)
fact("relaxed_completeness_proteins_added", sum(.added),
     "proteins gained by requiring >= 2 of 3 replicates per group instead of 12/12")
fact("relaxed_completeness_added_median_unique_pep",
     sprintf("%.0f", median(.add_pep, na.rm = TRUE)),
     sprintf("vs %.0f for the tested set; the recovered proteins are thinly supported",
             median(anno$n_unique[keep], na.rm = TRUE)))
fact("relaxed_completeness_added_median_log2_abundance",
     sprintf("%.1f", median(.add_abun, na.rm = TRUE)),
     sprintf("vs %.1f for the tested set", median(log2(apply(mat_raw[keep, , drop = FALSE],
             1, median, na.rm = TRUE)), na.rm = TRUE)))
fact("relaxed_completeness_net_gain_pct", sprintf("%.1f", 100 * .add_pass / sum(keep)),
     sprintf("%d of %d added proteins also pass the unique-peptide filter",
             .add_pass, sum(.added)))
## Register the two comparators explicitly, on the same definitions used above,
## so the methods sentence quotes both sides of the comparison from this file
## rather than one side from here and the other from a note field.
fact("tested_median_unique_peptides", median(anno$n_unique[keep], na.rm = TRUE),
     "comparator for relaxed_completeness_added_median_unique_pep; same definition")
fact("tested_median_log2_abundance",
     sprintf("%.1f", median(log2(apply(mat_raw[keep, , drop = FALSE], 1, median,
                                       na.rm = TRUE)), na.rm = TRUE)),
     "comparator for relaxed_completeness_added_median_log2_abundance; per-protein median across samples, pre-normalization")

mat_f  <- mat_raw[keep, , drop = FALSE]
anno_f <- anno[keep, , drop = FALSE]
cnt_f  <- prot_cnt$count_for_DEqMS[match(anno_f$Accession, prot_cnt$Accession)]
names(cnt_f) <- anno_f$Accession
stopifnot(!any(is.na(cnt_f)), all(cnt_f >= 1))

## --- 7b. Confidence and match-between-runs WITHIN the tested set ------------
## Section 3 and Section 12 of the methods quote these; the export-wide
## crosstab above does not answer them, because the question is what survives
## the two filters. Register them so the text is quoting the run, not a
## hand calculation.
.conf_f <- anno_f$fdr_conf
.q_f    <- anno_f$exp_qvalue
fact("tested_confidence_all_high", sum(.conf_f == "High", na.rm = TRUE),
     sprintf("of %d tested; %d Medium, %d Low", nrow(anno_f),
             sum(.conf_f == "Medium", na.rm = TRUE), sum(.conf_f == "Low", na.rm = TRUE)))
fact("tested_max_exp_qvalue", sprintf("%.3f", max(.q_f, na.rm = TRUE)),
     sprintf("%d tested proteins have q >= 0.01", sum(.q_f >= 0.01, na.rm = TRUE)))
## The complete-case set BEFORE the peptide filter: shows the peptide filter
## subsumes a confidence filter, which is why no separate one is applied.
.conf_c <- anno$fdr_conf[keep_complete]
.q_c    <- anno$exp_qvalue[keep_complete]
fact("complete_case_medium_or_low", sum(.conf_c != "High", na.rm = TRUE),
     sprintf("of %d quantified in all 12; %d have q >= 0.01; all are removed by the unique-peptide filter",
             sum(keep_complete), sum(.q_c >= 0.01, na.rm = TRUE)))
.mbr_f  <- mbr_frac[match(anno_f$Accession, names(mbr_frac))]
fact("tested_peak_found_cells", sum(fis[keep, , drop = FALSE] == "Peak Found", na.rm = TRUE),
     sprintf("of %d tested protein-by-sample values", nrow(anno_f) * ncol(mat_f)))
fact("tested_proteins_mbr_half_or_more", sum(.mbr_f >= 0.5, na.rm = TRUE),
     "tested proteins whose quantification is >= 50% match-between-runs; flagged in the results tables")


## Quantify the divergence between the two quantities WITHIN the tested set, so
## the methods sentence can carry the real number instead of a hedge.
.mn_f  <- cnt_f - PSEUDOCOUNT              # per-sample minimum, tested proteins
.pd_f  <- anno_f$n_unique                  # PD protein-level count, tested proteins
.n_dif <- sum(.mn_f != .pd_f, na.rm = TRUE)
.n_lt2 <- sum(.mn_f < 2, na.rm = TRUE)
.n_zer <- sum(.mn_f == 0, na.rm = TRUE)
cat(sprintf(paste0("   two peptide quantities: differ for %d of %d tested proteins (%.1f%%);\n",
                   "   %d tested proteins have a per-sample minimum < 2, of which %d are 0\n"),
            .n_dif, length(.mn_f), 100 * .n_dif / length(.mn_f), .n_lt2, .n_zer))
fact("tested_proteins_two_peptide_counts_differ", .n_dif,
     sprintf("of %d tested; PD protein-level count vs per-sample minimum", length(.mn_f)))
fact("tested_proteins_min_unique_below_2", .n_lt2,
     "pass the inclusion filter on PD's protein-level count but enter DEqMS with min < 2")
fact("tested_proteins_min_unique_zero", .n_zer,
     "per-sample minimum of 0; enter DEqMS at count 1 after the pseudocount")

## The 51-ish proteins with a per-sample minimum of 0 sit at the extreme edge of
## the variance-count trend (log2(1) = 0 after the pseudocount). They are
## RETAINED - they satisfy the stated inclusion rule, which is PD's protein-level
## count - but they are flagged in every results table so a reader can exclude
## them. Note also that DEqMS's own zero-count guard would NOT have caught them:
## its source reads min(fit$count == 0), the minimum of a logical vector, rather
## than min(fit$count) == 0, so it fires only if EVERY protein has count zero.
## The pseudocount added upstream is what actually prevents log2(0) = -Inf here.
write.csv(data.frame(Accession = anno_f$Accession[.mn_f == 0],
                     Gene_Symbol = anno_f$Gene_Symbol[.mn_f == 0],
                     pd_unique_peptides = .pd_f[.mn_f == 0],
                     min_unique_across_12 = 0),
          file.path(DIR$qc, "tested_proteins_min_unique_zero.csv"), row.names = FALSE)

## The claim that matters most for the reviewer response: zero imputed values.
stopifnot(sum(is.na(mat_f)) == 0)
cat("   ASSERTED: 0 missing values in the tested matrix (no imputation)\n")
fact("imputed_values", 0, "complete-case analysis; no imputation performed")

## Confirms the >=2 unique peptide filter subsumes the identification-confidence
## filter, so no separate confidence filter is needed. Computed, not assumed.
cat(sprintf("   tested set: %.1f%% High confidence, %.1f%% Exp.q<0.01 (max q = %.4f)\n",
            100 * mean(anno_f$fdr_conf == "High", na.rm = TRUE),
            100 * mean(anno_f$exp_qvalue < 0.01, na.rm = TRUE),
            max(anno_f$exp_qvalue, na.rm = TRUE)))
fact("tested_pct_high_confidence",
     sprintf("%.1f", 100 * mean(anno_f$fdr_conf == "High", na.rm = TRUE)))


## =============================================================================
## 8. LOG2 TRANSFORM AND EQUAL-MEDIAN NORMALIZATION
## =============================================================================
## Same method as the existing script (DEqMS::equalMedianNormalization) - no
## change needed there. The substantive difference from the VENDOR analysis is
## ORDER: they normalized after imputation, which lets invented constants set the
## normalization target. Here there is nothing to impute.
## The vendor also used equal-MEAN; on this data mean and median differ by at most
## 0.031 log2 per sample, so the centre choice is immaterial and the ordering is
## the real issue.

cat("\n== 8. Normalization ==\n")
log_f  <- log2(mat_f)

## FINITE GUARD. log2() of a zero gives -Inf and of a negative gives NaN, and
## neither raises an error - they travel silently into lmFit, where they poison
## the residual variance of that protein and, through the pooled variance trend,
## the prior for every other protein. The current export happens to contain no
## zeros and no negatives (minimum positive abundance is ~1.0e4), so this guard
## costs nothing today. It exists for the re-export that does contain one.
if (!all(is.finite(log_f))) {
  bad <- which(!is.finite(log_f), arr.ind = TRUE)
  stop(sprintf(paste0("Non-finite values after log2 transform in %d cell(s); ",
                      "first is protein '%s', sample '%s' (linear value %s).\n",
                      "Zero or negative abundances cannot be log-transformed. ",
                      "Inspect the export before proceeding - do NOT floor them silently."),
               nrow(bad), rownames(log_f)[bad[1, "row"]], colnames(log_f)[bad[1, "col"]],
               format(mat_f[bad[1, "row"], bad[1, "col"]])), call. = FALSE)
}
cat("   ASSERTED: all log2 values finite (no zero or negative abundances)\n")

## EQUAL-MEDIAN NORMALIZATION, then restore the grand mean.
## DEqMS::equalMedianNormalization is column-median SUBTRACTION and nothing else:
##     sizefactor <- matrixStats::colMedians(dat); sweep(dat, 2, sizefactor)
## Applied alone it leaves every column median at exactly 0, which means the
## group-mean columns published in the results tables are median-centred
## DEVIATIONS, not log2 abundances - a reader who takes mean_log2 = -3.2 as an
## abundance is misled by our own column name. Adding the grand mean of the
## pre-normalization sample medians back restores the abundance scale while
## leaving the between-sample correction untouched.
## This is a rigid shift applied identically to all 12 columns, so it cancels
## exactly in every contrast: no log2FC, t, p or CI changes by adding it. It
## affects only the interpretability of the group-mean columns and the y-axis of
## the normalization boxplots. Asserted below rather than asserted in prose.
## NOTE ON THE NAME. grand_med is the MEAN of the twelve sample medians, not
## their median, and the register item it feeds is likewise named
## norm_grand_median_restored_log2. Both names are kept as they are so that the
## methods document's citations continue to resolve against the frozen register;
## the methods text states the quantity correctly as a grand mean. On this
## dataset the twelve medians are close enough that mean and median agree to
## three decimals, and the constant is common to every cell, so it cancels in
## every contrast either way. If these names are ever changed, the methods text
## and RUN_FACTS_FROZEN.csv must be regenerated together.
med_before <- apply(log_f, 2, median, na.rm = TRUE)
grand_med  <- mean(med_before)
mat_norm   <- equalMedianNormalization(log_f) + grand_med
med_after  <- apply(mat_norm, 2, median, na.rm = TRUE)

## The two properties the methods sentence claims, both checked on the real run.
stopifnot(abs(diff(range(med_after))) < 1e-8)        # all sample medians equal
stopifnot(abs(mean(med_after) - grand_med) < 1e-8)   # and equal to the grand mean

cat(sprintf("   median spread before %.3f log2 -> after %.3f log2; max shift %.3f\n",
            diff(range(med_before)), diff(range(med_after)),
            max(abs(med_before - grand_med))))
cat(sprintf("   common sample median after normalization: %.3f log2 (grand mean restored)\n",
            med_after[1]))
fact("norm_median_spread_before_log2", sprintf("%.3f", diff(range(med_before))),
     "filtered, tested matrix (4,042 proteins) - NOT the raw unfiltered matrix")
fact("norm_median_spread_after_log2", sprintf("%.3g", diff(range(med_after))),
     "zero by construction: equal-median normalization equalizes sample medians")
fact("norm_max_sample_shift_log2", sprintf("%.3f", max(abs(med_before - grand_med))),
     "largest per-sample shift applied, filtered matrix")
fact("norm_grand_median_restored_log2", sprintf("%.3f", grand_med),
     "grand mean of pre-normalization sample medians, added back to all columns")
## Section 6 claims the previous analysis's equal-MEAN scaling differs
## immaterially from equal-median on this dataset. Quantify rather than assert.
.mean_before <- colMeans(log_f, na.rm = TRUE)
.shift_med   <- med_before   - grand_med
.shift_mean  <- .mean_before - mean(.mean_before)
fact("norm_median_vs_mean_max_diff_log2",
     sprintf("%.3f", max(abs(.shift_med - .shift_mean))),
     sprintf("largest per-sample disagreement between equal-median and equal-mean scaling (%.1f%% of the median shift)",
             100 * max(abs(.shift_med - .shift_mean)) / max(abs(.shift_med))))

write.csv(data.frame(Accession = rownames(mat_norm), mat_norm, check.names = FALSE),
          file.path(DIR$processed, "normalized_log2_abundance.csv"), row.names = FALSE)

## --- 8a. Post-normalization sample-relationship figures ---------------------
pca <- prcomp(t(mat_norm), scale. = TRUE)
pv  <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                     Sample = colnames(mat_norm),
                     Group = GRP_DISPLAY[match(samp_group, GRP_INTERNAL)])
p_pca <- ggplot(pca_df, aes(PC1, PC2, colour = Group, label = Sample)) +
  geom_point(size = 3) + geom_text_repel(size = 2.8, show.legend = FALSE) +
  scale_colour_manual(values = setNames(GRP_COLS, GRP_DISPLAY)) +
  labs(title = "PCA, normalized log2 abundance",
       subtitle = "0-3 = 0.3 mg/ml collagen, 1-5 = 1.5 mg/ml; A = anoxia, N = normoxia",
       x = paste0("PC1 (", pv[1], "%)"), y = paste0("PC2 (", pv[2], "%)")) +
  theme_bw(base_size = 10)
ggsave(file.path(DIR$qc, "PCA_normalized.png"), p_pca, width = 6.5, height = 5, dpi = FIG_DPI)

png(file.path(DIR$qc, "sample_correlation_heatmap.png"), width = 8, height = 7, units = "in", res = FIG_DPI)
cm <- cor(mat_norm, method = "spearman")
pheatmap(cm, main = "Spearman correlation, normalized log2",
         color = colorRampPalette(c("red3", "white", "royalblue4"))(100),
         breaks = seq(min(cm), 1, length.out = 101), border_color = NA)
dev.off()

png(file.path(DIR$qc, "normalization_boxplots.png"), width = 12, height = 5, units = "in", res = FIG_DPI)
op <- par(mfrow = c(1, 2), mar = c(9, 5, 4, 2))
boxplot(log_f, outline = FALSE, las = 2, col = GRP_COLS[samp_group],
        main = "Before normalization", ylab = "log2 abundance")
boxplot(mat_norm, outline = FALSE, las = 2, col = GRP_COLS[samp_group],
        main = "After equal-median normalization", ylab = "log2 abundance")
par(op); dev.off()


## =============================================================================
## 9. ONE LINEAR MODEL ON ALL TWELVE SAMPLES
## =============================================================================
## THE most consequential change from the vendor analysis. They ran six separate
## two-group tests (4 residual df each). One pooled fit gives 8 residual df on the
## same data, and gives every protein ONE variance estimate instead of six
## mutually inconsistent ones.
##
## ASSUMPTION: homogeneous residual variance across the four groups. Diagnosed
## below rather than assumed silently; limma::arrayWeights is the escape hatch if
## the diagnostic fails badly, but it is not applied automatically.

cat("\n== 9. Linear model ==\n")
grp_factor <- factor(samp_group, levels = GRP_INTERNAL)
design <- model.matrix(~ 0 + grp_factor)
colnames(design) <- levels(grp_factor)
stopifnot(identical(colnames(design), GRP_INTERNAL))

fit <- lmFit(mat_norm, design)
cat(sprintf("   residual df = %d (vendor's per-pair design would give %d)\n",
            unique(fit$df.residual)[1], 2 * (N_REP - 1)))
fact("residual_df_pooled", unique(fit$df.residual)[1])
fact("residual_df_vendor_equivalent", 2 * (N_REP - 1))

## Per-group residual SD, so the pooling assumption is visible in the supplement.
grp_sd <- sapply(GRP_INTERNAL, function(g) {
  sub <- mat_norm[, samp_group == g, drop = FALSE]
  median(apply(sub, 1, sd), na.rm = TRUE)
})
var_diag <- data.frame(group_display = GRP_DISPLAY,
                       median_within_group_SD = round(grp_sd, 4),
                       ratio_to_smallest = round(grp_sd / min(grp_sd), 2))
print(var_diag)
write.csv(var_diag, file.path(DIR$qc, "per_group_variance_diagnostic.csv"), row.names = FALSE)
if (max(grp_sd) / min(grp_sd) > 2)
  warning("Per-group residual SD differs by more than 2-fold. The pooled-variance ",
          "assumption is questionable; consider limma::arrayWeights.", call. = FALSE)

## Contrasts. Built from the fixed strings in section 2b so the internal names are
## taken literally.
cont_mat <- makeContrasts(contrasts = unname(CONTRASTS), levels = design)
colnames(cont_mat) <- names(CONTRASTS)

## Interaction appended as an extra column so it shares the variance model with
## the six primary contrasts. It is kept OUT of names(CONTRASTS) throughout, so
## every downstream loop and the significance summary still see exactly six
## contrasts.
cont_int <- makeContrasts(contrasts = unname(CONTRAST_INTERACTION), levels = design)
colnames(cont_int) <- names(CONTRAST_INTERACTION)
cont_all <- cbind(cont_mat, cont_int)
IDX_INTERACTION <- ncol(cont_all)

## The interaction must be a proper contrast: coefficients +1,-1,-1,+1 summing to
## zero. If a future edit breaks the parenthesisation this catches it.
stopifnot(abs(sum(cont_int[, 1])) < 1e-12,
          sum(cont_int[, 1] != 0) == 4,
          cont_int["G0_3_A", 1] ==  1, cont_int["G0_3_N", 1] == -1,
          cont_int["G1_5_A", 1] == -1, cont_int["G1_5_N", 1] ==  1)
cat("   ASSERTED: interaction contrast is (+1,-1,-1,+1) and sums to zero\n")

## ASSERTION - guards the orientation. For each contrast the numerator group must
## carry coefficient +1 and the denominator -1. Without this, a future edit could
## silently flip every sign in every table and figure, and nothing would catch it.
for (cn in names(CONTRASTS)) {
  parts <- trimws(strsplit(CONTRASTS[[cn]], "-", fixed = TRUE)[[1]])
  stopifnot(length(parts) == 2,
            cont_mat[parts[1], cn] == 1,
            cont_mat[parts[2], cn] == -1)
}
cat("   ASSERTED: for all 6 contrasts the named numerator carries +1 and the\n")
cat("   named denominator -1, so positive log2FC = higher in the numerator.\n")

fit2 <- contrasts.fit(fit, cont_all)   # six primary + the interaction
fit2 <- eBayes(fit2)   # required: spectraCounteBayes reads fit2$sigma and df


## =============================================================================
## 10. DEqMS VARIANCE MODERATION
## =============================================================================
## Empirical-Bayes moderation in which the PRIOR CENTRE depends on peptide count
## instead of being a single constant. Everything else is limma: same design
## matrix, same shrinkage algebra, same moderated t. "Moderated" refers only to
## the variance in the DENOMINATOR of t being shrunk toward the prior; the log2FC
## numerator is untouched. It is NOT a multiple-testing correction (that is
## section 11).
##
## Why not plain eBayes: on this dataset residual variance depends strongly on
## peptide count (rho about -0.72), the observed spread of log residual variance
## (SD about 1.32) far exceeds the sampling-noise floor at 8 df (0.53), and
## limma's own prior fit therefore drives d0 toward zero - switching moderation
## off. Removing the count trend restores a well-determined prior. The diagnostic
## is recomputed below on the actual run so the justification is measured, not
## quoted.

cat("\n== 10. DEqMS variance moderation ==\n")
fit2$count <- cnt_f[rownames(fit2$coefficients)]
stopifnot(!any(is.na(fit2$count)), all(fit2$count >= 1))

## The diagnostic that justifies the count-dependent prior, on THIS run.
rho <- suppressWarnings(cor(log(fit2$sigma^2), log2(fit2$count), method = "spearman"))
sd_obs <- sd(log(fit2$sigma^2))
df_res <- unique(fit$df.residual)[1]
sd_noise <- sqrt(2 / df_res)   # asymptotic SD of log(s2) under a single true variance
cat(sprintf("   rho(log residual var, log2 count) = %.3f\n", rho))
cat(sprintf("   SD of log residual var: observed %.3f vs sampling floor %.3f\n",
            sd_obs, sd_noise))
fact("rho_variance_vs_peptide_count", sprintf("%.3f", rho))
fact("sd_log_residual_variance_observed", sprintf("%.3f", sd_obs))
fact("sd_log_residual_variance_sampling_floor", sprintf("%.2f", sd_noise),
     sprintf("sqrt(2/%d); the spread expected if every protein shared one true variance", df_res))
## Section 8 quotes the variance trend as a range across peptide-count bins.
## Bin on the same quantity the prior uses, so the sentence matches the model.
.vb <- data.frame(v = fit2$sigma^2, k = fit2$count)
.vb$bin <- cut(.vb$k, breaks = quantile(.vb$k, probs = seq(0, 1, 0.2), na.rm = TRUE),
               include.lowest = TRUE)
.vmed <- tapply(.vb$v, .vb$bin, median, na.rm = TRUE)
.vmed <- .vmed[!is.na(.vmed)]
fact("variance_median_lowest_peptide_bin", sprintf("%.3f", .vmed[1]))
fact("variance_median_highest_peptide_bin", sprintf("%.3f", .vmed[length(.vmed)]))
fact("variance_bin_fold_range", sprintf("%.0f", .vmed[1] / .vmed[length(.vmed)]),
     sprintf("median residual variance, lowest vs highest peptide-count quintile (%d bins)",
             length(.vmed)))
fact("variance_bins_monotonic", all(diff(.vmed) < 0),
     "TRUE means median variance declines in every successive peptide-count bin")

## Section 8 states that a per-contrast peptide count would differ from the
## twelve-sample minimum for a large share of proteins. Quantify it on the
## tested set rather than asserting it.
.per_smp  <- as.matrix(prot_cnt[match(anno_f$Accession, prot_cnt$Accession),
                                sample_sheet$sample_internal, drop = FALSE])
.min12    <- apply(.per_smp, 1, min, na.rm = TRUE)
.differs  <- vapply(CONTRASTS, function(cf) {
  gs <- unique(unlist(regmatches(cf, gregexpr("G[0-9]_[0-9]_[AN]", cf))))
  ix <- which(sample_sheet$group_internal %in% gs)
  sum(apply(.per_smp[, ix, drop = FALSE], 1, min, na.rm = TRUE) != .min12, na.rm = TRUE)
}, numeric(1))
fact("per_contrast_count_would_differ_max", max(.differs),
     sprintf("of %d tested proteins, in the contrast where it differs most; the pooled prior requires one count per protein",
             nrow(anno_f)))
if (rho > -0.2)
  warning("Variance-count dependence is weak on this run (rho = ", round(rho, 3),
          "). The DEqMS justification in the Methods assumes a strong trend; ",
          "re-check before quoting it.", call. = FALSE)

## coef_col PASSED EXPLICITLY - do not drop it back to spectraCounteBayes(fit2).
##
## coef_col has NO default in the current DEqMS source (verified against the
## installed build, not the vignette). Omitting it does NOT error, because the
## argument is used only inside a subscript - fit$coefficients[, coef_col] - and
## a missing argument reaching a subscript position is R's EMPTY INDEX, which
## selects all columns. So the omitted form happens to compute all six contrasts
## correctly, and the six p-value vectors it produces are genuinely distinct.
##
## It is still wrong to rely on. The correctness rests on an undocumented
## interaction between argument missingness and subscript semantics. A DEqMS
## release that adds a missing() check, a match.arg, or a default of 1 turns
## this from correct into silently wrong - every contrast would report contrast
## 1's statistics against its own fold changes, with no error raised. Passing
## the columns explicitly costs nothing and removes the dependency entirely.
fit3 <- spectraCounteBayes(fit2, coef_col = seq_len(ncol(cont_all)))

## Assert the moderated statistics really do span all six contrasts. This is the
## check that would catch the failure described above on a future DEqMS version.
stopifnot(ncol(as.matrix(fit3$sca.t)) == ncol(cont_all),
          ncol(as.matrix(fit3$sca.p)) == ncol(cont_all))
cat(sprintf("   ASSERTED: moderated statistics computed for all %d contrasts\n",
            ncol(cont_all)))

## Independent check that no contrast is silently carrying another contrast's
## statistics: within each contrast, the moderated t and the fold change must
## share a sign on every finite row. They are the same quantity up to a positive
## scale factor, so any disagreement means coef_col is misaligned.
for (j in seq_len(ncol(cont_all))) {
  fc <- fit3$coefficients[, j]
  tt <- as.matrix(fit3$sca.t)[, j]
  ok <- is.finite(fc) & is.finite(tt) & fc != 0
  if (!all(sign(fc[ok]) == sign(tt[ok])))
    stop(sprintf(paste0("Contrast %d ('%s'): moderated t and log2FC disagree in sign ",
                        "on %d of %d proteins. coef_col is misaligned."),
                 j, colnames(cont_all)[j], sum(sign(fc[ok]) != sign(tt[ok])), sum(ok)),
         call. = FALSE)
}
cat("   ASSERTED: moderated t and log2FC share sign in every contrast\n")
fact("deqms_contrasts_moderated", ncol(cont_all),
     "coef_col passed explicitly; sca.t column count and t/log2FC sign asserted")

## CRITICAL FIELD DISTINCTION - do not "simplify" this.
##   fit3$df.prior     is NOT overwritten by spectraCounteBayes. It is carried
##                     over unchanged from the upstream limma eBayes call, so
##                     reading it gives limma's single-prior value and silently
##                     understates the DEqMS prior.
##   fit3$sca.dfprior  is the count-dependent prior degrees of freedom that
##                     DEqMS actually fitted, and the value its own sca.P.Value
##                     is computed against.
## Using df.prior here produces p-values and confidence intervals that disagree
## with each other (verified: 35 proteins in one contrast reported p < 0.05 while
## their CI spanned zero). Always use sca.dfprior downstream.
d0_deqms <- fit3$sca.dfprior
cat(sprintf("   DEqMS prior df (d0) = %.3f  ->  total df = %.3f\n",
            d0_deqms, d0_deqms + df_res))
fact("deqms_prior_df_d0", sprintf("%.3f", d0_deqms))
fact("deqms_total_df", sprintf("%.3f", d0_deqms + df_res))

## For the record: the plain-limma prior, showing what the prior would have been
## WITHOUT the count covariate. This is the number that makes the choice
## defensible, so it is computed rather than cited. fit2$df.prior is the correct
## field HERE, because fit2 is the plain-limma fit.
cat(sprintf("   for comparison, plain limma prior df = %.3f\n", fit2$df.prior[1]))
fact("plain_limma_prior_df", sprintf("%.3f", fit2$df.prior[1]))
fact("prior_df_gain_from_count_covariate",
     sprintf("%.2fx", d0_deqms / fit2$df.prior[1]),
     "DEqMS count-dependent prior df relative to limma single global prior")

## Guard: if the count covariate ever stops helping on a re-run, the Methods
## justification for choosing DEqMS no longer holds and must not be quoted.
if (d0_deqms <= fit2$df.prior[1])
  warning("The count-dependent prior is no longer better determined than the ",
          "single global prior (DEqMS d0 = ", round(d0_deqms, 3),
          " vs limma d0 = ", round(fit2$df.prior[1], 3),
          "). Re-check the variance-moderation justification in the Methods.",
          call. = FALSE)

png(file.path(DIR$figures, "DEqMS_variance_diagnostics.png"),
    width = 11, height = 5, units = "in", res = FIG_DPI)
op <- par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))
VarianceBoxplot(fit3, n = 20, xlab = "unique peptides + 1",
                main = "Residual variance vs peptide count")
VarianceScatterplot(fit3, xlab = "unique peptides + 1",
                    main = "DEqMS vs limma prior")
abline(h = log(fit2$s2.prior), col = "green", lwd = 2)
legend("topright", legend = c("DEqMS prior (count-dependent)", "limma prior (single)"),
       col = c("red", "green"), lwd = 2, bty = "n", cex = 0.8)
par(op); dev.off()


## =============================================================================
## 11. CONTRAST RESULTS, BH CORRECTION, AND PER-CONTRAST TABLES
## =============================================================================
## BH is applied WITHIN each contrast (DEqMS::outputResult does this), which is
## the limma::topTable convention. Stated explicitly in the Methods so it is not
## mistaken for an oversight. The six contrasts are not independent - they share
## samples and the variance prior - but BH is valid under positive dependence.

cat("\n== 11. Contrast results ==\n")
## sca.dfprior, not df.prior - see the field-distinction note in section 10.
## This is the same df DEqMS uses for its own p-values, so the CI and the
## p-value agree with each other by construction.
df_total <- d0_deqms + df_res
tcrit    <- qt(0.975, df_total)     # for the moderated 95% CI

## Group mean log2 abundance per protein, appended to every table so a reader can
## see the underlying values behind a fold change without opening another file.
grp_means <- sapply(GRP_INTERNAL, function(g)
  rowMeans(mat_norm[, samp_group == g, drop = FALSE]))
colnames(grp_means) <- paste0("mean_log2_", GRP_SLUG)

res_list <- list(); summ <- list()
for (k in seq_along(CONTRASTS)) {
  cn   <- names(CONTRASTS)[k]
  gnum <- trimws(strsplit(CONTRASTS[[cn]], "-", fixed = TRUE)[[1]])[1]
  gden <- trimws(strsplit(CONTRASTS[[cn]], "-", fixed = TRUE)[[1]])[2]
  disp_num <- GRP_DISPLAY[match(gnum, GRP_INTERNAL)]
  disp_den <- GRP_DISPLAY[match(gden, GRP_INTERNAL)]

  ## outputResult returns limma's ordinary columns PLUS the DEqMS moderated ones
  ## (sca.t, sca.P.Value, sca.adj.pval). The sca.* columns are the ones to use;
  ## the plain t/P.Value are limma's unmoderated-prior values and are kept only
  ## for transparency.
  tt <- outputResult(fit3, coef_col = k) %>% rownames_to_column("Accession")

  ## Moderated 95% CI, computed from the DEqMS posterior variance. limma's
  ## confint=TRUE would give the ORDINARY-prior interval, which is inconsistent
  ## with the p-value being reported, hence the manual computation.
  se <- fit3$stdev.unscaled[tt$Accession, k] * sqrt(fit3$sca.postvar[tt$Accession])
  tt$CI_lower <- tt$logFC - tcrit * se
  tt$CI_upper <- tt$logFC + tcrit * se

  out <- tt %>%
    left_join(anno_f, by = "Accession") %>%
    mutate(
      contrast          = cn,
      numerator         = disp_num,
      denominator       = disp_den,
      ## Self-documenting direction column: removes any need to know the
      ## orientation convention in order to read the table correctly.
      higher_in_if_positive = disp_num,
      direction         = ifelse(logFC > 0, paste("higher in", disp_num),
                                            paste("higher in", disp_den)),
      ## Signed fold change for display. The predicate must be >= 0, not > 0:
      ## with > 0 a protein at exactly logFC = 0 (no change) is reported as
      ## -1, i.e. a one-fold DECREASE, which is both wrong and alarming to read.
      fold_change       = ifelse(logFC >= 0, 2^logFC, -1 / (2^logFC)),
      min_unique_pep    = fit3$count[Accession] - PSEUDOCOUNT,
      mbr_fraction      = round(mbr_frac[Accession], 3),
      mbr_dependent     = mbr_frac[Accession] >= 0.5,
      significant       = sca.adj.pval < ADJP_CUTOFF & abs(logFC) >= LFC_PRIMARY
    ) %>%
    dplyr::select(Accession, Gene_Symbol, Description, contrast, numerator, denominator,
                  higher_in_if_positive, log2FC = logFC, fold_change,
                  CI_lower, CI_upper, moderated_t = sca.t, p_value = sca.P.Value,
                  adj_p_BH = sca.adj.pval, limma_t = t, limma_p = P.Value,
                  direction, significant, min_unique_pep, n_unique_total = n_unique,
                  n_psms, mbr_fraction, mbr_dependent, fdr_conf, exp_qvalue) %>%
    bind_cols(as.data.frame(grp_means[.$Accession, , drop = FALSE])) %>%
    arrange(adj_p_BH, desc(abs(log2FC)))

  res_list[[cn]] <- out
  write.csv(out, file.path(DIR$tables, paste0(cn, "_results.csv")), row.names = FALSE)

  ## Threshold sensitivity grid, per contrast, for the supplement.
  gr <- sapply(LFC_GRID, function(L) sum(out$adj_p_BH < ADJP_CUTOFF & abs(out$log2FC) >= L))
  summ[[cn]] <- data.frame(
    contrast = cn, numerator = disp_num, denominator = disp_den,
    n_tested = nrow(out),
    n_raw_p_lt_0.05 = sum(out$p_value < 0.05),
    n_adjp_only = gr[1], n_adjp_1.5x = gr[2], n_adjp_2x = gr[3],
    n_up_primary   = sum(out$significant & out$log2FC > 0),
    n_down_primary = sum(out$significant & out$log2FC < 0),
    stringsAsFactors = FALSE
  )
  cat(sprintf("   %-16s adj p<0.05: %4d | +1.5x: %4d | +2x: %4d\n",
              cn, gr[1], gr[2], gr[3]))
}
summ_df <- bind_rows(summ)
write.csv(summ_df, file.path(DIR$results, "significance_summary_by_contrast.csv"),
          row.names = FALSE)

## Note on limma::treat. A fold-change-aware test was considered and is NOT used
## or computed here. The user's instruction stands: if a reviewer has not asked
## for a piece of machinery and omitting it does not make the analysis
## incorrect, it is not included. The reported rule is BH adj p < ADJP_CUTOFF
## followed by the effect-size filter, and the sensitivity of every count to
## that filter is already shown by the LFC_GRID columns (1x, 1.5x, 2x) in
## significance_summary_by_contrast.csv. Do not re-add a treat block.

## --- 11a. SUPPLEMENTARY: interaction contrast --------------------------------
## Deliberately built AFTER summ_df and kept out of res_list, so that
## significance_summary_by_contrast.csv, the enrichment loop and the volcano
## loop continue to see exactly the six primary contrasts.
## Written to its own file, with its own header row in RUN_SUMMARY.
cat("\n== 11a. Interaction contrast (supplementary) ==\n")
tt_i <- outputResult(fit3, coef_col = IDX_INTERACTION) %>% rownames_to_column("Accession")
se_i <- fit3$stdev.unscaled[tt_i$Accession, IDX_INTERACTION] *
        sqrt(fit3$sca.postvar[tt_i$Accession])
int_res <- tt_i %>%
  left_join(anno_f, by = "Accession") %>%
  mutate(
    contrast   = names(CONTRAST_INTERACTION),
    CI_lower   = logFC - tcrit * se_i,
    CI_upper   = logFC + tcrit * se_i,
    ## The estimate is a DIFFERENCE OF DIFFERENCES in log2 units, not a fold
    ## change of one group over another. No fold_change column is emitted,
    ## because 2^estimate has no single-comparison interpretation here and
    ## would certainly be quoted as one.
    interpretation = ifelse(logFC > 0,
                            "anoxia response larger at 0.3 mg/ml",
                            "anoxia response larger at 1.5 mg/ml"),
    min_unique_pep = fit3$count[Accession] - PSEUDOCOUNT,
    significant    = sca.adj.pval < ADJP_CUTOFF & abs(logFC) >= LFC_PRIMARY
  ) %>%
  dplyr::select(Accession, Gene_Symbol, Description, contrast,
                interaction_log2 = logFC, CI_lower, CI_upper,
                moderated_t = sca.t, p_value = sca.P.Value, adj_p_BH = sca.adj.pval,
                interpretation, significant, min_unique_pep,
                n_unique_total = n_unique) %>%
  bind_cols(as.data.frame(grp_means[.$Accession, , drop = FALSE])) %>%
  arrange(adj_p_BH, desc(abs(interaction_log2)))

write.csv(int_res, file.path(DIR$results, "interaction_oxygen_by_collagen.csv"),
          row.names = FALSE)
n_int_adjp <- sum(int_res$adj_p_BH < ADJP_CUTOFF)
n_int_prim <- sum(int_res$significant)
cat(sprintf("   %-30s adj p<0.05: %4d | with |effect| >= %.4f log2: %4d\n",
            names(CONTRAST_INTERACTION), n_int_adjp, LFC_PRIMARY, n_int_prim))
fact("interaction_n_adjp_only", n_int_adjp,
     "supplementary oxygen x collagen interaction, BH adj p < 0.05")
fact("interaction_n_primary", n_int_prim,
     sprintf("interaction, adj p < %.2f AND |effect| >= %.4f log2", ADJP_CUTOFF, LFC_PRIMARY))

## Interaction volcano, so the supplementary contrast is presented on the same
## footing as the six primary ones rather than as a bare table.
## Y-AXIS IS BH-ADJUSTED p, matching the six primary volcanoes. An earlier
## version plotted the unadjusted moderated p here while the primary panels used
## adjusted p. Two proteins meet the call rule on this contrast, so an
## unadjusted-p panel next to six adjusted-p panels makes a near-null result
## look populated purely through a change of axis. The dashed horizontal line is
## the same ADJP_CUTOFF used everywhere else.
##
## The vertical lines mark |effect| >= LFC_PRIMARY. Note that this is a
## difference of differences, so 2^estimate is NOT a fold change of one group
## over another; the threshold is applied here as a magnitude convention only.
## The interaction sheet deliberately carries no fold_change column and no
## legend should introduce one.
p_int <- ggplot(int_res, aes(interaction_log2, -log10(adj_p_BH))) +
  geom_point(aes(colour = significant), size = 0.7, alpha = 0.6) +
  scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "firebrick"),
                      guide = "none") +
  geom_hline(yintercept = -log10(ADJP_CUTOFF), linetype = 2, linewidth = 0.3) +
  geom_vline(xintercept = c(-LFC_PRIMARY, LFC_PRIMARY), linetype = 2, linewidth = 0.3) +
  labs(title = "Oxygen x collagen interaction (SUPPLEMENTARY)",
       subtitle = paste0("(0.3 A - 0.3 N) - (1.5 A - 1.5 N); positive = anoxia ",
                         "response larger at 0.3 mg/ml. y = BH adjusted p, as in ",
                         "the six primary volcanoes."),
       x = "interaction estimate (log2, difference of differences)",
       y = "-log10 BH adjusted p") +
  theme_bw(base_size = 10)
ggsave(file.path(DIR$volcano, "volcano_INTERACTION_oxygen_by_collagen.png"),
       p_int, width = 6.5, height = 5, dpi = FIG_DPI)


## =============================================================================
## 12. RESULTS WORKBOOK AND SUMMARY SHEETS
## =============================================================================
## One workbook holding everything a co-author or reviewer needs, with the method
## parameters and the group mapping embedded so the file is interpretable on its
## own. Sheet names are truncated to Excel's 31-character limit.

cat("\n== 12. Writing results workbook ==\n")
wb <- createWorkbook()

addWorksheet(wb, "README")
writeData(wb, "README", data.frame(field = c(
  "Analysis", "Input file", "Run date", "Orientation", "Positive log2FC means",
  "Significance", "Multiple testing", "Variance model", "Normalization",
  "Missing values", "Collagen", "Oxygen"),
  value = c(
  "PID430 proteomics re-analysis (limma + DEqMS)", basename(INPUT_XLSX),
  format(Sys.time(), "%Y-%m-%d %H:%M"),
  "per contrast; see the numerator column of each contrast sheet",
  "higher in THAT contrast's numerator group (no single global rule applies)",
  sprintf("BH adj p < %.2f AND |log2FC| >= %.7f (= log2(%.1f), %.1f-fold)",
          ADJP_CUTOFF, LFC_PRIMARY, FC_PRIMARY, FC_PRIMARY),
  "Benjamini-Hochberg, within each contrast separately",
  "DEqMS spectraCounteBayes, prior centre a function of unique peptide count",
  ## DEqMS::equalMedianNormalization subtracts each column's median and stops
  ## there; the grand mean of the sample medians is added back by this script.
  ## Naming only the DEqMS function credits it with a step it does not perform.
  "log2, then column-median subtraction (DEqMS equalMedianNormalization), then addition of the grand mean of the twelve sample medians",
  "none imputed; complete-case only",
  "0-3 = 0.3 mg/ml, 1-5 = 1.5 mg/ml", "A = anoxia, N = normoxia")))

addWorksheet(wb, "group_mapping")
writeData(wb, "group_mapping", data.frame(
  display = GRP_DISPLAY, collagen = GRP_COLLAGEN, oxygen = GRP_OXYGEN,
  vendor_column_label = GRP_VENDOR, internal = GRP_INTERNAL,
  note = c("unprefixed vendor label = SOFT 0.3 mg/ml", "unprefixed vendor label = SOFT 0.3 mg/ml",
           "", "")))

addWorksheet(wb, "run_facts")
writeData(wb, "run_facts", bind_rows(.FACTS))

addWorksheet(wb, "summary_by_contrast"); writeData(wb, "summary_by_contrast", summ_df)
addWorksheet(wb, "filter_cascade");      writeData(wb, "filter_cascade", cascade)
addWorksheet(wb, "sample_sheet");        writeData(wb, "sample_sheet", sample_sheet)
addWorksheet(wb, "variance_diagnostic"); writeData(wb, "variance_diagnostic", var_diag)
addWorksheet(wb, "peptide_definitions");  writeData(wb, "peptide_definitions", pep_defs)
## Supplementary interaction, placed BEFORE the six contrast sheets and named so
## it cannot be mistaken for a seventh pairwise comparison.
addWorksheet(wb, "SUPP_interaction");     writeData(wb, "SUPP_interaction", int_res)

for (cn in names(res_list)) {
  sn <- substr(cn, 1, 31)
  addWorksheet(wb, sn); writeData(wb, sn, res_list[[cn]])
}
saveWorkbook(wb, file.path(DIR$results, "PID430_DEqMS_all_contrasts.xlsx"), overwrite = TRUE)
cat("   PID430_DEqMS_all_contrasts.xlsx written\n")

## Union of significant proteins across contrasts - the list most often wanted.
union_sig <- bind_rows(res_list) %>% filter(significant) %>%
  arrange(Accession, contrast)
write.csv(union_sig, file.path(DIR$results, "significant_union_all_contrasts.csv"),
          row.names = FALSE)
cat(sprintf("   union of significant protein-contrast pairs: %d (%d unique proteins)\n",
            nrow(union_sig), length(unique(union_sig$Accession))))
fact("union_significant_pairs", nrow(union_sig))
fact("union_significant_proteins", length(unique(union_sig$Accession)))


## =============================================================================
## 13. FIGURES: VOLCANOES AND HEATMAP
## =============================================================================

cat("\n== 13. Figures ==\n")
for (cn in names(res_list)) {
  d <- res_list[[cn]]
  d$cls <- ifelse(d$adj_p_BH >= ADJP_CUTOFF, "not significant",
           ifelse(abs(d$log2FC) < LFC_PRIMARY, "adj p only",
           ifelse(d$log2FC > 0, paste("higher in", d$numerator[1]),
                                paste("higher in", d$denominator[1]))))
  lab <- d %>% filter(significant) %>% slice_head(n = 15)
  p <- ggplot(d, aes(log2FC, -log10(adj_p_BH), colour = cls)) +
    geom_point(size = 0.8, alpha = 0.75) +
    geom_hline(yintercept = -log10(ADJP_CUTOFF), linetype = 2, linewidth = 0.3) +
    geom_vline(xintercept = c(-LFC_PRIMARY, LFC_PRIMARY), linetype = 2, linewidth = 0.3) +
    geom_text_repel(data = lab, aes(label = Gene_Symbol), size = 2.3,
                    max.overlaps = 20, show.legend = FALSE) +
    scale_colour_manual(values = setNames(
      c("grey75", "grey45", "#E15759", "#4E79A7"),
      c("not significant", "adj p only",
        paste("higher in", d$numerator[1]), paste("higher in", d$denominator[1])))) +
    labs(title = gsub("_", " ", cn),
         ## State the fold change, not a rounded log2. "0.58" is not the value
         ## the code applied (log2(1.5) = 0.5849625) and a reader cannot
         ## reconstruct the rule from it.
         subtitle = sprintf("BH adj p < %.2f and fold change >= %.1fx; positive = higher in %s",
                            ADJP_CUTOFF, FC_PRIMARY, d$numerator[1]),
         x = "log2 fold change", y = "-log10 BH adjusted p", colour = NULL) +
    theme_bw(base_size = 9) + theme(legend.position = "bottom")
  ggsave(file.path(DIR$volcano, paste0("volcano_", cn, ".png")), p,
         width = 5.5, height = 5, dpi = FIG_DPI)
}
cat("   6 volcano plots written\n")

## Heatmap of the union of significant proteins, row-scaled (z-score per protein)
## so that patterns across conditions are visible rather than abundance level.
if (length(unique(union_sig$Accession)) >= 2) {
  hm <- mat_norm[unique(union_sig$Accession), , drop = FALSE]
  rn <- anno_f$Gene_Symbol[match(rownames(hm), anno_f$Accession)]
  rn[is.na(rn) | rn == "NA" | rn == ""] <- rownames(hm)[is.na(rn) | rn == "NA" | rn == ""]
  rownames(hm) <- make.unique(rn)
  ann_c <- data.frame(Group = GRP_DISPLAY[match(samp_group, GRP_INTERNAL)])
  rownames(ann_c) <- colnames(hm)
  png(file.path(DIR$figures, "heatmap_significant_union.png"),
      width = 8.5, height = max(7, 0.11 * nrow(hm)), units = "in", res = FIG_DPI)
  pheatmap(hm, scale = "row", show_rownames = nrow(hm) <= 80,
           annotation_col = ann_c,
           annotation_colors = list(Group = setNames(GRP_COLS, GRP_DISPLAY)),
           color = colorRampPalette(c("#4E79A7", "white", "#E15759"))(100),
           border_color = NA, cluster_cols = FALSE,
           main = sprintf("Union of significant proteins (n = %d), row z-score", nrow(hm)))
  dev.off()
  cat(sprintf("   heatmap written (%d proteins)\n", nrow(hm)))
} else {
  cat("   heatmap skipped: fewer than 2 significant proteins\n")
}


## =============================================================================
## 14. FUNCTIONAL ENRICHMENT
## =============================================================================
## Three deliberate changes from the vendor enrichment, each answering a specific
## element of the reviewer's criticism:
##
##  (1) BACKGROUND = the tested proteome (all proteins that entered statistics),
##      NOT the mouse genome. This is the change that actually removes the
##      artefact risk: only part of the proteome is detectable by LC-MS, and that
##      part is biased toward abundant proteins, so a genome background returns
##      ribosome/spliceosome/metabolism terms regardless of biology.
##
##  (2) INPUT LIST = the SAME list reported as differentially abundant
##      (d$significant: BH adj p < ADJP_CUTOFF AND |log2FC| >= LFC_PRIMARY).
##      The reviewer's criticism was that pathway results rested on an
##      insufficiently filtered foundational list, so pathways are tested on
##      exactly the proteins the tables and text name - one list, no split.
##      Cost, stated in the methods rather than worked around: smaller lists
##      give less over-representation power, and a direction subset with fewer
##      than 5 symbols is not tested at all.
##
##  (3) BH-ADJUSTED enrichment p-values and quantitative dot plots, replacing the
##      vendor's uncorrected p < 0.05 and chord diagrams.
##
##  (4) RESOURCES: GO BP and Reactome, the two resources of the vendor analysis,
##      are both retained; GO MF/CC and KEGG are added. Reactome comes from the
##      local, versioned reactome.db package; KEGG is fetched online at run time,
##      so its status and access date are written to RUN_SUMMARY.csv.

cat("\n== 14. Enrichment ==\n")

## Map to Entrez via org.Mm.eg.db. Symbols are mapped rather than trusting the
## vendor's Entrez column, which is sparsely populated.
## Symbol -> Entrez. Two distinct de-duplications happen here and the methods
## must describe the one that actually runs:
##   - distinct(SYMBOL): one symbol can map to several Entrez ids; the first is
##     kept. This is de-duplication of the MAPPING.
##   - the `rank_by` argument below: several proteins (distinct accessions, e.g.
##     isoforms or a protein group split) can carry the SAME gene symbol. This is
##     de-duplication of the INPUT LIST, and it is the one the methods sentence
##     "duplicates are collapsed to the most significant entry" refers to. It was
##     previously claimed in the methods but never implemented, so the collapse
##     was whatever order the table happened to arrive in. It is implemented now.
map_entrez <- function(symbols, rank_by = NULL) {
  keep <- !is.na(symbols) & symbols != "" & symbols != "NA"
  symbols <- symbols[keep]
  if (!is.null(rank_by)) {
    rank_by <- rank_by[keep]
    ## Most significant entry per symbol: smallest ranking statistic wins.
    ord     <- order(rank_by, na.last = TRUE)
    symbols <- symbols[ord]
  }
  s <- unique(symbols)               # order-preserving: first occurrence survives
  suppressMessages(
    AnnotationDbi::select(org.Mm.eg.db, keys = s, keytype = "SYMBOL", columns = "ENTREZID")
  ) %>% filter(!is.na(ENTREZID)) %>% distinct(SYMBOL, .keep_all = TRUE)
}

univ_map <- map_entrez(anno_f$Gene_Symbol)
universe <- univ_map$ENTREZID
cat(sprintf("   background: %d tested proteins -> %d unique symbols -> %d Entrez IDs\n",
            nrow(anno_f), length(unique(anno_f$Gene_Symbol)), length(universe)))
fact("enrichment_background_n", length(universe), "tested proteome, NOT genome")
## Section 11 quotes the symbol accounting for the tested set.
.sym_na  <- sum(is.na(anno_f$Gene_Symbol) | anno_f$Gene_Symbol == "")
.sym_ok  <- anno_f$Gene_Symbol[!is.na(anno_f$Gene_Symbol) & anno_f$Gene_Symbol != ""]
fact("tested_proteins_without_gene_symbol", .sym_na,
     sprintf("of %d tested proteins", nrow(anno_f)))
fact("tested_unique_gene_symbols", length(unique(.sym_ok)),
     sprintf("%d named proteins collapse to this many distinct symbols (%d duplicated)",
             length(.sym_ok), length(.sym_ok) - length(unique(.sym_ok))))
writeLines(sort(unique(anno_f$Gene_Symbol[!is.na(anno_f$Gene_Symbol)])),
           file.path(DIR$onto_gl, "background_tested_proteome_symbols.txt"))

## enrichGO/enrichKEGG return NULL or an empty result when nothing is enriched;
## both are normal outcomes and must not stop the pipeline.
safe_enrich <- function(expr, label) {
  out <- tryCatch(expr, error = function(e) {
    warning("Enrichment failed for ", label, ": ", conditionMessage(e), call. = FALSE)
    NULL
  })
  if (is.null(out)) return(NULL)
  df <- as.data.frame(out)
  if (nrow(df) == 0) return(NULL)
  list(obj = out, df = df)
}

## KEGG mapping cache for the fallback path (see the KEGG block below).
.KEGG_T2G <- NULL; .KEGG_T2N <- NULL; .KEGG_UNAVAILABLE <- FALSE
## KEGG provenance. KEGG is queried online and changes continuously, so the
## methods can only cite it by access date (and release, where the KEGG REST
## service reports one). Both are captured at the first KEGG query and written to
## RUN_SUMMARY.csv after the loop, together with whether KEGG actually ran.
.KEGG_ACCESS_DATE <- NA_character_; .KEGG_RELEASE <- NA_character_
## Direction subsets too small to test (< 5 symbols), recorded for the methods.
.ENRICH_SKIPPED <- character(0)

enrich_index <- list()
for (cn in names(res_list)) {
  d <- res_list[[cn]]
  cdir <- file.path(DIR$onto, cn)
  dir.create(cdir, recursive = TRUE, showWarnings = FALSE)

  ## Three gene sets per contrast: all significant, and the two directions.
  ## Direction matters here: "up" means higher in THIS contrast's numerator
  ## group, whichever group that is. The folder names are built from the actual
  ## numerator (higher_in_<numerator>), not from a global stiffness/oxygen rule -
  ## for 0-3_N_vs_1-5_A the numerator is the 0.3 mg/ml normoxic member.
  ## Enrichment input = d$significant, the identical rule used by the results
  ## tables, volcano calls and heatmap (adj p AND fold change). Project owner's
  ## decision, replacing an earlier adj-p-only input: pathways are tested on
  ## exactly the proteins reported as differentially abundant.
  ## Consequence to expect when reading outputs: n_sig_proteins in
  ## enrichment_index.csv equals n_adjp_1.5x (all) / n_up_primary / n_down_primary
  ## in the significance summary; n_input is a few lower because unmapped symbols
  ## and duplicates are dropped.
  ## If you change this line, change Methods section 11 note (b) in the same edit.
  sig <- d %>% filter(significant)
  ## Carry the adjusted p alongside each symbol so that where two accessions
  ## share a gene symbol, the collapse keeps the most significant one rather than
  ## an arbitrary one. Matters for the symbol reported, not for set membership.
  sets <- list(
    all  = sig$Gene_Symbol,
    up   = sig$Gene_Symbol[sig$log2FC > 0],
    down = sig$Gene_Symbol[sig$log2FC < 0]
  )
  setp <- list(
    all  = sig$adj_p_BH,
    up   = sig$adj_p_BH[sig$log2FC > 0],
    down = sig$adj_p_BH[sig$log2FC < 0]
  )
  names(sets) <- c("all_significant",
                   paste0("higher_in_", gsub(" ", "", d$numerator[1])),
                   paste0("higher_in_", gsub(" ", "", d$denominator[1])))

  for (si in seq_along(sets)) {
    sn <- names(sets)[si]
    raw_sym <- sets[[si]]; raw_p <- setp[[si]]
    ok      <- !is.na(raw_sym) & raw_sym != "" & raw_sym != "NA"
    raw_sym <- raw_sym[ok]; raw_p <- raw_p[ok]
    syms    <- unique(raw_sym)
    n_dup   <- length(raw_sym) - length(syms)
    writeLines(sort(syms), file.path(DIR$onto_gl, paste0(cn, "__", sn, "_symbols.txt")))
    if (length(syms) < 5) {
      cat(sprintf("   %-16s %-28s %3d symbols - skipped (too few)\n", cn, sn, length(syms)))
      .ENRICH_SKIPPED <- c(.ENRICH_SKIPPED, sprintf("%s/%s (%d)", cn, sn, length(syms)))
      next
    }
    ## Ranked collapse (see map_entrez): duplicate symbols resolve to the entry
    ## with the smallest adjusted p. Both counts reported, as the methods claim.
    emap    <- map_entrez(raw_sym, rank_by = raw_p)
    ids     <- emap$ENTREZID
    n_unmap <- length(syms) - nrow(emap)
    sdir <- file.path(cdir, sn); dir.create(sdir, showWarnings = FALSE)

    for (ont in GO_ONTOLOGIES) {
      r <- safe_enrich(enrichGO(gene = ids, universe = universe, OrgDb = org.Mm.eg.db,
                                ont = ont, keyType = "ENTREZID",
                                pAdjustMethod = "BH", pvalueCutoff = ENRICH_PADJ,
                                qvalueCutoff = 1,
                                minGSSize = ENRICH_MIN_GS, maxGSSize = ENRICH_MAX_GS,
                                readable = TRUE),
                       paste(cn, sn, "GO", ont))
      if (is.null(r)) next
      write.csv(r$df, file.path(sdir, paste0("GO_", ont, ".csv")), row.names = FALSE)
      p <- dotplot(r$obj, showCategory = 15) +
        labs(title = paste0("GO ", ont, ": ", gsub("_", " ", cn)), subtitle = sn) +
        theme_bw(base_size = 8)
      ggsave(file.path(sdir, paste0("GO_", ont, "_dotplot.png")), p,
             width = 7, height = 5.5, dpi = FIG_DPI)
      enrich_index[[length(enrich_index) + 1]] <- data.frame(
        contrast = cn, gene_set = sn, database = paste0("GO_", ont),
        n_sig_proteins = length(raw_sym), n_symbols = length(syms),
        n_duplicate_symbols = n_dup, n_unmapped_symbols = n_unmap,
        n_input = length(ids), n_terms = nrow(r$df),
        top_term = r$df$Description[1], top_padj = signif(r$df$p.adjust[1], 3))
    }

    if (RUN_KEGG) {
      ## enrichKEGG downloads two things: the pathway->gene mapping (needed) and
      ## kegg_category.rda (a cosmetic pathway-class lookup). On networks that
      ## block the second host, the failed download raises a hard error that
      ## discards the ENTIRE result even though the usable mapping was already
      ## fetched - silently producing zero KEGG output on an otherwise clean run.
      ## Fallback: fetch the mapping directly and run generic over-representation
      ## against it, which is the identical hypergeometric test without the
      ## category lookup. .KEGG_T2G caches so this hits the KEGG API once.
      if (is.na(.KEGG_ACCESS_DATE)) {
        .KEGG_ACCESS_DATE <- format(Sys.Date(), "%Y-%m-%d")
        ## "kegg  Release 11x.0+/MM-DD, Mon YY" -> "Release 11x.0+/MM-DD, Mon YY".
        ## Failure leaves NA and never stops the run; the date is still recorded.
        .KEGG_RELEASE <- tryCatch({
          .ki <- grep("Release", readLines("https://rest.kegg.jp/info/kegg", warn = FALSE),
                      value = TRUE)[1]
          if (is.na(.ki)) NA_character_ else trimws(sub("^\\S+\\s+", "", .ki))
        }, error = function(e) NA_character_, warning = function(w) NA_character_)
      }
      r <- safe_enrich(enrichKEGG(gene = ids, universe = universe,
                                  organism = KEGG_ORGANISM, keyType = "kegg",
                                  pAdjustMethod = "BH", pvalueCutoff = ENRICH_PADJ,
                                  qvalueCutoff = 1,
                                  minGSSize = ENRICH_MIN_GS, maxGSSize = ENRICH_MAX_GS),
                       paste(cn, sn, "KEGG"))
      if (is.null(r) && !isTRUE(.KEGG_UNAVAILABLE)) {
        if (is.null(.KEGG_T2G)) {
          kg <- tryCatch(clusterProfiler::download_KEGG(KEGG_ORGANISM),
                         error = function(e) NULL)
          if (is.null(kg)) {
            .KEGG_UNAVAILABLE <<- TRUE
            warning("KEGG pathway mapping could not be downloaded; KEGG ",
                    "enrichment skipped for this run. GO results are unaffected.",
                    call. = FALSE)
          } else {
            .KEGG_T2G  <<- kg$KEGGPATHID2EXTID
            .KEGG_T2N  <<- kg$KEGGPATHID2NAME
            cat("   (KEGG via direct mapping fallback; category lookup unavailable)\n")
          }
        }
        if (!is.null(.KEGG_T2G))
          r <- safe_enrich(enricher(gene = ids, universe = universe,
                                    TERM2GENE = .KEGG_T2G, TERM2NAME = .KEGG_T2N,
                                    pAdjustMethod = "BH", pvalueCutoff = ENRICH_PADJ,
                                    qvalueCutoff = 1,
                                    minGSSize = ENRICH_MIN_GS,
                                    maxGSSize = ENRICH_MAX_GS),
                           paste(cn, sn, "KEGG-fallback"))
      }
      if (!is.null(r)) {
        write.csv(r$df, file.path(sdir, "KEGG.csv"), row.names = FALSE)
        p <- dotplot(r$obj, showCategory = 15) +
          labs(title = paste0("KEGG: ", gsub("_", " ", cn)), subtitle = sn) +
          theme_bw(base_size = 8)
        ggsave(file.path(sdir, "KEGG_dotplot.png"), p, width = 7, height = 5.5, dpi = FIG_DPI)
        enrich_index[[length(enrich_index) + 1]] <- data.frame(
          contrast = cn, gene_set = sn, database = "KEGG",
          n_sig_proteins = length(raw_sym), n_symbols = length(syms),
          n_duplicate_symbols = n_dup, n_unmapped_symbols = n_unmap,
          n_input = length(ids), n_terms = nrow(r$df),
          top_term = r$df$Description[1], top_padj = signif(r$df$p.adjust[1], 3))
      }
    }

    ## Reactome: the pathway resource of the vendor analysis, retained here with
    ## the same background, BH correction, size bounds and cut-off as GO/KEGG.
    ## ReactomePA reads the local reactome.db package (no download), keyed on
    ## Entrez ids - the same `ids`/`universe` vectors GO and KEGG receive.
    ## qvalueCutoff = 1 for the same reason as above: selection is on BH p only.
    r <- safe_enrich(ReactomePA::enrichPathway(
                       gene = ids, universe = universe, organism = REACTOME_ORGANISM,
                       pAdjustMethod = "BH", pvalueCutoff = ENRICH_PADJ,
                       qvalueCutoff = 1,
                       minGSSize = ENRICH_MIN_GS, maxGSSize = ENRICH_MAX_GS),
                     paste(cn, sn, "Reactome"))
    if (!is.null(r)) {
      write.csv(r$df, file.path(sdir, "Reactome.csv"), row.names = FALSE)
      p <- dotplot(r$obj, showCategory = 15) +
        labs(title = paste0("Reactome: ", gsub("_", " ", cn)), subtitle = sn) +
        theme_bw(base_size = 8)
      ggsave(file.path(sdir, "Reactome_dotplot.png"), p, width = 7, height = 5.5, dpi = FIG_DPI)
      enrich_index[[length(enrich_index) + 1]] <- data.frame(
        contrast = cn, gene_set = sn, database = "Reactome",
        n_sig_proteins = length(raw_sym), n_symbols = length(syms),
        n_duplicate_symbols = n_dup, n_unmapped_symbols = n_unmap,
        n_input = length(ids), n_terms = nrow(r$df),
        top_term = r$df$Description[1], top_padj = signif(r$df$p.adjust[1], 3))
    }
    cat(sprintf("   %-16s %-28s %3d proteins -> %3d symbols (%d dup) -> %3d mapped\n",
                cn, sn, length(raw_sym), length(syms), n_dup, length(ids)))
  }
}

if (length(enrich_index)) {
  ei <- bind_rows(enrich_index)
  write.csv(ei, file.path(DIR$onto, "enrichment_index.csv"), row.names = FALSE)
  cat(sprintf("   enrichment index: %d result sets\n", nrow(ei)))
} else {
  cat("   NOTE: no enriched terms anywhere. Gene lists are still in 00_gene_lists/.\n")
}

## Enrichment provenance for the methods. KEGG status is registered explicitly so
## a run in which KEGG silently failed cannot be written up as if it had run.
fact("enrichment_input_rule", "significant (adj p AND |log2FC| >= LFC_PRIMARY)",
     "same list as the differential-abundance tables")
fact("enrichment_sets_skipped_lt5", length(.ENRICH_SKIPPED),
     paste(.ENRICH_SKIPPED, collapse = "; "))
fact("kegg_status", if (!RUN_KEGG) "not run (RUN_KEGG = FALSE)" else
       if (is.na(.KEGG_ACCESS_DATE)) "not queried (no gene set reached 5 symbols)" else
       if (isTRUE(.KEGG_UNAVAILABLE)) "FAILED - mapping not downloadable; no KEGG results" else
       if (!is.null(.KEGG_T2G)) "ran via direct-mapping fallback" else "ran via enrichKEGG",
     "cite KEGG in the methods only if this says 'ran'")
fact("kegg_access_date", .KEGG_ACCESS_DATE, "KEGG is queried online on this date")
fact("kegg_release", .KEGG_RELEASE, "from rest.kegg.jp/info/kegg; NA if not reported")
fact("reactome_n_result_sets",
     if (length(enrich_index)) sum(bind_rows(enrich_index)$database == "Reactome") else 0L)


## =============================================================================
## 15. PRESENCE/ABSENCE TABLE FOR PROTEINS EXCLUDED BY THE COMPLETENESS FILTER
## =============================================================================
## Converts the main criticism of complete-case filtering ("you threw away
## biology") into a supplementary table, WITHOUT letting invented numbers into any
## statistical test. These proteins are reported as detected / not detected only -
## never with a fold change or a p-value.

cat("\n== 15. Presence/absence ==\n")
excl <- which(!keep_complete)
if (length(excl)) {
  det_grp <- sapply(GRP_INTERNAL, function(g)
    rowSums(!is.na(mat_raw[excl, samp_group == g, drop = FALSE])))
  colnames(det_grp) <- paste0("n_detected_", GRP_SLUG)
  pa <- data.frame(anno[excl, c("Accession", "Gene_Symbol", "Description",
                                "n_unique", "n_peptides")],
                   det_grp, check.names = FALSE)
  ## "on/off" = present in every replicate of at least one group AND absent from
  ## every replicate of at least one other. These are the candidates worth
  ## looking at; anything else is sporadic detection.
  full <- rowSums(det_grp == N_REP); zero <- rowSums(det_grp == 0)
  pa$pattern <- ifelse(full >= 1 & zero >= 1, "on/off candidate", "sporadic")
  pa <- pa %>% arrange(pattern, desc(n_unique))
  write.csv(pa, file.path(DIR$qual, "excluded_proteins_presence_absence.csv"),
            row.names = FALSE)
  cat(sprintf("   %d excluded proteins: %d on/off candidates, %d sporadic\n",
              nrow(pa), sum(pa$pattern == "on/off candidate"),
              sum(pa$pattern == "sporadic")))
  fact("excluded_by_completeness", nrow(pa))
  fact("on_off_candidates", sum(pa$pattern == "on/off candidate"),
       "reported qualitatively only; never with a fold change")
}


## =============================================================================
## 16. SESSION INFO AND FINAL SUMMARY
## =============================================================================
## sessionInfo is written to disk so the Methods can quote exact package versions
## instead of remembered ones.

cat("\n== 16. Finalising ==\n")
writeLines(capture.output(sessionInfo()), file.path(DIR$logs, "session_info.txt"))

## PROVENANCE STAMP. The methods text must be written against ONE identified run,
## and "the numbers in RUN_SUMMARY.csv" is only a reproducible instruction if the
## file says which run it describes. These fields let any later reader confirm
## that a quoted number and the file in front of them came from the same run.
fact("run_output_folder", basename(DIR$root))
fact("run_input_md5", tryCatch(tools::md5sum(normalizePath(INPUT_XLSX))[[1]],
                               error = function(e) NA_character_),
     "md5 of the input workbook; changes if the vendor re-exports")
fact("run_r_version", R.version.string)
## Register EVERY package version the methods text quotes, not just the four
## analysis packages. The methods previously took readxl / openxlsx / ggplot2 /
## ggrepel / pheatmap / enrichplot from session_info.txt. That is correct only
## while the quoted session_info belongs to the run being cited; a version typed
## from a later session in which one of them moved would be undetectable. With
## them in the register, the same mechanical check that validates every other
## quoted number validates these too.
for (.pk in c("limma", "DEqMS", "clusterProfiler", "org.Mm.eg.db",
              "ReactomePA", "reactome.db",
              "readxl", "openxlsx", "ggplot2", "ggrepel", "pheatmap", "enrichplot"))
  fact(paste0("pkg_", .pk), as.character(packageVersion(.pk)))
## The Reactome release is not the package version: reactome.db records the
## Reactome data release it was built from in its metadata table.
fact("reactome_db_source", tryCatch({
  .md <- AnnotationDbi::metadata(reactome.db::reactome.db)
  .hit <- .md[grepl("release|version|date", .md[[1]], ignore.case = TRUE), , drop = FALSE]
  if (nrow(.hit)) paste(paste0(.hit[[1]], "=", .hit[[2]]), collapse = "; ") else NA_character_
}, error = function(e) NA_character_), "cite this Reactome release in the methods")
fact("figure_dpi", FIG_DPI)
fact("methods_text_written_against", basename(DIR$root),
     "PROTEOMICS_METHODS.txt quotes this run; re-check every number if re-run")

## Rewrite run_facts now that later sections have added to it.
wb2 <- loadWorkbook(file.path(DIR$results, "PID430_DEqMS_all_contrasts.xlsx"))
removeWorksheet(wb2, "run_facts"); addWorksheet(wb2, "run_facts")
writeData(wb2, "run_facts", bind_rows(.FACTS))
saveWorkbook(wb2, file.path(DIR$results, "PID430_DEqMS_all_contrasts.xlsx"), overwrite = TRUE)
write.csv(bind_rows(.FACTS), file.path(DIR$root, "RUN_SUMMARY.csv"), row.names = FALSE)

cat("\n=============================================================\n")
cat("PIPELINE COMPLETE\n")
cat("=============================================================\n")
cat("Output root:", DIR$root, "\n\n")
cat("  01_QC/                      sample sheet, QC figures, filter cascade\n")
cat("  02_processed_data/          normalized matrix, peptide counts\n")
cat("  03_differential_abundance/  results workbook, per-contrast tables, summary\n")
cat("  04_figures/                 volcanoes, heatmap, DEqMS diagnostics\n")
cat("  05_ontology/                GO, KEGG and Reactome per contrast and direction\n")
cat("  06_presence_absence/        qualitative table for excluded proteins\n")
cat("  07_logs/                    run log, session info\n\n")
print(summ_df[, c("contrast", "n_adjp_only", "n_adjp_1.5x", "n_adjp_2x")])
cat("\nREMINDER: positive log2FC = higher in THAT CONTRAST'S numerator group.\n")
cat("There is no single global rule: for 0-3_N_vs_1-5_A the numerator is the\n")
cat("0.3 mg/ml normoxic member. Read the numerator column of each table.\n")
cat("Quote counts from RUN_SUMMARY.csv, not from any earlier draft.\n")

## =============================================================================
## END
## =============================================================================
