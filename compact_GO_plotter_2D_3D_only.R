library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.Mm.eg.db)

base_dir <- "D:/ESoftBio Lab/Omics Analysis/Oxo-Mech 2D vs 3D transcriptome/fragment counts 010726/Subset Analysis 060726/3D samples only"
# 3D: change the last folder to 3D samples only

csv_dir <- file.path(base_dir, "CSV")
out_dir <- file.path(base_dir, "GO_compact") 
if (!dir.exists(base_dir)) stop("base_dir does not exist: ", base_dir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

csv_files <- list.files(csv_dir, pattern = "_results\\.csv$", full.names = TRUE)
if (length(csv_files) == 0) stop("No *_results.csv files in: ", csv_dir)

for (f in csv_files) {
  contrast_name <- sub("_results\\.csv$", "", basename(f))
  res_df <- read.csv(f, stringsAsFactors = FALSE)
  
  need <- c("gene_id", "padj", "regulation")
  if (!all(need %in% names(res_df))) {
    stop("Missing columns in ", basename(f), ": ", paste(setdiff(need, names(res_df)), collapse = ", "))
  }
  
  background_genes <- unique(na.omit(res_df$gene_id[!is.na(res_df$padj)]))
  
  for (dir_lab in c("Up", "Down")) {
    genes <- unique(na.omit(res_df$gene_id[res_df$regulation == dir_lab]))
    if (length(genes) < 10) next
    
    ego <- tryCatch(
      enrichGO(
        gene          = genes,
        universe      = background_genes,
        OrgDb         = org.Mm.eg.db,
        keyType       = "ENSEMBL",
        ont           = "BP",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.2,
        readable      = TRUE
      ),
      error = function(e) { message(contrast_name, " ", dir_lab, ": ", conditionMessage(e)); NULL }
    )
    if (is.null(ego) || nrow(as.data.frame(ego)) == 0) next
    
    p <- dotplot(ego, showCategory = 10) +
      labs(title = paste(contrast_name, "- GO", toupper(dir_lab))) +
      coord_cartesian(clip = "off") +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 8),
        plot.title  = element_text(size = 12, face = "bold"),
        plot.margin = margin(4, 10, 4, 4)
      )
    
    ggsave(file.path(out_dir, paste0(contrast_name, "_GO_", toupper(dir_lab), ".png")),
           plot = p, width = 6, height = 8, dpi = 300)
  }
}