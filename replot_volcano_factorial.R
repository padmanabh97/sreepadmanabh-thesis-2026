suppressPackageStartupMessages(library(ggplot2))

output_dir <- Sys.getenv("OUTPUT_DIR", "DESeq2_Factorial_cellmeans")  # same default as cellmeans
ALPHA   <- 0.05
LFC_CUT <- 1
csv_dir <- file.path(output_dir, "01_tables", "CSV")
vol_dir <- file.path(output_dir, "04_figures", "Volcano_rawLFC_format")
dir.create(vol_dir, recursive = TRUE, showWarnings = FALSE)

EFFECTS <- c("Oxygen_A_vs_N", "Stiffness_0.3_vs_1.5", "Dim_3D_vs_2D",
             "DimxOxygen", "DimxStiffness", "StiffnessxOxygen", "DimxStiffxOxygen")

for (eff in EFFECTS) {
  f <- file.path(csv_dir, paste0(eff, "_results.csv"))
  if (!file.exists(f)) stop("Missing: ", f)
  df  <- read.csv(f, stringsAsFactors = FALSE)
  vdf <- df[!is.na(df$padj), ]
  vdf$regulation <- factor(vdf$regulation, levels = c("Down","NS","Up"))
  p_v <- ggplot(vdf, aes(log2FoldChange, -log10(padj), colour = regulation)) +
    geom_point(alpha = 0.6, size = 1.3) +
    scale_colour_manual(values = c(Down = "#2166AC", NS = "grey75", Up = "#B2182B")) +
    geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = -log10(ALPHA),     linetype = "dashed", colour = "grey40") +
    labs(title = eff, x = "log2 fold change (raw)", y = "-log10(padj)") +
    theme_bw(base_size = 13)
  ggsave(file.path(vol_dir, paste0(eff, "_volcano.png")),
         plot = p_v, width = 7, height = 6, dpi = 300)
}