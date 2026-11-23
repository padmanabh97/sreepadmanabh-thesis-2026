# =============================================================================
# DESeq2 ANALYSIS - 2D SAMPLES ONLY (Raw LFC)
# =============================================================================

library(DESeq2)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(pheatmap)

# ====================== SETUP ======================
# Set working directory to Subset analyses folder
setwd("D:/ESoftBio Lab/Omics Analysis/Oxo-Mech 2D vs 3D transcriptome/fragment counts 010726/Subset Analysis 060726")

counts_file       <- "counts.txt"
sample_sheet_file <- "sample_sheet.tsv"
gene_list_dir     <- "D:/ESoftBio Lab/Omics Analysis/Oxo-Mech 2D vs 3D transcriptome/gene lists"

# Create output folder for 2D samples only
output_dir <- "2D samples only"
dir.create(output_dir, showWarnings = FALSE)

dir.create(file.path(output_dir, "CSV"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "Volcano"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "GO"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "Heatmaps"), recursive = TRUE, showWarnings = FALSE)

# ====================== 1. LOAD DATA ======================
countData <- read.delim(counts_file, row.names = 1, skip = 1, check.names = FALSE)

# Clean column names
colnames(countData) <- basename(colnames(countData))
colnames(countData) <- sub("\\.bam$", "", colnames(countData))

sample_info <- read.delim(sample_sheet_file, stringsAsFactors = FALSE)

# ====================== FILTER ONLY 2D SAMPLES ======================
sample_info_2D <- sample_info[sample_info$dim == "2D", ]
valid_samples <- sample_info_2D$sample
count_columns <- intersect(valid_samples, colnames(countData))

countData <- countData[, count_columns, drop = FALSE]
countData <- as.matrix(countData)
mode(countData) <- "numeric"

sampleTable <- sample_info_2D[sample_info_2D$sample %in% count_columns, ]
sampleTable$group <- paste(sampleTable$stiffness, sampleTable$oxygen, sep = "_")
sampleTable$group <- factor(sampleTable$group)
rownames(sampleTable) <- sampleTable$sample

# ====================== 2. DESeq2 ======================
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = sampleTable,
                              design = ~ group)

keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]

dds <- DESeq(dds)
normCounts <- counts(dds, normalized = TRUE)

# ====================== 3. PCA ======================
vsd <- vst(dds, blind = TRUE)
pcaData <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p_pca <- ggplot(pcaData, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3.5) +
  geom_text_repel(aes(label = name), size = 3.8, max.overlaps = 50,
                  show.legend = FALSE, box.padding = 0.4, point.padding = 0.3) +
  scale_color_brewer(palette = "Set2") +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA - 2D Samples Only") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(output_dir, "PCA_with_labels.pdf"), plot = p_pca, width = 12, height = 9, device = "pdf")

# ====================== 4. ALL PAIRWISE COMPARISONS ======================
groups <- levels(dds$group)
comparisons <- combn(groups, 2, simplify = FALSE)

summary_list <- list()

for (comp in comparisons) {
  g1 <- comp[1]
  g2 <- comp[2]
  contrast_name <- paste0(g1, "_vs_", g2)
  
  cat("Processing:", contrast_name, "\n")
  
  res <- results(dds, contrast = c("group", g1, g2), alpha = 0.05)
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  
  # Gene symbols
  res_df$gene_name <- mapIds(
    org.Mm.eg.db,
    keys      = res_df$gene_id,
    column    = "SYMBOL",
    keytype   = "ENSEMBL",
    multiVals = "first"
  )
  res_df$gene_name <- ifelse(is.na(res_df$gene_name), res_df$gene_id, res_df$gene_name)
  
  # Regulation
  res_df$regulation <- "NS"
  res_df$regulation[res_df$log2FoldChange >  1 & res_df$padj < 0.05] <- "Up"
  res_df$regulation[res_df$log2FoldChange < -1 & res_df$padj < 0.05] <- "Down"
  
  # Save CSV
  write.csv(res_df, file.path(output_dir, "CSV", paste0(contrast_name, "_results.csv")), row.names = FALSE)
  
  # Volcano
  res_df$neg_log10_p <- -log10(res_df$pvalue)
  p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = neg_log10_p)) +
    geom_point(aes(color = regulation), size = 1.2, alpha = 0.7) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey70")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.4) +
    labs(title = paste("Volcano -", contrast_name)) +
    theme_minimal(base_size = 13)
  
  ggsave(file.path(output_dir, "Volcano", paste0(contrast_name, "_volcano.png")),
         plot = p_volcano, width = 9, height = 7, dpi = 300)
  
  # ====================== GO ENRICHMENT ======================
  up_genes   <- res_df$gene_id[res_df$regulation == "Up"]
  down_genes <- res_df$gene_id[res_df$regulation == "Down"]
  background_genes <- res_df$gene_id[!is.na(res_df$padj)]
  
  if (length(up_genes) >= 10) {
    ego_up <- enrichGO(
      gene          = up_genes,
      universe      = background_genes,
      OrgDb         = org.Mm.eg.db,
      keyType       = "ENSEMBL",
      ont           = "BP",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2,
      readable      = TRUE
    )
    
    if (!is.null(ego_up) && nrow(as.data.frame(ego_up)) > 0) {
      dotplot(ego_up, showCategory = 15) + 
        labs(title = paste(contrast_name, "- GO UP")) + 
        theme_minimal()
      
      ggsave(file.path(output_dir, "GO", paste0(contrast_name, "_GO_UP.png")), 
             width = 10, height = 8, dpi = 300)
    }
  }
  
  if (length(down_genes) >= 10) {
    ego_down <- enrichGO(
      gene          = down_genes,
      universe      = background_genes,
      OrgDb         = org.Mm.eg.db,
      keyType       = "ENSEMBL",
      ont           = "BP",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2,
      readable      = TRUE
    )
    
    if (!is.null(ego_down) && nrow(as.data.frame(ego_down)) > 0) {
      dotplot(ego_down, showCategory = 15) + 
        labs(title = paste(contrast_name, "- GO DOWN")) + 
        theme_minimal()
      
      ggsave(file.path(output_dir, "GO", paste0(contrast_name, "_GO_DOWN.png")), 
             width = 10, height = 8, dpi = 300)
    }
  }
  
  # Summary
  summary_list[[contrast_name]] <- data.frame(
    Contrast = contrast_name,
    Upregulated = sum(res_df$regulation == "Up"),
    Downregulated = sum(res_df$regulation == "Down"),
    Total_DE = sum(res_df$regulation != "NS"),
    Total_Genes = nrow(res_df)
  )
}

# ====================== SUMMARY ======================
summary_df <- bind_rows(summary_list) %>% arrange(desc(Total_DE))
write.csv(summary_df, file.path(output_dir, "DE_Summary_2D_Only.csv"), row.names = FALSE)

cat("\n=== 2D-ONLY ANALYSIS COMPLETE ===\n")
cat("All results saved in:", output_dir, "\n")