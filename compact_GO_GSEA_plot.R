library(ggplot2)
library(dplyr)

root <- "D:/ESoftBio Lab/Omics Analysis/Oxo-Mech 2D vs 3D transcriptome/Factorial Design/DESeq2_Factorial_cellmeans"

ora_in   <- file.path(root, "02_GO_ORA")
gsea_in  <- file.path(root, "03_GSEA")
ora_out  <- file.path(root, "02_GO_ORA_new")
gsea_out <- file.path(root, "03_GSEA_new")
dir.create(ora_out,  recursive = TRUE, showWarnings = FALSE)
dir.create(gsea_out, recursive = TRUE, showWarnings = FALSE)

parse_ratio <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  vapply(as.character(x), function(s) {
    p <- strsplit(s, "/", fixed = TRUE)[[1]]
    if (length(p) == 2) as.numeric(p[1]) / as.numeric(p[2])
    else suppressWarnings(as.numeric(s))
  }, numeric(1))
}

ora_files  <- list.files(ora_in,  pattern = "_GO_BP\\.csv$",      full.names = TRUE)
gsea_files <- list.files(gsea_in, pattern = "_GSEA_GO_BP\\.csv$", full.names = TRUE)
message("ORA csvs: ", length(ora_files), " | GSEA csvs: ", length(gsea_files))
if (length(ora_files) == 0 && length(gsea_files) == 0)
  stop("No enrichment CSVs found. Check 'root'.")

# ---------- GO BP ORA ----------
for (f in ora_files) {
  tab <- read.csv(f, stringsAsFactors = FALSE)
  need <- c("Description", "GeneRatio", "p.adjust", "Count")
  if (!all(need %in% names(tab)) || nrow(tab) == 0) next
  
  tab <- tab %>%
    mutate(GeneRatio_num = parse_ratio(GeneRatio)) %>%
    filter(is.finite(GeneRatio_num), is.finite(p.adjust)) %>%
    arrange(p.adjust) %>%
    slice_head(n = 10) %>%
    mutate(Description = factor(Description, levels = rev(unique(Description))))
  
  if (nrow(tab) == 0) next
  
  tag <- sub("_GO_BP\\.csv$", "", basename(f))
  p <- ggplot(tab, aes(GeneRatio_num, Description, size = Count, colour = p.adjust)) +
    geom_point() +
    scale_colour_gradient(low = "#B2182B", high = "#2166AC") +
    scale_x_continuous(expand = expansion(mult = c(0.10, 0.12))) +
    coord_cartesian(clip = "off") +
    labs(title = paste(tag, "- GO BP ORA"), x = "GeneRatio", y = NULL,
         size = "Count", colour = "p.adjust") +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y = element_text(size = 8),
      plot.title  = element_text(size = 11, face = "bold"),
      plot.margin = margin(4, 10, 4, 4)
    )
  ggsave(file.path(ora_out, paste0(tag, "_GO_BP.png")),
         plot = p, width = 6, height = 8, dpi = 300)
}

# ---------- GSEA ----------
for (f in gsea_files) {
  tab <- read.csv(f, stringsAsFactors = FALSE)
  need <- c("Description", "NES", "p.adjust", "setSize", "core_enrichment")
  if (!all(need %in% names(tab)) || nrow(tab) == 0) next
  
  tab <- tab %>%
    filter(is.finite(NES), is.finite(p.adjust), is.finite(setSize), setSize > 0) %>%
    mutate(
      .sign = factor(ifelse(NES >= 0, "activated", "repressed"),
                     levels = c("activated", "repressed")),
      Count = lengths(strsplit(as.character(core_enrichment), "/", fixed = TRUE)),
      GeneRatio = Count / setSize
    ) %>%
    filter(is.finite(GeneRatio)) %>%
    group_by(.sign) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice_head(n = 10) %>%
    ungroup() %>%
    arrange(.sign, p.adjust) %>%
    mutate(Description = factor(Description, levels = rev(unique(as.character(Description)))))
  
  if (nrow(tab) == 0) next
  
  tag <- sub("_GSEA_GO_BP\\.csv$", "", basename(f))
  p <- ggplot(tab, aes(GeneRatio, Description, size = Count, colour = p.adjust)) +
    geom_point() +
    facet_grid(. ~ .sign) +
    scale_colour_gradient(low = "#B2182B", high = "#2166AC") +
    scale_x_continuous(expand = expansion(mult = c(0.10, 0.12))) +
    coord_cartesian(clip = "off") +
    ggtitle(paste(tag, "- GSEA GO BP")) +
    theme_bw() +
    theme(
      axis.text.y     = element_text(size = 8),
      strip.text      = element_text(size = 20),
      panel.spacing.x = unit(0.3, "lines"),
      plot.margin     = margin(4, 6, 4, 4)
    )
  ggsave(file.path(gsea_out, paste0(tag, "_GSEA_dotplot.png")),
         plot = p, width = 6, height = 8, dpi = 300)
}