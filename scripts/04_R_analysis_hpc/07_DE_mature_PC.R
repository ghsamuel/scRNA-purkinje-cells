#!/usr/bin/env Rscript
# Single-cell DE analysis: Mature Purkinje cells E16.5 vs E18.5
# Uses FindMarkers (Wilcoxon test) since only 1 E18.5 replicate

library(Seurat)
library(dplyr)

setwd("/core/projects/GAP/GDA/gsamuel/scRNAseq")

# Create output directory
dir.create("results/04_R_analysis/07_DE_mature_PC", recursive = TRUE, showWarnings = FALSE)

# Load Purkinje subtypes
purkinje <- readRDS("results/04_R_analysis/05_purkinje_subtypes/purkinje_subtypes.rds")

cat("Loaded", ncol(purkinje), "Purkinje cells\n")

################################################################################
# SUBSET TO MATURE PURKINJE ONLY
################################################################################

# Extract only mature PC (Pcp2+)
mature_pc <- subset(purkinje, cell_type == "PC")

cat("\nExtracted", ncol(mature_pc), "mature Purkinje cells\n")
cat("E16.5:", sum(mature_pc$stage == "E16.5"), "cells\n")
cat("E18.5:", sum(mature_pc$stage == "E18.5"), "cells\n\n")

################################################################################
# DIFFERENTIAL EXPRESSION E16.5 vs E18.5
################################################################################

cat("Running differential expression (Wilcoxon test)...\n")

# Set identity to stage
Idents(mature_pc) <- "stage"

# Find markers: E18.5 vs E16.5 (positive log2FC = upregulated in E18.5)
de_results <- FindMarkers(
  mature_pc,
  ident.1 = "E18.5",
  ident.2 = "E16.5",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.1,
  verbose = FALSE
)

# Add gene names and sort
de_results <- de_results %>%
  tibble::rownames_to_column("gene") %>%
  arrange(p_val_adj)

cat("Found", nrow(de_results), "genes tested\n")
cat("Significant genes (p_val_adj < 0.05):", sum(de_results$p_val_adj < 0.05), "\n")
cat("  Upregulated in E18.5:", sum(de_results$p_val_adj < 0.05 & de_results$avg_log2FC > 0), "\n")
cat("  Downregulated in E18.5:", sum(de_results$p_val_adj < 0.05 & de_results$avg_log2FC < 0), "\n\n")

# Save all results
write.csv(de_results, 
          "results/04_R_analysis/07_DE_mature_PC/mature_PC_E18.5_vs_E16.5_all_genes.csv",
          row.names = FALSE)

# Save significant genes only
de_sig <- de_results %>% filter(p_val_adj < 0.05)
write.csv(de_sig, 
          "results/04_R_analysis/07_DE_mature_PC/mature_PC_E18.5_vs_E16.5_significant.csv",
          row.names = FALSE)

################################################################################
# VOLCANO PLOT
################################################################################

cat("Creating volcano plot...\n")

# Add significance category
de_results <- de_results %>%
  mutate(
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0.5 ~ "Up in E18.5",
      p_val_adj < 0.05 & avg_log2FC < -0.5 ~ "Down in E18.5",
      p_val_adj < 0.05 ~ "Significant",
      TRUE ~ "NS"
    )
  )

library(ggplot2)

p_volcano <- ggplot(de_results, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "Up in E18.5" = "#E74C3C",
    "Down in E18.5" = "#3498DB",
    "Significant" = "gray50",
    "NS" = "gray80"
  )) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray30") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray30") +
  labs(
    x = "log2 Fold Change (E18.5 / E16.5)",
    y = "-log10(adjusted p-value)",
    title = "Differential Expression: Mature Purkinje Cells E18.5 vs E16.5",
    color = NULL
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "top"
  )

ggsave("results/04_R_analysis/07_DE_mature_PC/volcano_plot.pdf",
       p_volcano, width = 10, height = 8)

################################################################################
# TOP DE GENES HEATMAP
################################################################################

cat("Creating heatmap of top DE genes...\n")

# Get top 50 DE genes (25 up, 25 down)
top_up <- de_sig %>% 
  filter(avg_log2FC > 0) %>% 
  arrange(p_val_adj) %>% 
  head(25) %>% 
  pull(gene)

top_down <- de_sig %>% 
  filter(avg_log2FC < 0) %>% 
  arrange(p_val_adj) %>% 
  head(25) %>% 
  pull(gene)

top_genes <- c(top_down, top_up)

if (length(top_genes) > 0) {
  p_heatmap <- DoHeatmap(
    mature_pc,
    features = top_genes,
    group.by = "stage",
    disp.min = -2,
    disp.max = 2
  ) +
    scale_fill_gradientn(colors = c("blue", "white", "red")) +
    ggtitle("Top DE Genes: Mature PC E18.5 vs E16.5")
  
  ggsave("results/04_R_analysis/07_DE_mature_PC/heatmap_top_DE_genes.pdf",
         p_heatmap, width = 10, height = 12)
}

################################################################################
# FEATURE PLOTS FOR KEY GENES
################################################################################

cat("Creating feature plots for key DE genes...\n")

# Top 6 upregulated
if (length(top_up) >= 6) {
  p_up <- FeaturePlot(mature_pc, features = top_up[1:6], ncol = 3, split.by = "stage")
  ggsave("results/04_R_analysis/07_DE_mature_PC/featureplot_top_upregulated.pdf",
         p_up, width = 18, height = 12)
}

# Top 6 downregulated
if (length(top_down) >= 6) {
  p_down <- FeaturePlot(mature_pc, features = top_down[1:6], ncol = 3, split.by = "stage")
  ggsave("results/04_R_analysis/07_DE_mature_PC/featureplot_top_downregulated.pdf",
         p_down, width = 18, height = 12)
}

cat("\n========================================\n")
cat("Differential expression analysis complete\n")
cat("========================================\n")
cat("Mature PC cells analyzed:\n")
cat("  E16.5: n =", sum(mature_pc$stage == "E16.5"), "\n")
cat("  E18.5: n =", sum(mature_pc$stage == "E18.5"), "\n")
cat("Significant genes:", sum(de_results$p_val_adj < 0.05), "\n")
cat("Outputs saved to: results/04_R_analysis/07_DE_mature_PC/\n")
cat("========================================\n")
