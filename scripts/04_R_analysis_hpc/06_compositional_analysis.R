#!/usr/bin/env Rscript
# Compositional analysis of Purkinje subtypes
# Chi-square test and proportion barplots

library(Seurat)
library(dplyr)
library(ggplot2)

setwd("/core/projects/GAP/GDA/gsamuel/scRNAseq")

# Create output directory
dir.create("results/04_R_analysis/06_compositional_analysis", recursive = TRUE, showWarnings = FALSE)
dir.create("results/04_R_analysis/06_compositional_analysis/figures", recursive = TRUE, showWarnings = FALSE)

# Load Purkinje subtypes from script 05
purkinje <- readRDS("results/04_R_analysis/05_purkinje_subtypes/purkinje_subtypes.rds")

cat("Loaded", ncol(purkinje), "Purkinje cells\n")
cat("Subtypes:", length(unique(purkinje$seurat_clusters)), "\n")
cat("E16.5 cells:", sum(purkinje$stage == "E16.5"), "\n")
cat("E18.5 cells:", sum(purkinje$stage == "E18.5"), "\n\n")

################################################################################
# CHI-SQUARE TEST - Subtype proportions by stage
################################################################################

# Distribution across stages
stage_dist <- table(purkinje$seurat_clusters, purkinje$stage)
write.csv(stage_dist, "results/04_R_analysis/06_compositional_analysis/subtype_by_stage.csv")

# Chi-square test
chisq_result <- chisq.test(stage_dist)
cat("Chi-square test for subtype composition differences:\n")
cat("X-squared =", chisq_result$statistic, "\n")
cat("p-value =", chisq_result$p.value, "\n\n")

# Save test results
sink("results/04_R_analysis/06_compositional_analysis/chisq_test_results.txt")
cat("Chi-square test: Purkinje subtype proportions E16.5 vs E18.5\n")
cat("================================================\n\n")
cat("Sample sizes:\n")
cat("  E16.5: n =", sum(purkinje$stage == "E16.5"), "cells (2 replicates)\n")
cat("  E18.5: n =", sum(purkinje$stage == "E18.5"), "cells (1 replicate)\n\n")
print(chisq_result)
cat("\n\nContingency table:\n")
print(stage_dist)
sink()

################################################################################
# BARPLOT - Normalized by stage
################################################################################

# Calculate proportions within each stage
stage_df <- as.data.frame(stage_dist)
colnames(stage_df) <- c("Subtype", "Stage", "Count")

stage_df <- stage_df %>%
  group_by(Stage) %>%
  mutate(Proportion = Count / sum(Count))

# Side-by-side barplot
p_stage <- ggplot(stage_df, aes(x = Subtype, y = Proportion, fill = Stage)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Purkinje Subtype", 
       y = "Proportion within stage", 
       title = paste0("Purkinje Subtype Distribution by Stage\n",
                     "Chi-square p = ", format.pval(chisq_result$p.value, digits = 3))) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold")
  )

ggsave("results/04_R_analysis/06_compositional_analysis/figures/subtype_by_stage.pdf", 
       p_stage, width = 12, height = 6)

################################################################################
# BARPLOT - PC vs PC.immature composition
################################################################################

# Distribution by cell type
celltype_dist <- table(purkinje$seurat_clusters, purkinje$cell_type)
write.csv(celltype_dist, "results/04_R_analysis/06_compositional_analysis/subtype_by_celltype.csv")

# Stacked barplot by cell type
celltype_df <- as.data.frame(celltype_dist)
colnames(celltype_df) <- c("Subtype", "CellType", "Count")

p_celltype <- ggplot(celltype_df, aes(x = Subtype, y = Count, fill = CellType)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Purkinje Subtype", 
       y = "Proportion", 
       title = "PC vs PC.immature Composition by Subtype") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold")
  )

ggsave("results/04_R_analysis/06_compositional_analysis/figures/subtype_by_celltype.pdf", 
       p_celltype, width = 10, height = 6)

cat("\n========================================\n")
cat("Compositional analysis complete\n")
cat("========================================\n")
cat("Outputs saved to: results/04_R_analysis/06_compositional_analysis/\n")
cat("========================================\n")
