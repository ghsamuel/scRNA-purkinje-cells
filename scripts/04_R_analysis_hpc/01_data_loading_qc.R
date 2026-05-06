#!/usr/bin/env Rscript
# Load Cell Ranger outputs and apply QC filtering
# Jini Samuel

library(Seurat)
library(ggplot2)
library(patchwork)

setwd("/core/projects/GAP/GDA/gsamuel/scRNAseq")

# Create output directories
dir.create("results/04_R_analysis/01_qc", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/04_R_analysis/01_qc", recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# E16.5 REPLICATE 1
# ============================================================================

cat("Processing E16.5 replicate 1...\n")

# Load Cell Ranger output
e16_rep1_data <- Read10X("results/03_counts/01_cellranger/SRR28065512/outs/filtered_feature_bc_matrix")

# Create Seurat object
e16_rep1 <- CreateSeuratObject(counts = e16_rep1_data, 
                                min.cells = 3, 
                                min.features = 200,
                                project = "Purkinje")

# Add metadata
e16_rep1$stage <- "E16.5"
e16_rep1$sample <- "E16_rep1"

# Calculate mitochondrial percentage
e16_rep1[["percent.mt"]] <- PercentageFeatureSet(e16_rep1, pattern = "^mt-")

# QC plots BEFORE filtering
p1 <- VlnPlot(e16_rep1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("figures/04_R_analysis/01_qc/violin_before_qc_e16_rep1.png", p1, width = 12, height = 4)

cat("  Cells before QC:", ncol(e16_rep1), "\n")

# Apply QC filters
e16_rep1_filtered <- subset(e16_rep1, 
                             subset = nFeature_RNA > 200 & 
                                      nFeature_RNA < 5000 & 
                                      nCount_RNA > 500 & 
                                      percent.mt < 10)

cat("  Cells after QC:", ncol(e16_rep1_filtered), "\n")

# QC plots AFTER filtering
p2 <- VlnPlot(e16_rep1_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("figures/04_R_analysis/01_qc/violin_after_qc_e16_rep1.png", p2, width = 12, height = 4)

# Save
saveRDS(e16_rep1_filtered, "results/04_R_analysis/01_qc/e16_rep1_filtered.rds")

# ============================================================================
# E16.5 REPLICATE 2
# ============================================================================

cat("\nProcessing E16.5 replicate 2...\n")

e16_rep2_data <- Read10X("results/03_counts/01_cellranger/SRR28065511/outs/filtered_feature_bc_matrix")

e16_rep2 <- CreateSeuratObject(counts = e16_rep2_data, 
                                min.cells = 3, 
                                min.features = 200,
                                project = "Purkinje")

e16_rep2$stage <- "E16.5"
e16_rep2$sample <- "E16_rep2"
e16_rep2[["percent.mt"]] <- PercentageFeatureSet(e16_rep2, pattern = "^mt-")

# QC plots BEFORE filtering
p3 <- VlnPlot(e16_rep2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("figures/04_R_analysis/01_qc/violin_before_qc_e16_rep2.png", p3, width = 12, height = 4)

cat("  Cells before QC:", ncol(e16_rep2), "\n")

e16_rep2_filtered <- subset(e16_rep2, 
                             subset = nFeature_RNA > 200 & 
                                      nFeature_RNA < 5000 & 
                                      nCount_RNA > 500 & 
                                      percent.mt < 10)

cat("  Cells after QC:", ncol(e16_rep2_filtered), "\n")

# QC plots AFTER filtering
p4 <- VlnPlot(e16_rep2_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("figures/04_R_analysis/01_qc/violin_after_qc_e16_rep2.png", p4, width = 12, height = 4)

saveRDS(e16_rep2_filtered, "results/04_R_analysis/01_qc/e16_rep2_filtered.rds")

# ============================================================================
# E18.5
# ============================================================================

cat("\nProcessing E18.5...\n")

e18_data <- Read10X("results/03_counts/01_cellranger/SRR28065510/outs/filtered_feature_bc_matrix")

e18 <- CreateSeuratObject(counts = e18_data, 
                          min.cells = 3, 
                          min.features = 200,
                          project = "Purkinje")

e18$stage <- "E18.5"
e18$sample <- "E18"
e18[["percent.mt"]] <- PercentageFeatureSet(e18, pattern = "^mt-")

# QC plots BEFORE filtering
p5 <- VlnPlot(e18, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("figures/04_R_analysis/01_qc/violin_before_qc_e18.png", p5, width = 12, height = 4)

cat("  Cells before QC:", ncol(e18), "\n")

e18_filtered <- subset(e18, 
                       subset = nFeature_RNA > 200 & 
                                nFeature_RNA < 5000 & 
                                nCount_RNA > 500 & 
                                percent.mt < 10)

cat("  Cells after QC:", ncol(e18_filtered), "\n")

# QC plots AFTER filtering
p6 <- VlnPlot(e18_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("figures/04_R_analysis/01_qc/violin_after_qc_e18.png", p6, width = 12, height = 4)

saveRDS(e18_filtered, "results/04_R_analysis/01_qc/e18_filtered.rds")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n=== QC SUMMARY ===\n")
cat("E16.5 rep1:", ncol(e16_rep1_filtered), "cells retained\n")
cat("E16.5 rep2:", ncol(e16_rep2_filtered), "cells retained\n")
cat("E18.5:", ncol(e18_filtered), "cells retained\n")
cat("Total:", ncol(e16_rep1_filtered) + ncol(e16_rep2_filtered) + ncol(e18_filtered), "cells\n")

# Write summary table
summary_df <- data.frame(
  Sample = c("E16_rep1", "E16_rep2", "E18"),
  Cells_Before = c(ncol(e16_rep1), ncol(e16_rep2), ncol(e18)),
  Cells_After = c(ncol(e16_rep1_filtered), ncol(e16_rep2_filtered), ncol(e18_filtered)),
  Percent_Retained = round(c(ncol(e16_rep1_filtered)/ncol(e16_rep1), 
                               ncol(e16_rep2_filtered)/ncol(e16_rep2), 
                               ncol(e18_filtered)/ncol(e18)) * 100, 1)
)

write.csv(summary_df, "results/04_R_analysis/01_qc/qc_summary.csv", row.names = FALSE)

cat("\nQC complete! Files saved to results/04_R_analysis/01_qc/ and figures/04_R_analysis/01_qc/\n")
