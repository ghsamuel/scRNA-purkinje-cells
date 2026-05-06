#!/usr/bin/env Rscript
# Load and QC filter Purkinje scRNA-seq data
# Data: GSE256438 - E16.5 (2 reps) and E18.5 cerebellum
# Jini Samuel

library(Seurat)
library(dplyr)
library(ggplot2)

setwd("/core/projects/GAP/GDA/gsamuel/scRNAseq")

# Paths
data_dir <- "results/03_counts/01_cellranger"
dir.create("results/01_qc", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/01_qc", recursive = TRUE, showWarnings = FALSE)

# Load Cell Ranger outputs
e16_rep1_data <- Read10X(file.path(data_dir, "SRR28065512/outs/filtered_feature_bc_matrix"))
e16_rep2_data <- Read10X(file.path(data_dir, "SRR28065511/outs/filtered_feature_bc_matrix"))
e18_data <- Read10X(file.path(data_dir, "SRR28065510/outs/filtered_feature_bc_matrix"))

# Create Seurat objects
e16_rep1 <- CreateSeuratObject(counts = e16_rep1_data, project = "Purkinje",
                                min.cells = 3, min.features = 200)
e16_rep1$stage <- "E16.5"
e16_rep1$sample <- "E16_rep1"
e16_rep1[["percent.mt"]] <- PercentageFeatureSet(e16_rep1, pattern = "^mt-")

e16_rep2 <- CreateSeuratObject(counts = e16_rep2_data, project = "Purkinje",
                                min.cells = 3, min.features = 200)
e16_rep2$stage <- "E16.5"
e16_rep2$sample <- "E16_rep2"
e16_rep2[["percent.mt"]] <- PercentageFeatureSet(e16_rep2, pattern = "^mt-")

e18 <- CreateSeuratObject(counts = e18_data, project = "Purkinje",
                          min.cells = 3, min.features = 200)
e18$stage <- "E18.5"
e18$sample <- "E18"
e18[["percent.mt"]] <- PercentageFeatureSet(e18, pattern = "^mt-")

# QC filtering: 200-5000 genes, >500 UMIs, <10% mito
e16_rep1_filtered <- subset(e16_rep1, 
                             subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & 
                                      nCount_RNA > 500 & percent.mt < 10)

e16_rep2_filtered <- subset(e16_rep2, 
                             subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & 
                                      nCount_RNA > 500 & percent.mt < 10)

e18_filtered <- subset(e18, 
                       subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & 
                                nCount_RNA > 500 & percent.mt < 10)

# Summary
cat("Cells retained:\n")
cat("E16 rep1:", ncol(e16_rep1_filtered), "/", ncol(e16_rep1), "\n")
cat("E16 rep2:", ncol(e16_rep2_filtered), "/", ncol(e16_rep2), "\n")
cat("E18:", ncol(e18_filtered), "/", ncol(e18), "\n")

# QC plots
p1 <- VlnPlot(e16_rep1_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              ncol = 3, pt.size = 0)
p2 <- VlnPlot(e16_rep2_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              ncol = 3, pt.size = 0)
p3 <- VlnPlot(e18_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              ncol = 3, pt.size = 0)

ggsave("figures/01_qc/qc_violin_e16_rep1.png", p1, width = 12, height = 6)
ggsave("figures/01_qc/qc_violin_e16_rep2.png", p2, width = 12, height = 6)
ggsave("figures/01_qc/qc_violin_e18.png", p3, width = 12, height = 6)

# Save
saveRDS(e16_rep1_filtered, "results/01_qc/e16_rep1_filtered.rds")
saveRDS(e16_rep2_filtered, "results/01_qc/e16_rep2_filtered.rds")
saveRDS(e18_filtered, "results/01_qc/e18_filtered.rds")

cat("\nQC complete\n")
