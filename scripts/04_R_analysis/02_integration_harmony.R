#!/usr/bin/env Rscript
# Merge samples and integrate with Harmony
# Fixes batch effect between E16.5 replicates


library(Seurat)
library(harmony)
library(ggplot2)

setwd("/core/projects/GAP/GDA/gsamuel/scRNAseq")

dir.create("results/02_integration", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/02_integration", recursive = TRUE, showWarnings = FALSE)

# Load filtered objects
e16_rep1 <- readRDS("results/01_qc/e16_rep1_filtered.rds")
e16_rep2 <- readRDS("results/01_qc/e16_rep2_filtered.rds")
e18 <- readRDS("results/01_qc/e18_filtered.rds")

# Merge
merged <- merge(e16_rep1, y = c(e16_rep2, e18),
                add.cell.ids = c("E16_rep1", "E16_rep2", "E18"))
merged$sample <- sub("_[^_]+$", "", colnames(merged))

# Standard workflow (pre-integration)
merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged, nfeatures = 2000)
merged <- ScaleData(merged)
merged <- RunPCA(merged, npcs = 30, verbose = FALSE)
merged <- FindNeighbors(merged, dims = 1:30)
merged <- FindClusters(merged, resolution = 0.5, verbose = FALSE)
merged <- RunUMAP(merged, dims = 1:30, verbose = FALSE)

# Check batch effect
p1 <- DimPlot(merged, group.by = "sample") + ggtitle("Pre-integration")
ggsave("figures/02_integration/umap_pre_integration_by_sample.png", p1, width = 10, height = 8)

# Harmony integration - corrects for sample-level batch effects
integrated <- RunHarmony(merged, group.by.vars = "sample", 
                         reduction = "pca", dims = 1:30, verbose = FALSE)

# Re-cluster on integrated embeddings
integrated <- FindNeighbors(integrated, reduction = "harmony", dims = 1:30)
integrated <- FindClusters(integrated, resolution = 0.5, verbose = FALSE)
integrated <- RunUMAP(integrated, reduction = "harmony", dims = 1:30, verbose = FALSE)

# Check integration quality
p2 <- DimPlot(integrated, group.by = "sample") + ggtitle("Post-integration")
p3 <- DimPlot(integrated, group.by = "stage") + ggtitle("By stage")
p4 <- DimPlot(integrated, label = TRUE) + ggtitle("Clusters")

ggsave("figures/02_integration/umap_post_integration_by_sample.png", p2, width = 10, height = 8)
ggsave("figures/02_integration/umap_post_integration_by_stage.png", p3, width = 10, height = 8)
ggsave("figures/02_integration/umap_post_integration_clusters.png", p4, width = 10, height = 8)

# Check for sample-specific clusters
cluster_dist <- table(integrated$seurat_clusters, integrated$sample)
write.csv(cluster_dist, "results/02_integration/cluster_sample_distribution.csv")

saveRDS(integrated, "results/02_integration/purkinje_integrated.rds")

cat("Integration complete -", ncol(integrated), "cells,", 
    length(unique(integrated$seurat_clusters)), "clusters\n")
