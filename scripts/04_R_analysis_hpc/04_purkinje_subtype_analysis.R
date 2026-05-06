#!/usr/bin/env Rscript
# High-resolution Purkinje cell subtyping
# Looking for Foxp1/Foxp2/Foxp4 expression patterns
# Jini Samuel

library(Seurat)
library(dplyr)
library(ggplot2)

setwd("/core/projects/GAP/GDA/gsamuel/scRNAseq")

# Create output directories
dir.create("results/04_purkinje_subtypes", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/04_purkinje_subtypes", recursive = TRUE, showWarnings = FALSE)

# Load and extract Purkinje cells
obj <- readRDS("results/03_annotation/purkinje_annotated.rds")
purkinje <- subset(obj, cell_type_integrated %in% c("Purkinje_cells", "Purkinje_subtype"))

cat("Extracted", ncol(purkinje), "Purkinje cells\n")

# Re-cluster at high resolution to find subtypes
purkinje <- NormalizeData(purkinje)
purkinje <- FindVariableFeatures(purkinje, nfeatures = 2000)
purkinje <- ScaleData(purkinje)
purkinje <- RunPCA(purkinje, npcs = 20, verbose = FALSE)
purkinje <- FindNeighbors(purkinje, dims = 1:15)
purkinje <- FindClusters(purkinje, resolution = 1.0, verbose = FALSE)
purkinje <- RunUMAP(purkinje, dims = 1:15, verbose = FALSE)

cat("Found", length(unique(purkinje$seurat_clusters)), "subtypes\n")

# Plot subtypes
p1 <- DimPlot(purkinje, label = TRUE, label.size = 4) + ggtitle("Purkinje Subtypes")
ggsave("figures/04_purkinje_subtypes/purkinje_subtypes_umap.png", p1, width = 10, height = 8)

# Foxp expression - key subtype markers from paper
foxp_genes <- c("Foxp1", "Foxp2", "Foxp4")

p2 <- DotPlot(purkinje, features = foxp_genes) + coord_flip()
ggsave("figures/04_purkinje_subtypes/foxp_dotplot.png", p2, width = 10, height = 5)

p3 <- FeaturePlot(purkinje, features = foxp_genes, ncol = 3, order = TRUE)
ggsave("figures/04_purkinje_subtypes/foxp_featureplots.png", p3, width = 18, height = 6)

# Find subtype markers
markers <- FindAllMarkers(purkinje, only.pos = TRUE, min.pct = 0.25, 
                          logfc.threshold = 0.25, verbose = FALSE)
write.csv(markers, "results/04_purkinje_subtypes/purkinje_subtype_markers.csv", 
          row.names = FALSE)

# Distribution across stages
stage_dist <- table(purkinje$seurat_clusters, purkinje$stage)
write.csv(stage_dist, "results/04_purkinje_subtypes/subtype_by_stage.csv")

saveRDS(purkinje, "results/04_purkinje_subtypes/purkinje_subtypes.rds")

cat("\nPurkinje subtyping complete\n")
