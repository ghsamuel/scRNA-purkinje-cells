#!/usr/bin/env Rscript
# High-resolution Purkinje cell subtyping
# Looking for Foxp1/Foxp2/Foxp4 expression patterns

library(Seurat)
library(dplyr)
library(ggplot2)

setwd("/core/projects/GAP/GDA/gsamuel/scRNAseq")

# Create output directories
dir.create("results/04_R_analysis/05_purkinje_subtypes", recursive = TRUE, showWarnings = FALSE)
dir.create("results/04_R_analysis/05_purkinje_subtypes/figures", recursive = TRUE, showWarnings = FALSE)

# Load annotated cerebellum
obj <- readRDS("results/04_R_analysis/04_annotate_clusters/cerebellum_annotated.rds")

# Foxp genes across whole cerebellum
foxp_genes <- c("Foxp1", "Foxp2", "Foxp4")

p_cerebellum <- FeaturePlot(obj, features = foxp_genes, ncol = 3, order = TRUE, pt.size = 0.5)
ggsave("results/04_R_analysis/05_purkinje_subtypes/figures/foxp_expression_cerebellum.pdf", 
       p_cerebellum, width = 18, height = 6)

# Extract Purkinje cells (clusters 9 + 11)
purkinje <- subset(obj, cell_type %in% c("PC", "PC.immature"))

cat("Extracted", ncol(purkinje), "Purkinje cells\n")

# Re-cluster at high resolution to find subtypes
purkinje <- NormalizeData(purkinje)
purkinje <- FindVariableFeatures(purkinje, nfeatures = 2000)
purkinje <- ScaleData(purkinje)
purkinje <- RunPCA(purkinje, npcs = 20, verbose = FALSE)
purkinje <- FindNeighbors(purkinje, dims = 1:15)
purkinje <- FindClusters(purkinje, resolution = 0.5, verbose = FALSE)
purkinje <- RunUMAP(purkinje, dims = 1:15, verbose = FALSE)

cat("Found", length(unique(purkinje$seurat_clusters)), "subtypes\n")

# Plot subtypes with improved styling - publication ready
p1 <- DimPlot(
  purkinje,
  label = TRUE,
  repel = TRUE,
  pt.size = 2.5,
  label.size = 10,
  raster = FALSE
) +
  ggtitle("Purkinje Cell Subtypes") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 28, face = "bold"),
    legend.text = element_text(size = 24, face = "bold"),
    axis.text = element_text(size = 24, face = "bold"),
    axis.title = element_text(size = 26, face = "bold"),
    axis.line = element_line(linewidth = 1.5, color = "black"),
    axis.ticks = element_line(linewidth = 1.2, color = "black")
  )

ggsave(
  "results/04_R_analysis/05_purkinje_subtypes/figures/purkinje_subtypes_umap.pdf",
  p1,
  width = 12,
  height = 8
)

# UMAP split by stage, colored by cell type (PC vs PC.immature)
p_split <- DimPlot(
  purkinje,
  group.by = "cell_type",
  split.by = "stage",
  label = TRUE,
  repel = TRUE,
  pt.size = 2.5,
  label.size = 10,
  raster = FALSE,
  ncol = 2
) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 28, face = "bold"),
    legend.text = element_text(size = 24, face = "bold"),
    axis.text = element_text(size = 24, face = "bold"),
    axis.title = element_text(size = 26, face = "bold"),
    axis.line = element_line(linewidth = 1.5, color = "black"),
    axis.ticks = element_line(linewidth = 1.2, color = "black"),
    strip.text = element_text(size = 24, face = "bold")
  )

ggsave(
  "results/04_R_analysis/05_purkinje_subtypes/figures/purkinje_celltype_by_stage.pdf",
  p_split,
  width = 20,
  height = 8
)

# Foxp expression in Purkinje subset
p2 <- DotPlot(purkinje, features = foxp_genes) + coord_flip()
ggsave("results/04_R_analysis/05_purkinje_subtypes/figures/foxp_dotplot.pdf", p2, width = 10, height = 5)

p3 <- FeaturePlot(purkinje, features = foxp_genes, ncol = 3, order = TRUE, pt.size = 0.8)
ggsave("results/04_R_analysis/05_purkinje_subtypes/figures/foxp_featureplots.pdf", p3, width = 18, height = 6)

# Find subtype markers
markers <- FindAllMarkers(purkinje, only.pos = TRUE, min.pct = 0.25, 
                          logfc.threshold = 0.25, verbose = FALSE)
write.csv(markers, "results/04_R_analysis/05_purkinje_subtypes/purkinje_subtype_markers.csv", 
          row.names = FALSE)

# Get top 5 markers per subtype for dotplot (like paper Figure 2b)
top5 <- markers %>%
  group_by(cluster) %>%
  top_n(5, avg_log2FC) %>%
  pull(gene) %>%
  unique()

# Dotplot of top markers (like Figure 2b from paper)
p_markers <- DotPlot(purkinje, features = top5) + 
  coord_flip() +
  ggtitle("Top Purkinje Subtype Markers")

ggsave("results/04_R_analysis/05_purkinje_subtypes/figures/subtype_markers_dotplot.pdf", 
       p_markers, width = 12, height = 10)

# Distribution across stages
stage_dist <- table(purkinje$seurat_clusters, purkinje$stage)
write.csv(stage_dist, "results/04_R_analysis/05_purkinje_subtypes/subtype_by_stage.csv")

# Stacked barplot of stage distribution
stage_df <- as.data.frame(stage_dist)
colnames(stage_df) <- c("Subtype", "Stage", "Count")

p_stage <- ggplot(stage_df, aes(x = Subtype, y = Count, fill = Stage)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Purkinje Subtype", y = "Proportion", 
       title = "Purkinje Subtype Distribution by Stage") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.text = element_text(size = 14, face = "bold")
  )

ggsave("results/04_R_analysis/05_purkinje_subtypes/figures/subtype_by_stage.pdf", 
       p_stage, width = 10, height = 6)

# Distribution by cell type (PC vs PC.immature)
celltype_dist <- table(purkinje$seurat_clusters, purkinje$cell_type)
write.csv(celltype_dist, "results/04_R_analysis/05_purkinje_subtypes/subtype_by_celltype.csv")

# Stacked barplot by cell type
celltype_df <- as.data.frame(celltype_dist)
colnames(celltype_df) <- c("Subtype", "CellType", "Count")

p_celltype_bar <- ggplot(celltype_df, aes(x = Subtype, y = Count, fill = CellType)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Purkinje Subtype", y = "Proportion", 
       title = "PC vs PC.immature by Subtype") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.text = element_text(size = 14, face = "bold")
  )

ggsave("results/04_R_analysis/05_purkinje_subtypes/figures/subtype_by_celltype.pdf", 
       p_celltype_bar, width = 10, height = 6)

saveRDS(purkinje, "results/04_R_analysis/05_purkinje_subtypes/purkinje_subtypes.rds")

cat("\nPurkinje subtyping complete\n")
cat("Generated", length(unique(purkinje$seurat_clusters)), "subtypes\n")
cat("E16.5 cells:", sum(purkinje$stage == "E16.5"), "\n")
cat("E18.5 cells:", sum(purkinje$stage == "E18.5"), "\n")
