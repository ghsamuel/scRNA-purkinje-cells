#!/usr/bin/env Rscript
# Annotate cell types using marker genes
# Jini Samuel

library(Seurat)
library(dplyr)
library(ggplot2)

setwd("/core/projects/GAP/GDA/gsamuel/scRNAseq")

# Create output directories
dir.create("results/03_annotation", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/03_annotation", recursive = TRUE, showWarnings = FALSE)

# Load integrated object
obj <- readRDS("results/02_integration/purkinje_integrated.rds")

# Find markers - need to join layers first in Seurat v5
obj <- JoinLayers(obj)
markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, 
                          logfc.threshold = 0.25, verbose = FALSE)

write.csv(markers, "results/03_annotation/all_cluster_markers.csv", row.names = FALSE)

# Top 10 per cluster
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
write.csv(top10, "results/03_annotation/top10_markers_per_cluster.csv", row.names = FALSE)

# Known markers
markers_list <- c("Pcp2", "Calb1", "Car8",  # Purkinje
                  "Slc1a3", "Aldoc", "Fabp7",  # Bergmann glia
                  "Neurod1", "Gabra6",  # Granule
                  "Pvalb", "Gad1",  # Interneurons
                  "Ascl1", "Ptf1a",  # Progenitors
                  "Mki67", "Top2a")  # Cycling

p1 <- DotPlot(obj, features = markers_list) + coord_flip()
ggsave("figures/03_annotation/marker_dotplot.png", p1, width = 14, height = 10)

# Annotate - based on marker expression from dotplot
annotations <- c(
  "0" = "Bergmann_glia",
  "1" = "Granule_progenitors",
  "2" = "Interneuron_progenitors",
  "3" = "Interneurons",
  "4" = "Granule_neurons",
  "5" = "Astrocytes",
  "6" = "Cycling_progenitors",
  "7" = "Purkinje_cells",
  "8" = "Oligodendrocyte_progenitors",
  "9" = "Microglia",
  "10" = "Endothelial",
  "11" = "Purkinje_subtype",
  "12" = "Unknown_1",
  "13" = "Unknown_2",
  "14" = "Pericytes"
)

obj$cell_type_integrated <- annotations[as.character(obj$seurat_clusters)]

# Plot annotated UMAP
p2 <- DimPlot(obj, group.by = "cell_type_integrated", label = TRUE, repel = TRUE)
ggsave("figures/03_annotation/umap_annotated.png", p2, width = 12, height = 10)

# Cell type counts
ct_counts <- table(obj$cell_type_integrated, obj$stage)
write.csv(ct_counts, "results/03_annotation/celltype_by_stage.csv")

saveRDS(obj, "results/03_annotation/purkinje_annotated.rds")

cat("Annotation complete -", length(unique(obj$cell_type_integrated)), "cell types identified\n")
