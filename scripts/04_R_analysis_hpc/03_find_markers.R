#!/usr/bin/env Rscript
# Find cluster markers - exploration step


library(Seurat)
library(dplyr)
library(ggplot2)

setwd("/core/projects/GAP/GDA/gsamuel/scRNAseq")

dir.create("results/04_R_analysis/03_find_markers", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/04_R_analysis/03_find_markers", recursive = TRUE, showWarnings = FALSE)

# Load integrated cerebellum object
# Note: file is currently named purkinje_integrated.rds, but it represents the cerebellum object
obj <- readRDS("results/04_R_analysis/02_integration/purkinje_integrated.rds")

# Join layers for Seurat v5 before marker finding
obj <- JoinLayers(obj)

# Find cluster markers
markers <- FindAllMarkers(
  obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  verbose = FALSE
)

write.csv(
  markers,
  "results/04_R_analysis/03_find_markers/all_cluster_markers.csv",
  row.names = FALSE
)

# Top 10 markers per cluster
top10 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

write.csv(
  top10,
  "results/04_R_analysis/03_find_markers/top10_markers_per_cluster.csv",
  row.names = FALSE
)

# Known cerebellum markers
markers_list <- c(
  "Pcp2", "Calb1", "Car8",        # Purkinje
  "Aldoc", "Slc1a3", "Fabp7",    # Bergmann glia / astroglia
  "Neurod1", "Gabra6",           # Granule neurons
  "Pvalb", "Gad1",               # Interneurons
  "Ascl1", "Ptf1a",              # Progenitors
  "Mki67", "Top2a"               # Cycling cells
)

markers_present <- markers_list[markers_list %in% rownames(obj)]

print("Markers present in object:")
print(markers_present)

# Marker dotplot
if (length(markers_present) > 0) {
  p1 <- DotPlot(obj, features = markers_present) + coord_flip()

  ggsave(
    "figures/04_R_analysis/03_find_markers/marker_dotplot.png",
    p1,
    width = 14,
    height = 10
  )
} else {
  warning("None of the marker genes were found in the object.")
}

# UMAP clusters before annotation
p2 <- DimPlot(obj, label = TRUE, repel = TRUE) +
  ggtitle("Clusters before annotation")

ggsave(
  "figures/04_R_analysis/03_find_markers/umap_clusters.png",
  p2,
  width = 12,
  height = 10
)

# Cluster distribution across stages
cluster_stage <- as.data.frame(table(obj$seurat_clusters, obj$stage))

write.csv(
  cluster_stage,
  "results/04_R_analysis/03_find_markers/cluster_by_stage.csv",
  row.names = FALSE
)

# Save checked object for annotation script
saveRDS(
  obj,
  "results/04_R_analysis/03_find_markers/cerebellum_markers_checked.rds"
)

cat("Marker finding complete - inspect dotplot and markers before annotating\n")
