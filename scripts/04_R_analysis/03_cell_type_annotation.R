#!/usr/bin/env Rscript
# Annotate cell types using marker genes
# Jini Samuel

library(Seurat)
library(dplyr)
library(ggplot2)

setwd("/core/projects/GAP/GDA/gsamuel/scRNAseq")

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
  "0" = "Bergmann_glia", "1" = "Cycling_cells", "2" = "Granule_neurons",
  "3" = "Granule_neurons", "4" = "Granule_neurons", "5" = "Granule_neurons",
  "6" = "Bergmann_glia", "7" = "Purkinje_cells", "8" = "Interneurons",
  "9" = "Interneurons", "10" = "Granule_neurons", "11" = "Purkinje_subtype",
  "12" = "GABAergic_neurons", "13" = "Progenitors", "14" = "Progenitors",
  "15" = "Progenitors", "16" = "RBCs", "17" = "OPCs", "18" = "Interneurons",
  "19" = "Bergmann_glia", "20" = "RBCs", "21" = "Pericytes",
  "22" = "Macrophages", "23" = "Immune_cells", "24" = "Epithelial",
  "25" = "Endothelial"
)

obj$cell_type_integrated <- annotations[as.character(obj$seurat_clusters)]

# Plot
p2 <- DimPlot(obj, group.by = "cell_type_integrated", label = TRUE, repel = TRUE)
ggsave("figures/03_annotation/umap_annotated.png", p2, width = 12, height = 10)

# Summary
cell_counts <- table(obj$cell_type_integrated)
summary_df <- data.frame(
  CellType = names(cell_counts),
  Count = as.numeric(cell_counts),
  Percentage = round(100 * as.numeric(cell_counts) / ncol(obj), 2)
) %>% arrange(desc(Count))

write.csv(summary_df, "results/03_annotation/cell_type_summary.csv", row.names = FALSE)
print(summary_df)

saveRDS(obj, "results/03_annotation/purkinje_annotated.rds")

cat("\nAnnotation complete\n")
