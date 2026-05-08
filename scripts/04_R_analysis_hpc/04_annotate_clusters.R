#!/usr/bin/env Rscript
# 04_annotate_clusters.R
# Annotate cerebellum clusters based on marker inspection


library(Seurat)
library(ggplot2)

setwd("/core/projects/GAP/GDA/gsamuel/scRNAseq")

dir.create("results/04_R_analysis/04_annotate_clusters", recursive = TRUE, showWarnings = FALSE)
dir.create("results/04_R_analysis/04_annotate_clusters/figures", recursive = TRUE, showWarnings = FALSE)

obj <- readRDS("results/04_R_analysis/03_find_markers/cerebellum_markers_checked.rds")

annotations <- c(
  "0" = "BG",
  "1" = "GCP.cyc",
  "2" = "IN",
  "3" = "GCP",
  "4" = "GC",
  "5" = "IN",
  "6" = "GABA.Pre",
  "7" = "BG",
  "8" = "GC",
  "9" = "PC",
  "10" = "OPC",
  "11" = "PC.immature",
  "12" = "UBC",
  "13" = "CPE",
  "14" = "GABA.diff",
  "15" = "BG.cycling",
  "16" = "Erythrocyte",
  "17" = "Oligodendrocytes",
  "18" = "CN",
  "19" = "Erythrocyte",
  "20" = "Endothelial",
  "21" = "Microglia",
  "22" = "Microglia",
  "23" = "Pericyte",
  "24" = "CPE"
)

cell_type_vec <- unname(annotations[as.character(obj$seurat_clusters)])

obj <- AddMetaData(
  obj,
  metadata = cell_type_vec,
  col.name = "cell_type"
)

# UMAP styled similar to integrated plot
p1 <- DimPlot(
  obj,
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE,
  pt.size = 1.2,
  label.size = 7,
  raster = FALSE
) +
  ggtitle("Annotated cerebellum cell types") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    legend.text = element_text(size = 20),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 22)
  )

# High-resolution PNG
ggsave(
  "results/04_R_analysis/04_annotate_clusters/figures/umap_annotated_celltypes.png",
  p1,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)

# Publication-quality PDF
ggsave(
  "results/04_R_analysis/04_annotate_clusters/figures/umap_annotated_celltypes.pdf",
  p1,
  width = 12,
  height = 8
)

# Cell type counts by stage
ct_counts <- as.data.frame(table(obj$cell_type, obj$stage))
colnames(ct_counts) <- c("cell_type", "stage", "count")

write.csv(
  ct_counts,
  "results/04_R_analysis/04_annotate_clusters/celltype_by_stage.csv",
  row.names = FALSE
)

# Cluster summary table
cluster_counts <- as.data.frame(table(obj$seurat_clusters))
colnames(cluster_counts) <- c("cluster", "n_cells")

summary_df <- data.frame(
  cluster = names(annotations),
  cell_type = unname(annotations)
)

summary_df <- merge(summary_df, cluster_counts, by = "cluster", all.x = TRUE)

write.csv(
  summary_df,
  "results/04_R_analysis/04_annotate_clusters/annotation_summary.csv",
  row.names = FALSE
)

saveRDS(
  obj,
  "results/04_R_analysis/04_annotate_clusters/cerebellum_annotated.rds"
)

cat("Annotation complete!\n")
cat("Total cells:", ncol(obj), "\n")
cat("Cell types identified:", length(unique(obj$cell_type)), "\n")
cat("High-resolution PNG and PDF saved.\n")
