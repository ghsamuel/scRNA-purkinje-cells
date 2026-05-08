#!/usr/bin/env Rscript
# CellChat - compare communication between E16.5 and E18.5


library(Seurat)
library(CellChat)
library(ggplot2)

setwd("/core/projects/GAP/GDA/gsamuel/scRNAseq")

dir.create("results/05_cellchat", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/05_cellchat", recursive = TRUE, showWarnings = FALSE)

# Load annotated object and split by stage
obj <- readRDS("results/03_annotation/purkinje_annotated.rds")
e16 <- subset(obj, stage == "E16.5")
e18 <- subset(obj, stage == "E18.5")

# Create CellChat objects
data_e16 <- GetAssayData(e16, layer = "data")
cellchat_e16 <- createCellChat(object = data_e16, meta = e16@meta.data, 
                                group.by = "cell_type_integrated")
cellchat_e16@DB <- CellChatDB.mouse

data_e18 <- GetAssayData(e18, layer = "data")
cellchat_e18 <- createCellChat(object = data_e18, meta = e18@meta.data, 
                                group.by = "cell_type_integrated")
cellchat_e18@DB <- CellChatDB.mouse

# Infer communication - E16.5
cellchat_e16 <- subsetData(cellchat_e16)
cellchat_e16 <- identifyOverExpressedGenes(cellchat_e16)
cellchat_e16 <- identifyOverExpressedInteractions(cellchat_e16)
cellchat_e16 <- computeCommunProb(cellchat_e16, type = "triMean")
cellchat_e16 <- filterCommunication(cellchat_e16, min.cells = 10)
cellchat_e16 <- computeCommunProbPathway(cellchat_e16)
cellchat_e16 <- aggregateNet(cellchat_e16)

# Infer communication - E18.5 (lower min.cells for smaller groups)
cellchat_e18 <- subsetData(cellchat_e18)
cellchat_e18 <- identifyOverExpressedGenes(cellchat_e18)
cellchat_e18 <- identifyOverExpressedInteractions(cellchat_e18)
cellchat_e18 <- computeCommunProb(cellchat_e18, type = "triMean")
cellchat_e18 <- filterCommunication(cellchat_e18, min.cells = 5)
cellchat_e18 <- computeCommunProbPathway(cellchat_e18)
cellchat_e18 <- aggregateNet(cellchat_e18)

saveRDS(cellchat_e16, "results/05_cellchat/cellchat_e16.rds")
saveRDS(cellchat_e18, "results/05_cellchat/cellchat_e18.rds")

# Compare
cat("\nCommunication summary:\n")
cat("E16.5:", sum(cellchat_e16@net$count), "interactions\n")
cat("E18.5:", sum(cellchat_e18@net$count), "interactions\n")

# Circle plots
png("figures/05_cellchat/circle_plot_e16.png", width = 800, height = 800)
netVisual_circle(cellchat_e16@net$count, weight.scale = TRUE, 
                 label.edge = FALSE, title.name = "E16.5")
dev.off()

png("figures/05_cellchat/circle_plot_e18.png", width = 800, height = 800)
netVisual_circle(cellchat_e18@net$count, weight.scale = TRUE, 
                 label.edge = FALSE, title.name = "E18.5")
dev.off()

# Heatmaps
png("figures/05_cellchat/heatmap_e16.png", width = 800, height = 800)
netVisual_heatmap(cellchat_e16, color.heatmap = "Reds")
dev.off()

png("figures/05_cellchat/heatmap_e18.png", width = 800, height = 800)
netVisual_heatmap(cellchat_e18, color.heatmap = "Reds")
dev.off()

# Differential analysis
cellchat_e16 <- netAnalysis_computeCentrality(cellchat_e16)
cellchat_e18 <- netAnalysis_computeCentrality(cellchat_e18)

merged <- mergeCellChat(list(E16.5 = cellchat_e16, E18.5 = cellchat_e18), 
                        add.names = c("E16.5", "E18.5"))
saveRDS(merged, "results/05_cellchat/cellchat_merged.rds")

# Difference heatmap
png("figures/05_cellchat/difference_heatmap.png", width = 1000, height = 900)
netVisual_diffInteraction(merged, weight.scale = TRUE, measure = "count")
dev.off()

cat("\nCellChat complete\n")
