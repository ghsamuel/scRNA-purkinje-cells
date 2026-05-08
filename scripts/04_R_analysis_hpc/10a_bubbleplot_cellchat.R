library(CellChat)
library(ggplot2)

setwd("/core/projects/GAP/GDA/gsamuel/scRNAseq")

# Load CellChat objects
cellchat_e16 <- readRDS("results/04_R_analysis/10_cellchat_analysis/cellchat_e16.rds")
cellchat_e18 <- readRDS("results/04_R_analysis/10_cellchat_analysis/cellchat_e18.rds")

# Merge for comparison
object.list <- list(E16.5 = cellchat_e16, E18.5 = cellchat_e18)
cellchat_merged <- mergeCellChat(object.list, add.names = names(object.list))

# Key cell types
key_cells <- c("PC", "PC.immature", "BG", "GC", "GCP", "GABA.diff", "IN")

# Get top 20 pathways based on overall communication strength
pathway_strength <- rowSums(cellchat_e16@netP$prob[key_cells, key_cells, ], dims = 2)
top_pathways <- names(sort(pathway_strength, decreasing = TRUE))[1:20]

cat("Top 20 pathways:\n")
print(top_pathways)

# Create bubble plot with only top 20 pathways
pdf("results/04_R_analysis/10_cellchat_analysis/figures/overall/pathway_bubble_top20.pdf", 
    width=14, height=10)

gg <- netVisual_bubble(cellchat_merged, 
                       sources.use = key_cells, 
                       targets.use = key_cells,
                       comparison = c(1, 2), 
                       angle.x = 45,
                       signaling = top_pathways,
                       remove.isolate = FALSE)
print(gg)

dev.off()

cat("\nSaved: pathway_bubble_top20.pdf\n")
