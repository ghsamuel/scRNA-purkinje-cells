library(CellChat)
library(ggplot2)

setwd("/core/projects/GAP/GDA/gsamuel/scRNAseq")

# Load CellChat objects
cellchat_e16 <- readRDS("results/04_R_analysis/10_cellchat_analysis/cellchat_e16.rds")
cellchat_e18 <- readRDS("results/04_R_analysis/10_cellchat_analysis/cellchat_e18.rds")

# Merge for comparison
object.list <- list(E16.5 = cellchat_e16, E18.5 = cellchat_e18)
cellchat_merged <- mergeCellChat(object.list, add.names = names(object.list))

# Define cell types
purkinje <- c("PC", "PC.immature")
partners <- c("BG", "GC", "GCP", "GABA.diff", "IN")

# ============================================================================
# PLOT 1: Pathways FROM Purkinje (what Purkinje SENDS)
# ============================================================================

cat("\n=== Analyzing pathways FROM Purkinje ===\n")

# Calculate pathway strength FROM Purkinje to partners
pathway_from_pc <- rowSums(cellchat_e16@netP$prob[purkinje, partners, ], dims = 2)
top_from_pc <- names(sort(pathway_from_pc, decreasing = TRUE))[1:10]

cat("Top 10 pathways FROM Purkinje:\n")
print(top_from_pc)

# Create bubble plot
png("results/04_R_analysis/10_cellchat_analysis/figures/purkinje/pathways_from_purkinje_bubble.png", 
    width=3600, height=2400, res=300)

gg <- netVisual_bubble(cellchat_merged, 
                       sources.use = purkinje, 
                       targets.use = partners,
                       comparison = c(1, 2), 
                       angle.x = 45,
                       signaling = top_from_pc,
                       remove.isolate = FALSE)
print(gg + ggtitle("Top 10 Pathways FROM Purkinje Cells"))

dev.off()

cat("Saved: pathways_from_purkinje_bubble.png\n")

# ============================================================================
# PLOT 2: Pathways TO Purkinje (what Purkinje RECEIVES)
# ============================================================================

cat("\n=== Analyzing pathways TO Purkinje ===\n")

# Calculate pathway strength TO Purkinje from partners
pathway_to_pc <- rowSums(cellchat_e16@netP$prob[partners, purkinje, ], dims = 2)
top_to_pc <- names(sort(pathway_to_pc, decreasing = TRUE))[1:10]

cat("Top 10 pathways TO Purkinje:\n")
print(top_to_pc)

# Create bubble plot
png("results/04_R_analysis/10_cellchat_analysis/figures/purkinje/pathways_to_purkinje_bubble.png", 
    width=3600, height=2400, res=300)

gg <- netVisual_bubble(cellchat_merged, 
                       sources.use = partners, 
                       targets.use = purkinje,
                       comparison = c(1, 2), 
                       angle.x = 45,
                       signaling = top_to_pc,
                       remove.isolate = FALSE)
print(gg + ggtitle("Top 10 Pathways TO Purkinje Cells"))

dev.off()

cat("Saved: pathways_to_purkinje_bubble.png\n")

cat("\n=== Complete! ===\n")
cat("Files saved in: results/04_R_analysis/10_cellchat_analysis/figures/purkinje/\n")
