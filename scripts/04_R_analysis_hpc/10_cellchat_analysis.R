library(Seurat)
library(CellChat)
library(dplyr)
library(ggplot2)
library(circlize)
library(patchwork)

setwd("/core/projects/GAP/GDA/gsamuel/scRNAseq")

# Create directories
dir.create("results/04_R_analysis/10_cellchat_analysis/figures/overall", recursive = TRUE)
dir.create("results/04_R_analysis/10_cellchat_analysis/figures/purkinje", recursive = TRUE)
dir.create("results/04_R_analysis/10_cellchat_analysis/tables", recursive = TRUE)

# Load data
cat("Loading data...\n")
obj <- readRDS("results/04_R_analysis/04_annotate_clusters/cerebellum_annotated.rds")

# Split by stage
e16 <- subset(obj, stage == "E16.5")
e18 <- subset(obj, stage == "E18.5")

# ============================================================================
# Create CellChat objects
# ============================================================================

cat("\nCreating CellChat objects...\n")

# E16.5
data_e16 <- GetAssayData(e16, layer = "data")
cellchat_e16 <- createCellChat(object = data_e16, meta = e16@meta.data, 
                               group.by = "cell_type")
cellchat_e16@DB <- CellChatDB.mouse

# E18.5
data_e18 <- GetAssayData(e18, layer = "data")
cellchat_e18 <- createCellChat(object = data_e18, meta = e18@meta.data,
                               group.by = "cell_type")
cellchat_e18@DB <- CellChatDB.mouse

# ============================================================================
# Preprocess
# ============================================================================

cat("Preprocessing E16.5...\n")
cellchat_e16 <- subsetData(cellchat_e16)
cellchat_e16 <- identifyOverExpressedGenes(cellchat_e16)
cellchat_e16 <- identifyOverExpressedInteractions(cellchat_e16)
cellchat_e16 <- computeCommunProb(cellchat_e16, type = "triMean")
cellchat_e16 <- filterCommunication(cellchat_e16, min.cells = 10)
cellchat_e16 <- computeCommunProbPathway(cellchat_e16)
cellchat_e16 <- aggregateNet(cellchat_e16)
cellchat_e16 <- netAnalysis_computeCentrality(cellchat_e16)

cat("Preprocessing E18.5...\n")
cellchat_e18 <- subsetData(cellchat_e18)
cellchat_e18 <- identifyOverExpressedGenes(cellchat_e18)
cellchat_e18 <- identifyOverExpressedInteractions(cellchat_e18)
cellchat_e18 <- computeCommunProb(cellchat_e18, type = "triMean")
cellchat_e18 <- filterCommunication(cellchat_e18, min.cells = 5)
cellchat_e18 <- computeCommunProbPathway(cellchat_e18)
cellchat_e18 <- aggregateNet(cellchat_e18)
cellchat_e18 <- netAnalysis_computeCentrality(cellchat_e18)

# Save
saveRDS(cellchat_e16, "results/04_R_analysis/10_cellchat_analysis/cellchat_e16.rds")
saveRDS(cellchat_e18, "results/04_R_analysis/10_cellchat_analysis/cellchat_e18.rds")

# ============================================================================
# SECTION 1: Overall Communication
# ============================================================================

cat("\n=== Overall Communication Analysis ===\n")

# Summary stats
total_e16 <- sum(cellchat_e16@net$count)
total_e18 <- sum(cellchat_e18@net$count)
strength_e16 <- sum(cellchat_e16@net$weight)
strength_e18 <- sum(cellchat_e18@net$weight)

cat("E16.5:", total_e16, "interactions, strength:", round(strength_e16, 1), "\n")
cat("E18.5:", total_e18, "interactions, strength:", round(strength_e18, 1), "\n")
cat("Change:", round((total_e18 - total_e16)/total_e16 * 100, 1), "%\n")

# 1. Stage comparison bar plots
comparison_df <- data.frame(
  Stage = c("E16.5", "E18.5"),
  Interactions = c(total_e16, total_e18),
  Strength = c(strength_e16, strength_e18)
)

pdf("results/04_R_analysis/10_cellchat_analysis/figures/overall/stage_comparison.pdf", 
    width=12, height=5)
p1 <- ggplot(comparison_df, aes(x=Stage, y=Interactions, fill=Stage)) +
  geom_bar(stat="identity", width=0.6) +
  theme_minimal(base_size=14) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
  ggtitle("Total Interactions") +
  theme(legend.position="none") +
  geom_text(aes(label=Interactions), vjust=-0.5, size=5)

p2 <- ggplot(comparison_df, aes(x=Stage, y=Strength, fill=Stage)) +
  geom_bar(stat="identity", width=0.6) +
  theme_minimal(base_size=14) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
  ggtitle("Interaction Strength") +
  theme(legend.position="none") +
  geom_text(aes(label=round(Strength, 1)), vjust=-0.5, size=5)

print(p1 + p2)
dev.off()

# 2. Heatmaps with marginal bars (shows signaling hubs)
pdf("results/04_R_analysis/10_cellchat_analysis/figures/overall/heatmap_hubs_e16.pdf", 
    width=12, height=10)
netVisual_heatmap(cellchat_e16, color.heatmap = "Reds", 
                  title.name = "E16.5 Interaction Count")
dev.off()

pdf("results/04_R_analysis/10_cellchat_analysis/figures/overall/heatmap_hubs_e18.pdf", 
    width=12, height=10)
netVisual_heatmap(cellchat_e18, color.heatmap = "Reds",
                  title.name = "E18.5 Interaction Count")
dev.off()

# 3. Chord diagrams (separate for each stage)
# Focus on key neural cells
key_cells <- c("PC", "PC.immature", "BG", "GC", "GCP", "GABA.diff", "IN")

# Manually subset to key cells only
net_e16_subset <- cellchat_e16@net$count[key_cells, key_cells]
net_e18_subset <- cellchat_e18@net$count[key_cells, key_cells]

# E16.5 chord
circos.clear()
pdf("results/04_R_analysis/10_cellchat_analysis/figures/overall/chord_e16.pdf", 
    width=12, height=12)
chordDiagram(net_e16_subset, 
             grid.col = rainbow(length(key_cells)),
             transparency = 0.5,
             annotationTrack = c("grid", "name"),
             annotationTrackHeight = c(0.03, 0.1))
circos.track(track.index = 2, panel.fun = function(x, y) {
  sector.name = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), ylim[1], sector.name, 
              facing = "bending.inside", 
              niceFacing = TRUE, 
              adj = c(0.5, 0),
              cex = 1.3, font = 2)
}, bg.border = NA)
title("E16.5 Cell-Cell Communication", cex.main = 2)
dev.off()
circos.clear()

# E18.5 chord
circos.clear()
pdf("results/04_R_analysis/10_cellchat_analysis/figures/overall/chord_e18.pdf", 
    width=12, height=12)
chordDiagram(net_e18_subset, 
             grid.col = rainbow(length(key_cells)),
             transparency = 0.5,
             annotationTrack = c("grid", "name"),
             annotationTrackHeight = c(0.03, 0.1))
circos.track(track.index = 2, panel.fun = function(x, y) {
  sector.name = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), ylim[1], sector.name, 
              facing = "bending.inside", 
              niceFacing = TRUE, 
              adj = c(0.5, 0),
              cex = 1.3, font = 2)
}, bg.border = NA)
title("E18.5 Cell-Cell Communication", cex.main = 2)
dev.off()
circos.clear()

# 4. Difference heatmap
diff_matrix <- cellchat_e16@net$count - cellchat_e18@net$count
diff_subset <- diff_matrix[key_cells, key_cells]

pdf("results/04_R_analysis/10_cellchat_analysis/figures/overall/difference_heatmap.pdf", 
    width=12, height=10)
heatmap(diff_subset, 
        Rowv = NA, Colv = NA,
        col = colorRampPalette(c("blue", "white", "red"))(100),
        scale = "none",
        margins = c(14, 14),
        cexRow = 1.2, cexCol = 1.2,
        main = "Change in interactions (E16.5 - E18.5)\nRed = more E16.5, Blue = more E18.5")
dev.off()

# 5. Bubble plot comparing pathways
object.list <- list(E16.5 = cellchat_e16, E18.5 = cellchat_e18)
cellchat_merged <- mergeCellChat(object.list, add.names = names(object.list))

pdf("results/04_R_analysis/10_cellchat_analysis/figures/overall/pathway_bubble.pdf", 
    width=10, height=8)
gg <- netVisual_bubble(cellchat_merged, sources.use = key_cells, targets.use = key_cells,
                       comparison = c(1, 2), angle.x = 45)
print(gg)
dev.off()

# ============================================================================
# SECTION 2: Purkinje-Focused Analysis
# ============================================================================

cat("\n=== Purkinje-Focused Analysis ===\n")

purkinje_cells <- c("PC", "PC.immature")

# 1. Incoming signals to Purkinje (who talks TO Purkinje)
incoming_e16 <- cellchat_e16@net$count[, purkinje_cells]
incoming_e18 <- cellchat_e18@net$count[, purkinje_cells]

# Sum across both PC types
incoming_total_e16 <- rowSums(incoming_e16)
incoming_total_e18 <- rowSums(incoming_e18)

incoming_df <- data.frame(
  Cell_Type = names(incoming_total_e16),
  E16.5 = incoming_total_e16,
  E18.5 = incoming_total_e18,
  Change = incoming_total_e18 - incoming_total_e16
) %>% arrange(desc(E16.5))

write.csv(incoming_df, 
          "results/04_R_analysis/10_cellchat_analysis/tables/purkinje_incoming.csv",
          row.names = FALSE)

pdf("results/04_R_analysis/10_cellchat_analysis/figures/purkinje/incoming_signals.pdf", 
    width=12, height=6)
incoming_long <- reshape2::melt(incoming_df[1:10,], id.vars = "Cell_Type", 
                                measure.vars = c("E16.5", "E18.5"))
p <- ggplot(incoming_long, aes(x=reorder(Cell_Type, -value), y=value, fill=variable)) +
  geom_bar(stat="identity", position="dodge") +
  theme_minimal(base_size=12) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
  labs(x="", y="Interaction Count", fill="Stage",
       title="Incoming Signals to Purkinje Cells (Top 10 Sources)")
print(p)
dev.off()

# 2. Outgoing signals from Purkinje (who Purkinje talks TO)
outgoing_e16 <- cellchat_e16@net$count[purkinje_cells, ]
outgoing_e18 <- cellchat_e18@net$count[purkinje_cells, ]

outgoing_total_e16 <- colSums(outgoing_e16)
outgoing_total_e18 <- colSums(outgoing_e18)

outgoing_df <- data.frame(
  Cell_Type = names(outgoing_total_e16),
  E16.5 = outgoing_total_e16,
  E18.5 = outgoing_total_e18,
  Change = outgoing_total_e18 - outgoing_total_e16
) %>% arrange(desc(E16.5))

write.csv(outgoing_df, 
          "results/04_R_analysis/10_cellchat_analysis/tables/purkinje_outgoing.csv",
          row.names = FALSE)

pdf("results/04_R_analysis/10_cellchat_analysis/figures/purkinje/outgoing_signals.pdf", 
    width=12, height=6)
outgoing_long <- reshape2::melt(outgoing_df[1:10,], id.vars = "Cell_Type",
                                measure.vars = c("E16.5", "E18.5"))
p <- ggplot(outgoing_long, aes(x=reorder(Cell_Type, -value), y=value, fill=variable)) +
  geom_bar(stat="identity", position="dodge") +
  theme_minimal(base_size=12) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
  labs(x="", y="Interaction Count", fill="Stage",
       title="Outgoing Signals from Purkinje Cells (Top 10 Targets)")
print(p)
dev.off()

# 3. PC ↔ PC.immature communication
pc_to_pc <- data.frame(
  Interaction = c("PC → PC", "PC → PC.immature", 
                  "PC.immature → PC", "PC.immature → PC.immature"),
  E16.5 = c(
    cellchat_e16@net$count["PC", "PC"],
    cellchat_e16@net$count["PC", "PC.immature"],
    cellchat_e16@net$count["PC.immature", "PC"],
    cellchat_e16@net$count["PC.immature", "PC.immature"]
  ),
  E18.5 = c(
    cellchat_e18@net$count["PC", "PC"],
    cellchat_e18@net$count["PC", "PC.immature"],
    cellchat_e18@net$count["PC.immature", "PC"],
    cellchat_e18@net$count["PC.immature", "PC.immature"]
  )
)

write.csv(pc_to_pc, 
          "results/04_R_analysis/10_cellchat_analysis/tables/purkinje_self_communication.csv",
          row.names = FALSE)

pdf("results/04_R_analysis/10_cellchat_analysis/figures/purkinje/pc_to_pc.pdf", 
    width=10, height=6)
pc_long <- reshape2::melt(pc_to_pc, id.vars = "Interaction")
p <- ggplot(pc_long, aes(x=Interaction, y=value, fill=variable)) +
  geom_bar(stat="identity", position="dodge") +
  theme_minimal(base_size=12) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
  labs(x="", y="Interaction Count", fill="Stage",
       title="Purkinje Cell Self-Communication")
print(p)
dev.off()

# 4. Top pathways TO Purkinje
pathways_to_pc_e16 <- cellchat_e16@netP$prob[, purkinje_cells, ]
pathways_to_pc_e18 <- cellchat_e18@netP$prob[, purkinje_cells, ]

# Sum across pathways and PC types
pathway_sums_e16 <- apply(pathways_to_pc_e16, 3, sum)
pathway_sums_e18 <- apply(pathways_to_pc_e18, 3, sum)

top_pathways_to_pc <- data.frame(
  Pathway = names(pathway_sums_e16),
  E16.5 = pathway_sums_e16,
  E18.5 = pathway_sums_e18,
  Change = pathway_sums_e18 - pathway_sums_e16
) %>% 
  filter(E16.5 > 0 | E18.5 > 0) %>%
  arrange(desc(E16.5)) %>%
  head(15)

write.csv(top_pathways_to_pc, 
          "results/04_R_analysis/10_cellchat_analysis/tables/top_pathways_to_purkinje.csv",
          row.names = FALSE)

pdf("results/04_R_analysis/10_cellchat_analysis/figures/purkinje/pathways_to_pc.pdf", 
    width=12, height=8)
pathway_long <- reshape2::melt(top_pathways_to_pc, id.vars = "Pathway",
                               measure.vars = c("E16.5", "E18.5"))
p <- ggplot(pathway_long, aes(x=reorder(Pathway, -value), y=value, fill=variable)) +
  geom_bar(stat="identity", position="dodge") +
  theme_minimal(base_size=12) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
  labs(x="", y="Pathway Probability", fill="Stage",
       title="Top Signaling Pathways TO Purkinje Cells")
print(p)
dev.off()

# 5. Top pathways FROM Purkinje
pathways_from_pc_e16 <- cellchat_e16@netP$prob[purkinje_cells, , ]
pathways_from_pc_e18 <- cellchat_e18@netP$prob[purkinje_cells, , ]

pathway_sums_from_e16 <- apply(pathways_from_pc_e16, 3, sum)
pathway_sums_from_e18 <- apply(pathways_from_pc_e18, 3, sum)

top_pathways_from_pc <- data.frame(
  Pathway = names(pathway_sums_from_e16),
  E16.5 = pathway_sums_from_e16,
  E18.5 = pathway_sums_from_e18,
  Change = pathway_sums_from_e18 - pathway_sums_from_e16
) %>% 
  filter(E16.5 > 0 | E18.5 > 0) %>%
  arrange(desc(E16.5)) %>%
  head(15)

write.csv(top_pathways_from_pc, 
          "results/04_R_analysis/10_cellchat_analysis/tables/top_pathways_from_purkinje.csv",
          row.names = FALSE)

pdf("results/04_R_analysis/10_cellchat_analysis/figures/purkinje/pathways_from_pc.pdf", 
    width=12, height=8)
pathway_long <- reshape2::melt(top_pathways_from_pc, id.vars = "Pathway",
                               measure.vars = c("E16.5", "E18.5"))
p <- ggplot(pathway_long, aes(x=reorder(Pathway, -value), y=value, fill=variable)) +
  geom_bar(stat="identity", position="dodge") +
  theme_minimal(base_size=12) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
  labs(x="", y="Pathway Probability", fill="Stage",
       title="Top Signaling Pathways FROM Purkinje Cells")
print(p)
dev.off()

# 6. Purkinje-specific pathway bubble plot
pdf("results/04_R_analysis/10_cellchat_analysis/figures/purkinje/pathway_bubble.pdf", 
    width=12, height=8)
gg <- netVisual_bubble(cellchat_merged, 
                       sources.use = purkinje_cells, 
                       targets.use = key_cells,
                       comparison = c(1, 2), 
                       angle.x = 45,
                       remove.isolate = FALSE)
print(gg)
dev.off()

# ============================================================================
# Summary table
# ============================================================================

summary_table <- data.frame(
  Metric = c("Total Cells", "Cell Types", "Total Interactions", 
             "Interaction Strength"),
  E16.5 = c(
    ncol(cellchat_e16@data.signaling),
    length(unique(cellchat_e16@idents)),
    total_e16,
    round(strength_e16, 1)
  ),
  E18.5 = c(
    ncol(cellchat_e18@data.signaling),
    length(unique(cellchat_e18@idents)),
    total_e18,
    round(strength_e18, 1)
  ),
  Percent_Change = c(
    "",
    "",
    paste0(round((total_e18 - total_e16)/total_e16 * 100, 1), "%"),
    paste0(round((strength_e18 - strength_e16)/strength_e16 * 100, 1), "%")
  )
)

write.csv(summary_table, 
          "results/04_R_analysis/10_cellchat_analysis/tables/overall_summary.csv",
          row.names = FALSE)

print(summary_table)

cat("\n=== Analysis Complete ===\n")
cat("Overall figures: results/04_R_analysis/10_cellchat_analysis/figures/overall/\n")
cat("Purkinje figures: results/04_R_analysis/10_cellchat_analysis/figures/purkinje/\n")
cat("Tables: results/04_R_analysis/10_cellchat_analysis/tables/\n")
