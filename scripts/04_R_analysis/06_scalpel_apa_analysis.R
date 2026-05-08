#!/usr/bin/env Rscript
# SCALPEL APA analysis - 3'UTR length changes


library(Seurat)
library(scalpelR)
library(rtracklayer)
library(dplyr)
library(ggplot2)

setwd("/core/projects/GAP/GDA/gsamuel/scRNAseq")

out_dir <- "results/06_SCALPEL_analysis"
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

# Load SCALPEL-integrated object (has isoform counts added)
obj <- readRDS(file.path(out_dir, "idge_cerebellum_final.rds"))

# Load annotation for UTR extraction
gtf <- import("data/references/mm39/GRCm39.gtf")

# Calculate weighted 3'UTR lengths by developmental stage
utr_lengths <- weightedUTR_cell(
  seurat_obj = obj,
  annotation_gr = gtf,
  group.by = "stage",
  assay = "RNA"
)

saveRDS(utr_lengths, file.path(out_dir, "weighted_3UTR_by_stage_full.rds"))
write.csv(utr_lengths, file.path(out_dir, "weighted_3UTR_by_stage_full.csv"), 
          row.names = FALSE)

# Find genes with biggest UTR changes
utr_wide <- utr_lengths %>%
  tidyr::pivot_wider(names_from = stage, values_from = weighted_utr_length, 
                      id_cols = gene)

if ("E16.5" %in% colnames(utr_wide) && "E18.5" %in% colnames(utr_wide)) {
  utr_wide$delta <- utr_wide$`E16.5` - utr_wide$`E18.5`
  
  top_apa <- utr_wide %>%
    filter(!is.na(delta)) %>%
    arrange(desc(abs(delta))) %>%
    head(130)
  
  write.csv(top_apa, file.path(out_dir, "top_APA_genes.csv"), row.names = FALSE)
  
  cat("Top APA genes:\n")
  cat("  Lengthening:", sum(top_apa$delta > 0, na.rm = TRUE), "\n")
  cat("  Shortening:", sum(top_apa$delta < 0, na.rm = TRUE), "\n")
}

# Plots
p1 <- ggplot(utr_lengths, aes(x = stage, y = weighted_utr_length)) +
  geom_violin(fill = "steelblue", alpha = 0.5) +
  geom_boxplot(width = 0.1) +
  labs(title = "3'UTR Length by Stage", x = "Stage", y = "UTR Length (nt)") +
  theme_minimal()

ggsave(file.path(out_dir, "figures/utr_length_by_stage.png"), p1, 
       width = 8, height = 6)

if (exists("top_apa")) {
  p2 <- ggplot(top_apa, aes(x = `E16.5`, y = `E18.5`)) +
    geom_point(alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(title = "UTR Length Changes", subtitle = "Top 130 genes",
         x = "E16.5 (nt)", y = "E18.5 (nt)") +
    theme_minimal()
  
  ggsave(file.path(out_dir, "figures/utr_scatter_e16_vs_e18.png"), p2, 
         width = 8, height = 6)
}

cat("\nSCALPEL APA analysis complete\n")
