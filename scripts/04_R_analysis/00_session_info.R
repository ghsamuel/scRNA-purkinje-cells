#!/usr/bin/env Rscript
# Document software versions for reproducibility

library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(patchwork)

# Try loading CellChat
tryCatch({
  library(CellChat)
}, error = function(e) {
  message("CellChat not installed - run: devtools::install_github('sqjin/CellChat')")
})

# Print to screen
cat("\n=== Software Versions ===\n")
cat("R:", as.character(getRversion()), "\n")
cat("Seurat:", as.character(packageVersion("Seurat")), "\n")
cat("harmony:", as.character(packageVersion("harmony")), "\n")

if ("CellChat" %in% loadedNamespaces()) {
  cat("CellChat:", as.character(packageVersion("CellChat")), "\n")
}

cat("\nCell Ranger: v9.0.1\n")
cat("Reference: GRCm39\n")
cat("Platform: UConn HPC\n")

# Save detailed session info
dir.create("results", showWarnings = FALSE)
sink("results/session_info.txt")
cat("Purkinje Cell Analysis - Session Information\n")
cat("Generated:", format(Sys.time()), "\n\n")
sessionInfo()
sink()

cat("\nSession info saved to results/session_info.txt\n")
