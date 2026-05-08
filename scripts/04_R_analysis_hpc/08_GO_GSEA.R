#!/usr/bin/env Rscript
# 08_GO_GSEA.R
# GO enrichment and GSEA for mature PC DE genes

library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(ggplot2)
library(enrichplot)
library(stringr)

setwd("/core/projects/GAP/GDA/gsamuel/scRNAseq")

outdir <- "results/04_R_analysis/08_GO_GSEA"
figdir <- "results/04_R_analysis/08_GO_GSEA/figures"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

de <- read.csv("results/04_R_analysis/07_DE_mature_PC/mature_PC_E18.5_vs_E16.5_all_genes.csv")

cat("Loaded genes:", nrow(de), "\n")
print(colnames(de))

# Background = all genes tested in DE
background <- bitr(
  de$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

# Split significant genes by direction
up_genes <- de %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0) %>%
  pull(gene)

down_genes <- de %>%
  filter(p_val_adj < 0.05, avg_log2FC < 0) %>%
  pull(gene)

cat("E18.5 up genes:", length(up_genes), "\n")
cat("E18.5 down genes:", length(down_genes), "\n")

run_go <- function(genes, name) {
  if (length(genes) < 10) {
    cat("Skipping", name, "- too few genes\n")
    return(NULL)
  }

  gene_ids <- bitr(
    genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Mm.eg.db
  )

  go <- enrichGO(
    gene = gene_ids$ENTREZID,
    universe = background$ENTREZID,
    OrgDb = org.Mm.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )

  if (!is.null(go) && nrow(go@result) > 0) {
    write.csv(
      as.data.frame(go),
      file.path(outdir, paste0("GO_BP_", name, ".csv")),
      row.names = FALSE
    )

    go_simple <- simplify(go, cutoff = 0.7, by = "p.adjust", select_fun = min)

    p <- dotplot(go_simple, showCategory = 12) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
      ggtitle(paste("GO BP", name)) +
      theme_bw(base_size = 14)

    ggsave(
      file.path(figdir, paste0("GO_BP_", name, "_dotplot.png")),
      p,
      width = 11,
      height = 8,
      dpi = 300,
      bg = "white"
    )

    ggsave(
      file.path(figdir, paste0("GO_BP_", name, "_dotplot.pdf")),
      p,
      width = 11,
      height = 8
    )
  }

  return(go)
}

go_up <- run_go(up_genes, "E18_up")
go_down <- run_go(down_genes, "E18_down")

# GSEA ranking
# If there is no Wald/stat column, use direction + significance
de <- de %>%
  filter(!is.na(gene), !is.na(avg_log2FC), !is.na(p_val)) %>%
  mutate(rank_metric = avg_log2FC * -log10(p_val + 1e-300))

gene_map <- bitr(
  de$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

rank_df <- de %>%
  select(gene, rank_metric) %>%
  inner_join(gene_map, by = c("gene" = "SYMBOL")) %>%
  group_by(ENTREZID) %>%
  slice_max(order_by = abs(rank_metric), n = 1) %>%
  ungroup()

gene_list <- rank_df$rank_metric
names(gene_list) <- rank_df$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

cat("Genes used for GSEA:", length(gene_list), "\n")

gsea_go <- gseGO(
  geneList = gene_list,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE
)

if (!is.null(gsea_go) && nrow(gsea_go@result) > 0) {
  gsea_df <- as.data.frame(gsea_go)

  write.csv(
    gsea_df,
    file.path(outdir, "GSEA_GO_BP.csv"),
    row.names = FALSE
  )

  # Dotplot split by direction
  p1 <- dotplot(gsea_go, showCategory = 12, split = ".sign") +
    facet_grid(. ~ .sign) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
    ggtitle("GSEA GO BP") +
    theme_bw(base_size = 14)

  ggsave(
    file.path(figdir, "GSEA_GO_BP_dotplot.png"),
    p1,
    width = 13,
    height = 8,
    dpi = 300,
    bg = "white"
  )

  ggsave(
    file.path(figdir, "GSEA_GO_BP_dotplot.pdf"),
    p1,
    width = 13,
    height = 8
  )

  # Bubble plot: red = E18.5, blue = E16.5
  bubble_df <- gsea_df %>%
    arrange(p.adjust) %>%
    head(20) %>%
    mutate(
      Description = str_wrap(Description, width = 40),
      Description = factor(Description, levels = rev(Description))
    )

  p2 <- ggplot(bubble_df, aes(x = NES, y = Description, size = setSize, color = NES)) +
    geom_point(alpha = 0.9) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(
      title = "GSEA GO BP Bubble Plot",
      x = "NES",
      y = "",
      size = "Gene set size",
      color = "NES"
    ) +
    theme_bw(base_size = 14)

  ggsave(
    file.path(figdir, "GSEA_GO_BP_bubble.png"),
    p2,
    width = 12,
    height = 8,
    dpi = 300,
    bg = "white"
  )

  ggsave(
    file.path(figdir, "GSEA_GO_BP_bubble.pdf"),
    p2,
    width = 12,
    height = 8
  )
}

cat("GO/GSEA analysis complete\n")
cat("Results saved to:", outdir, "\n")
