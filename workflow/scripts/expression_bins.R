#!/usr/bin/env Rscript

# Aggregate single-cell expression by cell type, normalize to CPM, and
# join to gene coordinates to produce a BED file for bedtools map.
# CPM values are taken from the returned object's "data" layer, which
# Seurat populates with RelativeCounts (scale.factor = 1e6). Requires
# Seurat >= 5.3.0: earlier versions silently drop scale.factor here and
# would produce columns summing to 1e4 instead of 1e6.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(tidyr)
})

stopifnot(packageVersion("Seurat") >= "5.3.0")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: expression_bins.R <seurat.rds> <genes.bed> <out.bed> <group_by>")
}

rds_path <- args[1]
genes_path <- args[2]
out_path <- args[3]
group_by <- args[4]

seurat_obj <- readRDS(rds_path)

pseudo <- AggregateExpression(
  seurat_obj,
  group.by = c(group_by),
  assays = "RNA",
  normalization.method = "RC",
  scale.factor = 1e6,
  return.seurat = TRUE
)

cpm <- LayerData(pseudo, assay = "RNA", layer = "data")

stopifnot(all(abs(colSums(cpm) - 1e6) < 0.01))

cpm <- as.data.frame(cpm)
cpm$gene <- rownames(cpm)

genes <- read_tsv(genes_path, col_names = c("chrom", "start", "end", "gene"))

genes <- genes %>%
  left_join(cpm, by = "gene") %>%
  drop_na() %>%
  arrange(chrom, start)

write.table(genes, out_path, quote = FALSE, row.names = FALSE, sep = "\t")
