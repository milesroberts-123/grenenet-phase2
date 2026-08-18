#!/usr/bin/env Rscript

# Aggregate single-cell expression by cell type, normalize to CPM, and
# join to gene coordinates to produce a BED file for bedtools map.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(tidyr)
})

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
  return.seurat = TRUE
)

counts <- LayerData(pseudo, assay = "RNA", layer = "counts")

calculate_cpm <- function(counts) {
  lib_sizes <- colSums(counts)
  cpm <- sweep(counts, 2, lib_sizes, "/") * 1e6
  return(cpm)
}

cpm <- calculate_cpm(counts)

#stopifnot(all(colSums(cpm) == 1e6))

cpm <- as.data.frame(cpm)
cpm$gene <- rownames(cpm)

genes <- read_tsv(genes_path, col_names = c("chrom", "start", "end", "gene"))

genes <- genes %>%
  left_join(cpm, by = "gene") %>%
  drop_na() %>%
  arrange(chrom, start)

write.table(genes, out_path, quote = FALSE, row.names = FALSE, sep = "\t")
