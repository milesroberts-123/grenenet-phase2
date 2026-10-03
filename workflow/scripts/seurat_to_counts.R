#!/usr/bin/env Rscript

# Export one Seurat object to per-cell matrices for entropy (ember-style)
# specificity metrics.
#
# Writes to <out_dir>:
#   norm.mtx   - genes x cells LogNormalize matrix (scale.factor = 1e4),
#                the input x for the entropy computation
#   counts.mtx - genes x cells raw counts (detection accounting)
#   cells.tsv  - barcode, CellType, block ("{stage}:{CellType}")
#   genes.tsv  - one gene name per line, row order of both matrices
#
# Usage: seurat_to_counts.R <seurat.rds> <out_dir>
# Stage is parsed from the rds filename between "GSE226097_" and the
# trailing date/asset tags, e.g. GSE226097_rosette_21d_230221.rds ->
# "rosette_21d".

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: seurat_to_counts.R <seurat.rds> <out_dir>")
}

rds_path <- args[1]
out_dir <- args[2]

stage <- sub("^GSE226097_", "", tools::file_path_sans_ext(basename(rds_path)))
stage <- sub("_[0-9]{6}$", "", stage)
stopifnot(nzchar(stage))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

seurat_obj <- readRDS(rds_path)

counts <- LayerData(seurat_obj, assay = "RNA", layer = "counts")
if (is.null(counts)) stop("no counts layer found in RNA assay")
stopifnot(all(counts@x >= 0))

cell_type <- as.character(seurat_obj@meta.data$CellType)
if (is.null(cell_type)) stop("no CellType column in meta.data")
cell_type[is.na(cell_type) | cell_type %in% c("", "NA")] <- "unknown"
stopifnot(length(unique(cell_type)) >= 5)

norm <- NormalizeData(
  seurat_obj,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 1e4,
  verbose = FALSE
)
x <- LayerData(norm, assay = "RNA", layer = "data")
stopifnot(all(x@x >= 0), all(is.finite(x@x)))

stopifnot(ncol(counts) == ncol(x), nrow(counts) == nrow(x))

Matrix::writeMM(counts, file.path(out_dir, "counts.mtx"))
Matrix::writeMM(x, file.path(out_dir, "norm.mtx"))

write.table(
  data.frame(
    barcode = colnames(seurat_obj),
    CellType = cell_type,
    block = paste(stage, cell_type, sep = ":")
  ),
  file.path(out_dir, "cells.tsv"),
  quote = FALSE, row.names = FALSE, sep = "\t"
)

writeLines(rownames(seurat_obj), file.path(out_dir, "genes.tsv"))

message(sprintf(
  "%s: %d genes x %d cells, %d blocks",
  stage, nrow(x), ncol(x), length(unique(paste(stage, cell_type, sep = ":")))
))