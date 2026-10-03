#!/usr/bin/env Rscript

# Entropy-based gene specificity metrics (Psi, Psi_block, Zeta) from
# per-cell expression, following the ember construction (Swarna et al.
# 2026) with blocks = stage:CellType.
#
#   p_i   = x_i / sum(x)                  per gene across all cells
#   E_T   = -sum(p_i log2 p_i)            total entropy
#   E_C   = -sum(q_jC log2 q_jC)          entropy within block C
#   p_C   = sum(x in C) / sum(x)
#   E_W   = sum_C p_C * E_C               within-block entropy
#   Psi   = E_W / E_T
#   psi_C = p_C * E_C / E_W               block specificity, sums to 1
#   Zeta  = 1 - H(psi_C) / log2(r)        specificity to the partition
#
# x is the per-sample LogNormalize matrix from seurat_to_counts.R. Samples
# are streamed one at a time; only per-gene sufficient statistics are
# kept in memory. Genes detected in fewer than <min_cells> cells (count
# > 0 in the raw counts matrix across all samples) are dropped.
#
# Usage: ember_entropy.R <cells_dirs...> <min_cells> <out_metrics_csv> <out_psi_block_csv>
#   cells_dirs each contain counts.mtx, norm.mtx, cells.tsv, genes.tsv

suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: ember_entropy.R <cells_dirs...> <min_cells> <out_metrics_csv> <out_psi_block_csv>")
}

cells_dirs <- args[1:(length(args) - 3)]
min_cells <- as.integer(args[length(args) - 2])
metrics_path <- args[length(args) - 1]
psi_block_path <- args[length(args)]
stopifnot(length(cells_dirs) >= 1, min_cells >= 1)

x_log2 <- function(x) if_else(x > 0, x * log2(x), 0)

# Per-gene entropy from a pooled sum S = sum(x) and moment m = sum(x*log2(x)):
#   E = log2(S) - (m/S). Note m must be the sum of x*log2(x) contributions,
#   computed at cell level and pooled; pooling S and m is exact because both
#   are additive over cells. Requires x >= 0 (LogNormalize guarantees this).

suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: ember_entropy.R <cells_dirs...> <min_cells> <out_metrics_csv> <out_psi_block_csv>")
}

cells_dirs <- args[1:(length(args) - 3)]
min_cells <- as.integer(args[length(args) - 2])
metrics_path <- args[length(args) - 1]
psi_block_path <- args[length(args)]
stopifnot(length(cells_dirs) >= 1, min_cells >= 1)

entropy_from_moments <- function(S, m) {
  e <- log2(S) - m / S
  e[!is.finite(e)] <- 0  # zero-expression genes (S = 0) or single-cell blocks (NaN)
  e
}

# accumulate per-gene sufficient statistics (matrices are genes x cells)
S_g <- NULL      # sum of x across all cells
m_g <- NULL      # sum of x * log2(x) at cell level
n_cells_expressed <- NULL
blocks <- list() # per-block list of (S_C, m_C)
genes_ref <- NULL

for (d in cells_dirs) {
  cells <- read.delim(file.path(d, "cells.tsv"), stringsAsFactors = FALSE)
  genes <- readLines(file.path(d, "genes.tsv"))
  x <- readMM(file.path(d, "norm.mtx"))
  counts <- readMM(file.path(d, "counts.mtx"))
  stopifnot(nrow(x) == nrow(counts), ncol(x) == ncol(counts), nrow(x) == length(genes))
  stopifnot(nrow(cells) == ncol(x))
  stopifnot(!any(duplicated(cells$barcode)))

  if (is.null(genes_ref)) {
    genes_ref <- genes
  } else {
    stopifnot(identical(genes, genes_ref))
  }

  # m = sum over cells of x*log2(x), per gene; aggregate sparse entries by gene
  # row index. Column-oriented format stores nonzeros by column (cell), with
  # row indices in the i slot, so group by @i + 1.
  x_csc <- as(x, "CsparseMatrix")
  m_sample <- rowsum(x_csc@x * log2(x_csc@x), group = x_csc@i + 1L, reorder = FALSE)
  m_sample_vec <- double(length(genes_ref))
  m_sample_vec[as.integer(rownames(m_sample))] <- m_sample[, 1]

  S_g <- rowSums(x) + if (is.null(S_g)) 0 else S_g
  m_g <- m_sample_vec + if (is.null(m_g)) 0 else m_g
  n_cells_expressed <- as.integer(rowSums(counts > 0)) + if (is.null(n_cells_expressed)) 0 else n_cells_expressed

  # per-block statistics (x is genes x cells)
  for (b in unique(cells$block)) {
    idx <- which(cells$block == b)
    xb <- x[, idx, drop = FALSE]
    S_C <- rowSums(xb)
    xb_csc <- as(xb, "CsparseMatrix")
    m_block <- rowsum(xb_csc@x * log2(xb_csc@x), group = xb_csc@i + 1L, reorder = FALSE)
    m_C <- double(length(genes_ref))
    m_C[as.integer(rownames(m_block))] <- m_block[, 1]
    if (b %in% names(blocks)) {
      blocks[[b]]$S <- blocks[[b]]$S + S_C
      blocks[[b]]$m <- blocks[[b]]$m + m_C
    } else {
      blocks[[b]] <- list(S = S_C, m = m_C)
    }
  }
  message(sprintf("processed %s: %d cells, %d genes", d, nrow(cells), length(genes)))
}

E_T <- entropy_from_moments(S_g, m_g)

genes <- genes_ref

keep <- n_cells_expressed >= min_cells & S_g > 0
message(sprintf("%d of %d genes pass the %d-cell filter", sum(keep), length(genes), min_cells))

block_names <- names(blocks)
r <- length(block_names)

psi_block_num <- matrix(0, nrow = length(genes), ncol = r, dimnames = list(genes, block_names))
E_W <- double(length(genes))
for (b in block_names) {
  S_C <- blocks[[b]]$S
  m_C <- blocks[[b]]$m
  E_C <- entropy_from_moments(S_C, m_C)
  p_C <- ifelse(S_g > 0, S_C / S_g, 0)
  contrib <- p_C * E_C
  E_W <- E_W + contrib
  psi_block_num[, b] <- contrib
}

Psi <- ifelse(E_T > 0, E_W / E_T, NA_real_)

psi_block <- sweep(psi_block_num, 1, E_W, "/")
psi_block[!is.finite(psi_block)] <- 0

# Zeta = 1 - H(psi_blocks) / log2(r) with H in bits (log2), matching ember
H_psi <- -apply(psi_block, 1, function(p) {
  pl <- p[p > 0]
  if (length(pl) <= 1) return(0)
  sum(pl * log2(pl))
})
Zeta <- 1 - H_psi / log2(r)
Zeta[!is.finite(Zeta)] <- NA_real_

metrics <- tibble::tibble(
  gene = genes,
  E_T = E_T,
  E_W = E_W,
  Psi = Psi,
  Zeta = Zeta,
  n_cells_expressed = n_cells_expressed
)

write.table(metrics[keep, ], metrics_path, quote = FALSE, row.names = FALSE, sep = "\t")
write.table(data.frame(gene = genes[keep], psi_block[keep, , drop = FALSE], check.names = FALSE),
            psi_block_path, quote = FALSE, row.names = FALSE, sep = "\t")

message(sprintf("wrote %s and %s (%d genes, %d blocks)",
                metrics_path, psi_block_path, sum(keep), r))