#!/usr/bin/env Rscript

# LD scalar s for nested panel sizes, one-generation windows only (stage 6a).
#
# The 100k chr3 panel is a contiguous window, so its first n loci are a nested
# sub-panel: s(10k), s(50k) and s(100k) describe the same window at three
# densities, which is what the window-size comparison needs. Only ngen = 1 is
# computed -- the analysis window set is one generation, where nR^(2*ngen) is
# exact.
#
# Usage: Rscript 30_vw_subsets_s.R <founder_cache.tsv> <panel_csv> <out_csv>

suppressPackageStartupMessages({
  library(data.table)
  library(Vw)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) stop("usage: <founder_cache.tsv> <panel_csv> <out_csv>")
cache_file <- args[[1]]
panel_csv  <- args[[2]]
out_csv    <- args[[3]]

rrate   <- 8.06452e-10
hap_len <- 30427671
sizes   <- c(10000L, 50000L, 100000L)

ts <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

C <- fread(cache_file, header = FALSE, showProgress = FALSE)
cache_pos <- C[[1]]
if (ncol(C) - 1L != 462L) stop("founder cache has ", ncol(C) - 1L, " haplotypes, expected 462")

panel <- fread(panel_csv, select = c("chrom", "pos"))
idx <- match(panel$pos, cache_pos)
if (anyNA(idx)) stop(sum(is.na(idx)), " panel loci are absent from the founder cache")
ts("panel ", nrow(panel), " loci, all present in the cache")

# s = mean over pairs of cor^2 * N^2, N the one-generation non-recombination
# probability; blockwise so the n x n matrix is never formed
w1 <- (1 - 0.5 * (1 - exp(-2 * (0:hap_len) * rrate)))^2

s_blockwise <- function(H, pos, block = 500L) {
  nh <- nrow(H); n <- ncol(H)
  Hc <- H - rep(colMeans(H), each = nh)
  vv <- colSums(Hc^2) / (nh - 1)
  acc <- 0
  for (i0 in seq(1, n, by = block)) {
    i1 <- min(i0 + block - 1L, n); k <- i0:i1
    Lb <- crossprod(Hc[, k, drop = FALSE], Hc) / (nh - 1)
    r2 <- (Lb^2) / outer(vv[k], vv)
    if (i1 < n) {
      cols <- (i1 + 1L):n
      D <- abs(outer(pos[k], pos[cols], "-"))
      acc <- acc + sum(r2[, cols, drop = FALSE] * w1[round(D) + 1L])
    }
    if (length(k) > 1L) {
      sub <- r2[, k, drop = FALSE]; ut <- upper.tri(sub)
      Ds <- abs(outer(pos[k], pos[k], "-"))
      acc <- acc + sum(sub[ut] * w1[round(Ds[ut]) + 1L])
    }
  }
  2 * acc / (n * (n - 1))
}

out <- rbindlist(lapply(sizes, function(n) {
  if (n > nrow(panel)) return(NULL)
  ts("s for the first ", n, " loci")
  H <- t(as.matrix(C[idx[seq_len(n)], -1, with = FALSE]))
  storage.mode(H) <- "double"
  data.table(n_loci = n, ngen = 1L, s = s_blockwise(H, panel$pos[seq_len(n)]))
}))

fwrite(out, out_csv)
print(out)
ts("DONE ", out_csv)
