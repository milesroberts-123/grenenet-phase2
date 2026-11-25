#!/usr/bin/env Rscript

# Stage 2 of the Vw analysis: compare two locus panels before any per-site work.
#
# Vw's estimator (est_Va_bc) consumes the founder LD structure only through a
# single scalar, computed in averageLD():
#
#   s = ( sum(cov2cor(L)^2 * nR^(2*ngen)) - n ) / ( n * (n - 1) )
#
# so "which loci" reduces to "which loci estimate s". This script computes s for
# three nested contiguous windows (10k, 50k, 100k), so differences between them
# are purely the extra loci and reports the pair-distance composition of each.
#
# L is never materialised at 50k (50000^2 doubles = 20 GB). Because nR depends
# only on distance, s is accumulated block-wise:
#   s = 2 * sum_{i<j} cor2_ij * w(d_ij) / (n*(n-1))
# and the result is validated against Vw::averageLD on the 10k window, where the
# full matrices fit.
#
# Usage: Rscript 26_vw_panel_compare.R <genotype_matrix.tsv> <out_dir> <chrom>

suppressPackageStartupMessages({
  library(data.table)
  library(Vw)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) stop("usage: <genotype_matrix.tsv> <out_dir> <chrom>")
mat_file <- args[[1]]
out_dir <- args[[2]]
chrom <- args[[3]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Recombination: stdpopsim AraTha (TAIR10, Salome et al. 2012 map), already
# corrected for selfing inside stdpopsim (sigma = 0.97, Nordborg 2000):
# 8.06452e-10 per bp per generation -- the same value the simulation grid uses.
# ---------------------------------------------------------------------------
rrate   <- 8.06452e-10
hap_len <- 30427671

ts <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))
t0 <- proc.time()

ts("reading ", mat_file)
M <- fread(mat_file, header = FALSE, showProgress = FALSE)
pos <- M[[1]]
H <- t(as.matrix(M[, -1, with = FALSE]))   # haplotypes x loci
storage.mode(H) <- "double"
rm(M); invisible(gc())
nloci <- ncol(H)
# 231 founder accessions, each contributing its two phased haplotypes
if (nrow(H) != 462) stop("expected 462 haplotypes (231 x 2), got ", nrow(H))
ts("chromosome ", chrom, ": ", nrow(H), " haplotypes x ", nloci,
   " loci (MAF >= 0.05, fully called)")

# ---------------------------------------------------------------------------
# nR(d) exactly as Vw::form_nR builds it, evaluated through a distance lookup
# ---------------------------------------------------------------------------
nR_raw <- function(d, atleast) {
  if (!atleast) return(1 - 0.5 * (1 - exp(-2 * d * rrate)))
  maxn <- actuar::qztpois(1 - 1e-5, rrate * hap_len)
  out <- numeric(length(d))
  for (nt in 0:maxn) {
    pnt <- actuar::dztpois(nt, rrate * hap_len)
    if (pnt == 0) next
    acc <- numeric(length(d))
    for (ne in seq(0, maxn, 2)) {
      acc <- acc + dbinom(ne, prob = d / hap_len, size = nt) * pnt
    }
    out <- out + acc
  }
  out / actuar::pztpois(maxn, rrate * hap_len)
}

# weights: nR(d)^(2*ngen) on an exact 1-bp distance grid, so pairs use the same
# value Vw::form_nR would give them (distances are integer bp)
w_cache <- new.env(parent = emptyenv())
w_table <- function(ngen, atleast) {
  key <- paste(ngen, atleast)
  if (!is.null(w_cache[[key]])) return(w_cache[[key]])
  ts("building weight table ngen=", ngen, " atleast=", atleast)
  w <- nR_raw(0:hap_len, atleast)^(2 * ngen)
  w_cache[[key]] <- w
  w
}

# ---------------------------------------------------------------------------
# block-wise s
# ---------------------------------------------------------------------------
s_banded <- function(H, pos, ngen, atleast, block = 500L) {
  nh <- nrow(H)
  n  <- ncol(H)
  Hc <- H - rep(colMeans(H), each = nh)
  vv <- colSums(Hc^2) / (nh - 1)              # diagonal of L
  wtab <- w_table(ngen, atleast)              # exact weight per integer distance

  weighted <- 0        # sum of cor2 * w over pairs i<j
  plain    <- 0        # sum of cor2 over pairs i<j (diagnostic)
  band_pairs <- numeric(0)

  for (i0 in seq(1, n, by = block)) {
    i1  <- min(i0 + block - 1L, n)
    idx <- i0:i1
    Lb  <- crossprod(Hc[, idx, drop = FALSE], Hc) / (nh - 1)   # block x n
    r2  <- (Lb^2) / outer(vv[idx], vv)

    # pairs between this block and all later loci (no mask needed)
    if (i1 < n) {
      cols <- (i1 + 1L):n
      D <- abs(outer(pos[idx], pos[cols], "-"))
      r2s <- r2[, cols, drop = FALSE]
      weighted <- weighted + sum(r2s * wtab[round(D) + 1L])
      plain    <- plain + sum(r2s)
      rm(D, r2s)
    }
    # pairs inside the block (upper triangle of a 500 x 500 tile)
    if (length(idx) > 1L) {
      sub <- r2[, idx, drop = FALSE]
      ut  <- upper.tri(sub)
      Ds  <- abs(outer(pos[idx], pos[idx], "-"))
      weighted <- weighted + sum(sub[ut] * wtab[round(Ds[ut]) + 1L])
      plain    <- plain + sum(sub[ut])
    }
  }
  list(s = 2 * weighted / (n * (n - 1)),
       mean_cor2 = 2 * plain / (n * (n - 1)))
}

# ---------------------------------------------------------------------------
# window selection: one random start, 50k loci from there, 10k = first 10k
# ---------------------------------------------------------------------------
sizes <- c(w10k = 10000L, w50k = 50000L, w100k = 100000L)
need  <- max(sizes)
if (nloci < need) stop("need at least ", need, " loci, have ", nloci)
set.seed(20261124)                     # recorded, so the windows are reproducible
start_i <- sample(seq_len(nloci - need + 1L), 1L)
wmax <- start_i:(start_i + need - 1L)  # all three windows are nested in this one

res_file <- file.path(out_dir, "vw_panel_compare.csv")
prev <- if (file.exists(res_file)) tryCatch(fread(res_file), error = function(e) NULL) else NULL
skip_compute <- !is.null(prev) && nrow(prev) == 4L * length(sizes) &&
                all(prev$chrom == chrom)
if (skip_compute) ts("reusing existing results: ", res_file)

res <- list(); diag_bands <- NULL
for (nm in if (skip_compute) character(0) else names(sizes)) {
  idx <- wmax[seq_len(sizes[[nm]])]
  Hi  <- H[, idx, drop = FALSE]
  pi_ <- pos[idx]
  for (atleast in c(FALSE, TRUE)) {
    for (ngen in c(1L, 5L)) {
      ts(sprintf("%s  atleast=%s ngen=%d", nm, atleast, ngen))
      out <- s_banded(Hi, pi_, ngen = ngen, atleast = atleast)
      res[[length(res) + 1L]] <- data.frame(
        chrom = chrom, window = nm, n_loci = length(idx),
        start_pos = pi_[1], end_pos = pi_[length(pi_)],
        atleast_one_recomb = atleast, ngen = ngen,
        s = out$s, mean_cor2 = out$mean_cor2
      )
      if (nm == "w10k" && !atleast && ngen == 1L) diag_bands <- NULL
    }
  }
}

if (skip_compute) {
  res <- as.data.table(prev)
} else {
  res <- rbindlist(res)
  res[, va_scale := 1 / s]              # Va is proportional to 1/s
}
# the s values are the product, so they are written before the guard runs
fwrite(res, file.path(out_dir, "vw_panel_compare.csv"))

# ---------------------------------------------------------------------------
# validation: full-matrix averageLD on the 10k window must match the blockwise s
# ---------------------------------------------------------------------------
ts("validating against Vw::averageLD on the 10k window")
Hv <- H[, wmax[seq_len(sizes[["w10k"]])], drop = FALSE]
nh <- nrow(Hv)
Hc <- Hv - rep(colMeans(Hv), each = nh)
L  <- cov(Hc) * (nh - 1) / nh          # cov of the loci (columns)
pos10 <- pos[wmax[seq_len(sizes[["w10k"]])]]
nR_false <- form_nR(pos10, rrate, hap_len, FALSE)
nR_true  <- form_nR(pos10, rrate, hap_len, TRUE)

v_atleast <- c(FALSE, FALSE, TRUE)
v_ngen    <- c(1L, 5L, 1L)
v_full <- c(
  as.numeric(averageLD(L, nR_false, 1, ngen = 1L)),
  as.numeric(averageLD(L, nR_false, 1, ngen = 5L)),
  as.numeric(averageLD(L, nR_true,  1, ngen = 1L))
)
w10 <- res[res$window == "w10k", ]
v_banded <- w10$s[match(paste(v_atleast, v_ngen),
                        paste(w10$atleast_one_recomb, w10$ngen))]
chk <- data.frame(atleast = v_atleast, ngen = v_ngen,
                  s_full = v_full, s_banded = v_banded)
chk$rel_err <- abs(chk$s_banded - chk$s_full) / abs(chk$s_full)
rm(L, nR_false, nR_true); invisible(gc())
print(chk)
fwrite(chk, file.path(out_dir, "vw_panel_validation.csv"))
if (any(chk$rel_err > 1e-8)) {
  stop("banded s does not match averageLD: max rel err ", max(chk$rel_err))
}

print(res[, .(window, atleast_one_recomb, ngen, s, mean_cor2, va_scale)])
ts("DONE elapsed ", round((proc.time() - t0)["elapsed"] / 60, 1), " min")
