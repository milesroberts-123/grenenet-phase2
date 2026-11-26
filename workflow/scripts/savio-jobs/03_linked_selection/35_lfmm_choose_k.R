#!/usr/bin/env Rscript

# LFMM stage 1c: turn the scree and the lambda sweep into the K actually used.
#
# Broken-stick counts components above an independence null, and linked records
# make it overshoot. The calibration sweep measures the genomic inflation factor
# per K; the criterion here is the smallest K whose lambda lands in a narrow band
# around 1, which is exactly the reference pipeline's rule ("the well-calibrated K
# is the smallest K with GIF ~ 1 and a flat null p-value histogram").
#
# Generations with no sweep fall back to their broken-stick K, so the output is a
# complete one-K-per-generation table either way. Every row records where its K
# came from, so a later reader can tell a calibrated value from a scree value.
#
# Usage: Rscript 35_lfmm_choose_k.R <k_selection_csv> <out_csv> [sweep_csv ...]

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("usage: <k_selection_csv> <out_csv> [sweep_csv ...]")
k_csv    <- args[[1]]
out_csv  <- args[[2]]
sweeps   <- args[-(1:2)]

lam_lo <- 0.95; lam_hi <- 1.05
ts <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

ks <- fread(k_csv)
base <- unique(ks[, .(generation, n_pools, n_records, K_broken_stick)])

sweep <- if (length(sweeps)) {
  rbindlist(lapply(sweeps, function(f) {
    if (!file.exists(f)) { ts("note: ", f, " does not exist, skipping"); return(NULL) }
    d <- fread(f)
    if (!"gif" %in% names(d)) { ts("note: ", f, " has no gif column, skipping"); return(NULL) }
    d
  }), fill = TRUE)
} else NULL

res <- base[, .(generation, n_pools, n_records,
                K_chosen = K_broken_stick, criterion = "broken-stick (no sweep)")]

if (!is.null(sweep) && nrow(sweep)) {
  for (g in unique(sweep$generation)) {
    d <- sweep[generation == g][order(K)]
    ok <- d[gif >= lam_lo & gif <= lam_hi]
    if (nrow(ok)) {
      pick <- ok[1]                                   # smallest calibrated K
      crit <- sprintf("calibrated: lambda %.3f in [%.2f, %.2f]", pick$gif, lam_lo, lam_hi)
    } else {
      pick <- d[which.min(abs(gif - 1))]
      crit <- sprintf("closest to lambda 1: %.3f at K=%d (none inside the band)",
                      pick$gif, pick$K)
    }
    ts("generation ", g, ": ", crit)
    res[generation == g, `:=`(K_chosen = as.integer(pick$K), criterion = crit)]
  }
}

res[, lambda_at_chosen := NA_real_]
if (!is.null(sweep) && nrow(sweep)) {
  m <- sweep[, .(generation, K, gif)][res[, .(generation, K = K_chosen)],
                                      on = .(generation, K)]
  res$lambda_at_chosen <- m$gif
}

fwrite(res, out_csv)
print(res)
ts("wrote ", out_csv)
