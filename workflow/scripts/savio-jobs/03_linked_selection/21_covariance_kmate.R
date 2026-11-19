#!/usr/bin/env Rscript

# Variance-covariance matrix of the full kmate increment matrix, uncorrected and
# corrected for finite sample size.
#
# The full matrix is the pairwise-complete covariance over all 893 interval
# columns, so it carries every within-series, cross-replicate and cross-site
# covariance. The sample-size
# correction is then applied only to the 352 within-(site, replicate) blocks: the
# sampling noise of two different replicates (or sites) is independent, so their
# covariances are unbiased as they stand and correcting them would remove real
# signal. Each block is corrected by rventures::correct_covmat_for_n() with the
# seedmix as the interval-start frequencies and total_flowers * 2 chromosomes per
# pool.
#
# Outputs (in data/):
#   ..._increment_covariances_kmate.tsv            (uncorrected)
#   ..._increment_covariances_kmate_corrected.tsv  (within-series blocks only)
#   ..._increment_covariances_kmate_qc.csv         (per-series correction record)

suppressMessages({
  library(data.table)
  library(rventures)
})

ts_message <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))
}

data_dir <- "/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/data/"
af_prefix <- paste0(data_dir, "PHASE_1_2_COMBINED_kmate_arch3v2_allele_frequency")

increments_csv <- paste0(af_prefix, "_increment_matrix.csv")
plan_csv       <- paste0(af_prefix, "_increment_plan.csv")
trays_csv      <- paste0(af_prefix, "_trays_randswitch_fixrm_biallelic.csv")
seedmix_csv    <- paste0(af_prefix, "_seedmix_randswitch_fixrm_biallelic.csv")
flower_csv     <- paste0(af_prefix, "_flower_counts_sitesfixed.csv")

out_uncorrected <- paste0(af_prefix, "_increment_covariances_kmate.tsv")
out_corrected   <- paste0(af_prefix, "_increment_covariances_kmate_corrected.tsv")
out_qc          <- paste0(af_prefix, "_increment_covariances_kmate_qc.csv")

for (f in c(increments_csv, plan_csv, trays_csv, seedmix_csv, flower_csv)) {
  if (!file.exists(f)) stop("missing input: ", f)
}

# sample size of the seedmix pool, matching the value used by the G(t) jobs
SEEDMIX_N <- 10000

## ---- uncorrected matrix ------------------------------------------------
ts_message("reading increments ...")
pdiffs <- fread(increments_csv)
ts_message("increments: ", format(nrow(pdiffs), big.mark = ","), " rows x ",
           ncol(pdiffs), " interval columns")

ts_message("computing covariances (pairwise complete, chunked) ...")

# stats::cov() cannot be used here: 6,207,205 x 893 doubles is a long vector
# (> 2^31 elements) and cov.c has no long-vector support ("long vectors not
# supported yet"). BLAS is int-indexed as well, so crossprod() over the whole
# matrix fails the same way. Accumulate the pairwise-complete covariance over
# row chunks instead:
#
#   cov_ij = (sum_ij x_i x_j - (sum_ij x_i)(sum_ij x_j)/n_ij) / (n_ij - 1)
#
# where sum_ij and n_ij are restricted to the rows where both columns are
# present, which is what use = "pairwise.complete.obs" computes.
p_int  <- ncol(pdiffs)
ints   <- colnames(pdiffs)
acc_xy <- matrix(0, p_int, p_int)   # sum of x_i x_j over shared rows
acc_x  <- matrix(0, p_int, p_int)   # sum of x_i over rows where x_j is present
acc_n  <- matrix(0, p_int, p_int)   # shared-row counts
chunk  <- 200000L
starts <- seq.int(1L, nrow(pdiffs), by = chunk)

for (k in seq_along(starts)) {
  idx <- starts[k]:min(starts[k] + chunk - 1L, nrow(pdiffs))
  # pdiffs is a data.table (fread), so a single index selects rows; a plain
  # matrix would be linear-indexed here, so fail loudly instead of silently
  # building a one-column chunk
  x <- as.matrix(pdiffs[idx])
  if (ncol(x) != p_int) {
    stop("chunk extraction returned ", ncol(x), " of ", p_int,
         " columns; pdiffs must be a data.table")
  }
  present <- !is.na(x)
  x[!present] <- 0                 # NA -> 0 makes every sum "both present"
  storage.mode(present) <- "double"

  acc_xy <- acc_xy + crossprod(x)
  acc_x  <- acc_x + crossprod(x, present)
  acc_n  <- acc_n + crossprod(present)

  if (k %% 5L == 0L) {
    ts_message("  chunk ", k, " of ", length(starts))
    gc(verbose = FALSE)
  }
}

covmat <- (acc_xy - acc_x * t(acc_x) / acc_n) / (acc_n - 1)
dimnames(covmat) <- list(ints, ints)
rm(acc_xy, acc_x, acc_n); gc(verbose = FALSE)
rm(pdiffs); gc(verbose = FALSE)

# written with the interval names as an explicit first column so a later
# fread() cannot mistake the row names for data
fwrite(data.table(interval = rownames(covmat), as.data.frame(covmat)),
       out_uncorrected, sep = "\t")
ts_message("wrote ", basename(out_uncorrected), " (", nrow(covmat), " x ",
           ncol(covmat), ")")

## ---- correction of the within-series blocks ----------------------------
ts_message("reading the frequency table for the correction ...")
gen0 <- fread(seedmix_csv, select = "X0_0_0")
trays <- fread(trays_csv, showProgress = FALSE)
flower_counts <- as.data.frame(fread(flower_csv))
plan <- fread(plan_csv)

if (!all(plan$interval %in% colnames(covmat))) {
  stop("interval plan and covariance matrix describe different columns")
}

corrected <- covmat
series <- unique(plan[, .(site, rep)])[order(site, rep)]
qc <- vector("list", nrow(series))
t0 <- Sys.time()

for (k in seq_len(nrow(series))) {
  s <- series$site[k]
  r <- series$rep[k]
  p <- plan[site == s & rep == r][order(to_gen)]

  idx <- match(p$interval, colnames(covmat))
  block <- covmat[idx, idx, drop = FALSE]

  # frequencies at each time point of this series, in interval order:
  # gen 0 from the seedmix, then the end generation of every interval
  pmat <- cbind(gen0$X0_0_0,
                as.matrix(trays[, p$to_col, with = FALSE]))
  colnames(pmat) <- c(paste0(s, "_0_0"), p$to_col)

  n_pool <- flower_counts$total_flowers[match(p$to_col, flower_counts$pool)] * 2
  if (anyNA(n_pool)) {
    stop("no sample size for site ", s, " rep ", r, ": ",
         paste(p$to_col[is.na(n_pool)], collapse = ", "))
  }
  n_vec <- c(SEEDMIX_N, n_pool)
  if (length(n_vec) != ncol(pmat)) stop("sample-size vector length mismatch")

  block_corrected <- rventures::correct_covmat_for_n(
    pmat = pmat, covmat = block, n = n_vec,
    input_asin_trans = FALSE, shared_sample_noise = TRUE
  )

  corrected[idx, idx] <- block_corrected

  n_clamped <- attr(block_corrected, "n_var_clamped")
  qc[[k]] <- data.table(
    site = s, rep = r,
    n_intervals = length(idx),
    n_clamped = if (is.null(n_clamped)) 0L else as.integer(n_clamped),
    min_var = min(diag(block_corrected), na.rm = TRUE),
    max_var = max(diag(block_corrected), na.rm = TRUE),
    n_na = sum(is.na(block_corrected))
  )

  rm(block, block_corrected, pmat)
  if (k %% 50L == 0L) {
    ts_message("corrected ", k, " of ", nrow(series), " series | elapsed ",
               round(as.numeric(difftime(Sys.time(), t0, units = "min")), 1), " min")
    gc(verbose = FALSE)
  }
}

fwrite(data.table(interval = rownames(corrected), as.data.frame(corrected)),
       out_corrected, sep = "\t")
fwrite(rbindlist(qc), out_qc)

ts_message("wrote ", basename(out_corrected))
ts_message("series corrected: ", nrow(series),
           " | clamped variances: ", sum(rbindlist(qc)$n_clamped))
ts_message("done")
