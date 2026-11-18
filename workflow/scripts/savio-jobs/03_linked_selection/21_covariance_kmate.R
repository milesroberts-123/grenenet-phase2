#!/usr/bin/env Rscript

# Variance-covariance matrix of the full kmate increment matrix, uncorrected and
# corrected for finite sample size.
#
# The full matrix is stats::cov() over all 541 interval columns, so it carries
# every within-series, cross-replicate and cross-site covariance. The sample-size
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

ts_message("computing covariances (pairwise complete) ...")
covmat <- stats::cov(pdiffs, use = "pairwise.complete.obs")
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
