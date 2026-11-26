#!/usr/bin/env Rscript

# How large an effect could the pipeline have missed?
#
# The positive control shows a planted bio1 effect of size 0.12 is recovered every
# time. That alone leaves "no signal" open-ended: it says nothing about effects
# smaller than the one planted. This walks the effect size down and reports the
# detection rate, turning a null into a bound.
#
# Fixture is the same shape as the control (200 pools, 20k records, two latent
# factors), with 50 planted records per effect size. Detection is p < 1e-5, the
# threshold used when reading the real sweep.
#
# Needs lfmm; run it where lfmm is installed (Savio).
# Run: Rscript lfmm_power_curve.R [out_csv]

suppressPackageStartupMessages({library(lfmm); library(data.table)})

args <- commandArgs(trailingOnly = TRUE)
out_csv <- if (length(args) >= 1L) args[[1]] else "lfmm_power_curve.csv"

set.seed(22)
n_pool <- 200L; n_rec <- 20000L; n_sig <- 50L
effects <- c(0, 0.01, 0.02, 0.04, 0.08, 0.12)

X <- rnorm(n_pool)
U <- matrix(rnorm(n_pool * 2L), n_pool, 2L)
V <- matrix(rnorm(2L * n_rec, 0, 0.05), 2L, n_rec)
base <- U %*% V + matrix(rnorm(n_pool * n_rec, 0, 0.05), n_pool, n_rec)
sig <- sample.int(n_rec, n_sig)

res <- rbindlist(lapply(effects, function(b) {
  Y <- base
  if (b > 0) Y[, sig] <- Y[, sig] + outer(X, rnorm(n_sig, b, 0))
  mod <- lfmm_ridge(Y = Y, X = X, K = 2L)
  pv  <- lfmm_test(Y = Y, X = X, lfmm = mod, calibrate = "gif")
  p <- as.numeric(pv$pvalue)
  det <- mean(p[sig] < 1e-5)
  null_rate <- mean(p[-sig] < 1e-5)
  cat(sprintf("effect %.3f  lambda %.3f  detected %3.0f%%  null rate %.2e\n",
              b, pv$gif, 100 * det, null_rate))
  data.table(effect = b, lambda = as.numeric(pv$gif),
             detection_rate = det, null_rate = null_rate, n_sig = n_sig)
}))

fwrite(res, out_csv)
cat("\nwrote ", out_csv, "\n")
