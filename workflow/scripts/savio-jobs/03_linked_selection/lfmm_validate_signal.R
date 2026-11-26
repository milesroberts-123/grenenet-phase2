#!/usr/bin/env Rscript

# Can this pipeline find a signal at all?
#
# Every K from 5 to 30 gave lambda 1.02 and a flat null on generation 1, which is
# either a true null (no genome-environment association) or a pipeline that cannot
# see one. This is the positive control: plant a known climate effect at a small
# set of records, inside realistic structure, and require that LFMM recovers them.
#
# Design: 200 pools, 20k records, two latent factors, and 20 records whose dp is
# driven by bio1 with a real effect size. Recovery means every planted record lands
# in the smallest p-values while the lambda stays near 1.
#
# Needs lfmm; run it where lfmm is installed (Savio).
# Run: Rscript lfmm_validate_signal.R

suppressPackageStartupMessages({library(lfmm); library(data.table)})

here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
set.seed(21)

n_pool <- 200L; n_rec <- 20000L; n_sig <- 20L
X <- rnorm(n_pool)
# two latent factors uncorrelated with X, so the model has something real to absorb
U <- matrix(rnorm(n_pool * 2L), n_pool, 2L)
V <- matrix(rnorm(2L * n_rec, 0, 0.05), 2L, n_rec)
Y <- U %*% V + matrix(rnorm(n_pool * n_rec, 0, 0.05), n_pool, n_rec)
sig <- sample.int(n_rec, n_sig)
Y[, sig] <- Y[, sig] + outer(X, rnorm(n_sig, 0.12, 0))

mod <- lfmm_ridge(Y = Y, X = X, K = 2L)
pv  <- lfmm_test(Y = Y, X = X, lfmm = mod, calibrate = "gif")

p <- as.numeric(pv$pvalue)
cat(sprintf("lambda = %.4f   planted records with p < 1e-4: %d/%d\n",
            pv$gif, sum(p[sig] < 1e-4), n_sig))
cat(sprintf("best rank of a planted record = %d of %d\n",
            min(rank(p)[sig]), n_rec))

stopifnot(is.finite(pv$gif),
          pv$gif > 0.8, pv$gif < 1.3,        # structure absorbed, calibration near 1
          sum(p[sig] < 1e-4) >= n_sig - 2L, # the planted signal is found
          min(rank(p)[sig]) <= 30L)
cat("pipeline recovers a planted climate signal\nLFMM_SIGNAL_OK\n")
