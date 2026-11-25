#!/usr/bin/env Rscript

# Does the locus-permutation null actually collapse a known signal?
#
# Fixture: 8 replicates, 5000 loci, a shared response across replicates (the
# biological signal) plus independent noise. Va should be clearly positive, and
# the permutation null should sit at ~0. If the null does not collapse here, it
# cannot be used on real data.
#
# Run: Rscript vw_validate_perm_null.R   (Vw only; no asreml fit)

suppressPackageStartupMessages(library(Vw))

set.seed(11)
nl <- 5000L; nr <- 8L
p1 <- matrix(runif(nr * nl, 0.2, 0.8), nr)
shared <- rnorm(nl, 0, 0.05)
p2 <- p1 + rep(shared, each = nr) + matrix(rnorm(nr * nl, 0, 0.01), nr)
p2 <- pmin(pmax(p2, 0.001), 0.999)

s_test <- 0.02   # stand-in for the LD scalar

va_from <- function(a, b, s) {
  Qm <- Q(a, b)
  q  <- Qm[upper.tri(Qm, diag = TRUE)]
  aa <- rep(0.5 * s, length(q))
  bb <- diag(nrow(Qm))[upper.tri(Qm, diag = TRUE)]
  unname(coef(lm(q ~ aa + bb - 1))["aa"])
}

va_null <- function(a, b, s) {
  perm <- t(vapply(seq_len(nrow(b)), function(i) sample.int(ncol(b)), integer(ncol(b))))
  bn <- b
  an <- a
  for (i in seq_len(nrow(b))) {
    bn[i, ] <- b[i, perm[i, ]]
    an[i, ] <- a[i, perm[i, ]]
  }
  va_from(an, bn, s)
}

obs <- va_from(p1, p2, s_test)
nul <- vapply(seq_len(100L), function(r) va_null(p1, p2, s_test), numeric(1))

cat(sprintf("observed Va      = %.6f\n", obs))
cat(sprintf("null mean        = %.6f (+/- %.6f)\n", mean(nul), sd(nul)))
cat(sprintf("null 97.5%% quant = %.6f\n", unname(quantile(nul, 0.975))))
cat(sprintf("null >= observed : %d of %d\n", sum(nul >= obs), length(nul)))

stopifnot(obs > 0, mean(nul) < 0.2 * obs)
cat("PERM_NULL_OK\n")
