#!/usr/bin/env Rscript

# Checks that the per-site shortcut in 28_vw_site_va.R reproduces Vw's own
# estimator. That script skips L and nR (they never fit at 100k loci) and uses
# only the scalar s plus the three-line regression est_Va_bc runs internally, so
# the equivalence has to be demonstrated on a case small enough for the real
# matrices.
#
# Run: Rscript vw_validate_mimic.R      (needs the asreml seat)

suppressPackageStartupMessages({
  library(Vw)
  library(data.table)
})

rrate <- 8.06452e-10
hap_len <- 30427671

set.seed(5)
nloci <- 40L
nind  <- 30L
nrep  <- 5L

# founder haplotypes -> individual genotypes, as analyse_parents() does
H <- matrix(rbinom(nind * 2L * nloci, 1, 0.3), nind * 2L, nloci)
pat <- seq(1, nind * 2L, 2L)
mat <- pat + 1L
c0  <- (H[pat, ] + H[mat, ]) / 2
L   <- cov(c0) * (nind - 1L) / nind
snps <- sort(sample.int(3e7, nloci))
nR  <- form_nR(snps, rrate, hap_len, FALSE)

# two time points across replicate populations
pbar1 <- matrix(runif(nrep * nloci, 0.1, 0.9), nrep, nloci)
pbar2 <- pbar1 + matrix(rnorm(nrep * nloci, 0, 0.05), nrep, nloci)
pbar2[pbar2 < 0.01] <- 0.01
pbar2[pbar2 > 0.99] <- 0.99

ref <- est_Va_bc(pbar1, pbar2, L = L, nR = nR, ngen = 1L)

# the shortcut used per site
s  <- averageLD(L, nR, nrep = nrep, ngen = 1L)
Qm <- Q(pbar1, pbar2)
q  <- Qm[upper.tri(Qm, diag = TRUE)]
a  <- rep(0.5 * s, length(q))
b  <- diag(nrep)[upper.tri(Qm, diag = TRUE)]
fit <- lm(q ~ a + b - 1)

va_mine <- unname(coef(fit)["a"])
ne_mine <- 0.5 / unname(coef(fit)["b"])
rel_va  <- abs(va_mine - ref$vA_BC) / abs(ref$vA_BC)
rel_ne  <- abs(ne_mine - ref$Ne_BC) / abs(ref$Ne_BC)

cat(sprintf("reference  Va=%.12g  Ne=%.12g\n", ref$vA_BC, ref$Ne_BC))
cat(sprintf("shortcut   Va=%.12g  Ne=%.12g\n", va_mine, ne_mine))
cat(sprintf("relative error  Va=%.2e  Ne=%.2e\n", rel_va, rel_ne))

stopifnot(rel_va < 1e-10, rel_ne < 1e-10)
cat("MIMIC_OK\n")
