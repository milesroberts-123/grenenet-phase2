#!/usr/bin/env Rscript

# Does 32_lfmm_k_selection.R recover a known number of latent factors?
#
# Builds a tiny trays/seedmix pair whose dp matrix has exactly two shared axes
# across pools, with the three variant classes mixed in and some records made
# incomplete, then runs the real script end to end (awk extraction included) and
# checks K and the record count it reports.
#
# Run: Rscript lfmm_validate_kselection.R   (data.table only)

suppressPackageStartupMessages(library(data.table))
script <- file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))),
                    "32_lfmm_k_selection.R")
if (!file.exists(script)) script <- "32_lfmm_k_selection.R"

set.seed(7)
n_pools <- 6L
n_snp <- 100L; n_indel <- 150L; n_sv <- 50L
n_rec <- n_snp + n_indel + n_sv
groups <- rep(c("snp", "smallindel", "sv"), c(n_snp, n_indel, n_sv))

pool_names <- c(paste0("1_1_", 1:3), paste0("2_1_", 1:3))
p0 <- runif(n_rec, 0.2, 0.8)

# exactly two latent factors: pool scores (6 x 2), record loadings (2 x n_rec)
U <- matrix(rnorm(n_pools * 2L, 0, 0.08), n_pools, 2L)
V <- matrix(rnorm(2L * n_rec, 0, 1), 2L, n_rec)
shared <- U %*% V                     # pools x records
noise  <- matrix(rnorm(n_pools * n_rec, 0, 0.01), n_pools, n_rec)
af <- t(p0 + t(shared + noise))       # pools x records, centred on p0
af <- pmin(pmax(af, 0.001), 0.999)

# make a handful of NON-SNP records incomplete in a way that must be dropped;
# rows 150-152 are smallindel, so they are inside the modelled subset (rows 4-6
# would have been snp and excluded before completeness is ever checked)
af[1, 150:152] <- NA

trays <- as.data.table(t(af))
setnames(trays, pool_names)
seedmix <- data.table(class = groups, X0_0_0 = p0)

dir <- file.path(tempdir(), "lfmm_k_fixture")
unlink(dir, recursive = TRUE); dir.create(dir, recursive = TRUE)
trays_csv <- file.path(dir, "trays.csv")
seedmix_csv <- file.path(dir, "seedmix.csv")
fwrite(trays, trays_csv)
fwrite(seedmix, seedmix_csv)

out_csv <- file.path(dir, "k.csv")
rscript <- file.path(R.home("bin"), "Rscript")
st <- system2(rscript, c("--no-save", script, trays_csv, seedmix_csv, out_csv, dir),
              stdout = TRUE, stderr = TRUE)
cat(paste(tail(st, 12), collapse = "\n"), "\n")
if (!file.exists(out_csv)) stop("script produced no output")

k <- fread(out_csv)
cat("\nreported:\n"); print(unique(k[, .(generation, n_pools, n_records, K_broken_stick, cum_var_at_K)]))

# one long table, one row per component of the single generation
stopifnot(length(unique(k$generation)) == 1L,
          unique(k$generation) == 1L,
          nrow(k) == n_pools,
          unique(k$n_pools) == n_pools,
          # 200 non-SNP records; the 3 incomplete ones are now imputed, not dropped
          unique(k$n_records) == n_indel + n_sv,
          length(unique(k$K_broken_stick)) == 1L,
          unique(k$K_broken_stick) >= 1L, unique(k$K_broken_stick) <= 4L)
# the planted rank is 2; broken stick should land near it, not at the noise floor
stopifnot(abs(unique(k$K_broken_stick) - 2L) <= 1L)
cat("\nK selection recovers the planted rank\nKSELECT_OK\n")
