#!/usr/bin/env Rscript

# Does 33_lfmm_kmate.R build valid LFMM inputs and join the results back?
#
# Fixture: 8 pools at generation 2 (two sites x four plots), 120 records across
# the three classes, two incomplete non-SNP records, and a climate table with
# bio1 for both sites. Runs the real runner against the vendored LFMM script and
# checks the row count, the joined columns and the recorded lambda.
#
# Needs the lfmm package, so run it where lfmm is installed (Savio).
# Run: Rscript lfmm_validate_kmate.R

suppressPackageStartupMessages(library(data.table))
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
runner <- file.path(here, "33_lfmm_kmate.R")
lfmm_script <- file.path(here, "run_lfmm_nogif.R")
stopifnot(file.exists(runner), file.exists(lfmm_script))

set.seed(3)
n_pools <- 8L
groups <- rep(c("snp", "smallindel", "sv"), c(40L, 100L, 30L))
n_rec <- length(groups)
pool_names <- c(paste0("1_2_", 1:4), paste0("2_2_", 1:4))
p0 <- runif(n_rec, 0.2, 0.8)
# pools x records, as in the k-selection fixture; the CSV is records x pools
af <- p0 + matrix(rnorm(n_pools * n_rec, 0, 0.05), n_pools, n_rec)
af <- pmin(pmax(af, 0.001), 0.999)
af[1, 45:46] <- NA          # records 41-140 are smallindel, so these are inside the model

trays <- as.data.table(t(af)); setnames(trays, pool_names)
seedmix <- data.table(chrom = "Chr3", pos = seq_along(groups), class = groups, X0_0_0 = p0)
climate <- data.table(site = rep(1:2, each = 5), year = rep(2018:2022, 2),
                      bio1 = c(10, 11, 12, 13, 14, 20, 21, 22, 23, 24))

dir <- file.path(tempdir(), "lfmm_kmate_fixture")
unlink(dir, recursive = TRUE); dir.create(dir, recursive = TRUE)
trays_csv <- file.path(dir, "trays.csv"); seedmix_csv <- file.path(dir, "seedmix.csv")
climate_rda <- file.path(dir, "climate.rda"); k_csv <- file.path(dir, "k.csv")
fwrite(trays, trays_csv); fwrite(seedmix, seedmix_csv)
climate_sites_yearly <- climate; save(climate_sites_yearly, file = climate_rda)
fwrite(data.table(generation = 2L, n_pools = n_pools, n_records = 130L,
                  component = 1L, eig_prop = 0.5, broken_stick = 0.2,
                  cum_var = 0.5, K_broken_stick = 2L, cum_var_at_K = 0.5), k_csv)

out_csv <- file.path(dir, "lfmm_out.csv")
st <- system2(file.path(R.home("bin"), "Rscript"),
              c("--no-save", runner, trays_csv, seedmix_csv, climate_rda, k_csv,
                "2", out_csv, dir, lfmm_script),
              stdout = TRUE, stderr = TRUE)
cat(paste(tail(st, 12), collapse = "\n"), "\n")
if (!file.exists(out_csv)) stop("runner produced no output")

r <- fread(out_csv)
print(dim(r)); print(names(r)); print(head(r[, .(chrom, pos, class, pval, z, K, lambda)], 3))

stopifnot(nrow(r) == 130L,
          all(c("chrom", "pos", "class", "pval", "B", "z", "generation", "K", "lambda") %in% names(r)),
          unique(r$generation) == 2L,
          unique(r$K) == 2L,
          is.finite(unique(r$lambda)[1]),
          all(is.finite(r$pval)),
          all(r$pval >= 0 & r$pval <= 1))
cat("\nplumbing OK\nLFMM_RUNNER_OK\n")
