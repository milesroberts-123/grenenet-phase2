#!/usr/bin/env Rscript

# LFMM stage 1b: K by calibration, not by scree alone.
#
# The broken-stick criterion (stage 1) counts components above an independence
# null, and linked records make many small components look like signal: on this
# data it returns K = 13-38, against the reference's K = 16 for a similar number
# of pools. The reference pipeline also carries a sweep for exactly this reason
# (kMate run_lfmm_ksweep.R): fit a range of K and keep the smallest K whose
# genomic inflation factor is ~1 with a flat null p-value histogram.
#
# This script builds the generation's matrix with lfmm_input.R -- the same code
# the production fit uses -- and calls the vendored sweep on it.
#
# Usage: Rscript 34_lfmm_ksweep.R <trays> <seedmix> <climate_rda> <generation> \
#                                 <out_csv> <workdir> <ksweep_script> <K1,K2,...>

suppressPackageStartupMessages(library(data.table))
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(here, "lfmm_input.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8L) stop("usage: <trays> <seedmix> <climate_rda> <generation> <out_csv> <workdir> <ksweep_script> <Ks>")
trays_csv     <- args[[1]]
seedmix_csv   <- args[[2]]
climate_rda   <- args[[3]]
gen           <- as.integer(args[[4]])
out_csv       <- args[[5]]
workdir       <- args[[6]]
ksweep_script <- args[[7]]
Ks            <- args[[8]]

ts <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))
stem_info <- build_lfmm_stem(trays_csv, seedmix_csv, climate_rda, gen, workdir)

st <- system2(file.path(R.home("bin"), "Rscript"),
              c("--no-save", ksweep_script, stem_info$stem, out_csv, Ks),
              stdout = TRUE, stderr = TRUE)
cat(paste(tail(st, 25), collapse = "\n"), "\n")
if (!file.exists(out_csv)) stop("K sweep produced no output")

res <- fread(out_csv)
res[, `:=`(generation = gen, n_pools = stem_info$n_pools, n_records = stem_info$nr)]
fwrite(res, out_csv)
print(res)
ts("wrote ", out_csv)
