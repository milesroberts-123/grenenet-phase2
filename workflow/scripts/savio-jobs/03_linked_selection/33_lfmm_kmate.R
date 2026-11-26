#!/usr/bin/env Rscript

# LFMM stage 2: fit the reference's LFMM for one generation of the non-SNP
# allele-frequency-change matrix.
#
# The matrix is built by lfmm_input.R, the same code the K calibration sweep uses,
# so K and the fit cannot describe different matrices. The fit itself is the
# vendored run_lfmm_nogif.R, unmodified; this script only joins the variant
# identity back onto its output and records lambda (never applied).
#
# Usage: Rscript 33_lfmm_kmate.R <trays_csv> <seedmix_csv> <climate_rda> <k_csv> \
#                               <generation> <out_csv> <workdir> <lfmm_script> [K_override]

suppressPackageStartupMessages(library(data.table))
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(here, "lfmm_input.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8L) stop("usage: <trays> <seedmix> <climate_rda> <k_csv> <generation> <out_csv> <workdir> <lfmm_script> [K]")
trays_csv   <- args[[1]]
seedmix_csv <- args[[2]]
climate_rda <- args[[3]]
k_csv       <- args[[4]]
gen         <- as.integer(args[[5]])
out_csv     <- args[[6]]
workdir     <- args[[7]]
lfmm_script <- args[[8]]
k_override  <- if (length(args) >= 9L) as.integer(args[[9]]) else NA_integer_

ts <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))
t0 <- proc.time()

stem_info <- build_lfmm_stem(trays_csv, seedmix_csv, climate_rda, gen, workdir)

# K resolution order: an explicit numeric argument, then a hand-chosen table
# written beside the broken-stick one (35_lfmm_choose_k.R), then broken stick.
# Whichever supplies K, if the table also records the record count it must match
# the matrix just built -- otherwise K describes a different matrix.
K <- k_override
if (is.na(K)) {
  chosen_csv <- file.path(dirname(k_csv), "lfmm_k_chosen.csv")
  src <- if (file.exists(chosen_csv)) chosen_csv else k_csv
  col <- if (src == chosen_csv) "K_chosen" else "K_broken_stick"
  kt <- fread(src)
  hit <- unique(kt[generation == gen, get(col)])
  if (length(hit) != 1L || is.na(hit)) {
    stop("no single ", col, " for generation ", gen, " in ", src)
  }
  K <- as.integer(hit)
  if ("n_records" %in% names(kt)) {
    n_k <- unique(kt[generation == gen, n_records])
    if (length(n_k) == 1L && !is.na(n_k) && n_k != stem_info$nr) {
      stop("the K table was built on ", n_k, " records but this matrix holds ",
           stem_info$nr, " -- rebuild the K table before fitting")
    }
  }
}
if (K < 1L) stop("K resolved to ", K, "; refusing to fit")

lam_file <- file.path(workdir, sprintf("lfmm_nonsnp_gen%d_gif.txt", gen))
raw_csv  <- file.path(workdir, sprintf("lfmm_nonsnp_gen%d_raw.csv", gen))
st <- system2(file.path(R.home("bin"), "Rscript"),
              c("--no-save", lfmm_script, stem_info$stem, K, raw_csv, lam_file),
              stdout = TRUE, stderr = TRUE)
cat(paste(tail(st, 15), collapse = "\n"), "\n")
if (!file.exists(raw_csv)) stop("LFMM produced no output for generation ", gen)

res <- fread(raw_csv)
if (nrow(res) != stem_info$nr) {
  stop("LFMM returned ", nrow(res), " rows for ", stem_info$nr, " records")
}
res <- cbind(stem_info$rec, res)
res[, `:=`(generation = gen, K = K, n_pools = stem_info$n_pools,
           n_imputed = stem_info$n_imputed,
           lambda = if (file.exists(lam_file)) as.numeric(readLines(lam_file)) else NA_real_)]
fwrite(res, out_csv)
ts("wrote ", out_csv, " (", nrow(res), " records, K = ", K, ", lambda = ",
   round(unique(res$lambda)[1], 4), ")  elapsed ",
   round((proc.time() - t0)["elapsed"] / 60, 1), " min")
