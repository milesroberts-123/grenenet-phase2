#!/usr/bin/env Rscript

# LFMM stage 2: fit the reference's LFMM for one generation of the non-SNP
# allele-frequency-change matrix.
#
# This script only builds the three input files the vendored
# run_lfmm_nogif.R expects, calls it, and joins the variant identity back on:
#
#   <stem>_Y.f64     dp = AF(pool) - p0, pools x records, float64, C-order
#   <stem>_dims.txt  "n_pools n_records"
#   <stem>_env.csv   one column, bio1 for the pool's site in that generation's year
#
# dp, the record set and the MAF rule are the reference's (see METHODS_lfmm.md in
# kMate): MAF on the contemporary mean AF across the generation's pools, records
# complete across those pools, K from 32_lfmm_k_selection.R. Lambda is written
# next to the p-values and never applied.
#
# Usage: Rscript 33_lfmm_kmate.R <trays_csv> <seedmix_csv> <climate_rda> <k_csv> \
#                               <generation> <out_csv> <workdir> <lfmm_script> [K_override]

suppressPackageStartupMessages(library(data.table))

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

maf_min <- 0.05
year_of <- function(g) 2017L + g      # generation 1 -> 2018, matching the notebook
dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
ts <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))
t0 <- proc.time()

# --------------------------------------------------------------------------
# record set and p0
# --------------------------------------------------------------------------
ts("reading seedmix")
sm <- fread(seedmix_csv, select = c("chrom", "pos", "class", "X0_0_0"))
keep_rec <- which(sm$class != "snp")
p0_all <- sm$X0_0_0
ts("non-SNP records: ", length(keep_rec))

# one awk pass for the non-SNP rows; reused when an earlier stage left the file
sub_csv <- file.path(workdir, "lfmm_nonsnp_rows.csv")
if (!file.exists(sub_csv)) {
  tmp_idx <- tempfile(fileext = ".txt")
  fwrite(data.table(keep_rec + 1L), tmp_idx, col.names = FALSE)
  ts("extracting non-SNP rows from the trays table")
  awk <- sprintf("awk -F, 'NR==FNR{k[$1]=1; next} FNR==1{print; next} (FNR in k){print}' %s %s > %s",
                 shQuote(tmp_idx), shQuote(trays_csv), shQuote(sub_csv))
  if (system(awk) != 0L) stop("awk extraction failed")
}
D <- fread(sub_csv, showProgress = FALSE)
if (nrow(D) != length(keep_rec)) {
  stop("subset holds ", nrow(D), " rows, expected ", length(keep_rec),
       " -- trays and seedmix are not row-aligned")
}

# --------------------------------------------------------------------------
# this generation's pools, and their environment
# --------------------------------------------------------------------------
gens <- as.integer(sub("^[0-9]+_([0-9]+)_.*$", "\\1", names(D)))
pool_names <- names(D)
cols <- which(gens == gen)
if (length(cols) < 3L) stop("generation ", gen, " has only ", length(cols), " pools")
sites <- as.integer(sub("_.*", "", pool_names[cols]))

e <- new.env(); load(climate_rda, envir = e)
clim <- get(ls(e)[1], envir = e)
clim <- as.data.table(clim)[year == year_of(gen), .(site = as.integer(site), bio1 = as.numeric(bio1))]
env <- clim[match(sites, site)]
if (anyNA(env$bio1)) {
  drop <- is.na(env$bio1)
  ts("dropping ", sum(drop), " pool(s) with no bio1 for ", year_of(gen))
  cols <- cols[!drop]; sites <- sites[!drop]; env <- env[!drop]
}
if (length(cols) < 3L) stop("too few pools with climate at generation ", gen)
X <- (env$bio1 - mean(env$bio1)) / sqrt(mean((env$bio1 - mean(env$bio1))^2))   # z, ddof = 0
ts("generation ", gen, ": ", length(cols), " pools, year ", year_of(gen),
   ", bio1 ", round(min(env$bio1), 2), " to ", round(max(env$bio1), 2))

# --------------------------------------------------------------------------
# response matrix
# --------------------------------------------------------------------------
Y <- as.matrix(D[, ..cols])
storage.mode(Y) <- "double"
maf <- rowMeans(Y, na.rm = TRUE); maf <- pmin(maf, 1 - maf)
complete <- rowSums(is.na(Y)) == 0L
keep <- is.finite(maf) & maf >= maf_min & complete
ts("records kept: ", sum(keep), " of ", nrow(Y), " (",
   sum(!complete), " incomplete, ", sum(!complete & is.finite(maf) & maf >= maf_min),
   " of those above MAF)")
if (sum(keep) < 100L) stop("too few records at generation ", gen)
Y <- Y[keep, , drop = FALSE]
rec <- sm[keep_rec[keep], .(chrom, pos, class)]
p0 <- p0_all[keep_rec][keep]
dp <- Y - p0
rm(Y, D); invisible(gc())

# pools x records, row-major on disk: t(dp) is records x pools, and its
# column-major layout is dp's row-major layout
Yt <- t(dp)
np <- ncol(Yt); nr <- nrow(Yt)

# --------------------------------------------------------------------------
# K: this generation's broken-stick value, unless overridden
# --------------------------------------------------------------------------
K <- k_override
if (is.na(K)) {
  kt <- fread(k_csv)
  hit <- unique(kt[generation == gen, K_broken_stick])
  if (length(hit) != 1L) stop("no single K for generation ", gen, " in ", k_csv)
  K <- as.integer(hit)
}
if (K < 1L) stop("K resolved to ", K, "; refusing to fit")

# --------------------------------------------------------------------------
# write the reference's inputs and call it
# --------------------------------------------------------------------------
stem <- file.path(workdir, sprintf("lfmm_nonsnp_gen%d", gen))
writeBin(as.numeric(Yt), paste0(stem, "_Y.f64"), size = 8L, endian = "little")
writeLines(paste(np, nr), paste0(stem, "_dims.txt"))
fwrite(data.table(X = X), paste0(stem, "_env.csv"))
ts("wrote ", stem, "_{Y.f64,dims.txt,env.csv}  Y = ", np, " x ", nr, " (",
   round(file.size(paste0(stem, "_Y.f64")) / 2^30, 2), " GB)")

lam_file <- file.path(workdir, sprintf("lfmm_nonsnp_gen%d_gif.txt", gen))
raw_csv  <- file.path(workdir, sprintf("lfmm_nonsnp_gen%d_raw.csv", gen))
st <- system2(file.path(R.home("bin"), "Rscript"),
              c("--no-save", lfmm_script, stem, K, raw_csv, lam_file),
              stdout = TRUE, stderr = TRUE)
cat(paste(tail(st, 15), collapse = "\n"), "\n")
if (!file.exists(raw_csv)) stop("LFMM produced no output for generation ", gen)

res <- fread(raw_csv)
if (nrow(res) != nr) stop("LFMM returned ", nrow(res), " rows for ", nr, " records")
res <- cbind(rec, res)
res[, `:=`(generation = gen, K = K,
           lambda = if (file.exists(lam_file)) as.numeric(readLines(lam_file)) else NA_real_)]
fwrite(res, out_csv)
ts("wrote ", out_csv, " (", nrow(res), " records, lambda = ",
   round(unique(res$lambda)[1], 4), ")  elapsed ",
   round((proc.time() - t0)["elapsed"] / 60, 1), " min")
