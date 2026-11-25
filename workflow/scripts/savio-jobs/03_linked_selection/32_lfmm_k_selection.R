#!/usr/bin/env Rscript

# LFMM stage 1: choose K, the number of latent factors, for the non-SNP
# allele-frequency-change matrix of each generation.
#
# Method is the one the reference pipeline uses (kMate
# analysis/grenenet_selection/r2_gea_nonsnp/build_lfmm_scree.py), which is also
# the lfmm tutorial's advice: PCA of Y, then take the number of components that
# stand above the broken-stick null.
#
#   Yc = Y - (mean over pools, per record)      centred per record
#   C  = Yc' Yc / n_records                     [pools x pools]
#   K  = first component whose eigenvalue falls below the broken-stick
#
# Y here is dp = AF(pool) - p0, the change from the founding seed mix, which is
# the reference's response variable. Records are non-SNP only (sv + smallindel),
# because that is the matrix LFMM is fitted to; K must describe the matrix it
# actually sees. MAF is defined on the contemporary mean AF across pools, the
# reference convention, not on p0.
#
# No dependency on the lfmm package: this is an eigendecomposition.
#
# Usage: Rscript 32_lfmm_k_selection.R <trays_csv> <seedmix_csv> <out_csv> [workdir]

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) stop("usage: <trays_csv> <seedmix_csv> <out_csv> [workdir]")
trays_csv   <- args[[1]]
seedmix_csv <- args[[2]]
out_csv     <- args[[3]]
workdir     <- if (length(args) >= 4L) args[[4]] else dirname(out_csv)
dir.create(workdir, recursive = TRUE, showWarnings = FALSE)

maf_min <- 0.05
ts <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))
t0 <- proc.time()

# --------------------------------------------------------------------------
# 1. p0 and the record set
# --------------------------------------------------------------------------
ts("reading seedmix")
sm <- fread(seedmix_csv, select = c("class", "X0_0_0"))
p0_all <- sm$X0_0_0
keep_rec <- which(sm$class != "snp")
ts("records: ", nrow(sm), "  non-SNP: ", length(keep_rec),
   "  (smallindel ", sum(sm$class == "smallindel"), ", sv ", sum(sm$class == "sv"), ")")
if (!length(keep_rec)) stop("no non-SNP records in ", seedmix_csv)
rm(sm); invisible(gc())

# --------------------------------------------------------------------------
# 2. one pass over the trays table, keeping only non-SNP rows
#
# The trays table has no chrom/pos/class columns: it is row-aligned with the
# seedmix table, which is the same build. The row set is selected by absolute
# line number (data line i is file line i+1), exactly as 27_vw_panel_extract.R
# does, so no parsing of the 42 GB table is needed to find the rows.
# --------------------------------------------------------------------------
tmp_idx <- tempfile(fileext = ".txt")
fwrite(data.table(keep_rec + 1L), tmp_idx, col.names = FALSE)
sub_csv <- file.path(workdir, "lfmm_nonsnp_rows.csv")
ts("extracting non-SNP rows from the trays table (single awk pass)")
# the trays header is line 1 and is not in the index, so it is reprinted
# explicitly; without it fread would take the first data row as the header
awk <- sprintf("awk -F, 'NR==FNR{k[$1]=1; next} FNR==1{print; next} (FNR in k){print}' %s %s > %s",
               shQuote(tmp_idx), shQuote(trays_csv), shQuote(sub_csv))
st <- system(awk)
if (st != 0L) stop("awk extraction failed with status ", st)
n_sub <- as.integer(trimws(system(paste("wc -l <", shQuote(sub_csv)), intern = TRUE))) - 1L
ts("extracted rows: ", n_sub)
if (n_sub != length(keep_rec)) {
  stop("extracted ", n_sub, " rows but expected ", length(keep_rec),
       " -- the trays and seedmix tables are not row-aligned")
}

ts("reading the non-SNP matrix")
D <- fread(sub_csv, showProgress = FALSE)
if (nrow(D) != n_sub) stop("fread read ", nrow(D), " of ", n_sub, " rows")
p0 <- p0_all[keep_rec]
rm(p0_all); invisible(gc())

pool_cols <- names(D)
gens <- as.integer(sub("^[0-9]+_([0-9]+)_.*$", "\\1", pool_cols))
if (anyNA(gens)) stop("could not parse generations from the trays header")
ts("matrix: ", nrow(D), " records x ", ncol(D), " pools; generations ",
   paste(sort(unique(gens)), collapse = ","))

# --------------------------------------------------------------------------
# 3. per generation: covariance of the centred dp matrix, then broken stick
# --------------------------------------------------------------------------
broken_stick <- function(p) {
  vapply(seq_len(p), function(k) sum(1 / seq(k, p)) / p, numeric(1))
}

res <- list()
for (g in sort(unique(gens))) {
  cols <- which(gens == g)
  if (length(cols) < 3L) {
    ts("generation ", g, ": only ", length(cols), " pools -- skipping")
    next
  }
  Y <- as.matrix(D[, ..cols])
  storage.mode(Y) <- "double"

  # contemporary MAF across this generation's pools, as in the reference
  maf <- rowMeans(Y, na.rm = TRUE)
  maf <- pmin(maf, 1 - maf)
  # a record must also be complete across this generation's pools: eigen() of a
  # covariance with NA fails, and the reference's imputation was inert because
  # their matrix had no missingness. Counted separately so the loss is visible.
  complete <- rowSums(is.na(Y)) == 0L
  keep <- is.finite(maf) & maf >= maf_min & complete
  ts("generation ", g, ": ", length(cols), " pools, ", sum(keep), " of ",
     nrow(Y), " records kept (", sum(!complete), " incomplete, ",
     sum(!complete & is.finite(maf) & maf >= maf_min), " of those would pass MAF)")
  if (sum(keep) < 100L) next
  Y <- Y[keep, , drop = FALSE]
  p0g <- p0[keep]

  # dp against the founding seed mix, then centre each record across pools
  dp <- Y - p0g
  dp <- dp - rowMeans(dp, na.rm = TRUE)

  nr <- nrow(dp)
  C  <- crossprod(dp) / nr
  ev <- eigen(C, symmetric = TRUE, only.values = TRUE)$values
  ev <- pmax(ev, 0)
  ev <- sort(ev, decreasing = TRUE)
  prop <- ev / sum(ev)
  bs   <- broken_stick(length(ev))
  cum  <- cumsum(prop)
  above <- prop > bs
  K_bs <- if (any(!above)) which(!above)[1] - 1L else length(ev)

  ts("generation ", g, ": K(broken-stick) = ", K_bs,
     "  cum var at that K = ", round(cum[max(K_bs, 1L)], 3))

  res[[length(res) + 1L]] <- data.table(
    generation = g, n_pools = length(cols), n_records = nr,
    component = seq_along(ev), eig_prop = prop, broken_stick = bs,
    cum_var = cum, K_broken_stick = K_bs,
    cum_var_at_K = cum[max(K_bs, 1L)]
  )
  rm(Y, dp, C); invisible(gc())
}

if (!length(res)) stop("no generation produced a covariance")
out <- rbindlist(res)
fwrite(out, out_csv)
print(out[, .(generation, n_pools, n_records, K_broken_stick, cum_var_at_K)])
ts("wrote ", out_csv, "  elapsed ", round((proc.time() - t0)["elapsed"] / 60, 1), " min")
