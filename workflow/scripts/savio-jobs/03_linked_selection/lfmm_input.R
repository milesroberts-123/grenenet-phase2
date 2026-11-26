#!/usr/bin/env Rscript

# Shared input construction for the LFMM stages.
#
# Both the K calibration sweep (34) and the production fit (33) must see exactly
# the same matrix, or K describes something other than what is fitted. This file
# holds that single construction path; the stages only differ in what they do
# with the result.
#
# build_lfmm_stem() writes the three files the vendored run_lfmm_nogif.R expects,
# following the reference pipeline's rules (kMate METHODS_lfmm.md):
#
#   Y = dp = AF(pool) - p0, pools x records, float64, written row-major
#   record set: non-SNP (sv + smallindel), MAF >= 0.05 on the contemporary mean
#               AF across this generation's pools
#   missing cells: filled with the record's mean dp across pools
#   X: bio1 for the pool's site in that generation's year, z-scored (ddof = 0)

suppressPackageStartupMessages(library(data.table))

lfmm_year_of <- function(g) 2017L + g      # generation 1 -> 2018, as in the notebook

lfmm_load_seedmix <- function(seedmix_csv) {
  sm <- fread(seedmix_csv, select = c("chrom", "pos", "class", "X0_0_0"))
  list(sm = sm, keep_rec = which(sm$class != "snp"))
}

lfmm_extract_nonsnp <- function(trays_csv, keep_rec, workdir) {
  sub_csv <- file.path(workdir, "lfmm_nonsnp_rows.csv")
  if (!file.exists(sub_csv)) {
    # The array tasks of 33 share this path, so a plain write could have several
    # processes filling the same 7 GB file at once. Each writes its own temp name
    # and renames, which is atomic on the same filesystem; the payload is identical
    # either way, so the last writer winning is harmless.
    tmp_idx <- tempfile(fileext = ".txt")
    fwrite(data.table(keep_rec + 1L), tmp_idx, col.names = FALSE)
    tmp_out <- paste0(sub_csv, ".", Sys.getpid(), ".tmp")
    on.exit(unlink(tmp_out), add = TRUE)
    awk <- sprintf("awk -F, 'NR==FNR{k[$1]=1; next} FNR==1{print; next} (FNR in k){print}' %s %s > %s",
                   shQuote(tmp_idx), shQuote(trays_csv), shQuote(tmp_out))
    if (system(awk) != 0L) stop("awk extraction failed")
    if (!file.rename(tmp_out, sub_csv)) stop("could not publish ", sub_csv)
  }
  D <- fread(sub_csv, showProgress = FALSE)
  if (nrow(D) != length(keep_rec)) {
    stop("subset holds ", nrow(D), " rows, expected ", length(keep_rec),
         " -- trays and seedmix are not row-aligned")
  }
  D
}

build_lfmm_stem <- function(trays_csv, seedmix_csv, climate_rda, gen, workdir,
                            maf_min = 0.05) {
  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
  ts <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

  smx <- lfmm_load_seedmix(seedmix_csv)
  D <- lfmm_extract_nonsnp(trays_csv, smx$keep_rec, workdir)

  pool_names <- names(D)
  gens <- as.integer(sub("^[0-9]+_([0-9]+)_.*$", "\\1", pool_names))
  cols <- which(gens == gen)
  if (length(cols) < 3L) stop("generation ", gen, " has only ", length(cols), " pools")
  sites <- as.integer(sub("_.*", "", pool_names[cols]))

  e <- new.env(); load(climate_rda, envir = e)
  clim <- as.data.table(get(ls(e)[1], envir = e))
  clim <- clim[year == lfmm_year_of(gen), .(site = as.integer(site), bio1 = as.numeric(bio1))]
  env <- clim[match(sites, site)]
  if (anyNA(env$bio1)) {
    drop <- is.na(env$bio1)
    ts("dropping ", sum(drop), " pool(s) with no bio1 for ", lfmm_year_of(gen))
    cols <- cols[!drop]; sites <- sites[!drop]; env <- env[!drop]
  }
  if (length(cols) < 3L) stop("too few pools with climate at generation ", gen)
  b <- env$bio1
  X <- (b - mean(b)) / sqrt(mean((b - mean(b))^2))

  Y <- as.matrix(D[, ..cols]); storage.mode(Y) <- "double"
  maf <- rowMeans(Y, na.rm = TRUE); maf <- pmin(maf, 1 - maf)
  keep <- is.finite(maf) & maf >= maf_min
  if (sum(keep) < 100L) stop("too few records at generation ", gen, " (", sum(keep), ")")
  Y <- Y[keep, , drop = FALSE]
  rec <- smx$sm[smx$keep_rec[keep], .(chrom, pos, class)]
  p0 <- smx$sm$X0_0_0[smx$keep_rec[keep]]

  dp <- Y - p0
  n_imp <- sum(!is.finite(dp))
  if (n_imp > 0L) {
    cmean <- rowMeans(dp, na.rm = TRUE)
    if (anyNA(cmean)) cmean[is.na(cmean)] <- 0
    idx <- which(!is.finite(dp), arr.ind = TRUE)
    dp[idx] <- cmean[idx[, 1L]]
  }
  rm(Y, D); invisible(gc())
  ts("generation ", gen, ": ", length(cols), " pools, ", nrow(dp), " records, ",
     n_imp, " imputed cells (", round(100 * n_imp / length(dp), 2), "%)")

  stem <- file.path(workdir, sprintf("lfmm_nonsnp_gen%d", gen))
  Yt <- t(dp)                                  # pools x records
  # dims are "n_pools n_records", and Yt is pools x records, so pools is the ROW
  # count; the column-major write of Yt is then pools x records row-major
  np <- nrow(Yt); nr <- ncol(Yt)
  writeBin(as.numeric(Yt), paste0(stem, "_Y.f64"), size = 8L, endian = "little")
  writeLines(paste(np, nr), paste0(stem, "_dims.txt"))
  fwrite(data.table(X = X), paste0(stem, "_env.csv"))
  ts("wrote ", stem, "_{Y.f64,dims.txt,env.csv}  Y = ", np, " pools x ", nr, " records (",
     round(file.size(paste0(stem, "_Y.f64")) / 2^30, 2), " GB)")

  list(stem = stem, np = np, nr = nr, rec = rec, n_imputed = n_imp,
       n_pools = length(cols), X = X, years = lfmm_year_of(gen))
}
