#!/usr/bin/env Rscript

# Generated from resources/03_linked_selection.Rmd chunk `merge_kmate_arch3v2`.
# Chunk body copied verbatim; only this header is added.
#
# One-time rebuild of the merged kMate arch3v2 allele frequency tables from
# resources/kmate_allele_frequencies/ (2,608 samples -> 893 pools, 8,536,388
# records). float16 .npy is decoded with numpy through reticulate; numpy lives in
# the miniforge base env, which the chunk pins via RETICULATE_PYTHON.
#
# Peak RSS is generation 1's pool accumulator: 325 x 6,269,335 doubles = 16.3 GB,
# plus 5.9 GB for non-SNP, so this needs a compute node (16 CPUs -> 64 GB).
#
# Writes, all in data/:
#   PHASE_1_2_COMBINED_kmate_arch3v2_allele_frequency_seedmix.csv
#   PHASE_1_2_COMBINED_kmate_arch3v2_allele_frequency_trays.csv
#   PHASE_1_2_COMBINED_kmate_arch3v2_allele_frequency_variant_info.csv
# plus a _trays.csv.done sentinel that makes reruns skip the rebuild. Delete the
# sentinel to force a rebuild.
#
# Assumed working directory: .../grenenet-phase2/resources

suppressMessages({library(data.table)})
options(warn = 1)

ts_message <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))
t0 <- proc.time()

kmate_dir <- "/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/resources/kmate_allele_frequencies/"
af_prefix <- "/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/data/PHASE_1_2_COMBINED_kmate_arch3v2_allele_frequency"

af_tmp  <- paste0(dirname(af_prefix), "/kmate_arch3v2_merge_tmp")
af_done <- paste0(af_prefix, "_trays.csv.done")

# float16 .npy cannot be read by any R package (RcppCNPy handles double/integer
# only; CRAN's `float` is 32-bit only), so decode through numpy, which has native
# IEEE 754 binary16 support. numpy is memory-mapped, so only one sample row is
# materialised at a time (12.5 MB float16 -> 25 MB float32 -> 50 MB in R).
Sys.setenv(RETICULATE_PYTHON = "/global/home/users/milesroberts/miniforge3/bin/python3")
np <- reticulate::import("numpy", convert = FALSE)

np_open <- function(path) {
  a <- np$load(path, mmap_mode = "r")
  # fully qualified: the job header attaches only data.table, so an unqualified
  # py_to_r() is out of scope there even though it resolves in the notebook
  shp <- as.integer(reticulate::py_to_r(a$shape))
  list(arr = a, nrow = shp[1], ncol = shp[2],
       dtype = as.character(reticulate::py_to_r(a$dtype)))
}

# i is 1-based; numpy indexing is 0-based
np_row <- function(h, i) {
  stopifnot(h$dtype == "float16", i >= 1L, i <= h$nrow)
  reticulate::py_to_r(np$asarray(h$arr[as.integer(i - 1L), ], dtype = "float32"))
}

if (file.exists(af_done)) {

  message("merged tables already built: ", af_prefix, "_*.csv")

} else {

  t_start <- Sys.time()
  dir.create(af_tmp, recursive = TRUE, showWarnings = FALSE)

  ## ---- record index: natural order per matrix, output order genomic -------
  idx_snp <- fread(paste0(kmate_dir, "years1-3/gen_matrices/snp_column_index.csv"))
  idx_ns  <- fread(paste0(kmate_dir, "years1-3/gen_matrices/nonsnp_column_index.csv"))
  stopifnot(nrow(idx_snp) == 6269335L, nrow(idx_ns) == 2267053L)

  rec <- rbind(
    data.table(src = "snp", j = idx_snp$col + 1L,
               chrom = idx_snp$chrom, pos = idx_snp$pos,
               ref_len = idx_snp$ref_len, alt_len = idx_snp$alt_len,
               class = "snp", n_called = idx_snp$n_called),
    data.table(src = "nonsnp", j = idx_ns$col + 1L,
               chrom = idx_ns$chrom, pos = idx_ns$pos,
               ref_len = idx_ns$ref_len, alt_len = idx_ns$alt_len,
               class = fifelse(abs(idx_ns$alt_len - idx_ns$ref_len) > 50L,
                               "sv", "smallindel"),
               n_called = idx_ns$n_called)
  )
  setorder(rec, chrom, pos, class, src, j)
  print(rec[, .N, by = class])

  ## ---- seedmix: mean alt_freq over the 8 founding replicates --------------
  p0_snp <- numeric(nrow(idx_snp))
  p0_ns  <- numeric(nrow(idx_ns))
  for (s in 1:8) {
    d <- fread(paste0(kmate_dir, "seedmix/SEEDMIX_S", s, ".tsv"),
               select = c("chrom", "pos", "ref_len", "alt_len", "alt_freq"))
    is_snp <- d$ref_len == 1L & d$alt_len == 1L
    if (s == 1L) {
      # positional alignment is the whole contract of this merge: verify it once
      stopifnot(sum(is_snp) == nrow(idx_snp),
                identical(d$chrom[is_snp], idx_snp$chrom),
                identical(d$pos[is_snp], idx_snp$pos),
                identical(d$chrom[!is_snp], idx_ns$chrom),
                identical(d$pos[!is_snp], idx_ns$pos))
    }
    p0_snp <- p0_snp + d$alt_freq[is_snp]
    p0_ns  <- p0_ns  + d$alt_freq[!is_snp]
    rm(d); gc(FALSE)
    message(sprintf("  seedmix replicate %d/8 (%.0f s)", s,
                    as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
  }
  p0_snp <- p0_snp / 8
  p0_ns  <- p0_ns / 8
  message("seedmix averaged over 8 replicates")

  ## ---- pools: flower-weighted aggregation, one generation at a time -------
  block_records <- 200000L
  blk   <- (seq_len(nrow(rec)) - 1L) %/% block_records + 1L
  nblk  <- max(blk)
  pool_names <- vector("list", 5L)
  npools     <- integer(5L)

  for (g in 1:5) {
    era  <- if (g <= 3L) "years1-3" else "years4-5"
    gm   <- paste0(kmate_dir, era, "/gen_matrices/")
    meta <- fread(paste0(gm, "gen", g, ".rowmeta.csv"))
    meta[, pool := paste(site, generation, plot, sep = "_")]
    pools <- unique(meta$pool)
    pidx  <- match(meta$pool, pools)
    tot   <- as.numeric(tapply(meta$flowerscollected, pidx, sum))
    w     <- meta$flowerscollected / tot[pidx]
    pool_names[[g]] <- pools
    npools[g] <- length(pools)
    message(sprintf("gen %d: %d samples -> %d pools", g, nrow(meta), length(pools)))

    acc_snp <- matrix(0, length(pools), nrow(idx_snp))
    acc_ns  <- matrix(0, length(pools), nrow(idx_ns))

    for (cl in c("snp", "nonsnp")) {
      f <- paste0(gm, "gen", g, "_", cl, "_af.npy")
      h <- np_open(f)
      stopifnot(h$nrow == nrow(meta), h$ncol == if (cl == "snp") nrow(idx_snp) else nrow(idx_ns))
      t_cls <- Sys.time()
      for (i in seq_len(h$nrow)) {
        v <- np_row(h, i)
        if (cl == "snp") acc_snp[pidx[i], ] <- acc_snp[pidx[i], ] + w[i] * v
        else             acc_ns[pidx[i], ]  <- acc_ns[pidx[i], ]  + w[i] * v
        # gen 1 has 892 rows, so a silent pass can run several minutes
        if (i %% 50L == 0L || i == h$nrow) {
          el <- as.numeric(difftime(Sys.time(), t_cls, units = "secs"))
          message(sprintf("    %s %d/%d (%.0f%%), %.0f s elapsed, ~%.0f s left",
                          cl, i, h$nrow, 100 * i / h$nrow, el,
                          el / i * (h$nrow - i)))
        }
      }
      rm(v); gc(FALSE)
      message("  ", cl, " aggregated (", format(Sys.time(), "%H:%M:%S"), ")")
    }

    ## write this generation's pool matrix in output order
    for (b in seq_len(nblk)) {
      ids <- which(blk == b)
      m   <- matrix(NA_real_, length(ids), npools[g])
      si  <- ids[rec$src[ids] == "snp"]
      ni  <- ids[rec$src[ids] == "nonsnp"]
      if (length(si)) m[match(si, ids), ] <- t(acc_snp[, rec$j[si], drop = FALSE])
      if (length(ni)) m[match(ni, ids), ] <- t(acc_ns[,  rec$j[ni], drop = FALSE])
      writeBin(as.numeric(m), file.path(af_tmp, sprintf("gen%d_block%03d.bin", g, b)),
               size = 4L)
      rm(m)
    }
    rm(acc_snp, acc_ns); gc(FALSE)
  }

  ## ---- assemble: all generations side by side, one row block at a time ----
  tray_cols <- unlist(pool_names)
  stopifnot(length(tray_cols) == sum(npools), !anyDuplicated(tray_cols))

  for (b in seq_len(nblk)) {
    ids   <- which(blk == b)
    parts <- vector("list", 5L)
    for (g in 1:5) {
      f <- file.path(af_tmp, sprintf("gen%d_block%03d.bin", g, b))
      parts[[g]] <- matrix(readBin(f, "numeric", n = length(ids) * npools[g], size = 4L),
                           nrow = length(ids), ncol = npools[g])
    }
    m <- do.call(cbind, parts); rm(parts)
    colnames(m) <- tray_cols
    m[m == 0 | m == 1] <- NA      # absorbing boundary -> NA
    m <- round(m, 5)              # float16 source carries ~3 significant digits
    fwrite(as.data.table(m), paste0(af_prefix, "_trays.csv"),
           append = b > 1L, col.names = b == 1L, quote = FALSE, na = "NA")
    rm(m); gc(FALSE)
    if (b %% 10L == 0L) message("  wrote block ", b, "/", nblk)
  }

  ## ---- seedmix and variant info, in the same record order -----------------
  p0_rec <- numeric(nrow(rec))
  si     <- rec$src == "snp"
  p0_rec[si]  <- p0_snp[rec$j[si]]
  p0_rec[!si] <- p0_ns[rec$j[!si]]
  p0_rec[p0_rec == 0 | p0_rec == 1] <- NA

  fwrite(data.table(chrom = rec$chrom, pos = rec$pos, class = rec$class,
                    X0_0_0 = p0_rec),
         paste0(af_prefix, "_seedmix.csv"), quote = FALSE, na = "NA")

  fwrite(rec[, .(chrom, pos, class, ref_len, alt_len, n_called)],
         paste0(af_prefix, "_variant_info.csv"), quote = FALSE)

  file.create(af_done)
  unlink(af_tmp, recursive = TRUE)
  message("merge complete in ",
          round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), " min")
}

ts_message("elapsed: ", round((proc.time() - t0)[["elapsed"]] / 60, 1), " min")
