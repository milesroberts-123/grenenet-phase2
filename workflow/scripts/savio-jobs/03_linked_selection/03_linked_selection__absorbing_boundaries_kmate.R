#!/usr/bin/env Rscript

# Generated from resources/03_linked_selection.Rmd chunk `absorbing_boundaries_kmate`.
#
# ONE-TIME: drop allele frequencies near absorbing boundaries.
#
# Two operations, both reading only the randswitch CSVs (no parquet cache):
#
#  1. seedmix: drop records whose founding frequency sits within 1/231 of 0 or 1
#     (or is NA, which is where the merge put exact 0/1).
#  2. pools: for each (site, plot), truncate the trajectory at its first NA and
#     leave everything after it NA. The merge turns exact 0/1 into NA, so this
#     is what stops a fixed or lost record from contributing observations after
#     the fixation event. Columns are ordered by generation within the group
#     before truncation.
#
# Reads : ..._seedmix_randswitch.csv, ..._trays_randswitch.csv   (~60 GB)
# Writes: ..._seedmix_randswitch_fixrm.csv
#         ..._trays_randswitch_fixrm.csv                         (~60 GB)
#
# Column order in the output is the input order. The old version rebuilt the
# table site-by-site, which reordered columns; nothing downstream addresses
# columns positionally, so preserving the input order is simpler and cheaper.

suppressMessages({
  library(data.table)
  library(readr)
})

ts_message <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))
t0 <- proc.time()

af_prefix <- "/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/data/PHASE_1_2_COMBINED_kmate_arch3v2_allele_frequency"

trays_in  <- paste0(af_prefix, "_trays_randswitch.csv")
trays_out <- paste0(af_prefix, "_trays_randswitch_fixrm.csv")

## ---- 1. seedmix ---------------------------------------------------------
ts_message("filtering seedmix")
gen0 <- fread(paste0(af_prefix, "_seedmix_randswitch.csv"))
keep <- !is.na(gen0$X0_0_0) &
  gen0$X0_0_0 >= 1/231 & gen0$X0_0_0 <= 1 - 1/231
ts_message("  keeping ", sum(keep), " of ", nrow(gen0), " records (",
           round(100 * mean(keep), 2), "%)")
fwrite(gen0[keep], paste0(af_prefix, "_seedmix_randswitch_fixrm.csv"), quote = FALSE)

## ---- 2. column groups: one per (site, plot), ordered by generation ------
hdr <- names(fread(trays_in, nrows = 0L))
parts <- tstrsplit(hdr, "_", fixed = TRUE)
site <- as.integer(parts[[1]]); gen <- as.integer(parts[[2]]); plot <- as.integer(parts[[3]])
grp <- paste(site, plot, sep = "_")
groups <- split(seq_along(hdr), grp)
groups <- lapply(groups, function(i) i[order(gen[i])])
groups <- groups[lengths(groups) > 1L]
ts_message("  ", length(groups), " (site, plot) trajectories to truncate")

## ---- 3. stream the pools ------------------------------------------------
if (file.exists(trays_out)) file.remove(trays_out)

n_in <- 0L; n_out <- 0L
cb <- function(x, pos) {
  rows <- pos:(pos + nrow(x) - 1L)
  m <- as.matrix(x)[keep[rows], , drop = FALSE]
  for (g in groups) {
    sub <- m[, g, drop = FALSE]
    if (!anyNA(sub)) next
    # cumulative any-NA along generations, then NA everything from the first hit
    ch <- is.na(sub)
    if (ncol(ch) > 1L) for (j in 2:ncol(ch)) ch[, j] <- ch[, j] | ch[, j - 1L]
    sub[ch] <- NA
    m[, g] <- sub
  }
  fwrite(as.data.table(m), trays_out, append = pos > 1, col.names = pos == 1,
         quote = FALSE)
  n_in  <<- n_in + nrow(x)
  n_out <<- n_out + nrow(m)
  if (n_in %% 500000L < nrow(x)) {
    ts_message(sprintf("  %s in, %s kept (%.0f%%)", format(n_in, big.mark = ","),
                       format(n_out, big.mark = ","), 100 * n_out / n_in))
  }
  rm(m, sub); gc(FALSE)
  NULL
}

ts_message("streaming trays -> fixrm")
read_csv_chunked(
  trays_in,
  callback = SideEffectChunkCallback$new(cb),
  chunk_size = 50000L,
  col_types = cols(.default = col_double())
)

na_frac <- NA_real_
ts_message("done: ", n_out, " records kept; elapsed ",
           round((proc.time() - t0)[["elapsed"]] / 60, 1), " min")
