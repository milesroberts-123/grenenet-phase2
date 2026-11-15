#!/usr/bin/env Rscript

# Generated from resources/03_linked_selection.Rmd chunk `biallelic_kmate`.
#
# ONE-TIME: keep only bi-allelic sites.
#
# The catalog is a set of bi-allelic *records*, but a position can carry several
# of them, so a "site" is only unambiguous when exactly one record sits there.
# 17.9% of the surviving records sit at multi-record positions. This keeps
# positions with exactly one record and drops the rest, together with every row
# of the pool matrix at those positions.
#
# The keep vector is derived from the seedmix file, which is written in the same
# row order as the pool matrix, rather than by matching on chrom/pos inside the
# stream. Rows are selected by absolute row counter, so no ordering assumption
# about the reader is involved.
#
# The row selection is a pure line filter: every kept line is copied verbatim,
# so nothing parses 893 numeric columns per row. Two earlier attempts failed
# here. read_csv_chunked (2026-11-14) reported EOF after ~6.9M of 7.56M rows and
# silently wrote a short table. An arrow RecordBatchReader then OOM-killed at
# 16 GB, because decoding and re-encoding every column grows memory without
# bound. Line filtering has neither failure mode, and the assertions below still
# fail loudly if the input and the seedmix disagree in length.
#
# Reads : ..._seedmix_randswitch_fixrm.csv, ..._trays_randswitch_fixrm.csv
# Writes: ..._seedmix_randswitch_fixrm_biallelic.csv
#         ..._trays_randswitch_fixrm_biallelic.csv

suppressMessages({
  library(data.table)
})

ts_message <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))
t0 <- proc.time()

af_prefix <- "/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/data/PHASE_1_2_COMBINED_kmate_arch3v2_allele_frequency"

trays_in  <- paste0(af_prefix, "_trays_randswitch_fixrm.csv")
trays_out <- paste0(af_prefix, "_trays_randswitch_fixrm_biallelic.csv")

## ---- keep vector, derived in the same row order as the pool matrix ------
ts_message("reading seedmix")
gen0 <- fread(paste0(af_prefix, "_seedmix_randswitch_fixrm.csv"))
gen0[, n_at_pos := .N, by = .(chrom, pos)]
keep <- gen0$n_at_pos == 1L

ts_message("  ", format(nrow(gen0), big.mark = ","), " records at ",
           format(uniqueN(paste(gen0$chrom, gen0$pos)), big.mark = ","), " positions")
ts_message("  keeping ", format(sum(keep), big.mark = ","), " (",
           round(100 * mean(keep), 2), "%), dropping ",
           format(sum(!keep), big.mark = ","), " at multi-record positions")
print(gen0[keep == TRUE, .N, by = class])
cat("dropped by class:\n"); print(gen0[keep == FALSE, .N, by = class])

fwrite(gen0[keep == TRUE, .(chrom, pos, class, X0_0_0)],
       paste0(af_prefix, "_seedmix_randswitch_fixrm_biallelic.csv"), quote = FALSE)

## ---- stream the pool matrix --------------------------------------------
# A pure line filter: every output row is an input row, byte for byte, so there
# is no reason to parse 893 numeric columns per row. Reading lines and writing
# the kept ones keeps memory flat (one chunk of text) and does not touch the
# values at all.
con_in  <- file(trays_in, "r")
con_out <- file(trays_out, "w")
on.exit(try(close(con_in), silent = TRUE), add = TRUE)
on.exit(try(close(con_out), silent = TRUE), add = TRUE)

writeLines(readLines(con_in, n = 1L), con_out)   # header

chunk_lines <- 50000L
n_in <- 0L; n_out <- 0L; n_chunk <- 0L
repeat {
  lines <- readLines(con_in, n = chunk_lines, warn = FALSE)
  if (!length(lines)) break
  n_chunk <- n_chunk + 1L
  rows <- (n_in + 1L):(n_in + length(lines))
  if (max(rows) > length(keep))
    stop("input has more rows than the seedmix: ", max(rows), " > ", length(keep))
  keep_rows <- lines[keep[rows]]
  writeLines(keep_rows, con_out)
  n_in  <- n_in + length(lines)
  n_out <- n_out + length(keep_rows)
  if (n_in %% 500000L < length(lines)) {
    ts_message(sprintf("  %s in -> %s kept (%.0f%%)", format(n_in, big.mark = ","),
                       format(n_out, big.mark = ","), 100 * n_out / n_in))
  }
  rm(lines, keep_rows)
}
close(con_in); close(con_out)

ts_message("  chunks: ", n_chunk, " | rows read: ", format(n_in, big.mark = ","))
stopifnot(
  n_in  == length(keep),      # every row of the input was seen
  n_out == sum(keep)          # every kept row was written
)

ts_message("done: ", format(n_out, big.mark = ","), " records; elapsed ",
           round((proc.time() - t0)[["elapsed"]] / 60, 1), " min")
