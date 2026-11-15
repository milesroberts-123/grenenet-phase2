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
# Reading is done with arrow (threads off) rather than readr: a 2026-11-14 run
# with read_csv_chunked reported EOF after ~6.9M of 7.56M rows and silently
# produced a short table. Arrow's RecordBatchReader lets the job assert both the
# rows read and the rows kept, so a short read fails loudly instead of writing a
# truncated file. Column types are pinned to double from the header: arrow's
# inference types an all-NA column as null, which would drop real values if the
# inference sample happened to be empty.
#
# Reads : ..._seedmix_randswitch_fixrm.csv, ..._trays_randswitch_fixrm.csv
# Writes: ..._seedmix_randswitch_fixrm_biallelic.csv
#         ..._trays_randswitch_fixrm_biallelic.csv

suppressMessages({
  library(data.table)
  library(arrow)
})

options(arrow.use_threads = FALSE)

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
hdr <- strsplit(readLines(trays_in, n = 1L), ",", fixed = TRUE)[[1]]
sch <- do.call(schema, setNames(rep(list(float64()), length(hdr)), hdr))

ds <- open_dataset(trays_in, format = "csv",
                   convert_options = csv_convert_options(col_types = sch),
                   read_options = csv_read_options(block_size = 64L * 1024L * 1024L))
if (file.exists(trays_out)) file.remove(trays_out)

n_in <- 0L; n_out <- 0L; n_batch <- 0L
reader <- ds$NewScan()$UseThreads(FALSE)$BatchSize(20000L)$Finish()$ToRecordBatchReader()
repeat {
  rb <- reader$read_next_batch()
  if (is.null(rb)) break
  n_batch <- n_batch + 1L
  rows <- (n_in + 1L):(n_in + rb$num_rows)
  m <- as.matrix(as.data.frame(rb))[keep[rows], , drop = FALSE]
  fwrite(as.data.table(m), trays_out, append = n_batch > 1L, col.names = n_batch == 1L,
         quote = FALSE)
  n_in  <- n_in + rb$num_rows
  n_out <- n_out + nrow(m)
  if (n_in %% 500000L < rb$num_rows) {
    ts_message(sprintf("  %s in -> %s kept (%.0f%%)", format(n_in, big.mark = ","),
                       format(n_out, big.mark = ","), 100 * n_out / n_in))
  }
  rm(m, rb); gc(FALSE)
  arrow::default_memory_pool()$ReleaseUnused()
}

ts_message("  batches: ", n_batch, " | rows read: ", format(n_in, big.mark = ","))
stopifnot(
  n_in  == length(keep),      # every row of the input was seen
  n_out == sum(keep)          # every kept row was written
)

ts_message("done: ", format(n_out, big.mark = ","), " records; elapsed ",
           round((proc.time() - t0)[["elapsed"]] / 60, 1), " min")
