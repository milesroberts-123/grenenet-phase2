#!/usr/bin/env Rscript

# Generated from resources/03_linked_selection.Rmd chunk `biallelic_kmate`.
#
# ONE-TIME: keep only bi-allelic sites.
#
# The catalog is a set of bi-allelic *records*, but a position can carry several
# of them, so a "site" is only unambiguous when exactly one record sits there.
# 8.19% of positions carry more than one record. This keeps positions with
# exactly one record and drops the rest, together with every row of the pool
# matrix at those positions.
#
# The keep vector is derived from the seedmix file, which is written in the same
# row order as the pool matrix, rather than by matching on chrom/pos inside the
# stream. Rows are selected by absolute row counter from the chunked reader, so
# no Arrow ordering assumption is involved.
#
# Reads : ..._seedmix_randswitch_fixrm.csv, ..._trays_randswitch_fixrm.csv
# Writes: ..._seedmix_randswitch_fixrm_biallelic.csv
#         ..._trays_randswitch_fixrm_biallelic.csv

suppressMessages({
  library(data.table)
  library(readr)
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
cat("positions with more than one record:", sum(gen0$n_at_pos > 1), "\n")

fwrite(gen0[keep == TRUE, .(chrom, pos, class, X0_0_0)],
       paste0(af_prefix, "_seedmix_randswitch_fixrm_biallelic.csv"), quote = FALSE)

## ---- stream the pool matrix --------------------------------------------
if (file.exists(trays_out)) file.remove(trays_out)

n_in <- 0L; n_out <- 0L
cb <- function(x, pos) {
  rows <- pos:(pos + nrow(x) - 1L)
  m <- as.matrix(x)[keep[rows], , drop = FALSE]
  fwrite(as.data.table(m), trays_out, append = pos > 1, col.names = pos == 1,
         quote = FALSE)
  n_in  <<- n_in + nrow(x)
  n_out <<- n_out + nrow(m)
  if (n_in %% 500000L < nrow(x)) {
    ts_message(sprintf("  %s in -> %s kept (%.0f%%)", format(n_in, big.mark = ","),
                       format(n_out, big.mark = ","), 100 * n_out / n_in))
  }
  rm(m); gc(FALSE)
  NULL
}

ts_message("streaming trays -> biallelic")
read_csv_chunked(
  trays_in,
  callback = SideEffectChunkCallback$new(cb),
  chunk_size = 50000L,
  col_types = cols(.default = col_double())
)

stopifnot(n_out == sum(keep))
ts_message("done: ", format(n_out, big.mark = ","), " records; elapsed ",
           round((proc.time() - t0)[["elapsed"]] / 60, 1), " min")
