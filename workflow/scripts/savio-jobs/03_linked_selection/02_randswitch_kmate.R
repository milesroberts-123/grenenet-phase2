#!/usr/bin/env Rscript

# Generated from resources/03_linked_selection.Rmd chunk `randswitch_kmate`.
#
# ONE-TIME: randomly switch which allele is tracked.
#
# The kMate matrices track the ALT allele throughout. For temporal covariance
# statistics the direction of the trajectory is arbitrary, so a random subset of
# records is flipped to track the ALT complement instead. The same draw is
# applied to the seedmix and to every pool column, so a record's whole
# trajectory is flipped together.
#
# Reads : ..._seedmix.csv, ..._trays.csv            (~60 GB)
# Writes: ..._seedmix_randswitch.csv
#         ..._trays_randswitch.csv                  (~60 GB)
#
# Chunked so memory stays bounded: 50k rows x 893 columns is ~360 MB as doubles.
# Output is written with fwrite; the previous version used write.table once per
# 1000-row chunk, which dominated the runtime.

suppressMessages({
  library(data.table)
  library(readr)
})

ts_message <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))
t0 <- proc.time()

af_prefix <- "/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/data/PHASE_1_2_COMBINED_kmate_arch3v2_allele_frequency"

trays_in  <- paste0(af_prefix, "_trays.csv")
trays_out <- paste0(af_prefix, "_trays_randswitch.csv")

ts_message("reading seedmix")
gen0 <- fread(paste0(af_prefix, "_seedmix.csv"))

# one draw per record, reused for every column of that record
set.seed(20261115)
idx <- sample(c(TRUE, FALSE), size = nrow(gen0), replace = TRUE)
ts_message("flipping ", sum(idx), " of ", length(idx), " records (",
           round(100 * mean(idx), 1), "%)")
gen0[idx, X0_0_0 := 1 - X0_0_0]
fwrite(gen0, paste0(af_prefix, "_seedmix_randswitch.csv"), quote = FALSE)

if (file.exists(trays_out)) file.remove(trays_out)

n_done <- 0L
cb <- function(x, pos) {
  rows <- pos:(pos + nrow(x) - 1L)
  flip <- which(idx[rows])
  m <- as.matrix(x)
  if (length(flip)) m[flip, ] <- 1 - m[flip, ]
  fwrite(as.data.table(m), trays_out, append = pos > 1, col.names = pos == 1,
         quote = FALSE)
  n_done <<- n_done + nrow(m)
  if (n_done %% 500000L < nrow(m)) {
    ts_message(sprintf("  %s records (%.0f%%)", format(n_done, big.mark = ","),
                       100 * n_done / nrow(gen0)))
  }
  rm(m); gc(FALSE)
  NULL
}

ts_message("streaming trays -> randswitch")
read_csv_chunked(
  trays_in,
  callback = SideEffectChunkCallback$new(cb),
  chunk_size = 50000L,
  col_types = cols(.default = col_double())
)

ts_message("done: ", n_done, " records; elapsed ",
           round((proc.time() - t0)[["elapsed"]] / 60, 1), " min")
