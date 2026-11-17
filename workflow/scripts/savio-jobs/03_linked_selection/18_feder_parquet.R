#!/usr/bin/env Rscript

# Rebuild data/feder_kmate_parquet from a per-variant Feder test table.
#
# 16_feder_kmate_summary aggregates each site's raw rows inside arrow, which is
# only fast if the source is parquet. The parquet directory is a cache, and
# 15_feder_kmate skips rebuilding it when the directory exists — which is how the
# 2026-11-16 summary ended up reading the 2026-11-06 cache (built from the old
# arch3 table) instead of the 247 GB table written that night. This job always
# rebuilds, so the cache cannot silently outlive the table it came from.
#
# Usage: Rscript 18_feder_parquet.R [date]
#   date: results folder holding feder_t_test_kmate.tsv (default: today)

suppressMessages({
  library(arrow)
  library(data.table)
})

ts_message <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))
}

# arrow's threaded scanner reassembles scan tasks out of order; the notebook sets
# this in its first chunk, which job scripts skip
options(arrow.use_threads = FALSE)

args <- commandArgs(trailingOnly = TRUE)
date <- if (length(args) >= 1L) args[[1]] else as.character(Sys.Date())

results_dir <- "/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/results/"
src  <- paste0(results_dir, date, "/feder_t_test_kmate.tsv")
dest <- "/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/data/feder_kmate_parquet"

if (!file.exists(src)) stop("source table not found: ", src)

ts_message("source: ", src, " (", round(file.size(src) / 1e9, 1), " GB)")

# rebuild from scratch: a stale directory is exactly the failure this job exists
# to prevent
if (dir.exists(dest)) {
  ts_message("removing previous parquet cache")
  unlink(dest, recursive = TRUE)
}
dir.create(dest, recursive = TRUE, showWarnings = FALSE)

ds <- open_dataset(src, format = "csv",
                   read_options = csv_read_options(block_size = 16L * 1024L * 1024L))

ts_message("converting to parquet ...")
write_dataset(ds, dest, format = "parquet", max_rows_per_file = 500000L)

files <- list.files(dest, pattern = "\\.parquet$")
ts_message("wrote ", length(files), " parquet files, ",
           round(sum(file.size(file.path(dest, files))) / 1e9, 1), " GB")

# sanity: the cache must carry the same site set as the table it replaced
chk <- open_dataset(dest)
site_n <- chk |> dplyr::distinct(site) |> dplyr::collect()
ts_message("sites in cache: ", nrow(site_n), " (",
           paste(sort(site_n$site), collapse = ","), ")")

ts_message("done")
