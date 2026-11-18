#!/usr/bin/env Rscript

# Split the per-variant Feder table into one CSV per site.
#
# Why this exists: the summary step needs each site's rows aggregated, and the
# direct route — arrow reading the 247 GB table — cannot be done safely here.
# Converting the whole table to parquet with a single write_dataset() grew the
# arrow memory pool to the 64 GB cgroup limit and was OOM-killed after 17 minutes
# (2026-11-17). Arrow 24 in R exposes no way to drain that pool.
#
# So the table is split once, streaming, into data/feder_sites/feder_site<site>.csv
# (~8 GB each). The per-site array then reads its own file with fread — no arrow,
# no pool, bounded memory — and every site is independent, as the array rule asks.
#
# The per-line bucketing is done by awk inside this script: R string operations
# over 1.57 billion lines take hours, awk does the same pass at C speed in
# minutes. This is pure data shuffling; everything analytical stays in R.
#
# Usage: Rscript 18_feder_split.R [date]
#   date: results folder holding feder_t_test_kmate.tsv (default: today)

ts_message <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))
}

args <- commandArgs(trailingOnly = TRUE)
date <- if (length(args) >= 1L) args[[1]] else as.character(Sys.Date())

results_dir <- "/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/results/"
src <- paste0(results_dir, date, "/feder_t_test_kmate.tsv")
out_dir <- "/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/data/feder_sites"

if (!file.exists(src)) stop("source table not found: ", src)
ts_message("source: ", src, " (", round(file.size(src) / 1e9, 1), " GB)")

# rebuild from scratch so a stale site file cannot survive into a new run
if (dir.exists(out_dir)) {
  ts_message("removing previous split")
  unlink(out_dir, recursive = TRUE)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# awk program: keep the header in every output file, bucket each row by its third
# comma-separated field (site), append as we go so memory stays flat
awk_prog <- file.path(tempdir(), "split_feder.awk")
writeLines(c(
  'BEGIN { FS = "," }',
  'NR == 1 { header = $0; next }',
  '{ f = sprintf("%s/feder_site%s.csv", dir, $3)',
  '  if (!($3 in seen)) { print header > f; seen[$3] = 1 }',
  '  print >> f }'
), awk_prog)

ts_message("splitting by site (single streaming pass) ...")
status <- system2(
  "awk",
  c("-F,", "-v", shQuote(paste0("dir=", out_dir)), "-f", shQuote(awk_prog), shQuote(src))
)
if (status != 0) stop("awk split failed with status ", status)

files <- list.files(out_dir, pattern = "^feder_site[0-9]+\\.csv$", full.names = TRUE)
if (!length(files)) stop("no per-site files were written")

sites <- as.integer(sub("^feder_site", "", sub("\\.csv$", "", basename(files))))
ord <- order(sites)
files <- files[ord]
sites <- sites[ord]

# wc -l over all files in one call; the trailing "total" line is dropped by
# indexing on the file paths
wc_out   <- system2("wc", c("-l", shQuote(files)), stdout = TRUE)
wc_lines <- as.integer(sub("^\\s*([0-9]+)\\s+.*$", "\\1", wc_out))
names(wc_lines) <- sub("^\\s*[0-9]+\\s+", "", wc_out)
wc_lines <- wc_lines[files]

manifest <- data.frame(
  site  = sites,
  file  = basename(files),
  bytes = file.size(files),
  lines = unname(wc_lines),
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(out_dir, "manifest.csv"), row.names = FALSE)

ts_message("sites: ", nrow(manifest),
           " | rows: ", format(sum(manifest$lines, na.rm = TRUE), big.mark = ","),
           " | biggest: site ", manifest$site[which.max(manifest$bytes)],
           " (", round(max(manifest$bytes) / 1e9, 1), " GB)")
ts_message("done")
