#!/usr/bin/env Rscript

# Concatenate the per-site resampling outputs into the single tables the
# notebook reads.
#
# 07_gt_curves_bootstrap and 08_gt_sign_permute run as job arrays, one site per
# task, because a single all-sites task needs ~100 h while the slowest site takes
# ~4 h. Each task writes its own file; this step stitches them back together and
# writes the canonical names:
#
#   gt_boots_kmate.tsv        <- gt_boots_kmate_site<site>.tsv
#   var_cov_boots_kmate.tsv   <- var_cov_boots_kmate_site<site>.tsv
#   gt_perms_kmate.tsv        <- gt_perms_kmate_site<site>.tsv
#   var_cov_perms_kmate.tsv   <- var_cov_perms_kmate_site<site>.tsv
#
# Run it after both arrays have finished. It refuses to write a combined table if
# any site file is missing, so a partial array cannot masquerade as a finished
# run.

suppressMessages(library(data.table))

ts_message <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

# Which result folders to read. The arrays stamp every task with Sys.Date() at
# start, so a task that begins after midnight writes into the next day's folder:
# the 2026-11-15 arrays put 23 sites in results/2026-11-15 and the 7 tasks that
# started on the 16th into results/2026-11-16. Usage:
#   Rscript 17_concat_gt_resampling.R <out_date> [<search_date> ...]
# The first date is where the combined tables are written; the remaining dates
# are read as well (default: read only the out date).
results_dir <- "/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/results/"

args <- commandArgs(trailingOnly = TRUE)
out_date <- if (length(args) >= 1L) args[[1]] else as.character(Sys.Date())
search_dates <- if (length(args) >= 2L) args[-1] else out_date
search_dates <- unique(search_dates)

dir_out <- file.path(results_dir, out_date)
dir_in <- file.path(results_dir, search_dates)
stopifnot(dir.exists(dir_out))
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

ts_message("writing to: ", out_date, " | reading: ", paste(search_dates, collapse = ", "))


# the site list the arrays iterate over; kept identical to the job scripts
af_prefix <- "/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/data/PHASE_1_2_COMBINED_kmate_arch3v2_allele_frequency"
trays <- paste0(af_prefix, "_trays_randswitch_fixrm_biallelic.csv")
hdr <- strsplit(readLines(trays, n = 1L), ",", fixed = TRUE)[[1]]
site_list <- sort(as.integer(unique(sub("_.*", "", hdr))))

ts_message("sites expected: ", length(site_list))

pairs <- list(
  list(combined = "gt_boots_kmate.tsv",      per_site = "gt_boots_kmate_site",      what = "bootstrap G(t)"),
  list(combined = "var_cov_boots_kmate.tsv", per_site = "var_cov_boots_kmate_site", what = "bootstrap var-cov"),
  list(combined = "gt_perms_kmate.tsv",      per_site = "gt_perms_kmate_site",      what = "permutation G(t)"),
  list(combined = "var_cov_perms_kmate.tsv", per_site = "var_cov_perms_kmate_site", what = "permutation var-cov")
)

for (p in pairs) {
  # one file per site, found across the searched dates
  files <- unlist(lapply(dir_in, function(d) {
    file.path(d, paste0(p$per_site, site_list, ".tsv"))
  }))
  files <- files[file.exists(files)]

  # a site must not appear twice, or the run is ambiguous about which copy is real
  site_of <- as.integer(sub(paste0(".*", p$per_site), "", sub("\\.tsv$", "", basename(files))))
  dup <- unique(site_of[duplicated(site_of)])
  if (length(dup)) {
    stop(p$what, ": site(s) ", paste(dup, collapse = ", "), " present in more than one ",
         "result folder among [", paste(search_dates, collapse = ", "), "]")
  }

  missing <- setdiff(site_list, site_of)
  if (length(missing)) {
    stop(p$what, ": ", length(missing), " of ", length(site_list),
         " site files missing (", paste(missing, collapse = ", "), "); not writing ",
         p$combined)
  }
  files <- files[match(site_list, site_of)]

  parts <- lapply(files, fread)
  out <- rbindlist(parts, use.names = TRUE, fill = TRUE)
  target <- file.path(dir_out, p$combined)
  if (file.exists(target)) file.remove(target)
  fwrite(out, target, sep = "\t")
  ts_message(sprintf("%-24s %6d site files -> %10d rows", p$combined, length(files), nrow(out)))
  rm(parts, out); gc()
}

ts_message("done")
