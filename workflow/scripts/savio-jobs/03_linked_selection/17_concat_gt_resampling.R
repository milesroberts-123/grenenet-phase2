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

today <- Sys.Date()
results_dir <- "/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/results/"
dir_today <- paste0(results_dir, today)

# the site list the arrays iterate over; kept identical to the job scripts
af_prefix <- "/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/data/PHASE_1_2_COMBINED_kmate_arch3v2_allele_frequency"
trays <- paste0(af_prefix, "_trays_randswitch_fixrm_biallelic.csv")
hdr <- strsplit(readLines(trays, n = 1L), ",", fixed = TRUE)[[1]]
site_list <- sort(as.integer(unique(sub("_.*", "", hdr))))

ts_message("today: ", today, " | sites expected: ", length(site_list))

pairs <- list(
  list(combined = "gt_boots_kmate.tsv",      per_site = "gt_boots_kmate_site",      what = "bootstrap G(t)"),
  list(combined = "var_cov_boots_kmate.tsv", per_site = "var_cov_boots_kmate_site", what = "bootstrap var-cov"),
  list(combined = "gt_perms_kmate.tsv",      per_site = "gt_perms_kmate_site",      what = "permutation G(t)"),
  list(combined = "var_cov_perms_kmate.tsv", per_site = "var_cov_perms_kmate_site", what = "permutation var-cov")
)

for (p in pairs) {
  files <- file.path(dir_today, paste0(p$per_site, site_list, ".tsv"))
  missing <- files[!file.exists(files)]
  if (length(missing)) {
    stop(p$what, ": ", length(missing), " of ", length(files),
         " site files missing (e.g. ", basename(missing[1]), "); not writing ",
         p$combined)
  }
  parts <- lapply(files, fread)
  out <- rbindlist(parts, use.names = TRUE, fill = TRUE)
  target <- file.path(dir_today, p$combined)
  if (file.exists(target)) file.remove(target)
  fwrite(out, target, sep = "\t")
  ts_message(sprintf("%-24s %6d site files -> %10d rows", p$combined, length(files), nrow(out)))
  rm(parts, out); gc()
}

ts_message("done")
