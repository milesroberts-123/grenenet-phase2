#!/usr/bin/env Rscript

# Manhattan plot of per-variant selection coefficients from one site's
# wright_s_estim table.
#
# 23_wright_s_estim_kmate writes one csv per site with one row per variant x
# segment. This plots each site as its own task: a variant can carry several
# segments, so the plot uses its longest one (ties broken by the largest
# denominator, i.e. the most informative), and the y axis is the signed
# coefficient on a linear scale via rventures::plot_manhattan(plot_effects = TRUE).
#
# The raw estimates have a long right tail: even with the starting-MAF filter a
# two-point segment can return |s| in the tens, so points beyond the site's 99th
# percentile of |s| are dropped and the cut is stated in the title.
#
# Usage: Rscript 24_wright_s_manhattan_kmate.R <s_date> <task_index>
#   s_date:     results folder holding wright_s_estim_site<site>.csv
#   task_index: 0-based, maps onto the sorted site list

suppressMessages({
  library(data.table)
  library(ggplot2)
  library(rventures)
})

ts_message <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("expected <s_date|latest> <task_index>")
s_date <- args[[1]]
task <- suppressWarnings(as.integer(args[[2]]))

results_root <- "/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/results/"

# "latest" resolves to the newest date folder that actually holds s tables. A job
# array stamps results with the date each task starts, and a long estimation run
# can start before midnight and finish after it, so a hard-coded or defaulted date
# silently points at an empty folder.
if (identical(s_date, "latest")) {
  candidates <- list.dirs(results_root, recursive = FALSE, full.names = TRUE)
  candidates <- sort(candidates[grepl("[0-9]{4}-[0-9]{2}-[0-9]{2}$", candidates)],
                     decreasing = TRUE)
  has_files <- vapply(candidates, function(d) {
    length(list.files(d, pattern = "^wright_s_estim_site[0-9]+\\.csv$")) > 0L
  }, logical(1))
  if (!any(has_files)) stop("no results folder under ", results_root, " holds s tables")
  s_dir <- candidates[which(has_files)[1]]
  s_dir <- paste0(s_dir, "/")
  ts_message("date 'latest' -> ", basename(s_dir))
} else {
  s_dir <- paste0(results_root, s_date, "/")
}

files <- list.files(s_dir, pattern = "^wright_s_estim_site[0-9]+\\.csv$",
                    full.names = TRUE)
if (!length(files)) stop("no per-site s tables in ", s_dir)

sites <- sort(as.integer(sub("^.*site([0-9]+)\\.csv$", "\\1", basename(files))))
files <- files[match(sites, as.integer(sub("^.*site([0-9]+)\\.csv$", "\\1", basename(files))))]

if (is.na(task) || task < 0L || task >= length(sites)) {
  stop("task index ", args[[2]], " outside 0..", length(sites) - 1L)
}
this_site <- sites[task + 1L]

ts_message("task ", task, " of ", length(sites) - 1L, " -> site ", this_site)

out_dir <- paste0(results_root, Sys.Date(), "/")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- read and reduce to one row per variant -----------------------------
s_site <- fread(files[task + 1L])
if (!nrow(s_site)) stop("site ", this_site, ": table is empty")

# longest segment; ties broken by the largest denominator (most information)
setorder(s_site, chrom, pos, -dt, -sum_den)
s_site <- s_site[, .SD[1L], by = .(chrom, pos)]
ts_message("variants: ", format(nrow(s_site), big.mark = ","))

finite_s <- s_site[is.finite(s)]
if (!nrow(finite_s)) stop("site ", this_site, ": no finite s")

clip <- as.numeric(stats::quantile(abs(finite_s$s), 0.99, names = FALSE))
s_plot <- finite_s[abs(s) <= clip]
setnames(s_plot, c("chrom", "pos"), c("CHROM", "POS"))
s_plot <- as.data.frame(s_plot)

ts_message("plotted ", format(nrow(s_plot), big.mark = ","),
           " variants | |s| <= ", signif(clip, 3))

## ---- plot ---------------------------------------------------------------
man <- plot_manhattan(
  s_plot,
  plot_effects = TRUE,
  effect_col = "s",
  title = paste0("Site ", this_site, ", |s| <= ", signif(clip, 2), " kept (99%)"),
  ylab = "selection coefficient (s)"
)

plot_name <- paste0(out_dir, "wright_s_manhattan_kmate_site", this_site)
ggsave(paste0(plot_name, ".pdf"), man, height = 4, width = 7)
ggsave(paste0(plot_name, ".png"), man, height = 4, width = 7)
ts_message("wrote ", basename(plot_name), ".pdf/.png")

ts_message("DONE site ", this_site)
