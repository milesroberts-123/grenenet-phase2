#!/usr/bin/env Rscript

# Manhattan plot from one site's kmate GLM coefficient table.
#
# 14_glm_kmate writes one table per site (~600 MB: one row per variant per
# block). Plotting all 30 in a knit takes 25-35 minutes, so each site is a task in
# this array instead.
#
# Only chrom/pos/class/p.value are selected before collect: four of the nine
# columns, which keeps a task's read to a few hundred MB rather than 600.
#
# Usage: Rscript 22_glm_manhattan_kmate.R <glm_date> <task_index>
#   glm_date:   results folder holding glm_quasibinomial_coef_site<site>.txt
#   task_index: 0-based, maps onto the sorted site list

suppressMessages({
  library(arrow)
  library(data.table)
  library(ggplot2)
  library(rventures)
})

ts_message <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))
}

# arrow's threaded CSV scanner reassembles scan tasks out of order; the notebook
# sets this in its first chunk, which job scripts skip
options(arrow.use_threads = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("expected <glm_date> <task_index>")
glm_date <- args[[1]]
task <- suppressWarnings(as.integer(args[[2]]))

results_root <- "/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/results/"
glm_dir <- paste0(results_root, glm_date, "/")

files <- list.files(glm_dir, pattern = "^glm_quasibinomial_coef_site[0-9]+\\.txt$",
                    full.names = TRUE)
if (!length(files)) stop("no per-site GLM tables in ", glm_dir)

sites <- sort(as.integer(sub("^.*site([0-9]+)\\.txt$", "\\1", basename(files))))
files <- files[match(sites, as.integer(sub("^.*site([0-9]+)\\.txt$", "\\1", basename(files))))]

if (is.na(task) || task < 0L || task >= length(sites)) {
  stop("task index ", args[[2]], " outside 0..", length(sites) - 1L)
}
this_site <- sites[task + 1L]

ts_message("task ", task, " of ", length(sites) - 1L, " -> site ", this_site)

out_dir <- paste0(results_root, Sys.Date(), "/")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- read and adjust ----------------------------------------------------
glm_site <- open_dataset(files[task + 1L], format = "csv") |>
  dplyr::select(chrom, pos, class, p.value) |>
  dplyr::collect()
setDT(glm_site)

if (!nrow(glm_site)) stop("site ", this_site, ": table is empty")
ts_message("variants: ", format(nrow(glm_site), big.mark = ","))

glm_site[, pvalue := p.adjust(p.value, method = "BH")]
glm_site[, highlight := fifelse(class == "sv", "yes", "no")]
setnames(glm_site, c("chrom", "pos"), c("CHROM", "POS"))

finite_p <- glm_site[is.finite(pvalue)]
if (!nrow(finite_p)) stop("site ", this_site, ": no finite p-values")

ts_message("BH < 0.05: ", format(finite_p[pvalue < 0.05, .N], big.mark = ","),
           " | min p: ", signif(finite_p[, min(pvalue)], 3))

## ---- plot ---------------------------------------------------------------
man_test <- plot_manhattan(glm_site, title = paste0("Site ", this_site))

# ggsave rather than rventures::save_plot: save_plot also writes symlinks under
# results/by-analysis, which this job does not need and which would fail if that
# folder did not exist yet
plot_name <- paste0(out_dir, "glm_manhattan_kmate_site", this_site)
ggsave(paste0(plot_name, ".pdf"), man_test, height = 4, width = 7)
ggsave(paste0(plot_name, ".png"), man_test, height = 4, width = 7)
ts_message("wrote ", basename(plot_name), ".pdf/.png")

ts_message("DONE site ", this_site)
