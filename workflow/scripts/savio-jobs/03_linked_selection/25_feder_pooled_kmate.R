#!/usr/bin/env Rscript

# Feder test pooled across all sites.
#
# The per-site test (16_feder_kmate_summary.R) groups a site's rows by variant and
# treats that site's replicates as the sample: ybar = mean(yi), s2 = var(yi),
# t = ybar / sqrt(s2/n), df = n - 1. Pooling collapses the grouping key to the
# variant alone, so every (site, rep) observation is a replicate and the same
# arithmetic answers a different question: is the mean allele-frequency change of
# this variant, over all sites and replicates, different from zero?
#
# Consequences worth keeping in view when reading the output:
#   * df rises from ~11 (reps within a site) to ~360, so much smaller |t| clears
#     significance than in the per-site test; hit counts are not comparable.
#   * a variant favoured at some sites and disfavoured at others cancels inside
#     ybar while inflating s2, which biases t toward zero. Heterogeneous selection
#     is therefore not detected by this test.
#
# Input : results/<date>/feder_t_test_kmate.tsv  (chrom,pos,site,rep,pool,yi)
# Output: feder_pooled_kmate_summary.csv
#         feder_pooled_kmate_manhattan.png / .pdf
#
# Usage: Rscript 25_feder_pooled_kmate.R [<results date> | latest]

ts_message <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))
}
t0 <- proc.time()
options(warn = 1)

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)     # the verbs arrow's lazy dataset API dispatches on
  library(data.table)
})

# overridable so the fixture test can run this exact script against a small
# results tree; unset on Savio, where the project root is fixed
root <- Sys.getenv(
  "GRENENET_ROOT",
  unset = "/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2"
)
args <- commandArgs(trailingOnly = TRUE)
want <- if (length(args) >= 1L) args[[1]] else "latest"
# second argument: "no" rebuilds the table without plotting, which lets the
# fixture exercise the statistics where rventures is not installed
do_plot <- !(length(args) >= 2L && args[[2]] %in% c("no", "false"))

# "latest" resolves to the results folder whose feder table was written last,
# not simply the newest folder: a later run may hold only plots or summaries
if (identical(want, "latest")) {
  dirs <- list.dirs(file.path(root, "results"), recursive = FALSE)
  tab <- file.path(dirs, "feder_t_test_kmate.tsv")
  keep <- file.exists(tab)
  if (!any(keep)) stop("no feder_t_test_kmate.tsv under ", file.path(root, "results"))
  src_dir <- dirs[keep][which.max(file.info(tab[keep])$mtime)]
} else {
  src_dir <- file.path(root, "results", want)
}

src <- file.path(src_dir, "feder_t_test_kmate.tsv")
if (!file.exists(src)) stop("missing: ", src)
ts_message("source: ", src, " (", round(file.size(src) / 2^30, 1), " GiB)")

# explicit schema: pos must stay integer (arrow's inferred int64 is fine as well,
# but pinning avoids a silent promotion of yi or a chr/pos mix-up), and the
# unused site/rep/pool columns have to be declared so they can be dropped
schema <- arrow::schema(
  chrom = string(), pos = int32(), site = int32(),
  rep = int32(), pool = string(), yi = float64()
)

ds <- open_dataset(
  src,
  format = "csv",
  schema = schema,
  # arrow requires column_names to match the schema names, and with a schema it
  # does not consume the header row, so skip_rows = 1 drops it explicitly
  read_options = csv_read_options(
    block_size = 64L * 1024L * 1024L,
    column_names = c("chrom", "pos", "site", "rep", "pool", "yi"),
    skip_rows = 1L
  )
)

# Arrow streams the table and aggregates on the fly, so the 235 GB never lands in
# R: only the ~6M per-variant moment rows are collected. Threading is safe here
# because a mean and a variance do not depend on row order (the ordering rules in
# the project notes concern positional alignment, which this job has none of).
ts_message("aggregating by variant (n, sum, sum of squares) ...")
agg <- ds |>
  select(chrom, pos, yi) |>
  filter(!is.na(yi)) |>
  group_by(chrom, pos) |>
  summarise(n = n(), s = sum(yi), ss = sum(yi * yi)) |>
  collect()
setDT(agg)
ts_message("variants aggregated: ", format(nrow(agg), big.mark = ","))

agg[, `:=`(n = as.numeric(n), s = as.numeric(s), ss = as.numeric(ss))]
agg <- agg[n > 1L]  # s2 needs at least two observations

# moments, then the same t statistic the per-site summary uses
agg[, ybar := s / n]
agg[, s2 := (ss - n * ybar^2) / (n - 1)]
agg[, `:=`(tvalue = ybar / sqrt(s2 / n), df = n - 1L)]
agg[, pvalue := 2 * pt(abs(tvalue), df = df, lower.tail = FALSE)]
agg[, p_adj := p.adjust(pvalue, method = "BH")]

# Variant class for the plot. The variant table has one row per record, and a
# position can hold several records (alleles) of different classes, so the merge
# key is not unique: keep one class per position, taking the most severe, or the
# merge would multiply variants.
vi <- fread(
  file.path(root, "data", "PHASE_1_2_COMBINED_kmate_arch3v2_allele_frequency_variant_info.csv"),
  select = c("chrom", "pos", "class")
)
class_rank <- c(sv = 1L, smallindel = 2L, snp = 3L)
vi[, rank := class_rank[class]]
setorder(vi, chrom, pos, rank)
vi <- unique(vi, by = c("chrom", "pos"))[, .(chrom, pos, class)]

n_before <- nrow(agg)
agg <- merge(agg, vi, by = c("chrom", "pos"), all.x = TRUE)
stopifnot(nrow(agg) == n_before)   # one class per position, so no row inflation
agg[, highlight := ifelse(class == "sv", "yes", "no")]

# the accumulator columns are only intermediates
agg[, c("s", "ss") := NULL]
setnames(agg, c("chrom", "pos"), c("CHROM", "POS"))
setcolorder(agg, c("CHROM", "POS", "n", "ybar", "s2", "tvalue", "df",
                   "pvalue", "p_adj", "class", "highlight"))

out_csv <- file.path(src_dir, "feder_pooled_kmate_summary.csv")
fwrite(agg, out_csv, row.names = FALSE, quote = FALSE)
ts_message("wrote ", out_csv)

if (do_plot) {
  suppressPackageStartupMessages(library(ggplot2))
  box::use(rventures[...])
  ts_message("plotting ...")
  man_test <- plot_manhattan(agg, title = "Feder test, all sites pooled")
  ggsave(file.path(src_dir, "feder_pooled_kmate_manhattan.png"),
         man_test, width = 9, height = 5)
  ggsave(file.path(src_dir, "feder_pooled_kmate_manhattan.pdf"),
         man_test, width = 9, height = 5, device = cairo_pdf)
}

ts_message(
  "DONE  variants=", format(nrow(agg), big.mark = ","),
  "  median n=", median(agg$n),
  "  p_adj<0.05: ", format(sum(agg$p_adj < 0.05, na.rm = TRUE), big.mark = ","),
  "  elapsed ", round((proc.time() - t0)["elapsed"] / 60, 1), " min"
)
