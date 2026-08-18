#!/usr/bin/env Rscript

# Plot binned expression along the genome, faceted by cell type and chromosome,
# plus a per-bin cell-type specificity (tau) track.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scico)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: plot_expression_bins.R <bins> <bed> <out.png> <out_tau.png> <width>")
}

bins_path <- args[1]
bed_path <- args[2]
out_path <- args[3]
tau_path <- args[4]
width <- as.numeric(args[5])

cell_types <- readLines(bed_path, n = 1) %>%
  strsplit("\t") %>%
  unlist() %>%
  tail(-4)

exp_bins <- read.table(
  bins_path,
  col.names = c("chrom", "start", "end", cell_types)
) %>%
  filter(!(chrom %in% c("M", "C")))

# Cell-type specificity (tau): 0 = ubiquitous, 1 = single cell type.
tau <- function(x) {
  n <- length(x)
  if (n <= 1 || max(x) == 0) return(NA_real_)
  x <- x / max(x)
  sum(1 - x) / (n - 1)
}

exp_bins$tau <- apply(exp_bins[, cell_types, drop = FALSE], 1, tau)

bins_long <- exp_bins %>%
  pivot_longer(
    cols = all_of(cell_types),
    names_to = "cell_type",
    values_to = "count"
  ) %>%
  mutate(bin_mid = (start + end) / 2,
         count = as.numeric(count))

p <- ggplot(bins_long, aes(x = bin_mid, y = count, fill = cell_type)) +
  geom_col() +
  facet_grid(cell_type ~ chrom, scales = "free") +
  scale_fill_scico_d(palette = "hawaii") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "None") +
  labs(x = "Window midpoint (bp)", y = "Expression (CPM)")

ggsave(out_path, p, width = width, height = width * 0.75, dpi = 150)

p_tau <- ggplot(exp_bins, aes(x = (start + end) / 2, y = tau)) +
  geom_col() +
  facet_grid(. ~ chrom, scales = "free") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  labs(x = "Window midpoint (bp)", y = "Cell-type Expression Specificity")

ggsave(tau_path, p_tau, width = width, height = width * 0.75, dpi = 150)
