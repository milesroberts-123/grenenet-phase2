#!/usr/bin/env Rscript

# Plot binned expression values along the genome, faceted by cell type and
# chromosome.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scico)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: plot_expression_bins.R <bins> <bed> <out.png> <width>")
}

bins_path <- args[1]
bed_path <- args[2]
out_path <- args[3]
width <- as.numeric(args[4])

cell_types <- readLines(bed_path, n = 1) %>%
  strsplit("\t") %>%
  unlist() %>%
  tail(-4)

exp_bins <- read.table(
  bins_path,
  col.names = c("chrom", "start", "end", cell_types)
)

exp_bins <- exp_bins %>%
  filter(!(chrom %in% c("M", "C")))

bins_long <- exp_bins %>%
  pivot_longer(
    cols = -c(chrom, start, end),
    names_to = "cell_type",
    values_to = "count"
  ) %>%
  mutate(bin_mid = (start + end) / 2,
         count = as.numeric(count))

p <- ggplot(bins_long, aes(x = bin_mid, y = count, color = cell_type)) +
  geom_col() +
  facet_grid(cell_type ~ chrom, scales = "free") +
  scale_color_scico_d(palette = "hawaii") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "None") +
  labs(x = "Window midpoint (bp)", y = "Expression (CPM)")

ggsave(out_path, p, width = width, height = width * 0.75, dpi = 150)
