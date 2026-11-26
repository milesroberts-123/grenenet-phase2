#!/usr/bin/env Rscript

# LFMM stage 3: assemble the per-generation fits and report calibration.
#
# Reads lfmm_nonsnp_gen<g>.csv (written by 33), concatenates them, and produces:
#
#   lfmm_nonsnp_all.csv      every record: chrom, pos, class, pval, B, z, generation,
#                            K, n_pools, n_imputed, lambda, q_BH
#   lfmm_lambda.csv          one row per generation: K, lambda, n_records, counts
#   lfmm_p_summary.csv       observed vs expected counts at fixed p thresholds, BH
#                            hits at 0.05, and hits per class
#   lfmm_qq_gen<g>.png       QQ of -log10 p against the exponential null
#   lfmm_manhattan_gen<g>.png
#
# No lfmm dependency. Lambda is reported, never applied (the reference's decision):
# where K over-corrects, deflating by lambda < 1 would inflate statistics.
#
# Usage: Rscript 36_lfmm_summary.R <results_dir> [generations]

suppressPackageStartupMessages({library(data.table); library(ggplot2)})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("usage: <results_dir> [generations]")
res_dir <- args[[1]]
gens <- if (length(args) >= 2L) as.integer(strsplit(args[[2]], ",")[[1]]) else 1:5

ts <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

files <- file.path(res_dir, sprintf("lfmm_nonsnp_gen%d.csv", gens))
files <- files[file.exists(files)]
if (!length(files)) stop("no lfmm_nonsnp_gen*.csv in ", res_dir)
ts("reading ", length(files), " generation file(s)")

all <- rbindlist(lapply(files, fread), fill = TRUE)
setorder(all, generation, chrom, pos)
ts("records: ", nrow(all))

# BH within each generation, since each is a separate scan
all[, q_BH := p.adjust(pval, method = "BH"), by = generation]

lam <- all[, .(n_records = .N,
               n_pools = unique(n_pools)[1],
               K = unique(K)[1],
               lambda = unique(lambda)[1],
               n_imputed = unique(n_imputed)[1],
               n_p_lt_1e5 = sum(pval < 1e-5),
               n_p_lt_1e3 = sum(pval < 1e-3),
               n_q_lt_05 = sum(q_BH < 0.05),
               min_p = min(pval)), by = generation]
lam[, expected_1e5 := n_records * 1e-5]
lam[, expected_1e3 := n_records * 1e-3]
print(lam)

# per class within generation
by_class <- all[, .(n = .N, n_p_lt_1e5 = sum(pval < 1e-5),
                    min_p = min(pval)), by = .(generation, class)]
print(by_class[order(generation, n_p_lt_1e5, decreasing = TRUE)][1:min(.N, 10)])

fwrite(all, file.path(res_dir, "lfmm_nonsnp_all.csv"))
fwrite(lam, file.path(res_dir, "lfmm_lambda.csv"))
fwrite(by_class, file.path(res_dir, "lfmm_p_summary.csv"))

# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------
for (g in unique(all$generation)) {
  d <- all[generation == g]
  n <- nrow(d)

  # QQ against the exponential null, thinned for plotting
  set.seed(1)
  idx <- if (n > 20000L) sort(sample.int(n, 20000L)) else seq_len(n)
  qq <- data.table(observed = -log10(sort(d$pval[idx])),
                   expected = -log10(ppoints(length(idx))))
  p_qq <- ggplot(qq, aes(expected, observed)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
    geom_point(size = 0.5, alpha = 0.5) +
    labs(x = "expected -log10 p", y = "observed -log10 p",
         title = sprintf("LFMM generation %d", g),
         subtitle = sprintf("%d records, K = %d, lambda = %.3f",
                            n, d$K[1], d$lambda[1])) +
    theme_minimal()
  ggsave(file.path(res_dir, sprintf("lfmm_qq_gen%d.png", g)), p_qq,
         width = 5, height = 4.5, dpi = 150)

  chr <- as.integer(gsub("[^0-9]", "", d$chrom))
  d[, `:=`(chr_num = chr)]
  # cumulative position, so chromosomes sit side by side
  offs <- d[, .(maxpos = max(pos)), by = chr_num][order(chr_num)]
  offs[, offset := cumsum(as.numeric(maxpos)) - maxpos]
  d <- merge(d, offs[, .(chr_num, offset)], by = "chr_num")
  d[, gpos := pos + offset]
  axis <- d[, .(mid = mean(range(gpos)), lab = chr_num), by = chr_num][order(chr_num)]

  p_man <- ggplot(d, aes(gpos, -log10(pval))) +
    geom_point(aes(colour = factor(chr_num %% 2 == 0)), size = 0.35, alpha = 0.6) +
    scale_colour_manual(values = c("grey45", "grey20"), guide = "none") +
    scale_x_continuous(breaks = axis$mid, labels = axis$lab) +
    geom_hline(yintercept = -log10(1e-5), linetype = 2, colour = "firebrick") +
    labs(x = "chromosome", y = "-log10 p",
         title = sprintf("LFMM bio1, generation %d", g),
         subtitle = sprintf("raw p, K = %d, lambda = %.3f (not applied)",
                            d$K[1], d$lambda[1])) +
    theme_minimal()
  ggsave(file.path(res_dir, sprintf("lfmm_manhattan_gen%d.png", g)), p_man,
         width = 11, height = 3.6, dpi = 150)
  ts("figures for generation ", g)
}

ts("DONE")
