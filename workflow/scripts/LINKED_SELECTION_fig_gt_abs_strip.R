## Summary-panel strip of Miles' Gabs(t) (genome-wide linked selection), one
## sparkline per garden, cold -> warm, with the permutation null in grey.
## Style: bottom panel of the PI's summary figure ("B genome-wide total
## selection": cold locations -- colour bar -- warm locations; thin per-garden
## curves; grey dashed null; tiny "1 ... 5 Years" axis under the first panel).
##
## Data: Miles Roberts, Drive folder 1wDf5Dt5YOzc8zBOsijI9_cVt98Yk5bBM, files
## dated 2026-10-19 ("latest G(t) attempt", Slack 2026-10-20), copied to
## data-intermediate/covariances_miles_2026-10-19/. SV class only (the perms
## exist only for sv; Miles: SV-only is within 1% of all variants).
##   gt_point_estimates_kmate.tsv  per site x rep (tray) x gen, gt_abs = Gabs(t)
##   gt_perms_kmate.tsv            250 permutations of the same
## gen 1 is t = 0 (Gabs = 0 by definition); gens 2-5 = years 1-4 of evolution.
##
## Per garden: thick line = mean Gabs(t) over replicate trays; thin lines = each
## tray; grey ribbon = 2.5-97.5% of the per-permutation garden mean, dotted =
## its median. Sites ordered and coloured by WorldClim bio1.
##
## Usage: Rscript scripts/LINKED_SELECTION_fig_gt_abs_strip.R [metric]
##        metric = gt_abs (default, Miles' Gabs(t)) | gt_absval (|signed G(t)|)
## Output: figs/LINKED_SELECTION_gt_abs_strip[_<metric>].{pdf,png}

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork); library(here) })
setwd(here::here())
args <- commandArgs(TRUE)
metric <- if (length(args) >= 1) args[1] else "gt_abs"
covdir <- "/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/results/2026-10-23"

g <- fread(file.path(covdir, "gt_point_estimates_kmate.tsv"))[class == "sv"]
p <- fread(file.path(covdir, "gt_perms_kmate.tsv"))[class %in% c("sv", "smallindel_sv_snp")]  # perm loop writes the pooled label; point tables use "sv"
for (d in list(g, p)) d[, gt_absval := abs(gt)]
clim <- fread("data/climate_sites_experiment.csv")[, .(site, bio1 = site_bio1)]

g[, year := gen]; p[, year := gen]                    # gen 1 = t0 -> "year 1" tick as in the reference
reps <- g[, .(y = mean(get(metric), na.rm = TRUE)), by = .(site, rep, year)]
obs  <- g[, .(y = mean(get(metric), na.rm = TRUE)), by = .(site, year)]
nul  <- p[, .(y = mean(get(metric), na.rm = TRUE)), by = .(site, year, perm)][
  is.finite(y)][
  , .(lo = quantile(y, .025), med = median(y), hi = quantile(y, .975)), by = .(site, year)]
ord <- clim[site %in% obs$site][order(bio1)]
ord[, ord := .I]
for (d in list(reps, obs, nul)) { d[clim, bio1 := i.bio1, on = "site"]; d[ord, ord := i.ord, on = "site"] }
cat(sprintf("%d gardens, bio1 %.1f to %.1f C; %d trays; max year per garden: %s\n",
            nrow(ord), min(ord$bio1), max(ord$bio1), uniqueN(reps[, .(site, rep)]),
            paste(table(obs[, max(year), by = site]$V1), collapse = "/")))

TLIM <- range(ord$bio1)
ymax <- max(c(obs$y, nul$hi)) * 1.05
## RdBu with the near-white centre replaced by a light grey so mid-temperature gardens stay visible
rdbu <- rev(RColorBrewer::brewer.pal(11, "RdBu")); rdbu[6] <- "#cfd8dc"
pal <- scale_colour_gradientn(colours = rdbu, limits = TLIM, guide = "none")

## tiny axis under the first panel only
ax <- data.table(ord = 1, x0 = 1, x1 = 5, y = -0.03 * ymax)

strip <- ggplot() +
  geom_ribbon(data = nul, aes(year, ymin = lo, ymax = hi, group = site), fill = "grey82", alpha = 0.9) +
  geom_line(data = nul, aes(year, med, group = site), colour = "grey50", linewidth = 0.35, linetype = "22") +
  geom_line(data = reps, aes(year, y, group = interaction(site, rep), colour = bio1), linewidth = 0.3, alpha = 0.35) +
  geom_line(data = obs, aes(year, y, group = site, colour = bio1), linewidth = 1.4, lineend = "round") +
  geom_segment(data = ax, aes(x = x0, xend = x1, y = y, yend = y), colour = "grey30", linewidth = 0.3) +
  geom_segment(data = ax, aes(x = x0 - 0.4, xend = x0 - 0.4, y = 0, yend = ymax), colour = "black", linewidth = 0.3) +
  geom_segment(data = ax, aes(x = x0, xend = x0, y = y, yend = y - 0.025 * ymax), colour = "grey30", linewidth = 0.3) +
  geom_segment(data = ax, aes(x = x1, xend = x1, y = y, yend = y - 0.025 * ymax), colour = "grey30", linewidth = 0.3) +
  geom_text(data = ax, aes(x = x0, y = y - 0.05 * ymax, label = "1"), size = 2.6, vjust = 1, colour = "grey30") +
  geom_text(data = ax, aes(x = x1, y = y - 0.05 * ymax, label = "5"), size = 2.6, vjust = 1, colour = "grey30") +
  geom_text(data = ax, aes(x = (x0 + x1) / 2, y = y - 0.17 * ymax, label = "Years"), size = 2.8, vjust = 1, colour = "grey30") +
  pal +
  facet_wrap(~ ord, nrow = 1) +
  scale_x_continuous(expand = expansion(add = 0.4)) +
  scale_y_continuous(limits = c(-0.25 * ymax, ymax), breaks = seq(0, 0.6, 0.2), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "Linked selection\n|G(t)|") +
  theme_classic(base_size = 10) +
  theme(strip.background = element_blank(), strip.text = element_blank(),
        panel.spacing.x = unit(1.5, "pt"),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        axis.line.y = element_blank(), axis.ticks.y = element_line(linewidth = 0.3),
        axis.title.y = element_text(size = 9, colour = "grey20"), axis.text.y = element_text(size = 7.5, colour = "grey30"),
        plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(2, 6, 4, 4))

## header: panel letter/title at left, colour bar "cold locations -- warm locations" centred
bar <- data.table(x = seq(TLIM[1], TLIM[2], length.out = 200))
header <- ggplot(bar) +
  geom_tile(aes(x = x, y = 0, fill = x), height = 0.6) +
  geom_segment(aes(x = TLIM[1] + c(0.25, 0.5, 0.75) * diff(TLIM), xend = TLIM[1] + c(0.25, 0.5, 0.75) * diff(TLIM), y = -0.45, yend = 0.45),
               data = data.table(), colour = "grey30", linewidth = 0.3) +
  annotate("text", x = TLIM[1] - 0.04 * diff(TLIM), y = 0, label = "cold locations", hjust = 1, size = 3.6, colour = "#2166ac") +
  annotate("text", x = TLIM[2] + 0.04 * diff(TLIM), y = 0, label = "warm locations", hjust = 0, size = 3.6, colour = "#b2182b") +
  scale_fill_gradientn(colours = rdbu, limits = TLIM, guide = "none") +
  scale_x_continuous(limits = TLIM + c(-2.2, 2.2) * diff(TLIM), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1, 1)) +
  theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA), plot.margin = margin(0, 0, 0, 0))

title <- ggplot() + theme_void() +
  annotate("text", 0, 0.6, label = "B", fontface = "bold", hjust = 0, size = 5) +
  annotate("text", 0.08, 0.6, label = "genome-wide\ntotal selection", hjust = 0, vjust = 0.5, size = 3.8, lineheight = 0.9) +
  scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(0, 1)) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

top <- (title | header | plot_spacer()) + plot_layout(widths = c(1.1, 3, 1.1))
fig <- top / strip + plot_layout(heights = c(0.22, 1))

suffix <- if (metric == "gt_abs") "" else paste0("_", metric)
ggsave(sprintf("figs/LINKED_SELECTION_gt_abs_strip%s.pdf", suffix), fig, width = 13, height = 3.6)
ggsave(sprintf("figs/LINKED_SELECTION_gt_abs_strip%s.png", suffix), fig, width = 13, height = 3.6, dpi = 250)
cat(sprintf("Saved figs/LINKED_SELECTION_gt_abs_strip%s.{pdf,png}\n", suffix))

## numbers for the caption
cat("\nmean over gardens of garden-mean Gabs(t) vs null median, by year:\n")
print(merge(obs[, .(obs = round(mean(y), 3), n_gardens = .N), by = year], nul[, .(null_med = round(mean(med), 3), null_hi = round(mean(hi), 3)), by = year], by = "year"))
cat("\ngardens whose mean exceeds the null 97.5% at their last year:\n")
last <- obs[, .SD[which.max(year)], by = site][nul, on = c("site", "year"), nomatch = 0]
cat(sprintf("  %d of %d\n", sum(last$y > last$hi), nrow(last)))
