#!/usr/bin/env Rscript

# Stage 4 of the Vw analysis: per-site estimate of Va from the frozen panel.
#
# For one site, each generation is a set of replicate populations (the plots).
# Vw's est_Va_bc regresses the among-replicate covariance of allele-frequency
# change on the expected covariance due to linked selection, which averageLD()
# collapses to a single scalar s:
#
#   q = Va * (0.5*s) + b * I(diagonal)          (q from Q(), diagonal = drift)
#   Va = coefficient on a,  Ne = 0.5 / coefficient on b
#
# s is computed once for the panel (27_vw_panel_extract.R) and passed in, so the
# full L and nR matrices never enter a per-site task. The regression below is the
# same three lines est_Va_bc runs; vw_validate_mimic.R checks it against
# est_Va_bc on a case small enough for the real matrices.
#
# pbar1 is the seed mix (generation 0) repeated over the site's plots, pbar2 is
# the site at generation g, so ngen = g.
#
# Usage: Rscript 28_vw_site_va.R <panel_csv> <s_csv> <task_index> <out_dir>

suppressPackageStartupMessages({
  library(data.table)
  library(Vw)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4L) stop("usage: <panel_csv> <s_csv> <task_index> <out_dir>")
panel_csv <- args[[1]]
s_csv     <- args[[2]]
task      <- as.integer(args[[3]])
out_dir   <- args[[4]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ts <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))
t0 <- proc.time()

panel <- fread(panel_csv)
s_tab <- fread(s_csv)
pool_cols <- setdiff(names(panel), c("chrom", "pos", "class", "seedmix"))
site_of <- sub("_.*", "", pool_cols)
sites <- sort(unique(as.integer(site_of)))
if (is.na(task) || task < 0L || task >= length(sites)) {
  stop("task ", task, " outside 0..", length(sites) - 1L)
}
this_site <- sites[task + 1L]
ts("task ", task, " -> site ", this_site, " (", length(sites), " sites)")

cols <- pool_cols[site_of == as.character(this_site)]
gens <- as.integer(sub("^[0-9]+_([0-9]+)_.*$", "\\1", cols))
plot_ids <- sub("^[0-9]+_[0-9]+_", "", cols)
if (anyNA(gens)) stop("could not parse generations from: ",
                      paste(head(cols[is.na(gens)]), collapse = ", "))

M    <- as.matrix(panel[, ..cols])       # loci x pools
seed <- panel$seedmix
ts("panel: ", nrow(M), " loci x ", ncol(M), " pools; generations ",
   paste(sort(unique(gens)), collapse = ","))

# seed mix is generation 0, so a site's generation g is a g-generation window
s_for <- function(g) {
  hit <- s_tab[ngen == g, s]
  if (length(hit) != 1L) stop("no s for ngen = ", g, " in ", s_csv)
  hit
}

# Windows to evaluate. "seedmix" windows run from the shared base, "consecutive"
# windows run between two measured generations (the paper's Box 1 design: change
# between generations t and tau after the base, not from the base itself). Only
# adjacent measured generations are used for the consecutive set, so ngen = 1 and
# the package's own nR^(2*ngen) multiplier is exactly the L_m the derivation wants.
mode <- if (length(args) >= 5L) args[[5]] else "both"
ug <- sort(unique(gens))
win_tab <- data.table(window = character(), t0 = integer(), tau = integer())
if (mode %in% c("seedmix", "both")) {
  win_tab <- rbind(win_tab, data.table(window = c("base_to_first", "base_to_last"),
                                       t0 = 0L, tau = c(min(ug), max(ug))))
}
if (mode %in% c("consecutive", "both")) {
  adj <- ug[which(diff(ug) == 1L)]
  if (length(adj)) {
    win_tab <- rbind(win_tab, data.table(window = paste0("gen", adj, "_to_gen", adj + 1L),
                                         t0 = adj, tau = adj + 1L))
  }
}
win_tab <- unique(win_tab)

res <- list()
for (w in seq_len(nrow(win_tab))) {
  win <- win_tab$window[w]
  t0  <- win_tab$t0[w]
  g   <- win_tab$tau[w]
  sel <- which(gens == g)
  if (length(sel) < 2L) {
    ts("site ", this_site, " ", win, ": only ", length(sel),
       " pool(s) at generation ", g, " -- skipping")
    next
  }
  base_sel <- if (t0 == 0L) NULL else which(gens == t0)
  if (t0 != 0L && length(base_sel) < 2L) {
    ts("site ", this_site, " ", win, ": only ", length(base_sel),
       " pool(s) at generation ", t0, " -- skipping")
    next
  }

  # Two measured generations must be paired by PLOT, not by position: the same
  # plot is the replicate, and a site need not have the same plots in both
  # generations. Only plots present at both times enter the comparison.
  if (t0 == 0L) {
    pbar2   <- t(M[, sel, drop = FALSE])
    base_mat <- matrix(rep(seed, each = length(sel)), nrow = length(sel))
  } else {
    p2_ids <- plot_ids[sel]
    p1_ids <- plot_ids[base_sel]
    common <- intersect(p1_ids, p2_ids)
    if (length(common) < 2L) {
      ts("site ", this_site, " ", win, ": only ", length(common),
         " plot(s) present at both generations -- skipping")
      next
    }
    pbar2    <- t(M[, sel[match(common, p2_ids)], drop = FALSE])
    base_mat <- t(M[, base_sel[match(common, p1_ids)], drop = FALSE])
  }
  sel <- seq_len(nrow(pbar2))

  # Vw::Q() fails on pre-existing NA (its `if (any(pbar2 == 0) | ...)` returns NA
  # and if() errors), so loci incomplete in either the seed mix or this site's
  # pools are dropped before the call. Frequencies of exactly 0/1 are handled
  # inside Q() itself and are left alone.
  # completeness per locus (columns), not per row: complete.cases() on a matrix
  # reports per-row completeness and would be recycled against seed, silently
  # keeping the wrong loci
  ok <- colSums(is.na(pbar2)) == 0L & colSums(is.na(base_mat)) == 0L
  if (sum(ok) < 100L) {
    ts("site ", this_site, " ", win, ": only ", sum(ok),
       " complete loci -- skipping")
    next
  }
  pbar2 <- pbar2[, ok, drop = FALSE]
  pbar1 <- base_mat[, ok, drop = FALSE]

  Qm <- Q(pbar1, pbar2)
  s  <- s_for(g - t0)
  q  <- Qm[upper.tri(Qm, diag = TRUE)]
  a  <- rep(0.5 * s, length(q))
  b  <- diag(nrow(Qm))[upper.tri(Qm, diag = TRUE)]
  fit <- lm(q ~ a + b - 1)

  res[[length(res) + 1L]] <- data.table(
    site = this_site, window = win, t0 = t0, ngen = g - t0,
    n_plots = length(sel), n_loci = ncol(pbar2),
    # M is loci x this site's pools, so the panel size is nrow(M)
    frac_loci_used = ncol(pbar2) / nrow(M),
    s = s,
    Va = unname(coef(fit)["a"]),
    Ne = 0.5 / unname(coef(fit)["b"]),
    n_q_offdiag = sum(upper.tri(Qm))
  )
  ts("site ", this_site, " ", win, ": ", t0, "->", g, "  Va=",
     signif(unname(coef(fit)["a"]), 4), "  Ne=",
     signif(0.5 / unname(coef(fit)["b"]), 4))
}

if (!length(res)) stop("site ", this_site, ": nothing to write")
out <- file.path(out_dir, sprintf("vw_site%02d.csv", this_site))
fwrite(rbindlist(res), out)
ts("wrote ", out, "  elapsed ", round((proc.time() - t0)["elapsed"] / 60, 1), " min")
