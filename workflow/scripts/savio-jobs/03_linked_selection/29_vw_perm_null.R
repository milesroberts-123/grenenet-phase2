#!/usr/bin/env Rscript

# Permutation null for the per-site Va estimate (stage 5).
#
# The estimator's signal is a covariance across LOCI between replicate pairs:
# Q() forms, for each pair of replicates (i,j), the covariance over loci of
# their allele-frequency changes. A biological/selection signal makes replicate
# i and replicate j move together at the SAME locus.
#
# A label permutation -- shuffling which replicate is which -- cannot test this:
# Q's off-diagonal entries are just relabelled, and the regression uses every
# upper-triangle entry with one shared coefficient, so the fit is unchanged.
# (Same trap as the sign-permutation null in the G(t) work.)
#
# The null that can fail is a LOCUS permutation: give every replicate its own
# random locus order (applied to that replicate's p1 and p2 rows identically, so
# the replicate's own change series is preserved), which destroys the
# locus-by-locus alignment between replicates while leaving each replicate's
# marginal distribution -- and therefore its diagonal drift/sampling variance --
# untouched. Under the null, Va should collapse toward zero.
#
# Usage: Rscript 29_vw_perm_null.R <panel_csv> <s_csv> <task_index> <out_dir> [n_perm]

suppressPackageStartupMessages({
  library(data.table)
  library(Vw)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4L) stop("usage: <panel_csv> <s_csv> <task_index> <out_dir> [n_perm]")
panel_csv <- args[[1]]
s_csv     <- args[[2]]
task      <- as.integer(args[[3]])
out_dir   <- args[[4]]
n_perm    <- if (length(args) >= 5L) as.integer(args[[5]]) else 200L
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ts <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))
t0 <- proc.time()

panel <- fread(panel_csv)
s_tab <- fread(s_csv)
pool_cols <- setdiff(names(panel), c("chrom", "pos", "class", "seedmix"))
site_of   <- sub("_.*", "", pool_cols)
sites     <- sort(unique(as.integer(site_of)))
if (task < 0L || task >= length(sites)) stop("task ", task, " outside 0..", length(sites) - 1L)
this_site <- sites[task + 1L]
set.seed(9000 + this_site)

cols     <- pool_cols[site_of == as.character(this_site)]
gens     <- as.integer(sub("^[0-9]+_([0-9]+)_.*$", "\\1", cols))
plot_ids <- sub("^[0-9]+_[0-9]+_", "", cols)
M    <- as.matrix(panel[, ..cols])
seed <- panel$seedmix
nl   <- nrow(M)
ts("site ", this_site, ": ", nl, " loci x ", ncol(M), " pools; gens ",
   paste(sort(unique(gens)), collapse = ","))

s_for <- function(g) {
  hit <- s_tab[ngen == g, s]
  if (length(hit) != 1L) stop("no s for ngen = ", g)
  hit
}

# Va from a (pbar1, pbar2) pair, both replicates x loci
va_from <- function(p1, p2, s) {
  Qm <- Q(p1, p2)
  q  <- Qm[upper.tri(Qm, diag = TRUE)]
  a  <- rep(0.5 * s, length(q))
  b  <- diag(nrow(Qm))[upper.tri(Qm, diag = TRUE)]
  unname(coef(lm(q ~ a + b - 1))["a"])
}

# null: permute loci within each replicate, jointly for its p1 and p2 rows
va_null <- function(p1, p2, s) {
  nr <- nrow(p2)
  perm <- t(vapply(seq_len(nr), function(i) sample.int(ncol(p2)), integer(ncol(p2))))
  p1n <- p1; p2n <- p2
  for (i in seq_len(nr)) {
    p1n[i, ] <- p1[i, perm[i, ]]
    p2n[i, ] <- p2[i, perm[i, ]]
  }
  va_from(p1n, p2n, s)
}

ug <- sort(unique(gens))
win_tab <- data.table(window = c("base_to_first", "base_to_last"),
                      t0 = c(0L, 0L), tau = c(min(ug), max(ug)))
adj <- ug[which(diff(ug) == 1L)]
if (length(adj)) {
  win_tab <- rbind(win_tab, data.table(window = paste0("gen", adj, "_to_gen", adj + 1L),
                                       t0 = adj, tau = adj + 1L))
}
win_tab <- unique(win_tab)

res <- list()
for (w in seq_len(nrow(win_tab))) {
  win <- win_tab$window[w]; tt <- win_tab$t0[w]; g <- win_tab$tau[w]
  sel <- which(gens == g)
  if (length(sel) < 2L) next
  if (tt == 0L) {
    pbar2 <- t(M[, sel, drop = FALSE])
    pbar1 <- matrix(rep(seed, each = length(sel)), nrow = length(sel))
  } else {
    bs <- which(gens == tt)
    common <- intersect(plot_ids[bs], plot_ids[sel])
    if (length(common) < 2L) next
    pbar1 <- t(M[, bs[match(common, plot_ids[bs])], drop = FALSE])
    pbar2 <- t(M[, sel[match(common, plot_ids[sel])], drop = FALSE])
  }

  # complete loci only; Vw::Q() aborts on NA
  ok <- colSums(is.na(pbar2)) == 0L & colSums(is.na(pbar1)) == 0L
  if (sum(ok) < 100L) next
  pbar1 <- pbar1[, ok, drop = FALSE]; pbar2 <- pbar2[, ok, drop = FALSE]
  s <- s_for(g - tt)

  va_obs <- va_from(pbar1, pbar2, s)
  va_nul <- vapply(seq_len(n_perm), function(r) va_null(pbar1, pbar2, s), numeric(1))

  res[[length(res) + 1L]] <- data.table(
    site = this_site, window = win, t0 = tt, ngen = g - tt,
    n_plots = nrow(pbar2), n_loci = ncol(pbar2), s = s,
    Va = va_obs,
    Va_null_mean = mean(va_nul), Va_null_sd = sd(va_nul),
    Va_null_q975 = unname(quantile(va_nul, 0.975)),
    p_perm = (1 + sum(va_nul >= va_obs)) / (n_perm + 1),
    n_perm = n_perm
  )
  ts(win, ": Va=", signif(va_obs, 4), "  null ", signif(mean(va_nul), 3),
     " +/- ", signif(sd(va_nul), 3), "  p=", signif((1 + sum(va_nul >= va_obs)) / (n_perm + 1), 3))
}

if (!length(res)) stop("site ", this_site, ": nothing computed")
out <- file.path(out_dir, sprintf("vw_perm_site%02d.csv", this_site))
fwrite(rbindlist(res), out)
ts("wrote ", out, "  elapsed ", round((proc.time() - t0)["elapsed"] / 60, 1), " min")
