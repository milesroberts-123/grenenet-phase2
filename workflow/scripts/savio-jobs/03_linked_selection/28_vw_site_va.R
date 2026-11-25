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

res <- list()
for (win in c("first", "last")) {
  g <- if (win == "first") min(gens) else max(gens)
  sel <- which(gens == g)
  if (length(sel) < 2L) {
    ts("site ", this_site, " ", win, ": only ", length(sel),
       " pool(s) at generation ", g, " -- skipping")
    next
  }

  pbar2 <- t(M[, sel, drop = FALSE])                              # plots x loci

  # Vw::Q() fails on pre-existing NA (its `if (any(pbar2 == 0) | ...)` returns NA
  # and if() errors), so loci incomplete in either the seed mix or this site's
  # pools are dropped before the call. Frequencies of exactly 0/1 are handled
  # inside Q() itself and are left alone.
  # completeness per locus (columns), not per row: complete.cases() on a matrix
  # reports per-row completeness and would be recycled against seed, silently
  # keeping the wrong loci
  ok <- colSums(is.na(pbar2)) == 0L & !is.na(seed)
  if (sum(ok) < 100L) {
    ts("site ", this_site, " ", win, ": only ", sum(ok),
       " complete loci -- skipping")
    next
  }
  pbar2 <- pbar2[, ok, drop = FALSE]
  pbar1 <- matrix(rep(seed[ok], each = length(sel)), nrow = length(sel))

  Qm <- Q(pbar1, pbar2)
  s  <- s_for(g)
  q  <- Qm[upper.tri(Qm, diag = TRUE)]
  a  <- rep(0.5 * s, length(q))
  b  <- diag(nrow(Qm))[upper.tri(Qm, diag = TRUE)]
  fit <- lm(q ~ a + b - 1)

  res[[length(res) + 1L]] <- data.table(
    site = this_site, window = win, ngen = g,
    n_plots = length(sel), n_loci = ncol(pbar2),
    # M is loci x this site's pools, so the panel size is nrow(M)
    frac_loci_used = ncol(pbar2) / nrow(M),
    s = s,
    Va = unname(coef(fit)["a"]),
    Ne = 0.5 / unname(coef(fit)["b"]),
    n_q_offdiag = sum(upper.tri(Qm))
  )
  ts("site ", this_site, " ", win, ": gen ", g, "  Va=",
     signif(unname(coef(fit)["a"]), 4), "  Ne=",
     signif(0.5 / unname(coef(fit)["b"]), 4))
}

if (!length(res)) stop("site ", this_site, ": nothing to write")
out <- file.path(out_dir, sprintf("vw_site%02d.csv", this_site))
fwrite(rbindlist(res), out)
ts("wrote ", out, "  elapsed ", round((proc.time() - t0)["elapsed"] / 60, 1), " min")
