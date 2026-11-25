#!/usr/bin/env Rscript

# Stage 3 of the Vw analysis: freeze one locus panel for the per-site Va runs.
#
# The panel is the first 100,000 founder-eligible loci of the seeded window that
# 26_vw_panel_compare.R used (seed 20261124), so its LD structure and s are the
# ones already reported. Everything per-site needs is written once:
#
#   data/vw_panel_<n>_chr<chrom>.csv   chrom,pos,class,seedmix,<893 pool columns>
#   data/vw_panel_<n>_chr<chrom>_s.csv s for ngen 1..5 (Haldane nR)
#
# Haldane (AtleastOneRecomb = FALSE) matches the rate: 8.06452e-10 is already
# selfing-corrected in stdpopsim (sigma 0.97, Nordborg 2000), so it counts
# effective crossovers, of which most meioses have none.
#
# Usage: Rscript 27_vw_panel_extract.R <founder_cache.tsv> <chrom> <n_loci> <out_dir>

suppressPackageStartupMessages({
  library(data.table)
  library(Vw)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4L) stop("usage: <founder_cache.tsv> <chrom> <n_loci> <out_dir>")
cache_file <- args[[1]]
chrom      <- args[[2]]
n_loci     <- as.integer(args[[3]])
root       <- args[[4]]

rrate   <- 8.06452e-10
hap_len <- 30427671
seedmix_csv <- file.path(root, "data",
  "PHASE_1_2_COMBINED_kmate_arch3v2_allele_frequency_seedmix_randswitch_fixrm_biallelic.csv")
trays_csv <- file.path(root, "data",
  "PHASE_1_2_COMBINED_kmate_arch3v2_allele_frequency_trays_randswitch_fixrm_biallelic.csv")
out_dir <- file.path(root, "results", format(Sys.Date()))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ts <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))
t0 <- proc.time()

# ---------------------------------------------------------------------------
# 1. the same window 26_vw_panel_compare.R drew
# ---------------------------------------------------------------------------
if (!file.exists(cache_file)) stop("founder cache missing: ", cache_file,
  "\n  rerun 26_vw_panel_compare.sbatch first")
C <- fread(cache_file, header = FALSE, showProgress = FALSE)
pos_all <- C[[1]]
if (length(pos_all) < n_loci) stop("cache holds ", length(pos_all), " loci, need ", n_loci)
set.seed(20261124)
start_i <- sample(seq_len(length(pos_all) - n_loci + 1L), 1L)
win <- start_i:(start_i + n_loci - 1L)
win_pos <- pos_all[win]
ts("panel: ", n_loci, " loci, chr", chrom, " ", min(win_pos), "-", max(win_pos))

# ---------------------------------------------------------------------------
# 2. s for ngen 1..5 on this panel (Haldane), from the founder haplotypes
# ---------------------------------------------------------------------------
w_table <- function(ngen) {
  d <- 0:hap_len
  (1 - 0.5 * (1 - exp(-2 * d * rrate)))^(2 * ngen)
}

H <- t(as.matrix(C[win, -1, with = FALSE]))
storage.mode(H) <- "double"
rm(C); invisible(gc())
if (nrow(H) != 462) stop("expected 462 haplotypes, got ", nrow(H))

s_blockwise <- function(H, pos, ngen, block = 500L) {
  nh <- nrow(H); n <- ncol(H)
  Hc <- H - rep(colMeans(H), each = nh)
  vv <- colSums(Hc^2) / (nh - 1)
  wtab <- w_table(ngen)
  weighted <- 0
  for (i0 in seq(1, n, by = block)) {
    i1 <- min(i0 + block - 1L, n); idx <- i0:i1
    Lb <- crossprod(Hc[, idx, drop = FALSE], Hc) / (nh - 1)
    r2 <- (Lb^2) / outer(vv[idx], vv)
    if (i1 < n) {
      cols <- (i1 + 1L):n
      D <- abs(outer(pos[idx], pos[cols], "-"))
      r2s <- r2[, cols, drop = FALSE]
      weighted <- weighted + sum(r2s * wtab[round(D) + 1L])
    }
    if (length(idx) > 1L) {
      sub <- r2[, idx, drop = FALSE]; ut <- upper.tri(sub)
      Ds <- abs(outer(pos[idx], pos[idx], "-"))
      weighted <- weighted + sum(sub[ut] * wtab[round(Ds[ut]) + 1L])
    }
  }
  2 * weighted / (n * (n - 1))
}

s_tab <- data.table(ngen = 1:5,
                    s = vapply(1:5, function(g) {
                      ts("s for ngen=", g)
                      s_blockwise(H, win_pos, g)
                    }, numeric(1)))
rm(H); invisible(gc())
fwrite(s_tab, file.path(root, "data",
                        sprintf("vw_panel_%dk_chr%s_s.csv", as.integer(n_loci / 1000), chrom)))
print(s_tab)

# ---------------------------------------------------------------------------
# 3. locate the panel rows in the seedmix table (row order defines the trays)
# ---------------------------------------------------------------------------
ts("reading seedmix")
sm <- fread(seedmix_csv, select = c("chrom", "pos", "class", "X0_0_0"))
if (nrow(sm) != 6207205L) ts("note: seedmix rows = ", nrow(sm))
idx <- which(sm$chrom == paste0("Chr", chrom) & sm$pos %in% win_pos)
ts("panel rows found in seedmix: ", length(idx), " of ", n_loci)
if (length(idx) < 0.95 * n_loci) stop("only ", length(idx), " panel loci present in the seedmix table")

# row numbers for awk (header excluded, so data line of seedmix row i is i+1)
tmp_idx <- tempfile(fileext = ".txt")
fwrite(data.table(idx + 1L), tmp_idx, col.names = FALSE)

# ---------------------------------------------------------------------------
# 4. stream the trays CSV once and keep only those rows
# ---------------------------------------------------------------------------
hdr <- strsplit(readLines(trays_csv, n = 1L), ",", fixed = TRUE)[[1]]
ts("trays columns: ", length(hdr))
panel_raw <- tempfile(fileext = ".csv")
awk <- sprintf(
  "awk -F, 'NR==FNR{k[$1]=1; next} (FNR in k){print}' %s %s > %s",
  shQuote(tmp_idx), shQuote(trays_csv), shQuote(panel_raw))
ts("extracting panel rows from the trays table (single pass)")
system(awk)
ts("extracted panel file: ", round(file.size(panel_raw) / 2^20), " MB")

P <- fread(panel_raw, header = FALSE, showProgress = FALSE)
if (nrow(P) != length(idx)) stop("extracted ", nrow(P), " rows, expected ", length(idx))
setnames(P, names(P), hdr)

panel <- cbind(sm[idx, .(chrom, pos, class, seedmix = X0_0_0)], P)
out_csv <- file.path(root, "data",
                     sprintf("vw_panel_%dk_chr%s.csv", as.integer(n_loci / 1000), chrom))
fwrite(panel, out_csv)
ts("wrote ", out_csv, " (", nrow(panel), " x ", ncol(panel), ")")

ts("DONE elapsed ", round((proc.time() - t0)["elapsed"] / 60, 1), " min")
