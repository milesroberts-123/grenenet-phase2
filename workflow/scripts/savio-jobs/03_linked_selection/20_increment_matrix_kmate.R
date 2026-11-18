#!/usr/bin/env Rscript

# Build the full kmate allele-frequency increment matrix.
#
# For every (site, replicate) series, take consecutive sampled generations and
# difference them, so each column of the output is one interval: allele frequency
# change from generation a to generation b within one replicate. Columns are
# named <site>_<a>-<b>_<rep> so a gap in the sampling schedule stays visible
# (e.g. 27_2-5_3) instead of being implied as contiguous.
#
# 893 pool columns with 352 (site, replicate) series give 541 intervals, and
# 6,207,205 variants -> ~27 GB in memory, ~40 GB as CSV. The source table is
# therefore streamed in row batches and appended to disk; the frequencies are
# never held in memory as a whole.
#
# Outputs (in data/):
#   PHASE_1_2_COMBINED_kmate_arch3v2_allele_frequency_increment_matrix.csv
#   PHASE_1_2_COMBINED_kmate_arch3v2_allele_frequency_increment_plan.csv

suppressMessages({
  library(arrow)
  library(data.table)
})

ts_message <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))
}

# arrow's threaded CSV scanner reassembles scan tasks out of order; the notebook
# sets this in its first chunk, which job scripts skip. Every cbind against the
# seedmix depends on file order.
options(arrow.use_threads = FALSE)

data_dir <- "/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/data/"
af_prefix <- paste0(data_dir, "PHASE_1_2_COMBINED_kmate_arch3v2_allele_frequency")
trays_csv <- paste0(af_prefix, "_trays_randswitch_fixrm_biallelic.csv")
seedmix_csv <- paste0(af_prefix, "_seedmix_randswitch_fixrm_biallelic.csv")

out_matrix <- paste0(af_prefix, "_increment_matrix.csv")
out_plan   <- paste0(af_prefix, "_increment_plan.csv")

if (!file.exists(trays_csv)) stop("trays table not found: ", trays_csv)

## ---- column map ---------------------------------------------------------
hdr <- strsplit(readLines(trays_csv, n = 1L), ",", fixed = TRUE)[[1]]
meta <- data.table(
  col  = hdr,
  site = as.integer(sub("_.*", "", hdr)),
  gen  = as.integer(sub("^[^_]*_([0-9]+)_.*$", "\\1", hdr)),
  rep  = sub("^.*_", "", hdr)
)
if (anyNA(meta$site) || anyNA(meta$gen)) stop("could not parse site/gen from the header")
ts_message("pool columns: ", nrow(meta), " | series: ",
           uniqueN(meta[, .(site, rep)]))

## ---- interval plan ------------------------------------------------------
# one row per interval: which column it is subtracted from, which from, and the
# generations involved. The first interval of every series starts at gen 0,
# whose frequencies come from the seedmix rather than from a trays column.
series <- unique(meta[, .(site, rep)])[order(site, rep)]
plan_list <- vector("list", nrow(series))
for (k in seq_len(nrow(series))) {
  d <- meta[site == series$site[k] & rep == series$rep[k]][order(gen)]
  if (!nrow(d)) next
  g <- d$gen
  # every sampled generation ends exactly one interval; the interval starts at
  # generation 0 (seedmix) or at the previous sampled generation
  from_gen <- c(0L, head(g, -1L))
  from_col <- c(NA_character_, head(d$col, -1L))
  if (length(from_gen) != length(g) || length(from_col) != length(g)) {
    stop("interval plan length mismatch for site ", series$site[k],
         " rep ", series$rep[k])
  }
  plan_list[[k]] <- data.table(
    site     = series$site[k],
    rep      = series$rep[k],
    from_gen = from_gen,
    to_gen   = g,
    from_col = from_col,
    to_col   = d$col
  )
}
plan <- rbindlist(plan_list)
plan[, interval := paste0(site, "_", from_gen, "-", to_gen, "_", rep)]
if (anyDuplicated(plan$interval)) stop("interval names are not unique")

# the plan must agree with the header: each to_col carries the generation it
# claims, and each from_col likewise. This is what catches an off-by-one in the
# interval construction, which the row count cannot.
if (!all(meta$gen[match(plan$to_col, meta$col)] == plan$to_gen)) {
  stop("plan/header mismatch in to_col generations")
}
has_from <- !is.na(plan$from_col)
if (!all(meta$gen[match(plan$from_col[has_from], meta$col)] == plan$from_gen[has_from])) {
  stop("plan/header mismatch in from_col generations")
}

ts_message("intervals: ", nrow(plan), " (from ", nrow(meta), " columns and ",
           nrow(series), " series)")

# expected number of variants, used as the streaming assertion
gen0 <- fread(seedmix_csv, select = "X0_0_0")
n_expected <- nrow(gen0)
ts_message("variants expected: ", format(n_expected, big.mark = ","))

## ---- stream the source table -------------------------------------------
ds <- open_dataset(
  trays_csv, format = "csv",
  read_options = csv_read_options(block_size = 128L * 1024L * 1024L)
)

to_idx   <- match(plan$to_col, hdr)
from_idx <- match(plan$from_col, hdr)   # NA for gen-0 intervals

if (file.exists(out_matrix)) unlink(out_matrix)

reader <- ds$NewScan()$UseThreads(FALSE)$BatchSize(100000L)$Finish()$ToRecordBatchReader()

n_done  <- 0L
n_batch <- 0L
t0 <- Sys.time()

while (TRUE) {
  batch <- reader$read_next_batch()
  if (is.null(batch)) break
  n_batch <- n_batch + 1L

  pm <- as.matrix(batch$to_data_frame())
  rows <- nrow(pm)

  # append the seedmix as an extra column so the gen-0 intervals can be taken
  # with the same subtraction as the rest
  src <- cbind(pm, `__gen0__` = gen0$X0_0_0[n_done + seq_len(rows)])
  from_idx2 <- ifelse(is.na(from_idx), ncol(src), from_idx)

  inc <- src[, to_idx, drop = FALSE] - src[, from_idx2, drop = FALSE]
  colnames(inc) <- plan$interval

  fwrite(as.data.frame(inc), out_matrix, append = n_batch > 1L,
         col.names = n_batch == 1L, quote = FALSE)

  n_done <- n_done + rows
  if (n_batch %% 50L == 0L) {
    ts_message("batches: ", n_batch, " | rows: ", format(n_done, big.mark = ","),
               " | elapsed: ", round(as.numeric(difftime(Sys.time(), t0, units = "min")), 1), " min")
  }
  rm(pm, src, inc); gc(verbose = FALSE)
}

if (n_done != n_expected) {
  stop("row count mismatch: wrote ", n_done, " but the seedmix has ", n_expected,
       " variants")
}
ts_message("wrote ", format(n_done, big.mark = ","), " rows x ", nrow(plan),
           " interval columns")

## ---- plan + spot check --------------------------------------------------
setcolorder(plan, c("site", "rep", "from_gen", "to_gen", "from_col", "to_col", "interval"))
fwrite(plan, out_plan, quote = FALSE)

# recompute row 1 of every interval directly from the raw table and compare with
# the first output row: catches a mis-mapped column, which the row count alone
# cannot detect
raw1 <- fread(trays_csv, nrows = 1L, header = TRUE, showProgress = FALSE)
to_val <- as.numeric(unlist(raw1[1, plan$to_col, with = FALSE], use.names = FALSE))
from_val <- rep(gen0$X0_0_0[1], nrow(plan))
has_from <- !is.na(plan$from_col)
from_val[has_from] <- as.numeric(
  unlist(raw1[1, plan$from_col[has_from], with = FALSE], use.names = FALSE)
)
expect <- to_val - from_val

got1 <- fread(out_matrix, nrows = 1L, showProgress = FALSE)
if (!identical(names(got1), plan$interval)) {
  stop("spot check: output column names do not match the interval plan")
}
got_val <- as.numeric(got1[1, ])

bad <- which(!(is.na(expect) & is.na(got_val)) &
               (is.na(expect) | is.na(got_val) | abs(expect - got_val) > 1e-10))
if (length(bad)) {
  stop("spot check failed for ", length(bad), " interval(s); first: ",
       plan$interval[bad[1]], " expected ", expect[bad[1]], " got ", got_val[bad[1]])
}
ts_message("spot check passed on row 1 across ", nrow(plan), " intervals")
ts_message("done")
