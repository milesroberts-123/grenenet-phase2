#!/usr/bin/env Rscript

# Per-variant selection coefficients from the kmate allele-frequency time series,
# one job-array task per site.
#
# Estimator: rventures::wright_s_estim(), the first-order approximation to the
# Mathieson & McVean (2013) MLE for a diploid Wright-Fisher population with
# additive selection (h = 0.5). It returns the numerator 2*(f_T - f_0) and the
# denominator dt*sum(f_t*(1-f_t)) when called with components = TRUE, which lets
# the replicates within a site be combined as a ratio of averages:
#
#   s_site = sum_over_reps(numerator) / sum_over_reps(denominator)
#
# rather than an average of ratios, so no replicate is over-weighted by a small
# denominator.
#
# Trajectories:
#   - one trajectory per (variant, replicate). Replicate = the pool/plot label,
#     the third field of a "<site>_<generation>_<replicate>" column name.
#   - the seedmix frequency is generation 0 and is prepended to any segment whose
#     first observed generation is 1.
#   - generations are used only while consecutive: a missing generation truncates
#     the trajectory, so no gap is ever bridged. dt is the span of the run.
#
# Usage: Rscript 23_wright_s_estim_kmate.R <task_index>   (0-based SLURM_ARRAY_TASK_ID)

ts_message <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))
t0 <- proc.time()
options(warn = 1)
# arrow's threaded CSV scanner reassembles scan tasks out of order, so each
# collect() of the large kmate trays table returns a different row permutation.
# The notebook sets this in its first chunk, which job scripts skip.
options(arrow.use_threads = FALSE)

suppressPackageStartupMessages({
    library(arrow)
    library(data.table)
    library(dplyr)
    library(stringr)
})

box::use(rventures[...])

ts_message("directories")
kmate_dir <- "/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/resources/kmate_allele_frequencies/"
af_prefix <- "/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/data/PHASE_1_2_COMBINED_kmate_arch3v2_allele_frequency"

today <- Sys.Date()
results_dir <- paste0("/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/results/", today)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

ts_message("open tables")
# seedmix: generation-0 frequencies (column X0_0_0) plus the variant key
gen0 <- fread(paste0(af_prefix, "_seedmix_randswitch_fixrm_biallelic.csv"))
hapfire <- open_dataset(paste0(af_prefix, "_trays_randswitch_fixrm_biallelic.csv"), format = "csv")


#' Estimate selection coefficients for every variant x segment at one site
#'
#' @param freq data.frame of allele-frequency columns for one site, named
#'   "<site>_<generation>_<replicate>"; rows are variants in catalogue order
#' @param seedmix numeric, the generation-0 frequency of each variant
#' @param min_start_maf variants whose trajectory starts at or below this minor
#'   allele frequency are ignored; applied symmetrically, so a start at or above
#'   1 - min_start_maf is ignored too
#' @return data.table with one row per variant x segment boundary pair
estimate_site <- function(freq, seedmix, min_start_maf = 0.05) {

    cols <- colnames(freq)
    # all columns belong to one site: "<site>_<generation>_<replicate>"
    site_prefix <- sub("_[0-9]+_[0-9]+$", "", cols[1L])
    parts <- strsplit(cols, "_", fixed = TRUE)
    gen <- as.integer(vapply(parts, `[`, character(1), 2L))
    rep <- as.integer(vapply(parts, `[`, character(1), 3L))

    # one entry per replicate, holding its segments as (variant rows) x (sums)
    pieces <- list()

    for (r in sort(unique(rep))) {
        gens_r <- sort(unique(gen[rep == r]))
        # split into runs of consecutive generations
        runs <- split(gens_r, cumsum(c(1L, diff(gens_r) > 1L)))

        for (run in runs) {
            # the seedmix is generation 0, so it belongs only to a run that
            # starts at generation 1; otherwise prepending it would bridge a gap
            use_seed <- min(run) == 1L
            seg_start <- if (use_seed) 0L else min(run)
            seg_end <- max(run)
            dt <- seg_end - seg_start
            if (dt < 1L) next

            seg_cols <- paste0(site_prefix, "_", run, "_", r)
            m <- as.matrix(freq[, seg_cols, drop = FALSE])
            if (use_seed) m <- cbind(seedmix, m)

            p_start <- m[, 1L]
            # variant filters: starting MAF, no missing values, no absorbing state
            keep <- !is.na(p_start) &
                p_start > min_start_maf & p_start < (1 - min_start_maf)
            keep <- keep & stats::complete.cases(m)
            # boundary check only where NAs are already excluded, otherwise
            # rowSums() propagates NA into keep
            if (any(keep)) {
                keep[keep] <- rowSums(m[keep, , drop = FALSE] <= 0 |
                                          m[keep, , drop = FALSE] >= 1) == 0L
            }

            n_drop_start <- sum(!is.na(p_start) &
                                    (p_start <= min_start_maf | p_start >= (1 - min_start_maf)))
            n_drop_na <- sum(!is.na(p_start) &
                                 p_start > min_start_maf & p_start < (1 - min_start_maf) &
                                 !stats::complete.cases(m))
            ts_message("rep ", r, " gens ", seg_start, "-", seg_end, ": kept ", sum(keep),
                       " of ", nrow(m), " variants (dropped ", n_drop_start,
                       " at the starting MAF, ", n_drop_na, " with NAs)")

            if (!any(keep)) next
            ms <- m[keep, , drop = FALSE]
            idx <- which(keep)

            comp <- lapply(seq_len(nrow(ms)), function(i) {
                suppressWarnings(
                    wright_s_estim(ms[i, ], dt = dt, components = TRUE)
                )
            })

            pieces[[length(pieces) + 1L]] <- data.table(
                idx = idx,
                seg_start = seg_start,
                seg_end = seg_end,
                dt = dt,
                rep = r,
                numerator = vapply(comp, `[[`, numeric(1), "numerator"),
                denominator = vapply(comp, `[[`, numeric(1), "denominator")
            )
        }
    }

    if (length(pieces) == 0L) return(NULL)
    per_rep <- rbindlist(pieces)

    # ratio of averages across replicates sharing the same segment boundaries
    out <- per_rep[, .(
        n_obs = .N,
        sum_num = sum(numerator),
        sum_den = sum(denominator)
    ), by = .(idx, seg_start, seg_end, dt)]
    out[, s := sum_num / sum_den]
    out[]
}


## ---- job array: one site per task ---------------------------------------
site_list <- sort(as.integer(unique(gsub("_.*", "", colnames(hapfire)))))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("expected a 0-based task index, e.g. Rscript 23_wright_s_estim_kmate.R 0")
task <- suppressWarnings(as.integer(args[[1]]))
if (is.na(task) || task < 0L || task >= length(site_list)) {
    stop("task index ", args[[1]], " outside 0..", length(site_list) - 1L)
}
site <- site_list[task + 1L]

output_s <- paste0(results_dir, "/wright_s_estim_site", site, ".csv")
check_file_exists(output_s)
ts_message("task ", task, " of ", length(site_list) - 1L, " -> site ", site,
           " -> ", basename(output_s))

ts_message("collect site ", site)
hap_site <- hapfire %>%
    dplyr::select(starts_with(paste0(site, "_"))) %>%
    collect()
hap_site <- as.data.frame(hap_site)

ts_message("estimate site ", site)
res <- estimate_site(freq = hap_site, seedmix = gen0$X0_0_0)

if (is.null(res)) {
    ts_message("no estimable segments at site ", site)
} else {
    # variant key from the catalogue order; every column of hap_site shares it
    res$chrom <- gen0$chrom[res$idx]
    res$pos   <- gen0$pos[res$idx]
    res$class <- gen0$class[res$idx]
    setcolorder(res, c("chrom", "pos", "class", "seg_start", "seg_end",
                       "dt", "n_obs", "s", "sum_num", "sum_den"))
    setorder(res, chrom, pos, seg_start, seg_end)
    res$idx <- NULL
    fwrite(res, output_s)
    ts_message("wrote ", nrow(res), " rows to ", basename(output_s))
}

ts_message("DONE site ", site, " elapsed ", round((proc.time() - t0)["elapsed"] / 60, 1), " minutes")
