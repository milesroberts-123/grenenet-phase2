#!/usr/bin/env Rscript

# One explicit correction step: relabel generation 4 in the biallelic pool table.
#
# The phase-2 sample sheet (years4-5/samples_data_phase2.csv, and therefore
# gen4.rowmeta.csv) labels two gardens as site 1 and site 2. They are site 5 and
# site 4. Generation 5 already uses the correct numbers, and site 4 only appears
# at generation 4, so the correction creates no collisions:
#
#   gen 4 before: 1, 2, 27, 42, 45, 46, 52, 55, 57
#   gen 4 after : 4, 5, 27, 42, 45, 46, 52, 55, 57
#
# Scope: this touches nothing upstream. Every source table keeps its own labels,
# which is why the correction lives here as its own stage rather than being
# folded back into gen4.rowmeta.csv. Note that flower_counts, which is derived
# from the rowmeta, therefore still carries the old labels at generation 4.
#
# Only the header is rewritten, in place. The substitution is single digit to
# single digit, so the line keeps its exact byte length and the 17.5 GB of data
# below it is never read or copied. The mapping is written alongside as a CSV.

suppressMessages(library(data.table))

ts_message <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

af_prefix <- "/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/data/PHASE_1_2_COMBINED_kmate_arch3v2_allele_frequency"
trays <- paste0(af_prefix, "_trays_randswitch_fixrm_biallelic.csv")

## ---- read the header line as raw bytes, terminator included -------------
read_header_raw <- function(path) {
  con <- file(path, "rb"); on.exit(close(con))
  buf <- raw(0)
  repeat {
    b <- readBin(con, "raw", 65536L)
    if (!length(b)) stop("no newline found in ", path)
    buf <- c(buf, b)
    i <- match(as.raw(10L), buf)
    if (!is.na(i)) return(buf[1:i])
  }
}

hdr_raw  <- read_header_raw(trays)
hdr_txt  <- rawToChar(hdr_raw)
eol      <- if (grepl("\r\n$", hdr_txt)) "\r\n" else "\n"
hdr_line <- sub("\r?\n$", "", hdr_txt)
old_names <- strsplit(hdr_line, ",", fixed = TRUE)[[1]]

## ---- the mapping --------------------------------------------------------
map <- data.table(
  from = c("^1_4_", "^2_4_"),
  to   = c("5_4_",   "4_4_")
)
new_names <- old_names
for (k in seq_len(nrow(map))) {
  hit <- grepl(map$from[k], old_names)
  new_names[hit] <- sub(map$from[k], map$to[k], old_names[hit])
  ts_message(sprintf("  %s -> %s : %d columns", map$from[k], map$to[k], sum(hit)))
}

if (identical(old_names, new_names)) {
  ts_message("nothing to rename: this table already carries the corrected labels")
  quit(save = "no")
}

## ---- guards before touching the file ------------------------------------
stopifnot(
  length(new_names) == length(old_names),
  !anyDuplicated(new_names),
  sum(grepl("^5_4_", new_names)) == 12L,
  sum(grepl("^4_4_", new_names)) == 11L,
  !any(grepl("^1_4_|^2_4_", new_names)),
  # untouched neighbours, to prove the anchors did not overreach
  identical(sort(old_names[grepl("^27_4_", old_names)]),
            sort(new_names[grepl("^27_4_", new_names)]))
)

new_line <- paste(new_names, collapse = ",")
stopifnot(nchar(hdr_line, type = "bytes") == nchar(new_line, type = "bytes"))

## ---- in-place rewrite of the header bytes only --------------------------
con <- file(trays, "r+b")
writeBin(charToRaw(paste0(new_line, eol)), con)
close(con)

## ---- verify by re-reading ----------------------------------------------
check <- strsplit(sub("\r?\n$", "", rawToChar(read_header_raw(trays))), ",", fixed = TRUE)[[1]]
stopifnot(identical(check, new_names))
ts_message("header verified: ", length(check), " columns, ",
           sum(grepl("^5_4_", check)), " at 5_4, ", sum(grepl("^4_4_", check)), " at 4_4")
# row count is untouched: only the header bytes are written

## ---- provenance ---------------------------------------------------------
mapping <- data.table(
  table = basename(trays),
  from = c("1_4", "2_4"), to = c("5_4", "4_4"),
  n_columns = c(sum(grepl("^5_4_", check)), sum(grepl("^4_4_", check)))
)
print(mapping)
fwrite(mapping, paste0(af_prefix, "_site_label_map.csv"), quote = FALSE)
ts_message("done")
