#!/usr/bin/env Rscript
# rMATS atlas summary -- the resource leg's numbers, without reading the resource.
#
# WHAT THIS IS FOR. `11-paper-figures.Rmd` section 1 quotes two things about the atlas: how many
# RNA-seq experiments it covers and how many splicing events it contains. Those are the numbers
# behind the paper's claim that the analysed cohort is a small slice of what was quantified -- the
# manuscript uses 0.8% of the PSI values, 5% of the events and 27% of the samples.
#
# WHY IT DOES NOT AGGREGATE ANYTHING. The obvious route is to count distinct uuids and IDs in the
# aggregated `<ET>.MATS.JC.txt.csv.gz` tables. That is 136 GB of gzip and takes hours -- measured:
# ~10 min/GB, and a first attempt was killed at the 90-minute wall having finished 2 of 6 files.
# It is also unnecessary. Both quantities are structural, not data-dependent:
#
#   uuids  -- the sample set is DECLARED, in `uuid_string` of the rMATS command table. rMATS was
#             invoked once per (paired/single x library type x read length) stratum, so the union
#             over the strata belonging to this transcript_filter IS the atlas cohort.
#   events -- rMATS derives its event set from the GTF, not from the reads, so `fromGTF.<ET>.txt`
#             row counts ARE the event counts. Verified identical across run directories (same GTF)
#             and verified against full decompression of the aggregated tables for four of the five
#             classes: SE 559,157 / RI 2,459 / A3SS 6,937 / A5SS 4,417, exact in every case.
#
# Cost: five annotation files (~22 MB) and a 42-row TSV, versus 136 GB. This is also what keeps
# `11` honest about its own rule that it must never read the rMATS tree -- that read happens once,
# here, on the small files.
#
# NOT the data release. Producing the full 1,522-sample long-format atlas is a separate job that
# must re-derive the per-sample counts from the raw `<ET>.MATS.JC.txt`; see the note in
# `revision/deferred-cosmetic-edits.md`. Nothing here reads or writes those tables.

suppressPackageStartupMessages(library(data.table))
setDTthreads(2L)

tf       <- Sys.getenv("TRANSCRIPT_FILTER", "biotype_filtered")
atlas_root <- Sys.getenv("RMATS_ATLAS_ROOT", "/nfs/data/IHEC/RNAseq/rmats")
cmd_file <- Sys.getenv("RMATS_COMMANDS", "data/rmats_split_post_commands.tsv")
out      <- Sys.getenv("ATLAS_OUT", sprintf("processed_data/atlas_summary_%s.csv", tf))

EVENT_TYPES <- c("SE", "RI", "A3SS", "A5SS", "MXE")

cmd <- fread(cmd_file)
# `out_folder` names the transcript_filter; the same physical strata are re-run per filter, so this
# grep is what restricts the atlas to the filter in question.
runs <- cmd[grepl(tf, out_folder, fixed = TRUE)]
stopifnot(nrow(runs) > 0L)

uuids <- unique(unlist(strsplit(runs$uuid_string, ",", fixed = TRUE)))
uuids <- uuids[nzchar(uuids)]
cat(sprintf("%s: %d rMATS runs (library strata), %d uuids DECLARED in the command table\n",
  tf, nrow(runs), length(uuids)))

# `N` and `uuid_string` are both DECLARATIONS from the command table, so checking one against the
# other checks the table against itself. The declaration is therefore also checked against the DATA:
# one row of one class per run gives the true field count. 14 small reads, nothing against 136 GB.
#
# This check once reported three runs as overstated (single_fr-unstranded_{42,36,62}, 2 declared vs
# 1 quantified) and the resource as 1,519. That was WRONG -- fread was reading the 2-sample value
# "10,11" as the number 10.11, taking the comma for a decimal separator. With colClasses=character
# the count is 1,522 and no run is short. The check is kept because a genuine gap would silently
# mislabel samples, but it is expected to find nothing.
n_declared <- sum(runs$N)
stopifnot(n_declared == sum(lengths(strsplit(runs$uuid_string, ",", fixed = TRUE))))

.probe <- rbindlist(lapply(seq_len(nrow(runs)), function(i) {
  rn <- runs$run_name[i]
  decl <- strsplit(runs$uuid_string[i], ",", fixed = TRUE)[[1]]
  f <- file.path(atlas_root, tf, rn, sprintf("%s.MATS.JC.txt", EVENT_TYPES[1]))
  if (!file.exists(f)) stop("missing ", f)
  # colClasses=character: fread otherwise reads a 2-sample value like "10,11" as the number 10.11,
  # treating the comma as a decimal separator and collapsing two samples into one. That misread is
  # what made this file report 1,519 quantified samples when the true count is 1,522.
  first <- fread(f, select = "IJC_SAMPLE_1", colClasses = c(IJC_SAMPLE_1 = "character"), nrows = 1L)
  nf <- length(strsplit(as.character(first[[1]][1]), ",", fixed = TRUE)[[1]])
  data.table(run_name = rn, declared = length(decl), quantified = nf,
    kept = list(decl[seq_len(min(nf, length(decl)))]))
}))

mism <- .probe[declared != quantified]
if (nrow(mism)) {
  cat("=== command table OVERSTATES these runs (declared > quantified) ===\n")
  print(mism[, .(run_name, declared, quantified)])
  cat("\n")
}

uuids_declared <- uuids
uuids <- unique(unlist(.probe$kept))
cat(sprintf("quantified uuids: %d (declared %d, %d never quantified)\n\n",
  length(uuids), length(uuids_declared), length(uuids_declared) - length(uuids)))

# Event counts from the GTF-derived annotation, read from ONE run directory and then verified to be
# identical in every other -- if rMATS were ever re-run against a different GTF for some strata,
# that assumption breaks silently, so it is checked rather than assumed.
run_dirs <- file.path(atlas_root, tf, runs$run_name)
run_dirs <- run_dirs[dir.exists(run_dirs)]
stopifnot(length(run_dirs) > 0L)

.count_events <- function(d, et) {
  f <- file.path(d, sprintf("fromGTF.%s.txt", et))
  if (!file.exists(f)) return(NA_integer_)
  # header line is not an event
  as.integer(length(readLines(f, warn = FALSE)) - 1L)
}

ev <- rbindlist(lapply(EVENT_TYPES, function(et) {
  per_run <- vapply(run_dirs, .count_events, integer(1), et = et)
  per_run <- per_run[!is.na(per_run)]
  stopifnot(length(per_run) > 0L)
  if (length(unique(per_run)) != 1L) {
    stop("fromGTF.", et, ".txt row count differs across run directories (",
      paste(sort(unique(per_run)), collapse = ", "),
      ") -- the strata were not run against a common GTF, so a single event count is not defined")
  }
  data.table(event_type = et, events = per_run[[1]], n_runs_checked = length(per_run))
}))

res <- data.table(
  transcript_filter = tf,
  event_type        = ev$event_type,
  events            = ev$events,
  # `uuids` is the QUANTIFIED count -- what rMATS actually emitted. `uuids_declared` is what the
  # command table asked for. Both are kept: the difference is a property of the resource worth
  # stating, and keeping only one invites the same overstatement to come back.
  uuids             = length(uuids),
  uuids_declared    = length(uuids_declared),
  runs              = nrow(runs),
  n_runs_checked    = ev$n_runs_checked,
  source            = "fromGTF.<event_type>.txt (events); MATS.JC.txt field counts (uuids)"
)
print(res)
cat(sprintf("\nTOTAL: %s events over %d classes, %d uuids, %d library strata\n",
  format(sum(res$events), big.mark = ","), nrow(res), length(uuids), nrow(runs)))

# The contrast the manuscript draws: how much of the quantified resource the analysis actually used.
mod <- sprintf("processed_data/cohort_counts_%s.csv", tf)
if (file.exists(mod)) {
  cc <- fread(mod)[scope == "modelled"]
  cat(sprintf("modelled cohort is %d of %d uuids (%.1f%%)\n",
    cc$uuids[1], length(uuids), 100 * cc$uuids[1] / length(uuids)))
}

dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
fwrite(res, out)
cat(sprintf("\nwritten: %s\n", out))
