## ---------------------------------------------------------------------------
## Tier-1 screen aggregator
##
## Collects the per-event screen rows written by 09s-ridge-screen.R, computes
## FDR q-values across all events (Storey qvalue → BH fallback), and writes:
##   processed_data/event_models/<tf>/screen_results.csv.gz   (all events + q)
##   processed_data/event_models/<tf>/tier1_hits_<tf>.txt      (hit event IDs)
##
## A "hit" = q < q_threshold AND effect > 0 (real screen_R2 above the matched
## feature-rotation null — event-specific epigenetic signal, not a cell-type
## fingerprint). The Tier-2 elastic-net (09zz via 09-1) runs on this hit list.
##
## CLI:
##   Rscript 09s-aggregate.R <cfg_rds> [q_threshold]
## ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
cfg_path <- args[1]
q_threshold <- if (length(args) >= 2L) as.numeric(args[2L]) else {
  getOption("EpiATLAS_AS_SCREEN_Q", 0.1)
}

cfg <- readRDS(cfg_path)
setwd(cfg$project_dir)
event_dir <- cfg$event_dir
tf <- basename(event_dir) # event_models/<tf>
screen_dir <- file.path(event_dir, "screen")

files <- list.files(screen_dir, pattern = "_screen\\.csv\\.gz$", full.names = TRUE)
if (length(files) == 0L) {
  stop("No per-event screen files in ", screen_dir)
}
# Completeness guard: every dispatched event should have written a file (even an
# NA row). Far fewer means the screen array is still running / partially failed —
# aggregating now would give incomplete q-values. Warn loudly but proceed (so a
# deliberate partial aggregation is still possible).
all_ids_file <- file.path("processed_data", sprintf("event_glmnet_all_ids_%s.txt", tf))
if (file.exists(all_ids_file)) {
  n_expected <- length(readLines(all_ids_file))
  if (length(files) < n_expected) {
    warning(sprintf(
      "Only %d/%d per-event screen files present — the screen array may still be ",
      length(files), n_expected
    ), "running; q-values will be computed on the partial set.", immediate. = TRUE)
  }
}
dt <- data.table::rbindlist(lapply(files, data.table::fread), fill = TRUE)

# LONG-REUSE (fallback): the previous omnibus screen (single ridge over ALL
# epigenetic features) IS the `long` feature set. 09s-ridge-screen.R's per-feature-set
# resume gate normally retains those rows in the per-event files themselves (tagging
# untagged legacy rows as "long"), so this hook only fires if the per-event long rows
# are gone — e.g. the screen dir was cleared but screen_long_results.csv.gz was kept.
long_reuse <- file.path(event_dir, "screen_long_results.csv.gz")
have_long <- "feature_set" %in% names(dt) && any(dt$feature_set == "long")
if (file.exists(long_reuse) && !have_long) {
  ldt <- data.table::fread(long_reuse)
  if (!"feature_set" %in% names(ldt)) ldt[, feature_set := "long"]
  if (!"transcript_filter" %in% names(ldt)) ldt[, transcript_filter := tf]
  if (!"Event Type" %in% names(ldt) || anyNA(ldt$`Event Type`)) {
    sess <- readRDS(cfg$session_rds)
    et_map <- sess$event_dt[, .(ID, ET = `Event Type`)]
    ldt[et_map, on = "ID", `Event Type` := i.ET]
    rm(sess)
  }
  dt <- dt[feature_set != "long"] # drop any freshly-computed long, prefer reused
  dt <- data.table::rbindlist(list(dt, ldt), fill = TRUE)
  message(sprintf("Reused %d long rows from %s", nrow(ldt), long_reuse))
}

# Feature-set / event-type keys for grouped FDR (older single-row screens lack them).
if (!"feature_set" %in% names(dt)) dt[, feature_set := "long"]
if (!"transcript_filter" %in% names(dt)) dt[, transcript_filter := tf]
message(sprintf("Aggregated %d rows (%d with finite p_emp) over %d feature set(s)",
                nrow(dt), sum(is.finite(dt$p_emp)),
                data.table::uniqueN(dt$feature_set)))

# --- FDR q-values, computed WITHIN each (transcript_filter, Event Type,
#     feature_set) group. Per-space + per-event-type FDR: a local-only signal is
#     tested against other local screens, not diluted by long/short. Storey qvalue
#     (estimates pi0 -> more power); BH fallback when a group is too small/degenerate.
qfun <- function(p) {
  if (length(p) == 0L) return(numeric(0))
  tryCatch(
    qvalue::qvalue(p)$qvalues,
    error = function(e) {
      message("qvalue failed (", conditionMessage(e), ") -> BH fallback for a group")
      stats::p.adjust(p, method = "BH")
    }
  )
}
dt[, q := NA_real_]
dt[
  is.finite(p_emp),
  q := qfun(p_emp),
  by = .(transcript_filter, `Event Type`, feature_set)
]

data.table::setorder(dt, `Event Type`, feature_set, q, -effect, na.last = TRUE)
out_results <- file.path(event_dir, "screen_results.csv.gz")
data.table::fwrite(dt, out_results)

# Hit = significant in ANY feature set (Tier-2 then decomposes long/short/local).
hit_rows <- dt[is.finite(q) & q < q_threshold & is.finite(effect) & effect > 0]
hits <- sort(unique(hit_rows$ID))
hits_file <- file.path(event_dir, sprintf("tier1_hits_%s.txt", tf))
writeLines(as.character(hits), hits_file)

message(sprintf(
  "q<%.3g & effect>0 -> %d unique hit events (%d feature-set rows). Wrote %s + %s",
  q_threshold, length(hits), nrow(hit_rows), out_results, hits_file
))
# Per-feature-set breakdown.
print(dt[is.finite(q), .(
  n_tested = .N,
  n_hits = sum(q < q_threshold & effect > 0, na.rm = TRUE)
), by = .(`Event Type`, feature_set)][order(`Event Type`, feature_set)])
