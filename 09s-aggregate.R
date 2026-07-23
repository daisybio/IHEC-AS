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
message(sprintf("Aggregated %d events (%d with finite p_emp)",
                nrow(dt), sum(is.finite(dt$p_emp))))

# --- FDR q-values across events with a finite empirical p -------------------
dt[, q := NA_real_]
ok <- is.finite(dt$p_emp)
if (sum(ok) > 0L) {
  p <- dt$p_emp[ok]
  qv <- tryCatch(
    {
      qobj <- qvalue::qvalue(p)
      qobj$qvalues
    },
    error = function(e) {
      message("qvalue failed (", conditionMessage(e), ") → BH fallback")
      stats::p.adjust(p, method = "BH")
    }
  )
  dt$q[ok] <- qv
}

data.table::setorder(dt, q, -effect, na.last = TRUE)
out_results <- file.path(event_dir, "screen_results.csv.gz")
data.table::fwrite(dt, out_results)

hits <- dt[is.finite(q) & q < q_threshold & is.finite(effect) & effect > 0, ID]
hits_file <- file.path(event_dir, sprintf("tier1_hits_%s.txt", tf))
writeLines(as.character(hits), hits_file)

message(sprintf(
  "q<%.3g & effect>0 → %d hits / %d tested (%.1f%%). Wrote %s + %s",
  q_threshold, length(hits), sum(ok), 100 * length(hits) / max(1L, sum(ok)),
  out_results, hits_file
))
