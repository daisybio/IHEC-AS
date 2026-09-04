## ---------------------------------------------------------------------------
## Tier-1 screen aggregator
##
## Collects the per-event screen rows written by 09s-ridge-screen.R, computes
## FDR q-values across all events (Storey qvalue → BH fallback), and writes:
##   processed_data/event_models/<tf>/screen_results.csv.gz   (all events + q)
##   processed_data/event_models/<tf>/tier1_hits_<tf>.txt      (hit event IDs)
##
## A "hit" = q < q_threshold. The one-sided exceedance p already encodes direction, so
## the former `AND effect > 0` clause was redundant; it was also outlier-dominated on the
## unbounded R2 (61% of real rows positive for numerical reasons), so it has been dropped
## (2026-08-01, see the note next to `min_effect` below). An event is a hit when it beats
## its matched feature-rotation null — event-specific epigenetic signal rather than a
## cell-type fingerprint. The Tier-2 elastic-net (09zz via 09-1) runs on this hit list.
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
# For SCREEN_STAT_VERSION — the provenance stamp the per-event outputs carry. The
# aggregator must know the current value to refuse stale rows (see the stamp gates
# below); it is defined next to the statistic itself in 09-ml-shared.R.
source("09-ml-shared.R")
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
# Parallel read: at full scale this is ~34k small gzipped files per feature set and
# a sequential lapply() over them dominated the whole rule (~45 min wall of a 45-min
# run). Same pattern as the null-sidecar read below. Falls back to lapply if
# pbmcapply is unavailable.
.read_one <- function(f) tryCatch(data.table::fread(f), error = function(e) NULL)
.cores_read <- max(1L, min(8L, data.table::getDTthreads()))
.file_parts <- if (requireNamespace("pbmcapply", quietly = TRUE)) {
  pbmcapply::pbmclapply(files, .read_one, mc.cores = .cores_read)
} else {
  lapply(files, .read_one)
}
.bad <- sum(!vapply(.file_parts, is.data.frame, logical(1)))
if (.bad > 0L) {
  warning(sprintf("%d of %d per-event screen files were unreadable and skipped",
                  .bad, length(files)), immediate. = TRUE)
}
dt <- data.table::rbindlist(
  .file_parts[vapply(.file_parts, is.data.frame, logical(1))], fill = TRUE
)
rm(.file_parts)

# LONG-REUSE (fallback): the previous omnibus screen (single ridge over ALL
# epigenetic features) IS the `long` feature set. 09s-ridge-screen.R's per-feature-set
# resume gate normally retains those rows in the per-event files themselves (tagging
# untagged legacy rows as "long"), so this hook only fires if the per-event long rows
# are gone — e.g. the screen dir was cleared but screen_long_results.csv.gz was kept.
# STAMP GATE: only a file produced by the CURRENT statistic may be folded back in.
# An unstamped file is version 1 (the pre-2026-07-29 leaky statistic) and must be
# refused -- this hook fires precisely when the per-event `long` rows are absent,
# i.e. right after someone deletes them to force a recompute, so without this check
# it would silently re-inject exactly the results they were trying to discard.
long_reuse <- file.path(event_dir, "screen_long_results.csv.gz")
have_long <- "feature_set" %in% names(dt) && any(dt$feature_set == "long")
if (file.exists(long_reuse) && !have_long) {
  .lv <- tryCatch(data.table::fread(long_reuse, nrows = 1L), error = function(e) NULL)
  .lver <- if (!is.null(.lv) && "stat_version" %in% names(.lv)) {
    suppressWarnings(as.integer(.lv$stat_version[1L]))
  } else {
    1L
  }
  if (is.na(.lver) || .lver != SCREEN_STAT_VERSION) {
    warning(sprintf(
      "IGNORING %s: stat_version %s != current %d (stale statistic). Delete it.",
      basename(long_reuse), if (is.na(.lver)) "NA/absent" else .lver,
      SCREEN_STAT_VERSION
    ), immediate. = TRUE)
    long_reuse <- ""
  }
}
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

# --- STAMP GATE: never mix statistic versions inside one FDR family ------------
# Pooling rows computed by different versions of the statistic into one FDR family
# is exactly the silent corruption the stamp exists to prevent: q-values would be
# computed across incomparable numbers. Unstamped rows are version 1 (pre-2026-07-29
# leaky statistic). Stale rows are DROPPED, not tolerated -- a partial screen is
# recoverable (re-run the missing events), a quietly-wrong q-value is not.
if (!"stat_version" %in% names(dt)) dt[, stat_version := NA_integer_]
dt[, stat_version := suppressWarnings(as.integer(stat_version))]
dt[is.na(stat_version), stat_version := 1L]
.vers <- sort(unique(dt$stat_version))
if (length(.vers) > 1L || .vers[1L] != SCREEN_STAT_VERSION) {
  .n_stale <- sum(dt$stat_version != SCREEN_STAT_VERSION)
  warning(sprintf(
    paste("screen rows carry stat_version %s but the current statistic is %d --",
          "DROPPING %d stale row(s). Re-run the screen for those events",
          "(deleting their per-event files is no longer necessary: the worker's",
          "resume gate now discards mismatched rows by itself)."),
    paste(.vers, collapse = "/"), SCREEN_STAT_VERSION, .n_stale
  ), immediate. = TRUE)
  dt <- dt[stat_version == SCREEN_STAT_VERSION]
  if (!nrow(dt)) {
    stop("No screen rows match the current stat_version (", SCREEN_STAT_VERSION,
         ") -- the screen must be re-run before it can be aggregated.")
  }
}

# Feature-set / event-type keys for grouped FDR.
#
# These must be filled ROW-WISE, not column-wise. The screen dir is routinely in a
# MIXED state — most per-event files predate the per-feature-set screen (no
# feature_set / Event Type / transcript_filter columns at all) while some carry all
# three. `rbindlist(fill = TRUE)` then creates the columns and leaves NA in the
# legacy rows, so a `if (!"feature_set" %in% names(dt))` guard sees the column
# present and fills nothing: every legacy row lands in one NA FDR family and the
# per-family FDR silently degenerates to the old ungrouped behaviour.
# (Real occurrence: 34,143 of 34,148 rows were NA on all three keys.)
if (!"feature_set" %in% names(dt)) dt[, feature_set := NA_character_]
if (!"transcript_filter" %in% names(dt)) dt[, transcript_filter := NA_character_]
if (!"Event Type" %in% names(dt)) dt[, `Event Type` := NA_character_]
dt[, feature_set := as.character(feature_set)]
dt[is.na(feature_set), feature_set := "long"] # untagged == the omnibus/long screen
dt[is.na(transcript_filter), transcript_filter := tf]
if (anyNA(dt$`Event Type`)) {
  .sess <- readRDS(cfg$session_rds)
  .et <- .sess$event_dt[, .(ID = as.integer(ID), k_et = `Event Type`)]
  rm(.sess)
  dt[, ID := as.integer(ID)]
  dt[.et, on = "ID", k_et := i.k_et]
  dt[is.na(`Event Type`), `Event Type` := k_et]
  dt[, k_et := NULL]
  if (anyNA(dt$`Event Type`)) {
    warning(sprintf(
      "%d rows still have no Event Type after the session join — they form their own FDR family",
      sum(is.na(dt$`Event Type`))
    ), immediate. = TRUE)
  }
}
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

## ---------------------------------------------------------------------------
## SENSITIVITY p: floor-free, family-stratified pooled z
##
## WHY. The exceedance p above is floored at 1/(R_used + 1) (= 0.00498 at R=200).
## BH at q<0.1 only passes a p if p <= k*q/m at its rank k, and measured on the
## real screen that is 6.6e-4 (SE, m=32368, k=215 at the floor) and 2.8e-4 (RI) --
## i.e. NO event can pass FDR at 200 rotations however strong its signal, so a
## "0 hits" result from the primary statistic alone is not interpretable. More
## rotations cannot fix it (the matched-control pool is smaller than the ~1500
## needed, and for RI smaller than the RI event count itself), and dropping
## low-variance events does not either (the BH threshold is scale-invariant).
##
## WHAT. z = (stat - mean(own matched null)) / sd(own matched null); the p is then
## read off a reference of leave-one-out-standardised NULL z. Only the
## dimensionless deviation SHAPE is borrowed across events -- each event's own R2
## level and scale stay local, so a cell-type-fingerprint event (high null MEAN)
## still gets a small z and is not rescued by this.
##
## CRITICAL: the reference is pooled WITHIN the FDR family
## (transcript_filter, Event Type, feature_set) -- the same grouping as q above.
## Pooling across Event Types is measurably wrong: leave-one-EVENT-out
## calibration gives obs/nominal 1.43 for SE (anti-conservative) and 0.41 for RI
## (conservative) at alpha=1e-4, which average to a deceptive 1.01 on a 95%-SE
## sample. Within-family it is 1.00 (SE) / 0.91-1.05 (RI). The tails genuinely
## differ (KS SE-z vs RI-z p=0.0072). Stratifying is free -- each family's own
## nulls give 3-4 orders more resolution than its own BH target.
## Do NOT substitute a median/MAD scale: it is 2.6-5.2x anti-conservative.
##
## This is a SENSITIVITY analysis. `tier1_hits` and `sig` stay on the primary
## exceedance p, which assumes only that an event's own controls are exchangeable
## with it; this one additionally assumes a family-common deviation shape.
## Validation: verification/check_pooled_z_strata.R
# Smallest usable reference (null z values left after removing the event's own).
# 10,000 resolves p down to ~1e-4, below the BH bound of every real family here.
min_ref_z <- getOption("EpiATLAS_AS_SCREEN_MIN_REF_Z", 10000L)
screen_dir_nulls <- list.files(
  screen_dir, pattern = "_screen_null\\.csv\\.gz$", full.names = TRUE
)
dt[, `:=`(z_screen = NA_real_, p_pooled_z = NA_real_, q_pooled_z = NA_real_)]
if (length(screen_dir_nulls) == 0L) {
  warning("No _screen_null sidecars found -- skipping the floor-free sensitivity p. ",
          "These sidecars are REQUIRED for it (they hold the per-control null ",
          "vectors); do not delete them after aggregation.", immediate. = TRUE)
} else {
  .read_null <- function(f) {
    d <- tryCatch(data.table::fread(f), error = function(e) NULL)
    if (is.null(d) || !nrow(d)) return(NULL)
    if (!"feature_set" %in% names(d)) d[, feature_set := "long"]
    # Stamp gate, same rule as the summary rows: a null vector from an older
    # statistic cannot be part of the reference distribution (unstamped == v1).
    .v <- if ("stat_version" %in% names(d)) {
      suppressWarnings(as.integer(d$stat_version))
    } else {
      rep(1L, nrow(d))
    }
    d <- d[!is.na(.v) & .v == SCREEN_STAT_VERSION]
    if (!nrow(d)) return(NULL)
    d <- d[is.finite(null_R2), .(ID = as.integer(ID),
                                 feature_set = as.character(feature_set), null_R2)]
    if (!nrow(d)) NULL else d
  }
  .cores <- max(1L, min(8L, data.table::getDTthreads()))
  .parts <- if (requireNamespace("pbmcapply", quietly = TRUE)) {
    pbmcapply::pbmclapply(screen_dir_nulls, .read_null, mc.cores = .cores)
  } else {
    lapply(screen_dir_nulls, .read_null)
  }
  nulls <- data.table::rbindlist(
    .parts[vapply(.parts, is.data.frame, logical(1))], fill = TRUE
  )
  rm(.parts)

  # Leave-one-out standardisation of each event's own null vector (closed form,
  # identical to the validated script). n>=3 needed for a LOO variance.
  nulls <- nulls[, if (.N >= 3L) {
    x <- null_R2
    n <- .N
    s <- sum(x)
    ss <- sum(x^2)
    m_i <- (s - x) / (n - 1)
    v_i <- (ss - x^2 - (n - 1) * m_i^2) / (n - 2)
    .(z = (x - m_i) / sqrt(pmax(v_i, .Machine$double.eps)))
  }, by = .(ID, feature_set)]
  nulls <- nulls[is.finite(z)]

  # Attach the family keys the reference is pooled within. NB `i.` cannot prefix a
  # backticked/spaced column name, so `Event Type` is aliased in the join table
  # (same gotcha as the long-reuse hook above).
  keys <- unique(dt[, .(ID = as.integer(ID), feature_set = as.character(feature_set),
                        k_tf = transcript_filter, k_et = `Event Type`)])
  nulls[keys, on = c("ID", "feature_set"),
        `:=`(transcript_filter = i.k_tf, `Event Type` = i.k_et)]
  nulls <- nulls[!is.na(`Event Type`)]

  dt[, z_screen := (screen_R2 - null_R2_mean) / null_R2_sd]
  dt[!is.finite(null_R2_sd) | null_R2_sd <= 0, z_screen := NA_real_]

  # Per family: p = (1 + #{reference z >= z_real}) / (1 + N_ref), with the event's
  # OWN null z removed from the reference (leave-one-event-out, matching how the
  # method was validated -- a value must never contribute to what scores it).
  fam <- unique(nulls[, .(transcript_filter, `Event Type`, feature_set)])
  data.table::setkey(nulls, transcript_filter, `Event Type`, feature_set, ID)
  for (i in seq_len(nrow(fam))) {
    f <- fam[i]
    ref <- nulls[.(f$transcript_filter, f$`Event Type`, f$feature_set)]
    if (!nrow(ref) || anyNA(ref$z)) ref <- ref[!is.na(z)]
    if (!nrow(ref)) next
    ref_sorted <- sort(ref$z)
    n_ref <- length(ref_sorted)
    idx <- dt[transcript_filter == f$transcript_filter &
                `Event Type` == f$`Event Type` &
                as.character(feature_set) == f$feature_set &
                is.finite(z_screen), which = TRUE]
    if (!length(idx)) next
    zq <- dt$z_screen[idx]
    ge_all <- n_ref - findInterval(zq - 1e-12, ref_sorted)
    # Own-event contribution, computed with a keyed grouped join rather than a
    # per-row lookup: `split()` + `[[<character>]]` inside a loop is an O(n) name
    # search per row, i.e. quadratic in the number of events (~1e9 ops and >23 min
    # at SE scale before this was vectorised).
    q_dt <- data.table::data.table(
      ID = as.integer(dt$ID[idx]), zq = zq, ord = seq_along(idx)
    )
    own <- ref[q_dt, on = "ID", by = .EACHI, .(
      ord = i.ord, own_ge = sum(z >= i.zq), own_tot = .N
    )]
    own_ge <- integer(length(idx))
    own_tot <- integer(length(idx))
    .ok <- !is.na(own$ord) & !is.na(own$own_ge)
    own_ge[own$ord[.ok]] <- own$own_ge[.ok]
    own_tot[own$ord[.ok]] <- own$own_tot[.ok]
    # MINIMUM REFERENCE SIZE. After removing the event's own nulls the reference must
    # still be large enough to resolve the alpha this family's FDR needs. A family
    # containing a single screened event leaves an EMPTY reference and would return
    # p = (1+0)/(1+0) = 1 for it — a meaningless value that looks like a real result
    # (seen for real: the one 3-set-tagged event in an otherwise legacy screen dir).
    # Require the floor 1/(n_ref_eff+1) to be below the family's own BH bound.
    n_ref_eff <- n_ref - own_tot
    too_small <- n_ref_eff < min_ref_z
    p_vals <- (1 + (ge_all - own_ge)) / (1 + n_ref_eff)
    p_vals[too_small] <- NA_real_
    if (any(too_small)) {
      message(sprintf(
        "  %s/%s/%s: %d of %d rows have an effective reference < %d null z values -> sensitivity p = NA",
        f$transcript_filter, f$`Event Type`, f$feature_set,
        sum(too_small), length(idx), min_ref_z
      ))
    }
    data.table::set(dt, i = idx, j = "p_pooled_z", value = p_vals)
  }
  dt[
    is.finite(p_pooled_z),
    q_pooled_z := qfun(p_pooled_z),
    by = .(transcript_filter, `Event Type`, feature_set)
  ]
  message(sprintf(
    "Sensitivity p: %d rows scored against family-stratified pooled-z references (%d null z values)",
    sum(is.finite(dt$p_pooled_z)), nrow(nulls)
  ))
  rm(nulls)
  gc()
}

data.table::setorder(dt, `Event Type`, feature_set, q, -effect, na.last = TRUE)
out_results <- file.path(event_dir, "screen_results.csv.gz")
data.table::fwrite(dt, out_results)

# Minimum effect size for a hit. Default 0 = the historical criterion, unchanged.
# Raise it to suppress variance-collapse artifacts: events whose z is large only
# because their null sd is proportionally tiny. Real numbers from the `long` screen
# show this is the right instrument and a null-sd cutoff is NOT -- null_R2_sd is
# never 0 (min 1.0e-5), while the top-z events have effects of only 0.005-0.013
# (e.g. ID 50526: screen_R2 0.0052, null_sd 3.0e-4, z 17.0). effect > 0.01 drops
# max z from 17.0 to 12.9; effect > 0.05 drops it to 6.7.
min_effect <- getOption("EpiATLAS_AS_SCREEN_MIN_EFFECT", 0)

# THE `effect > 0` CLAUSE HAS BEEN DROPPED FROM THE HIT RULE (2026-08-01).
# Two reasons:
#   1. It is redundant. p_emp is a ONE-SIDED exceedance probability,
#      (1 + #{null >= real})/(R + 1), so it already encodes direction: an event cannot
#      reach a small p while sitting below its own null distribution.
#   2. As computed on the raw statistic it was actively misleading. `effect` is
#      mean-based on an unbounded R2, so it is outlier-dominated: measured across 678
#      real event x feature-set combos, 61% of effects were positive for purely
#      numerical reasons, and the median effect flipped sign by |R2| stratum
#      (+0.31 / +0.93 / -11.9 / -222.6). A clause that passes ~3/4 of events for
#      arithmetic reasons is not a filter.
# `min_effect` survives as an OPT-IN instrument for a DIFFERENT job -- suppressing
# variance-collapse artifacts (see the note above) -- and now acts on `effect_bounded`,
# a difference of bounded quantities, rather than on the corrupted raw `effect`.
# Default 0 = no effect filter at all (NOT "effect > 0", which is what it used to mean).
use_effect_filter <- is.finite(min_effect) && min_effect > 0
if (use_effect_filter && !"effect_bounded" %in% names(dt)) {
  stop(
    "EpiATLAS_AS_SCREEN_MIN_EFFECT > 0 but `effect_bounded` is absent -- these screen ",
    "rows predate stat_version 3. Re-screen, or set the option back to 0."
  )
}

# Hit = significant in ANY feature set (Tier-2 then decomposes long/short/local).
# PRIMARY criterion only -- the sensitivity p never feeds the Tier-2 hit list.
hit_rows <- if (use_effect_filter) {
  dt[is.finite(q) & q < q_threshold &
    is.finite(effect_bounded) & effect_bounded > min_effect]
} else {
  dt[is.finite(q) & q < q_threshold]
}
hits <- sort(unique(hit_rows$ID))
hits_file <- file.path(event_dir, sprintf("tier1_hits_%s.txt", tf))
writeLines(as.character(hits), hits_file)

message(sprintf(
  "q<%.3g%s -> %d unique hit events (%d feature-set rows). Wrote %s + %s",
  q_threshold,
  if (use_effect_filter) sprintf(" & effect_bounded>%g", min_effect) else "",
  length(hits), nrow(hit_rows), out_results, hits_file
))
# Per-feature-set breakdown, primary vs sensitivity side by side.
print(dt[is.finite(q), .(
  n_tested = .N,
  n_hits = sum(
    q < q_threshold & (!use_effect_filter | effect_bounded > min_effect),
    na.rm = TRUE
  ),
  n_hits_pooled_z = sum(
    q_pooled_z < q_threshold & (!use_effect_filter | effect_bounded > min_effect),
    na.rm = TRUE
  ),
  min_p_emp = signif(min(p_emp, na.rm = TRUE), 3),
  min_p_pooled_z = if (any(is.finite(p_pooled_z))) {
    signif(min(p_pooled_z, na.rm = TRUE), 3)
  } else {
    NA_real_
  }
), by = .(`Event Type`, feature_set)][order(`Event Type`, feature_set)])

# --- p-FLOOR arithmetic: can the primary statistic pass FDR AT ALL? -----------
# The exceedance p cannot go below 1/(R_used+1). BH passes a p at rank k only if
# p <= k*q/m. If the floor exceeds that bound the family is UNRESOLVABLE by the
# primary statistic and a 0-hit result there says nothing on its own -- which is
# exactly why the sensitivity p above exists. Printed so this is never implicit.
floor_tab <- dt[is.finite(p_emp), {
  m <- .N
  fl <- 1 / (max(R_used, na.rm = TRUE) + 1)
  k <- sum(
    p_emp <= fl + 1e-12 & (!use_effect_filter | effect_bounded > min_effect),
    na.rm = TRUE
  )
  .(m_tests = m, R_used_max = max(R_used, na.rm = TRUE),
    p_floor = signif(fl, 3), k_at_floor = k,
    bh_allows_at_k = if (k > 0L) signif(k * q_threshold / m, 3) else NA_real_,
    floor_can_pass = k > 0L && fl <= k * q_threshold / m,
    k_needed_at_floor = ceiling(fl * m / q_threshold))
}, by = .(transcript_filter, `Event Type`, feature_set)]
message("\np-floor reachability of the PRIMARY exceedance p (per FDR family):")
print(floor_tab[order(`Event Type`, feature_set)])
if (any(!floor_tab$floor_can_pass)) {
  message(
    "NOTE: in the families above with floor_can_pass = FALSE, NO event could pass ",
    "FDR via the exceedance p regardless of signal strength (its floor is below ",
    "what BH allows). Read the pooled-z sensitivity columns before concluding ",
    "anything from a 0-hit result there."
  )
}
