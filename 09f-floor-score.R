#!/usr/bin/env Rscript
# Detection-floor scorer -- `rule floor_score` in the Snakefile.
# Builds the pooled-z reference, self-validates it, then scores every injected event.
# See revision/EXECUTION-PLAN.md item 1.
#
# WHY THIS EXISTS SEPARATELY. Injected events must be scored against the REAL screen's pooled-z
# reference, not one built from the injected sample's own nulls. Running 09s-aggregate.R on an
# injection directory would pool ~600 events x 200 controls = ~120k null z values, resolving p to
# ~8e-6; the real screen's reference holds 24,092,316 and resolves to ~4e-8. Three orders. It is
# also the correct QUESTION -- "would this effect have been detected in the actual screen?" -- not
# merely the cheaper one.
#
# This script builds and validates that reference NOW, before any injected data exists, by
# reproducing 09s-aggregate.R's construction and checking it reproduces the real screen's own
# p_pooled_z. If it does, the scoring step for the floor is trustworthy.
#
# Construction, mirroring 09s-aggregate.R exactly:
#   1. read every screen/*_screen_null.csv.gz  (the REAL screen only -- list.files is
#      non-recursive and screen_dir is pinned, so injected sidecars can never leak in)
#   2. per (ID, feature_set) with >= 3 controls, LEAVE-ONE-OUT standardise null_R2 -> z
#   3. pool z within the FDR family (Event Type x feature_set)
#   4. score any statistic as p = (1 + #{ref z >= z}) / (1 + N_ref)
#
# Reads production read-only. Writes only to .claude/scratch/.

suppressPackageStartupMessages(library(data.table))
setDTthreads(4L)

tf <- Sys.getenv("TRANSCRIPT_FILTER", "biotype_filtered")
event_dir <- file.path("processed_data", "event_models", tf)
screen_dir <- file.path(event_dir, "screen") # PINNED -- never point this at an injection dir
N_VALIDATE <- as.integer(Sys.getenv("N_VALIDATE", "400"))

scr <- fread(file.path(event_dir, "screen_results.csv.gz"))
keys <- unique(scr[, .(ID, feature_set, `Event Type`)])

nulls_files <- list.files(screen_dir, pattern = "_screen_null\\.csv\\.gz$", full.names = TRUE)
cat(sprintf("null sidecars: %s\n", format(length(nulls_files), big.mark = ",")))

.read_null <- function(f) {
  d <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(d) || !nrow(d)) return(NULL)
  # stat_version guard: an unstamped or mismatched sidecar is a DIFFERENT statistic and must
  # not enter the reference (unstamped == v1, per 09s-aggregate.R).
  if ("stat_version" %in% names(d)) d <- d[stat_version == max(scr$stat_version)]
  if (!nrow(d)) return(NULL)
  d[is.finite(null_R2), .(ID = as.integer(ID), feature_set = as.character(feature_set), null_R2)]
}

cores <- max(1L, min(8L, parallel::detectCores() - 1L))
parts <- if (requireNamespace("pbmcapply", quietly = TRUE)) {
  pbmcapply::pbmclapply(nulls_files, .read_null, mc.cores = cores)
} else {
  parallel::mclapply(nulls_files, .read_null, mc.cores = cores)
}
nulls <- rbindlist(Filter(function(x) is.data.table(x) && nrow(x), parts), fill = TRUE)
cat(sprintf("raw null R2 values: %s\n", format(nrow(nulls), big.mark = ",")))

# leave-one-out standardisation within (ID, feature_set) -- a control must not contribute to the
# mean/sd it is then measured against
nulls <- nulls[, if (.N >= 3L) {
  x <- null_R2
  n <- .N
  s <- sum(x)
  ss <- sum(x^2)
  mu_i <- (s - x) / (n - 1L)
  var_i <- (ss - x^2 - (n - 1L) * mu_i^2) / (n - 2L)
  # MUST be .Machine$double.eps, not 0. 09s-aggregate.R floors the LOO variance at
  # double.eps, so an event whose LOO variance goes non-positive still yields a huge but
  # FINITE z that stays in the reference. Flooring at 0 gives z = Inf, which is_finite()
  # then drops -- costing 367,104 values (1.5% of the reference) and biasing every p by
  # ~1.5% relatively, i.e. the ~2.4e-5 median discrepancy seen on the first two passes.
  .(z = (x - mu_i) / sqrt(pmax(var_i, .Machine$double.eps)))
} else {
  .(z = numeric(0))
}, by = .(ID, feature_set)]
nulls <- nulls[is.finite(z)]
# `i.` cannot prefix a backticked/spaced column name -- alias in the join table's j.
# CLAUDE.md documents this; it is its FIFTH occurrence in this repo.
nulls[keys[, .(ID, feature_set, k_et = `Event Type`)], on = c("ID", "feature_set"),
  `Event Type` := i.k_et]
nulls <- nulls[!is.na(`Event Type`)]

cat(sprintf("usable reference z values: %s\n\n", format(nrow(nulls), big.mark = ",")))
cat("=== reference size and resolution per FDR family ===\n")
ref_sz <- nulls[, .(n_ref = .N, floor_p = signif(1 / (.N + 1), 3),
  q999 = round(quantile(z, 0.999), 2), max_z = round(max(z), 1)),
  by = .(`Event Type`, feature_set)][order(`Event Type`, feature_set)]
print(ref_sz)

# ---- VALIDATION: reproduce the real screen's own p_pooled_z from this reference ----------
cat("\n=== validation: recompute p_pooled_z for a sample of REAL events ===\n")
set.seed(1L)
real <- scr[is.finite(z_screen) & is.finite(p_pooled_z)]
samp <- real[sample.int(.N, min(N_VALIDATE, .N))]
# LEAVE-ONE-EVENT-OUT. 09s-aggregate.R removes the event's OWN null z values from the
# reference before scoring it -- "a value must never contribute to what scores it". Omitting
# this was the entire validation failure on the first pass (median |diff| 2.6e-5, which is
# ~200/23.7M, exactly the size of one event's own contribution).
#
# NB for INJECTED events this term is zero: their nulls live in the injection directory and
# are not part of the real reference at all. So floor scoring is the simpler case -- but the
# validation must reproduce the real-event path to prove the reference itself is correct.
samp[, p_recomputed := {
  ref <- nulls[`Event Type` == .BY$`Event Type` & feature_set == .BY$feature_set,
    .(ID, z)]
  zs <- sort(ref$z)
  n_ref <- length(zs)
  ids <- ID
  vapply(seq_along(z_screen), function(k) {
    zz <- z_screen[k]
    ge_all <- n_ref - findInterval(zz - 1e-12, zs)
    own <- ref[ID == ids[k], z]
    (1 + ge_all - sum(own >= zz)) / (1 + n_ref - length(own))
  }, numeric(1))
}, by = .(`Event Type`, feature_set)]
samp[, d := abs(p_recomputed - p_pooled_z)]
cat(sprintf("  n=%d  max |diff| = %.3g  median |diff| = %.3g  identical(<1e-12): %.1f%%\n",
  nrow(samp), max(samp$d), median(samp$d), 100 * mean(samp$d < 1e-12)))
if (max(samp$d) > 1e-6) {
  cat("  !! reference does NOT reproduce the aggregator -- do not use it to score the floor\n")
  print(head(samp[order(-d), .(ID, `Event Type`, feature_set, z_screen, p_pooled_z, p_recomputed, d)], 5))
} else {
  cat("  reference reproduces the aggregator -- safe to score injected events against it\n")
}

if (max(samp$d) > 1e-6) stop("reference does not reproduce the aggregator -- refusing to score")

# ============================ SCORE THE INJECTED EVENTS ============================
floor_root <- file.path(event_dir, "floor")
inj_files <- list.files(floor_root, pattern = "_screen\\.csv\\.gz$",
  full.names = TRUE, recursive = TRUE)
cat(sprintf("\n=== injected rows: %d files under %s ===\n", length(inj_files), floor_root))
if (!length(inj_files)) stop("no injected outputs found under ", floor_root)

inj <- rbindlist(lapply(inj_files, function(f) {
  d <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(d) || !nrow(d)) return(NULL)
  d
}), fill = TRUE)
cat(sprintf("rows: %s over %d events x %d rho values\n",
  format(nrow(inj), big.mark = ","), uniqueN(inj$ID), uniqueN(inj$inject_rho)))

# z from the injected row's OWN null moments. NB no leave-one-event-out term here: an
# injected event's nulls live under floor/, never in the real reference, so its own
# contribution to that reference is exactly zero.
inj[, z_inj := (screen_R2 - null_R2_mean) / null_R2_sd]
inj[!is.finite(null_R2_sd) | null_R2_sd <= 0, z_inj := NA_real_]

inj[, p_pooled_z := {
  ref <- sort(nulls[`Event Type` == .BY$`Event Type` & feature_set == .BY$feature_set, z])
  n_ref <- length(ref)
  fifelse(is.na(z_inj), NA_real_,
    (1 + n_ref - findInterval(z_inj - 1e-12, ref)) / (1 + n_ref))
}, by = .(`Event Type`, feature_set)]

# DETECTION. BH admits a p at rank k iff p <= k*q/m. Splice the injected p into the REAL
# family's p-vector: its rank is 1 + #{real p < p_inj}, and m is the real family size.
# Exact up to the one real p this event would have displaced (shifts rank by <= 1).
QTHR <- 0.1
# Largest rho whose response still carries noise in a meaningful amount. rho enters as
# sqrt(rho)*signal + sqrt(1-rho)*noise, so rho = 1 is noiseless and destabilises the
# matched null (see the long comment at the anchors block). Anchors and the power curve
# are restricted to rho <= this; rho = 1 stays in the sweep as a diagnostic.
RHO_MAX_VALID <- 0.4
famp <- scr[is.finite(p_pooled_z), .(p = sort(p_pooled_z), m = .N),
  by = .(`Event Type`, feature_set)]
inj[, detected := {
  fp <- famp[`Event Type` == .BY$`Event Type` & feature_set == .BY$feature_set]
  ps <- fp$p
  m <- fp$m[1]
  k <- findInterval(p_pooled_z, ps) + 1L
  !is.na(p_pooled_z) & p_pooled_z <= k * QTHR / m
}, by = .(`Event Type`, feature_set)]

ev <- fread(sprintf("processed_data/event_annotations_dt_%s.csv.gz", tf),
  select = c("ID", "Variability"))
inj[ev[, .(ID, k_vb = Variability)], on = "ID", Variability := i.k_vb]

cat("\n=== FLOOR: detection rate vs achieved out-of-fold R2 ===\n")
cat("rho is a knob, NOT the reported quantity -- the rho->R2 map is nonlinear\n")
cat("(at rho=1: gaussian 1.000, logit 0.809, rankmap 0.518).\n\n")
print(inj[feature_set == "local", .(
  n = .N,
  med_R2 = round(median(screen_R2, na.rm = TRUE), 4),
  det_rate_pct = round(100 * mean(detected, na.rm = TRUE), 1)
), by = .(`Event Type`, Variability, inject_rho)][order(`Event Type`, Variability, inject_rho)])

# Blow-up status must come from the REAL screen (floor_events.tsv), NOT be recomputed
# from the injected R2. Recomputing it makes the strata drift across rho -- a strong
# injection stabilises a fit that blows up on real PSI, so the event silently changes
# stratum and the per-rho rows stop being comparable. Measured 2026-09-01: the post-hoc
# flag gave n = 184/177/181/176/183/184/180/228 across the rho grid for one stratum
# (it should be constant), and moved RI's rho=1 detection rate from 10.7% to 7.0%.
.ev_file <- Sys.getenv("FLOOR_EVENTS",
  file.path(event_dir, "floor", "floor_events.tsv"))
if (file.exists(.ev_file)) {
  .ev <- fread(.ev_file, select = c("ID", "blowup"))
  inj[.ev[, .(ID, k_bl = blowup)], on = "ID", blowup := i.k_bl]
} else {
  warning("floor_events.tsv not found; falling back to post-hoc blow-up status, ",
    "whose strata are NOT comparable across rho")
  inj[, blowup := abs(screen_R2) > 1]
}

cat("\n=== same, stratified by REAL blow-up status (power is NOT uniform) ===\n")
cat("blowup is the event's status in the REAL screen, fixed across rho.\n")
cat("Blown-up events are unrecoverable at ANY effect size -- 0% at every rho.\n\n")
print(inj[feature_set == "local", .(n = .N,
  med_R2  = round(median(screen_R2, na.rm = TRUE), 4),
  med_nsd = round(median(null_R2_sd, na.rm = TRUE), 5),
  med_z   = round(median(z_inj, na.rm = TRUE), 2),
  det_rate_pct = round(100 * mean(detected, na.rm = TRUE), 1)
), by = .(`Event Type`, blowup, inject_rho)][order(`Event Type`, blowup, inject_rho)])

cat("\n=== ANCHORS ===\n")
print(inj[feature_set == "local" & inject_rho %in% c(0, RHO_MAX_VALID), .(
  n = .N, med_R2 = round(median(screen_R2, na.rm = TRUE), 4),
  med_p = signif(median(p_pooled_z, na.rm = TRUE), 3),
  det_rate_pct = round(100 * mean(detected, na.rm = TRUE), 1)
), by = .(inject_rho)][order(inject_rho)])
cat("rho=0 must give a calibrated null: detection at or below the nominal rate.\n")
cat(sprintf(
  "Detection must rise MONOTONICALLY over the noisy range rho <= %g.\n", RHO_MAX_VALID))
cat("There is deliberately NO 'rho=1 must be detected' anchor -- see below.\n")

# WHY rho = 1 IS EXCLUDED (measured 2026-09-01, do not re-add it as an anchor).
# The worker builds the latent response as sqrt(rho)*signal + sqrt(1-rho)*noise, so at
# rho = 1 the noise coefficient is exactly 0: the target is a noiseless, deterministic
# function of s=10 columns. Kernel ridge on the CONTROL features fails against such a
# target wildly and erratically, which inflates the null's SPREAD without shifting its
# mean (null_R2_mean stays ~0; frac(null_R2_mean > 0) stays flat at 0.33-0.43, so the
# controls do not predict better -- they fail more variably). Since
# z = (R2 - null_mean)/null_sd, that sinks z even as R2 rises:
#   RI  rho 0.4 -> 1:  med R2 0.072 -> 0.718 (x10)  but med null_sd 0.056 -> 1.027 (x18)
#                      => med z 0.86 -> 0.72, detection 20.7% -> 10.7%
#   SE  rho 0.4 -> 1:  med null_sd 0.009 -> 0.238 (x26)
# rho = 1 is therefore a DEGENERATE setting, not a strong-signal setting, and its low
# detection rate says nothing about the screen's sensitivity. Keep it in the sweep as a
# diagnostic; never treat it as the ceiling.
cat("\n=== POWER CURVE: detection rate vs achieved R2 (clean stratum, valid rho) ===\n")
.pc <- inj[feature_set == "local" & blowup %in% FALSE & inject_rho <= RHO_MAX_VALID,
  .(n = .N,
    med_R2 = round(median(screen_R2, na.rm = TRUE), 4),
    power_pct = round(100 * mean(detected, na.rm = TRUE), 1)),
  by = .(`Event Type`, inject_rho)][order(`Event Type`, inject_rho)]
print(.pc)

# Report POWER at a stated effect size, not the smallest R2 that happened to be detected.
# min(screen_R2 | detected) is an outlier statistic and reads as a floor when it is not:
# detection is decided on z, so an event whose null variance has collapsed is detected at
# a trivial R2 while the median event at 10x that R2 is not. Quoting it as "the floor"
# would overstate the screen's sensitivity by an order of magnitude.
cat("\n=== POWER AT A TARGET: is there any effect size with >= 50% power? ===\n")
POWER_TARGET <- 0.5
.pk <- .pc[, .(
  peak_power_pct = max(power_pct),
  med_R2_at_peak = med_R2[which.max(power_pct)],
  rho_at_peak    = inject_rho[which.max(power_pct)],
  reaches_target = any(power_pct >= 100 * POWER_TARGET)
), by = .(`Event Type`)]
print(.pk)
if (!any(.pk$reaches_target)) {
  cat(sprintf(paste0(
    "\nNO effect size reaches %.0f%% power in any event type. The screen therefore has\n",
    "no detection FLOOR to quote -- report the power curve above instead, and state that\n",
    "blown-up events (~50%% of the real screen) have 0%% power at EVERY effect size.\n"),
    100 * POWER_TARGET))
}

cat("\n=== for reference only: smallest achieved R2 that happened to be detected ===\n")
cat("NOT a floor -- an outlier statistic. See the comment above. Do not quote it.\n")
print(inj[detected == TRUE & feature_set == "local",
  .(min_detected_R2 = round(min(screen_R2, na.rm = TRUE), 4),
    n_detected = .N), by = .(`Event Type`, Variability)][order(`Event Type`, Variability)])

f <- Sys.getenv("FLOOR_SCORED",
  file.path(event_dir, "floor", "floor_results.csv.gz"))
dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
fwrite(inj, f, compress = "gzip")
cat(sprintf("\nwritten: %s\n", f))
