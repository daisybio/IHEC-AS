# Paired-epirr protocol effect validation (reviewer §4.12g, §4.13h).
# See revision/file-changes/06c-paired-protocol.R.md for the full spec.
#
# Uses the epirrs with BOTH mRNA-Seq and total-RNA-Seq UUIDs to estimate the
# true within-sample protocol effect on PSI via a paired design, and validates
# the mixed-pool regression estimate (beta_protocol ~= +0.05 for RI, ~=0 for SE).

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
})

set.seed(getOption("EpiATLAS_AS_SEED", 42L))

this_transcript_filter <- Sys.getenv(
  "TRANSCRIPT_FILTER",
  getOption("EpiATLAS_AS_PRIMARY_FILTER", "biotype_filtered")
)

sanity_csv <- "qc/sanity_summary.csv"
append_sanity <- function(check_name, status, value, threshold, notes = "") {
  new_row <- data.table(
    check_name = check_name,
    status = status,
    value = value,
    threshold = threshold,
    notes = notes
  )
  if (file.exists(sanity_csv)) {
    # upsert by check_name: drop any prior row for this check before writing
    # the new one, so sanity_summary.csv doesn't accumulate stale/duplicate
    # rows across reruns.
    existing <- fread(sanity_csv)
    existing <- existing[check_name != new_row$check_name]
    fwrite(rbind(existing, new_row, fill = TRUE), sanity_csv)
  } else {
    fwrite(new_row, sanity_csv)
  }
}

message(sprintf("06c-paired-protocol: transcript_filter=%s", this_transcript_filter))

aggregated_dt <- fread(
  sprintf("processed_data/aggregated_dt_filtered_%s.csv.gz", this_transcript_filter),
  select = c("ID", "IHEC", "uuid", "protocol", "PSI", "Event Type")
)

## 1. Identify paired epirrs (both protocol UUIDs present) --------------------
protocol_counts <- aggregated_dt[, .(n_protocols = uniqueN(protocol)), by = IHEC]
paired_epirrs <- protocol_counts[n_protocols >= 2, IHEC]
message(sprintf(
  "Paired epirrs (both protocols): %d / %d total",
  length(paired_epirrs), nrow(protocol_counts)
))
append_sanity(
  sprintf("paired_protocol_n_paired_epirrs_%s", this_transcript_filter),
  if (length(paired_epirrs) > 0L) "PASS" else "FAIL",
  length(paired_epirrs), 1,
  "epirrs with both mRNA-Seq and total-RNA-Seq UUIDs; expected 36 per prior RNA-Seq layer count (441 UUIDs - 405 epirrs)"
)

## 2. Per-event paired PSI, wide by protocol -----------------------------------
paired_dt <- aggregated_dt[IHEC %in% paired_epirrs]
# collapse to one PSI per (ID, IHEC, protocol) in case of duplicate uuids under
# the same protocol for a given epirr (not expected, but a safe aggregation).
paired_wide <- dcast(
  paired_dt,
  ID + IHEC + `Event Type` ~ protocol,
  value.var = "PSI",
  fun.aggregate = mean
)
stopifnot(
  "Expected exactly the two known protocol levels as dcast columns" = all(
    c("mRNA-Seq", "total-RNA-Seq") %in% names(paired_wide)
  )
)
paired_wide <- paired_wide[!is.na(`mRNA-Seq`) & !is.na(`total-RNA-Seq`)]
paired_wide[, dPSI := `total-RNA-Seq` - `mRNA-Seq`]

message(sprintf(
  "Paired (event x epirr) observations with both protocols present: %d",
  nrow(paired_wide)
))
append_sanity(
  sprintf("paired_protocol_n_paired_obs_%s", this_transcript_filter),
  if (nrow(paired_wide) > 0L) "PASS" else "FAIL",
  nrow(paired_wide), 1,
  "event x paired-epirr rows with PSI present under both protocols"
)

## 3. Paired t-test / mixed model per Event Type -------------------------------
# Two estimators, deliberately redundant as a cross-check:
#   - mixed model dPSI ~ 1 + (1 | IHEC): accounts for multiple events sharing
#     the same epirr (non-independent), gives a fixed-effect intercept + a
#     normal-approximation Wald p-value (lmerTest is not installed, so the
#     p-value here is 2*pnorm(-abs(t)), a standard large-sample approximation,
#     not a Satterthwaite/KR-corrected one).
#   - one-sample t-test on each epirr's OWN mean dPSI (collapses the repeated
#     events within an epirr into one independent observation per epirr first,
#     sidestepping the clustering problem entirely) -- the spec's "equivalent"
#     alternative, and a dependency-free sanity check on the mixed model.
test_dPSI <- function(dt_sub) {
  mm <- tryCatch(
    suppressMessages(suppressWarnings(lmer(dPSI ~ 1 + (1 | IHEC), data = dt_sub))),
    error = function(e) NULL
  )
  mm_est <- NA_real_
  mm_se <- NA_real_
  mm_p <- NA_real_
  if (!is.null(mm)) {
    co <- summary(mm)$coefficients
    mm_est <- co[1, "Estimate"]
    mm_se <- co[1, "Std. Error"]
    mm_p <- 2 * pnorm(-abs(mm_est / mm_se))
  }

  per_epirr_mean <- dt_sub[, .(mean_dPSI = mean(dPSI)), by = IHEC]
  tt <- t.test(per_epirr_mean$mean_dPSI)

  list(
    n_obs = nrow(dt_sub),
    n_epirrs = uniqueN(dt_sub$IHEC),
    mixed_model_estimate = mm_est,
    mixed_model_se = mm_se,
    mixed_model_p_value = mm_p,
    paired_ttest_estimate = unname(tt$estimate),
    paired_ttest_ci_low = tt$conf.int[1],
    paired_ttest_ci_high = tt$conf.int[2],
    paired_ttest_p_value = tt$p.value
  )
}

dpsi_results <- paired_wide[, test_dPSI(.SD), by = `Event Type`]
fwrite(
  dpsi_results,
  sprintf("processed_data/paired_protocol_dPSI_summary_%s.csv.gz", this_transcript_filter)
)
print(dpsi_results)

for (i in seq_len(nrow(dpsi_results))) {
  et <- dpsi_results[i, `Event Type`]
  append_sanity(
    sprintf("paired_protocol_dPSI_%s_%s", et, this_transcript_filter),
    "INFO",
    dpsi_results[i, paired_ttest_estimate],
    0,
    sprintf(
      "paired-epirr mean dPSI (total-RNA-Seq - mRNA-Seq) for %s; mixed_model_est=%.4f (p=%.3g), paired_ttest_est=%.4f (p=%.3g, n_epirrs=%d); compare to mixed-pool regression beta_protocol (RI~+0.05, SE~0)",
      et,
      dpsi_results[i, mixed_model_estimate], dpsi_results[i, mixed_model_p_value],
      dpsi_results[i, paired_ttest_estimate], dpsi_results[i, paired_ttest_p_value],
      dpsi_results[i, n_epirrs]
    )
  )
}

## 4. Cross-protocol PSI concordance (§4.13h) ----------------------------------
concordance_dt <- paired_wide[, .(
  n = .N,
  pearson_r = suppressWarnings(cor(`mRNA-Seq`, `total-RNA-Seq`, method = "pearson")),
  spearman_r = suppressWarnings(cor(`mRNA-Seq`, `total-RNA-Seq`, method = "spearman"))
), by = `Event Type`]
fwrite(
  concordance_dt,
  sprintf("processed_data/paired_protocol_concordance_%s.csv.gz", this_transcript_filter)
)
print(concordance_dt)

for (i in seq_len(nrow(concordance_dt))) {
  et <- concordance_dt[i, `Event Type`]
  expected_min <- if (et == "SE") 0.85 else 0.70
  observed_r <- concordance_dt[i, pearson_r]
  append_sanity(
    sprintf("paired_protocol_concordance_%s_%s", et, this_transcript_filter),
    if (!is.na(observed_r) && observed_r >= expected_min) "PASS" else "WARN",
    observed_r, expected_min,
    sprintf(
      "cross-protocol Pearson r of paired PSI for %s (n=%d); expected >= %.2f per spec",
      et, concordance_dt[i, n], expected_min
    )
  )
}

## 5. Per-event concordance / protocol-sensitivity flag ------------------------
minimum_paired_obs_per_event <- 5
per_event_concordance <- paired_wide[, {
  if (.N >= minimum_paired_obs_per_event) {
    list(n_pairs = .N, r = suppressWarnings(cor(`mRNA-Seq`, `total-RNA-Seq`, method = "pearson")))
  } else {
    list(n_pairs = .N, r = NA_real_)
  }
}, by = .(ID, `Event Type`)]
per_event_concordance[, protocol_sensitive := !is.na(r) & r < 0.5]
fwrite(
  per_event_concordance,
  sprintf("processed_data/paired_protocol_per_event_concordance_%s.csv.gz", this_transcript_filter)
)

n_flagged <- per_event_concordance[, sum(protocol_sensitive, na.rm = TRUE)]
n_testable <- per_event_concordance[!is.na(r), .N]
message(sprintf(
  "Per-event concordance: %d/%d testable events (>= %d paired obs) flagged protocol-sensitive (r < 0.5)",
  n_flagged, n_testable, minimum_paired_obs_per_event
))
append_sanity(
  sprintf("paired_protocol_flagged_events_%s", this_transcript_filter),
  "INFO",
  n_flagged, n_testable,
  sprintf(
    "events with cross-protocol PSI r < 0.5 (>= %d paired obs required to test)",
    minimum_paired_obs_per_event
  )
)
