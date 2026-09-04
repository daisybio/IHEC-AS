#!/usr/bin/env Rscript
# Verify every number 09-2-ml-local-new.Rmd now states in prose, plus the 09-1
# session objects it loads, against the real on-disk Tier-1 outputs.
#
# Fit-free: reads screen_results.csv.gz + 09-1's session rds only. Safe on the
# login node. Prints EXPECT/GOT/verdict per claim; nonzero exit if any FAIL.

suppressPackageStartupMessages(library(data.table))
setDTthreads(2L)

tf <- "biotype_filtered"
event_dir <- file.path("processed_data", "event_models", tf)
fails <- 0L

chk <- function(label, expect, got, ok) {
  verdict <- if (isTRUE(ok)) "PASS" else {
    fails <<- fails + 1L
    "FAIL"
  }
  cat(sprintf("  [%s] %-52s expect %-22s got %s\n", verdict, label, expect, got))
}
near <- function(a, b, tol) is.finite(a) && is.finite(b) && abs(a - b) <= tol
inrange <- function(x, lo, hi) all(is.finite(x)) && all(x >= lo & x <= hi)

dt <- fread(file.path(event_dir, "screen_results.csv.gz"))
cat(sprintf(
  "\nscreen_results.csv.gz: %s rows x %d cols\n",
  format(nrow(dt), big.mark = ","), ncol(dt)
))

# ---------------------------------------------------------------- 1. integrity
cat("\n=== 1. File integrity / provenance ===\n")
chk("stat_version unique", "5", paste(sort(unique(dt$stat_version)), collapse = ","),
  identical(sort(unique(dt$stat_version)), 5L))
chk("transcript_filter unique", tf, paste(unique(dt$transcript_filter), collapse = ","),
  identical(as.character(unique(dt$transcript_filter)), tf))
n_ev <- uniqueN(dt$ID)
chk("events x 3 feature sets == rows", "34146 x 3 = 102438",
  sprintf("%d x %d = %d", n_ev, uniqueN(dt$feature_set), nrow(dt)),
  n_ev == 34146L && nrow(dt) == 102438L)
chk("feature_set levels", "local,long,short",
  paste(sort(unique(as.character(dt$feature_set))), collapse = ","),
  setequal(unique(as.character(dt$feature_set)), c("local", "long", "short")))
chk("no duplicate (ID, feature_set)", "0 dups",
  sprintf("%d dups", sum(duplicated(dt, by = c("ID", "feature_set")))),
  !any(duplicated(dt, by = c("ID", "feature_set"))))

# --------------------------------------------------------------- 2. rotations
cat("\n=== 2. Rotation counts (per-Event-Type, 818 RI / 200 SE) ===\n")
# R_used is 0, not NA, on the 59 refused events (note != ""), so filter on a finite p
# rather than on R_used itself — otherwise the observed range reads 0-200 / 0-818.
rot <- dt[is.finite(p_emp), .(
  req = paste(sort(unique(n_rotations_requested)), collapse = ","),
  R_min = min(R_used), R_max = max(R_used),
  elig_min = min(n_eligible_controls), elig_med = as.numeric(median(as.numeric(n_eligible_controls)))
), by = `Event Type`][order(`Event Type`)]
print(rot)
chk("RI requested 818", "818", rot[`Event Type` == "RI", req], rot[`Event Type` == "RI", req] == "818")
chk("SE requested 200", "200", rot[`Event Type` == "SE", req], rot[`Event Type` == "SE", req] == "200")
# R_used tops out AT the requested cap but can dip below it: a control whose feature table
# is unreadable, or whose rows do not cover this event's uuids, returns NA from
# control_fs_R2 and is dropped. Real spread: RI 794-818 (633 of 1772 rows below max),
# SE 198-200 (56 rows). So assert "reaches the cap, never exceeds it, stays close".
chk("RI R_used reaches cap 818, no more", "max == 818",
  sprintf("%g-%g", rot[`Event Type` == "RI", R_min], rot[`Event Type` == "RI", R_max]),
  rot[`Event Type` == "RI", R_max] == 818 && rot[`Event Type` == "RI", R_min] >= 780)
chk("SE R_used reaches cap 200, no more", "max == 200",
  sprintf("%g-%g", rot[`Event Type` == "SE", R_min], rot[`Event Type` == "SE", R_max]),
  rot[`Event Type` == "SE", R_max] == 200 && rot[`Event Type` == "SE", R_min] >= 190)

# ------------------------------------------------------------------- 3. hits
cat("\n=== 3. Zero hits, under 09-2's CORRECTED sig rule and the old one ===\n")
qthr <- 0.1
sig_new <- dt[, is.finite(q) & q < qthr]
sig_old <- dt[, is.finite(q) & q < qthr & is.finite(effect) & effect > 0]
chk("hits, q-only rule (09-2 now)", "0", sum(sig_new), sum(sig_new) == 0L)
chk("hits, dropped effect>0 rule", "0 (same at 0 hits)", sum(sig_old), sum(sig_old) == 0L)
# stat_version 5: the pooled-z sensitivity now disagrees with the primary rule -- 7
# candidates, not 0 (see .claude/scratch/investigate_pooled_z_candidates_v5.R). This is
# a tracked, real, open finding, not a bug in this script -- assert the count exactly so
# it is caught if the candidate set ever changes again.
chk("hits, pooled-z sensitivity == 7 (v5 open finding)", "7", sum(dt$q_pooled_z < qthr, na.rm = TRUE),
  sum(dt$q_pooled_z < qthr, na.rm = TRUE) == 7L)
hits_file <- file.path(event_dir, sprintf("tier1_hits_%s.txt", tf))
chk("tier1_hits file empty", "0 bytes", sprintf("%d bytes", file.info(hits_file)$size),
  file.info(hits_file)$size == 0)

# ------------------------------------------------------- 4. per-family claims
cat("\n=== 4. Per-family q minima and p-uniformity (09-2 prose) ===\n")
fam <- dt[is.finite(q), .(
  n_tested = .N,
  min_q = signif(min(q, na.rm = TRUE), 3),
  min_q_pz = signif(min(q_pooled_z, na.rm = TRUE), 3),
  min_p_pz = signif(min(p_pooled_z, na.rm = TRUE), 3),
  pct_p05 = round(100 * mean(p_emp <= 0.05, na.rm = TRUE), 2),
  pct_ppz05 = round(100 * mean(p_pooled_z <= 0.05, na.rm = TRUE), 2),
  pct_ppz001 = round(100 * mean(p_pooled_z <= 0.001, na.rm = TRUE), 3)
), by = .(`Event Type`, feature_set)][order(`Event Type`, feature_set)]
print(fam)
chk("min q over all families == 0.551", "0.551", min(fam$min_q), near(min(fam$min_q), 0.551, 0.005))
chk("min q_pooled_z == 0.0145", "0.0145", min(fam$min_q_pz), near(min(fam$min_q_pz), 0.0145, 0.0003))
chk("p_emp<=0.05 range 4.39-5.97%", "[4.39, 5.97]",
  sprintf("[%.2f, %.2f]", min(fam$pct_p05), max(fam$pct_p05)),
  near(min(fam$pct_p05), 4.39, 0.02) && near(max(fam$pct_p05), 5.97, 0.02))
chk("p_pooled_z<=0.05 range 3.77-6.20%", "[3.77, 6.20]",
  sprintf("[%.2f, %.2f]", min(fam$pct_ppz05), max(fam$pct_ppz05)),
  near(min(fam$pct_ppz05), 3.77, 0.02) && near(max(fam$pct_ppz05), 6.20, 0.02))
chk("p_pooled_z<=0.001 range 0.056-0.169%", "[0.056, 0.169]",
  sprintf("[%.3f, %.3f]", min(fam$pct_ppz001), max(fam$pct_ppz001)),
  near(min(fam$pct_ppz001), 0.056, 0.002) && near(max(fam$pct_ppz001), 0.169, 0.002))

# ---------------------------------------------- 5. effect corruption (Task B)
cat("\n=== 5. `effect` is outlier-corrupted (the claim Task B replaced) ===\n")
fin <- dt[is.finite(p_emp) & is.finite(screen_R2) & is.finite(effect)]
pct_eff_pos <- 100 * mean(fin$effect > 0)
pct_blow <- 100 * mean(abs(fin$screen_R2) > 1)
chk("effect>0 over all finite rows == 79.3%", "79.3%", sprintf("%.1f%%", pct_eff_pos),
  near(pct_eff_pos, 79.3, 0.15))
# Blow-up fraction roughly doubled under the LOOCV+1SE lambda change (29.6% -> 50.1%) --
# more conservative regularization did not reduce out-of-support extrapolation blow-ups.
# Real, unexplained, tracked as an open question, not a bug in this check.
chk("|screen_R2|>1 == 50.1% (v5, up from 29.6%)", "50.1%", sprintf("%.1f%%", pct_blow),
  near(pct_blow, 50.1, 0.15))
# NB two DIFFERENT stratifications live in 09-2's check 5, and mixing them up produces
# spurious mismatches (it did, on the first pass at this script):
#   .blow    = abs(screen_R2) > 10   -> the n_samples / lambda comparison (v5: 279/334 vs 318/1559;
#            v4 was 90/116 vs 328/647, and the v4 "blow-ups are small-n" story is DEAD)
#   .r2_mag  = abs(screen_R2) bins   -> the median_effect / pct_effect_pos magnitude table
#   r2_stratum (check 2, SIGNED R2, cuts at 0/0.1/0.2/0.5) -> the fingerprint table, which is where
#            09-2's prose median_p gradient comes from. The headline blow-up fraction (v5: 50.1%) is
#            the >1 threshold, NOT the .blow group.
blow <- fin[abs(screen_R2) > 10]
oth <- fin[abs(screen_R2) <= 10]
chk("blow-up median n_samples == 279", "279", median(blow$n_samples), near(median(blow$n_samples), 279, 2))
chk("blow-up median lambda == 333.8", "333.8", signif(median(blow$lambda), 4),
  near(median(blow$lambda), 333.8, 5))
chk("rest median n / lambda == 318 / 1559", "318 / 1559",
  sprintf("%g / %g", median(oth$n_samples), signif(median(oth$lambda), 4)),
  near(median(oth$n_samples), 318, 2) && near(median(oth$lambda), 1559, 15))

strat <- fin[, .(
  n = .N,
  median_effect = round(median(effect), 3),
  pct_effect_pos = round(100 * mean(effect > 0), 1),
  median_p = round(median(p_emp), 3)
), by = .(r2_mag = cut(abs(screen_R2), c(-Inf, 1, 10, 100, Inf),
  labels = c("<=1", "(1,10]", "(10,100]", ">100")))][order(r2_mag)]
cat("  (magnitude table, 09-2 check 5 -- NOT the source of the prose numbers)\n")
print(strat)

# --------------------------------- 6. the fingerprint table (Task B's numbers)
# THIS is 09-2 check 2's stratification: SIGNED screen_R2, cuts at 0/0.1/0.2/0.5,
# grouped per (Event Type x feature_set). Every figure Task B quotes comes from here.
cat("\n=== 6. Fingerprint table: signed-R2 strata per family (Task B's source) ===\n")
fin[, r2_stratum := cut(screen_R2, c(-Inf, 0, 0.1, 0.2, 0.5, Inf),
  labels = c("<=0", "(0,0.1]", "(0.1,0.2]", "(0.2,0.5]", ">0.5"))]
fp <- fin[, .(
  n = .N,
  median_R2 = round(median(screen_R2), 3),
  median_effect = round(median(effect), 4),
  pct_effect_pos = round(100 * mean(effect > 0), 1),
  median_p = round(median(p_emp), 3),
  min_q = signif(min(q, na.rm = TRUE), 3)
), by = .(`Event Type`, feature_set, r2_stratum)][order(`Event Type`, feature_set, r2_stratum)]
print(fp)

# NB v5's fingerprint table is far noisier than v4's, driven by RI's small event counts
# (only 1775 tested split across 5 strata; the top stratum holds as few as 1 event for
# RI/short and 0 for RI/long). Assertions below are written to tolerate that noise rather
# than paper over it -- see 09-2's "How to read a report with no hits" for the honest
# version of this table.
base_str <- fp[r2_stratum == "<=0"]
chk("median effect in the <=0 stratum: 0.004..22.76", "[0.004, 22.76]",
  sprintf("[%.3f, %.2f]", min(base_str$median_effect), max(base_str$median_effect)),
  near(min(base_str$median_effect), 0.004, 0.002) && near(max(base_str$median_effect), 22.76, 0.02))
pos <- fp[r2_stratum != "<=0"]
# One n=1 stratum (RI/short, >0.5) is a genuine outlier at 0% -- exclude n<5 strata for the
# "healthy" claim and separately flag the outlier.
pos_robust <- pos[n >= 5L]
chk("pct_effect_pos in strata above 0, n>=5: 80-100%", "[80, 100]",
  sprintf("[%.1f, %.1f]", min(pos_robust$pct_effect_pos), max(pos_robust$pct_effect_pos)),
  min(pos_robust$pct_effect_pos) >= 80 && max(pos_robust$pct_effect_pos) <= 100)
chk("exactly one n<5 outlier stratum below that range", "RI/short >0.5, n=1, 0%",
  sprintf("%d strata with n<5 and pct<80", pos[n < 5L & pct_effect_pos < 80, .N]),
  pos[n < 5L & pct_effect_pos < 80, .N] == 1L)
chk("median_effect not ~0 anywhere (all clearly nonzero)", "all |x| > 0.001",
  sprintf("min |x| = %.4f", min(abs(fp$median_effect))), min(abs(fp$median_effect)) > 0.001)

# The two families that stay strictly monotone under v5 -- SE/local and SE/long. RI/long
# and SE/short (v4's illustrative pair) are NOT monotone under v5; do not reuse them.
g_se_l <- fp[`Event Type` == "SE" & feature_set == "local"]
g_se_g <- fp[`Event Type` == "SE" & feature_set == "long"]
chk("SE/local median_p 0.716 -> 0.030", "0.716 -> 0.030",
  sprintf("%.3f -> %.3f", g_se_l[r2_stratum == "<=0", median_p], g_se_l[r2_stratum == ">0.5", median_p]),
  near(g_se_l[r2_stratum == "<=0", median_p], 0.716, 0.001) &&
    near(g_se_l[r2_stratum == ">0.5", median_p], 0.030, 0.001))
chk("SE/long median_p 0.522 -> 0.127", "0.522 -> 0.127",
  sprintf("%.3f -> %.3f", g_se_g[r2_stratum == "<=0", median_p], g_se_g[r2_stratum == ">0.5", median_p]),
  near(g_se_g[r2_stratum == "<=0", median_p], 0.522, 0.001) &&
    near(g_se_g[r2_stratum == ">0.5", median_p], 0.127, 0.001))
n_mono <- fp[, all(diff(median_p) <= 0), by = .(`Event Type`, feature_set)][, sum(V1)]
n_fam <- uniqueN(fp[, .(`Event Type`, feature_set)])
chk("median_p falls monotonically with signed R2 (only 2/6 now)", "2/6 families",
  sprintf("%d/%d families", n_mono, n_fam), n_mono == 2L)
top <- fp[r2_stratum == ">0.5"]
chk("top stratum holds 1-83 events (RI/long has 0, excluded)", "[1, 83]",
  sprintf("[%d, %d]", min(top$n), max(top$n)),
  min(top$n) >= 1 && max(top$n) <= 83)
chk("top stratum never significant", "min q >= 0.1",
  sprintf("min q = %g", min(top$min_q)), min(top$min_q) >= qthr)

# ------------------------------------------------------ 7. floor / k_at_floor
cat("\n=== 7. Exceedance floor and k_at_floor (Stage 2 mootness) ===\n")
m_fam <- dt[is.finite(q), .N, by = .(`Event Type`, feature_set)]
# `fl` MUST use max(R_used), matching 09s-aggregate.R's floor_tab. R_used varies within
# RI (794-818; 633 of 1772 rows sit below the max), so picking any other representative
# inflates the floor and over-counts k_at_floor — that is what produced a spurious
# RI/local k=2 against the real 1 on this script's first pass.
kf <- dt[is.finite(p_emp), {
  fl <- 1 / (max(R_used, na.rm = TRUE) + 1)
  .(R = max(R_used, na.rm = TRUE), floor = signif(fl, 3), m = .N,
    k_at_floor = sum(p_emp <= fl + 1e-12),
    k_needed = ceiling(fl * .N / qthr))
}, by = .(`Event Type`, feature_set)][order(`Event Type`, feature_set)]
print(kf)
chk("RI k_at_floor local/long/short == 2/2/2", "2,2,2",
  paste(kf[`Event Type` == "RI"][order(feature_set), k_at_floor], collapse = ","),
  identical(as.integer(kf[`Event Type` == "RI"][order(feature_set), k_at_floor]), c(2L, 2L, 2L)))
chk("SE k_at_floor local/long/short == 199/191/208", "199,191,208",
  paste(kf[`Event Type` == "SE"][order(feature_set), k_at_floor], collapse = ","),
  identical(as.integer(kf[`Event Type` == "SE"][order(feature_set), k_at_floor]), c(199L, 191L, 208L)))
chk("RI k_needed == 22", "22", paste(unique(kf[`Event Type` == "RI", k_needed]), collapse = ","),
  all(kf[`Event Type` == "RI", k_needed] == 22L))
chk("SE k_needed == 1608", "1608", paste(unique(kf[`Event Type` == "SE", k_needed]), collapse = ","),
  all(kf[`Event Type` == "SE", k_needed] == 1608L))
chk("no family can pass at the floor", "floor > k*q/m everywhere",
  sprintf("%d families can pass", kf[, sum(floor <= k_at_floor * qthr / m)]),
  kf[, all(floor > k_at_floor * qthr / m)])

# ------------------------------------------------------------------ 8. notes
cat("\n=== 8. Refusal notes ===\n")
nt <- dt[!is.na(note) & note != "", .N, by = note][order(-N)]
print(nt)
gn <- function(x) { v <- nt[note == x, N]; if (length(v)) v else 0L }
chk("fold_train_too_small:1 == 153", "153", gn("fold_train_too_small:1"), gn("fold_train_too_small:1") == 153L)
# spliceosome_* added as confounds (v5) -> more fold-level confound collinearity.
# Handled safely (aliased columns contribute 0, not a crash), but a real 45x rise.
chk("fold_rank_deficient == 3753 (v5, up from 84)", "3753", gn("fold_rank_deficient"),
  gn("fold_rank_deficient") == 3753L)
chk("collapsed_to_one_group == 6", "6", gn("collapsed_to_one_group"), gn("collapsed_to_one_group") == 6L)
n_excl <- uniqueN(dt[!is.finite(q), ID])
chk("events excluded from testing == 53", "53", n_excl, n_excl == 53L)

# --------------------------------------------------- 9. 09-1 session objects
cat("\n=== 9. 09-1 outputs that 09-2 loads ===\n")
sess_file <- sprintf("processed_data/session_09_1_ml_local_%s.rds", tf)
chk("09-1 session rds exists", "exists", file.exists(sess_file), file.exists(sess_file))
sess <- readRDS(sess_file)
cat(sprintf("  session objects: %s\n", paste(names(sess), collapse = ", ")))
ev <- as.data.table(sess$event_dt)
chk("event_dt has required cols", "ID, Event Type, seqnames, Variability",
  paste(intersect(names(ev), c("ID", "Event Type", "seqnames", "Variability")), collapse = ", "),
  all(c("ID", "Event Type", "seqnames", "Variability") %in% names(ev)))
chk("every screened ID present in event_dt", "0 missing",
  sprintf("%d missing", length(setdiff(unique(dt$ID), ev$ID))),
  length(setdiff(unique(dt$ID), ev$ID)) == 0L)
et_screen <- dt[, .N, by = .(ID, `Event Type`)][, .N, by = `Event Type`][order(`Event Type`)]
print(et_screen)
chk("screened events == m per Event Type", "RI 1776 + SE 32370 = 34146",
  sprintf("RI %d + SE %d = %d", et_screen[`Event Type` == "RI", N],
    et_screen[`Event Type` == "SE", N], sum(et_screen$N)),
  sum(et_screen$N) == 34146L)
# Event Type must agree between the screen rows and 09-1's event_dt
mrg <- merge(unique(dt[, .(ID, et_screen = as.character(`Event Type`))]),
  unique(ev[, .(ID, et_sess = as.character(`Event Type`))]), by = "ID")
chk("Event Type agrees screen vs 09-1", "0 mismatches",
  sprintf("%d mismatches", mrg[et_screen != et_sess, .N]), mrg[et_screen != et_sess, .N] == 0L)
# keep_rows_manual is an INTEGER VECTOR OF EVENT IDs (39,176 of them, range 2..65,124), NOT a logical
# mask. An earlier version of this line printed sum() of it as a "TRUE" count -- 1,300,517,697, five
# orders of magnitude above the row total -- and the impossibility went unnoticed. Assert the type.
if (!is.null(sess$keep_rows_manual)) {
  .kr <- sess$keep_rows_manual
  chk("keep_rows_manual is an ID vector, not a mask", "integer, not logical",
    sprintf("%s, n=%s", paste(class(.kr), collapse = "/"), format(length(.kr), big.mark = ",")),
    is.integer(.kr) && !is.logical(.kr) && length(.kr) == uniqueN(.kr))
  chk("keep_rows_manual IDs cover every screened event", "0 missing",
    sprintf("%d missing", length(setdiff(unique(dt$ID), .kr))),
    length(setdiff(unique(dt$ID), .kr)) == 0L)
}

# ------------------------------------------------- 10. Tier-2 (must be empty)
cat("\n=== 10. Tier-2 / en_available (must be FALSE at 0 hits) ===\n")
for (pat in c("_all_metrics", "_best_params", "_fit_metrics", "_event_summary",
              "_nonzero_coefs", "_fold_features")) {
  n <- length(list.files(event_dir, pattern = sprintf("^\\d+%s\\.csv\\.gz$", pat)))
  chk(sprintf("EN files %s", pat), "0", n, n == 0L)
}
n_skip <- length(list.files(event_dir, pattern = "^\\d+_skipped\\.csv\\.gz$"))
cat(sprintf("  (%d *_skipped.csv.gz stubs present — Tier-2 skip markers, not EN results)\n", n_skip))

cat(sprintf("\n==== %s ====\n", if (fails == 0L) "ALL CHECKS PASS" else
  sprintf("%d CHECK(S) FAILED", fails)))
quit(save = "no", status = if (fails == 0L) 0L else 1L)
