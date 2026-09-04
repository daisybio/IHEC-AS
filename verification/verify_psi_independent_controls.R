## ---------------------------------------------------------------------------
## Verification for the PSI-independent Tier-1 control change
## (FEATURE_TABLE_VERSION 2 / SCREEN_STAT_VERSION 4).
##
##   Rscript .claude/scratch/verify_psi_independent_controls.R [n_events]
##
## Exercises the REAL assembly (09-ml-shared.R's build_full_event_rows) rather than
## a stand-in. Fit-free apart from one small closed-form ridge, so it is safe on the
## login node.
##
## Covers §5 checks 1, 2, 3 and 4's precondition. Check 5 (byte-identical reruns of
## the real worker) needs the rebuilt v2 tables and is therefore only runnable after
## 05 + build_feature_tables have actually run.
##
## GRID VINTAGE. 05 has not re-run yet, so processed_data/aggregated_dt_{tf}.csv.gz
## does not exist and the only grid on disk is the legacy unsuffixed one, in RAW beta
## and dated 2026-07-17. That is the same vintage as aggregated_dt_filtered, so the
## two are mutually consistent and the join-key checks below are valid. The script
## applies 05's own M-value formula to the grid slice it uses, which is exactly what
## the new 05 will write -- so the assembly is exercised on correctly-scaled input
## without waiting for the rerun. It prefers the per-filter grid when that exists.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))
source("09-ml-shared.R")

tf <- "biotype_filtered"
n_events <- {
  a <- commandArgs(trailingOnly = TRUE)
  if (length(a) >= 1L) as.integer(a[1L]) else 12L
}
set.seed(42L)

pass <- 0L
fail <- 0L
ok <- function(cond, label, detail = "") {
  if (isTRUE(cond)) {
    pass <<- pass + 1L
    cat(sprintf("  PASS  %s%s\n", label, if (nzchar(detail)) paste0(" -- ", detail) else ""))
  } else {
    fail <<- fail + 1L
    cat(sprintf("  FAIL  %s%s\n", label, if (nzchar(detail)) paste0(" -- ", detail) else ""))
  }
}

grid_file_new <- sprintf("processed_data/aggregated_dt_%s.csv.gz", tf)
grid_file_old <- "processed_data/aggregated_dt.csv.gz"
grid_file <- if (file.exists(grid_file_new)) grid_file_new else grid_file_old
grid_is_current <- identical(grid_file, grid_file_new)
cat(sprintf(
  "Grid source: %s (%s)\n", grid_file,
  if (grid_is_current) "per-filter, written by the new 05" else "LEGACY unsuffixed, pre-05-rerun"
))

cat("\nLoading aggregated_dt (filtered) ...\n")
ag <- fread(
  sprintf("processed_data/aggregated_dt_filtered_%s.csv.gz", tf),
  stringsAsFactors = TRUE
)
cols <- classify_feature_columns(names(ag))
# self-consistent with `ag`, so an ID-space change from a newer 03 cannot skew this
all_uuids <- sort(unique(as.character(ag$uuid)))
cat(sprintf("  %d rows, %d cols, %d uuids\n", nrow(ag), ncol(ag), length(all_uuids)))

ids <- sample(ag[, unique(ID)], n_events)
cat(sprintf("Sampled %d events: %s\n", length(ids), paste(head(ids, 8L), collapse = ", ")))

cat("\nLoading grid slice ...\n")
grid <- fread(grid_file, stringsAsFactors = FALSE)
grid <- grid[ID %in% ids]
dnam_cols <- grep("^DNAm;", names(grid), value = TRUE)
grid_dnam_min <- suppressWarnings(min(vapply(
  dnam_cols, function(cn) min(grid[[cn]], na.rm = TRUE), numeric(1)
)))
grid_raw_beta <- is.finite(grid_dnam_min) && grid_dnam_min >= 0
cat(sprintf(
  "  grid slice: %d rows; DNAm min = %.4f -> %s\n", nrow(grid), grid_dnam_min,
  if (grid_raw_beta) "RAW beta" else "M-values"
))

## ===========================================================================
## CHECK 2 -- per-(ID, IHEC) agreement between the grid and the filtered table
## ===========================================================================
## The real assertion behind "shared-row identity": every event-proximal column must
## be a function of (ID, IHEC) alone, and the grid must agree with the filtered table
## on it. If this fails, the grid is the wrong source or the join key is wrong -- and
## since chromHMM aside these ARE the control features, a mismatch would corrupt every
## null silently.
cat("\n[Check 2] grid vs filtered on per-(ID, IHEC) columns\n")
cmp_cols <- base::setdiff(cols$per_epigenome, dnam_cols) # DNAm scale handled below
# The legacy grid predates 05's CpGs NA -> 0 mirror (that fill runs in
# expand-aggregated-dt, i.e. AFTER the old grid write), so on a legacy grid the CpGs
# columns are expected to differ from the filtered table by exactly NA-vs-0. Apply the
# fill here so this check tests what the NEW 05 writes; on a current grid it is a
# no-op and any remaining difference is a genuine failure.
cpgs_cols <- grep("^CpGs", names(grid), value = TRUE)
if (!grid_is_current && length(cpgs_cols)) {
  cat(sprintf(
    "  applying 05's CpGs NA->0 fill to the legacy grid (%d cols) before comparing\n",
    length(cpgs_cols)
  ))
  for (.i in cpgs_cols) {
    set(grid, i = which(is.na(grid[[.i]])), j = .i, value = 0)
  }
}
num_cmp <- cmp_cols[vapply(cmp_cols, function(cn) is.numeric(ag[[cn]]), logical(1))]
fac_cmp <- base::setdiff(cmp_cols, num_cmp)

ag_slice <- unique(ag[ID %in% ids, c("ID", "IHEC", cmp_cols), with = FALSE])
ok(
  nrow(ag_slice) == uniqueN(ag[ID %in% ids, .(ID, IHEC)]),
  "per-(ID,IHEC) columns are constant within (ID, IHEC) in the filtered table",
  sprintf("%d unique rows", nrow(ag_slice))
)

setkey(grid, ID, IHEC)
gk <- grid[, .(ID, IHEC = as.character(IHEC))]
ak <- ag_slice[, .(ID, IHEC = as.character(IHEC))]
i <- match(paste(ak$ID, ak$IHEC), paste(gk$ID, gk$IHEC))
ok(
  !anyNA(i), "every (ID, IHEC) in the filtered slice is present in the grid",
  sprintf("%d/%d matched", sum(!is.na(i)), length(i))
)

if (!anyNA(i)) {
  worst <- 0
  bad <- character(0)
  for (cn in num_cmp) {
    a <- as.numeric(ag_slice[[cn]])
    b <- as.numeric(grid[[cn]][i])
    d <- suppressWarnings(max(abs(a - b), na.rm = TRUE))
    if (!is.finite(d)) d <- 0
    if (!identical(is.na(a), is.na(b))) bad <- c(bad, paste0(cn, "(NA pattern)"))
    if (d > 1e-9) bad <- c(bad, sprintf("%s(%.3g)", cn, d))
    worst <- max(worst, d)
  }
  ok(
    length(bad) == 0L,
    sprintf("all %d numeric per-(ID,IHEC) columns agree exactly", length(num_cmp)),
    if (length(bad)) paste("mismatched:", paste(head(bad, 5L), collapse = ", "))
    else sprintf("max |diff| = %.3g", worst)
  )
  fbad <- character(0)
  for (cn in fac_cmp) {
    a <- as.character(ag_slice[[cn]])
    b <- as.character(grid[[cn]][i])
    if (!identical(a, b)) fbad <- c(fbad, cn)
  }
  ok(
    length(fbad) == 0L,
    sprintf("all %d categorical per-(ID,IHEC) columns agree exactly", length(fac_cmp)),
    if (length(fbad)) paste("mismatched:", paste(fbad, collapse = ", ")) else ""
  )
}

## DNAm: only checkable as an exact match once 05 has re-run. Until then verify the
## transform relationship instead, which is the property that actually matters.
cat("\n[Check 2b] DNAm scale\n")
if (grid_raw_beta) {
  b <- as.numeric(grid[[dnam_cols[1L]]][i])
  a <- as.numeric(ag_slice_dnam <- unique(
    ag[ID %in% ids, c("ID", "IHEC", dnam_cols[1L]), with = FALSE]
  )[[dnam_cols[1L]]])
  beta <- b / 100
  expect <- log2((beta + 0.001) / (1 - beta + 0.001))
  d <- suppressWarnings(max(abs(a - expect), na.rm = TRUE))
  ok(
    is.finite(d) && d < 1e-9,
    "legacy raw-beta grid reproduces the filtered table's M-values under 05's formula",
    sprintf("max |diff| = %.3g on %s", d, dnam_cols[1L])
  )
  cat("  NOTE  exact-match DNAm check deferred until 05 re-runs (grid is raw beta).\n")
} else {
  d <- 0
  for (cn in dnam_cols) {
    a <- as.numeric(unique(ag[ID %in% ids, c("ID", "IHEC", cn), with = FALSE])[[cn]])
    d <- max(d, suppressWarnings(max(abs(a - as.numeric(grid[[cn]][i])), na.rm = TRUE)))
  }
  ok(d < 1e-9, "DNAm agrees exactly between grid and filtered table",
     sprintf("max |diff| = %.3g", d))
}

## Put the grid on the M-value scale the new 05 will write, so the assembly below is
## exercised on correctly-scaled input regardless of vintage.
if (grid_raw_beta && length(dnam_cols)) {
  grid[, (dnam_cols) := lapply(.SD, function(x) {
    beta <- x / 100
    log2((beta + 0.001) / (1 - beta + 0.001))
  }), .SDcols = dnam_cols]
}
grid[, IHEC := factor(as.character(IHEC), levels = levels(ag$IHEC))]
setkey(grid, ID, IHEC)

## ===========================================================================
## Assembly inputs
## ===========================================================================
sample_cov <- unique(ag[, c("uuid", "IHEC", cols$per_uuid), with = FALSE])
stopifnot(!anyDuplicated(sample_cov$uuid))
setkey(sample_cov, uuid)

rbp_score_dt <- fread(sprintf("processed_data/rbp_score_dt_%s.csv.gz", tf))
rbp_score_dt <- rbp_score_dt[
  ID %in% ids,
  c("ID", "uuid", intersect(cols$per_event_sample, names(rbp_score_dt))),
  with = FALSE
]
setkey(rbp_score_dt, ID, uuid)

gene_expr_dt <- fread(sprintf("processed_data/gene_expression_normalised_%s.csv.gz", tf))
gene_expr_dt <- gene_expr_dt[, c("gene_id", "uuid", cols$per_gene), with = FALSE]
# version-suffix strip, mirroring 09-1 (gene_expression_normalised is versioned,
# aggregated_dt is not)
gene_expr_dt[, gene_id := sub("\\.\\d+$", "", gene_id)]
stopifnot(!anyDuplicated(gene_expr_dt, by = c("gene_id", "uuid")))
setkey(gene_expr_dt, gene_id, uuid)

## ===========================================================================
## CHECKS 1 + 3 -- run the real assembly per event
## ===========================================================================
cat("\n[Checks 1 + 3] real build_full_event_rows() per event\n")
res <- rbindlist(lapply(ids, function(id) {
  obs <- ag[ID == id]
  full <- build_full_event_rows(
    obs, id, all_uuids, sample_cov, cols, grid, rbp_score_dt, gene_expr_dt
  )
  # CHECK 1, structurally: the worker loads feature_data[!is.na(PSI)], so the focal
  # fit's input must be byte-identical to what v1 gave it.
  back <- full[!is.na(PSI)]
  setcolorder(back, names(obs))
  identical_focal <- isTRUE(all.equal(
    as.data.frame(back[order(uuid)]), as.data.frame(obs[order(uuid)]),
    check.attributes = FALSE
  ))
  epi_na <- mean(unlist(lapply(
    cols$per_epigenome, function(cn) is.na(full[[cn]])
  )))
  epi_na_added <- mean(unlist(lapply(
    cols$per_epigenome, function(cn) is.na(full[[cn]][is.na(full$PSI)])
  )))
  data.table(
    ID = id, n_obs = nrow(obs), n_full = nrow(full),
    cols_match = identical(names(full), names(obs)),
    focal_identical = identical_focal,
    types_match = all(vapply(names(obs), function(cn) {
      identical(class(full[[cn]]), class(obs[[cn]])) &&
        (!is.factor(obs[[cn]]) || identical(levels(full[[cn]]), levels(obs[[cn]])))
    }, logical(1))),
    epi_na_frac = epi_na, epi_na_frac_added = epi_na_added
  )
}))

ok(all(res$n_full == length(all_uuids)),
   sprintf("every rebuilt table has the full cohort row set (%d)", length(all_uuids)),
   sprintf("observed %d-%d -> %d-%d rows",
           min(res$n_obs), max(res$n_obs), min(res$n_full), max(res$n_full)))
ok(all(res$cols_match), "column set and order unchanged")
ok(all(res$types_match), "column types and factor levels preserved")
ok(all(res$focal_identical),
   "CHECK 1: focal input identical -- full[!is.na(PSI)] == aggregated_dt[ID == id]",
   sprintf("%d/%d events", sum(res$focal_identical), nrow(res)))
ok(max(res$epi_na_frac_added) < 0.5,
   "added rows are actually populated with epigenetic features",
   sprintf("NA fraction on added rows: %.3f-%.3f (all rows: %.3f-%.3f)",
           min(res$epi_na_frac_added), max(res$epi_na_frac_added),
           min(res$epi_na_frac), max(res$epi_na_frac)))
cat(sprintf(
  "  row growth: mean %.1f -> %d (x%.2f)\n",
  mean(res$n_obs), length(all_uuids), length(all_uuids) / mean(res$n_obs)
))

## ===========================================================================
## CHECK 4 precondition -- eligible control pool under the new criteria
## ===========================================================================
cat("\n[Check 4 precondition] matched-control pool size, per Event Type\n")
ea <- fread(sprintf("processed_data/event_annotations_dt_%s.csv.gz", tf))
# The screen runs on all_modelable_ids, NOT on every ID in the filtered table: 09-1
# additionally requires >= minimum_events samples and sd(PSI) > 0. Using the filtered
# table's raw ID set inflates both m and the pool (RI 1,846 vs the real 1,776), and m
# is what the BH bound divides by -- so reproduce 09-1's filters here. Verified to
# reproduce event_glmnet_all_ids exactly (34,146 ids, setequal TRUE).
.min_events <- 25L # .Rprofile's minimum_events
.n_by <- ag[, .N, by = ID]
.sd_by <- ag[, .(sd = sd(PSI, na.rm = TRUE)), by = ID]
modelable <- base::intersect(.n_by[N >= .min_events, ID], .sd_by[sd > 0, ID])
cat(sprintf(
  "  modelable events: %d of %d in the filtered table (minimum_events=%d, sd>0)\n",
  length(modelable), uniqueN(ag$ID), .min_events
))
e <- ea[ID %in% modelable, .(ID, ET = `Event Type`, seqnames, Variability,
                             tfc = transcript_filter)]
e[, grp_n := .N, by = .(ET, Variability, tfc)]
e[, same_chr := .N, by = .(ET, Variability, tfc, seqnames)]
e[, pool := grp_n - same_chr]
pool <- e[, .(events = .N, min_pool = min(pool),
              median_pool = as.numeric(median(pool)), max_pool = max(pool)),
          by = ET][order(ET)]
print(pool)

spec <- c(RI = 818L, SE = 200L)
for (et in pool$ET) {
  if (!et %in% names(spec)) next
  R <- resolve_screen_rotations(spec, et)
  m <- pool[ET == et, events]
  med <- pool[ET == et, median_pool]
  R_eff <- min(R, med)
  floor_p <- 1 / (R_eff + 1)
  k_needed <- ceiling(floor_p * m / 0.1)
  cat(sprintf(
    "  %s: m=%d  pool(median)=%.0f  cap=%d  R_eff=%.0f  p_floor=%.3g  k needed at q=0.1: %d\n",
    et, m, med, R, R_eff, floor_p, k_needed
  ))
}
ok(pool[ET == "RI", median_pool] > 700,
   "CHECK 4: RI's matched pool clears 700 (vs ~92 measured under the PSI gate)",
   sprintf("median %.0f", pool[ET == "RI", median_pool]))

## ===========================================================================
cat(sprintf("\n%d passed, %d failed\n", pass, fail))
if (!grid_is_current) {
  cat(paste0(
    "\nREMAINING, after 05 + build_feature_tables run:\n",
    "  * re-run this script (it will switch to the per-filter grid and check DNAm exactly)\n",
    "  * CHECK 4 proper: read R_used / n_eligible_controls / n_rotations_requested\n",
    "    out of the real screen rows and confirm RI's R_used actually rose\n",
    "  * CHECK 5: two fresh runs of one event via\n",
    "    bash .claude/scratch/run_screen_examples.sh fresh <id>\n",
    "    must give byte-identical _screen.csv.gz and _screen_null.csv.gz\n"
  ))
}
quit(save = "no", status = if (fail > 0L) 1L else 0L)
