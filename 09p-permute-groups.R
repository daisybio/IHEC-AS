#!/usr/bin/env Rscript
# Permuted-groups control for the Tier-1 screen (EXECUTION-PLAN item 1b).
#
# QUESTION: how much of Tier-1's performance comes from the test sample's group being
# represented in training?
#
# DESIGN: permute the supergroup LABELS across samples rather than splitting randomly. That
# preserves the fold count and the fold-size distribution exactly, so the only thing that
# changes is whether each fold corresponds to a real group. A random split would change
# sizes too and confound the two effects.
#
# WHAT IT MEASURES: the magnitude of exploitable group structure. NOT cell-type identity --
# a permuted fold exposes the test sample's batch, lab and protocol along with its ontology,
# so this measures whatever `ontology` proxies. See PATH-FORWARD sections 1 and 3 item 4c.
#
# WHAT IT ALSO PICKS UP: interpolation vs extrapolation. Held-out real groups sit far outside
# training support, so part of any gap is simply that permuted folds interpolate. That part is
# measurable -- max_abs_z / frac_z_gt10 are recorded for both arms so the gap can be
# stratified by out-of-support severity rather than reported as one pooled number.
#
# CHANGES NOTHING UPSTREAM: sources 09-ml-shared.R read-only and calls the production statistic
# functions (screen_partition_columns / resolve_supergroup_folds / prep_event / ridge_screen_stat),
# so it cannot drift from the screen it is a control for. No rotations -- p_emp under a permuted
# grouping answers nothing worth quoting, and skipping them is ~200x cheaper.
#
# STAGE 09p, and it CANNOT sit at 03b beside the variance decomposition despite both being
# "promoted controls": it samples its events FROM screen_results.csv.gz and reads the per-event
# feature tables, so it depends on screen_aggregate and build_feature_tables. Structurally it is
# the twin of 09f-floor-sample.R, and it is placed beside it for that reason.
# Promoted from an undeclared prototype 2026-09-01: 7 of paper_numbers' keys come from here, and Nature
# requires the code behind a quoted number to be available at submission.

suppressPackageStartupMessages(library(data.table))
setDTthreads(4L)
source("09-ml-shared.R")

tf <- "biotype_filtered"
# Production default is the size the reported result was measured at (n = 100, 2026-08-28), NOT
# the 12 the pilot used -- the stage re-runs at whatever is baked in here, so the pilot value
# would silently downgrade the published number. Snakemake passes it from
# config["permuted_groups_n_events"].
N_EVENTS <- as.integer(Sys.getenv("N_EVENTS", "100"))
GROUPING <- "ontology"
NFOLDS <- 5L

cfg <- readRDS(sprintf("processed_data/event_glmnet_cfg_%s.rds", tf))
ftdir <- cfg$feature_table_dir
seed_base <- if (!is.null(cfg$seed)) cfg$seed else 42L

# Sample events across Event Type x n_samples, from those that screened cleanly.
scr <- fread(sprintf("processed_data/event_models/%s/screen_results.csv.gz", tf),
  select = c("ID", "Event Type", "feature_set", "n_samples", "screen_R2", "note"))
scr <- scr[feature_set == "local" & (is.na(note) | note == "") & is.finite(screen_R2)]
set.seed(11)
pick <- rbindlist(lapply(c("RI", "SE"), function(et) {
  d <- scr[`Event Type` == et][order(n_samples)]
  d[unique(round(seq(1, .N, length.out = N_EVENTS %/% 2)))]
}))
cat(sprintf("events: %d (%s)\n\n", nrow(pick),
  paste(pick[, .N, by = `Event Type`][, sprintf("%s=%d", `Event Type`, N)], collapse = " ")))

one_event <- function(id) {
  f <- file.path(ftdir, sprintf("feature_table_%d.csv.gz", id))
  if (!file.exists(f)) return(NULL)
  fd <- fread(f)[!is.na(PSI)]
  if (nrow(fd) < 6L) return(NULL)
  groups <- resolve_supergroup_folds(fd[[GROUPING]], NFOLDS)
  if (length(unique(groups)) < 2L) return(NULL)
  # same refusal the worker applies: a 1-row train partition is not screenable
  if (length(groups) - max(table(groups)) < 2L) return(NULL)

  parts <- screen_partition_columns(names(fd), GROUPING)
  cdf <- as.data.frame(fd[, parts$confound_cols, with = FALSE])
  xa <- parts$x_cols
  chm <- xa[grepl("chromhmm", xa, fixed = TRUE)]
  spaces <- list(local = base::setdiff(xa, chm), long = xa)

  # permute the LABELS -> fold count and sizes identical, group identity destroyed
  set.seed(seed_base + id)
  groups_perm <- sample(groups)

  rbindlist(lapply(names(spaces), function(sp) {
    cols <- spaces[[sp]]
    if (!length(cols)) return(NULL)
    X <- as.matrix(fd[, cols, with = FALSE])
    out <- lapply(list(real = groups, perm = groups_perm), function(g) {
      p <- prep_event(fd[["PSI"]], cdf, g)
      ridge_screen_stat(p, X, g)
    })
    data.table(
      ID = id, feature_set = sp, n_samples = nrow(fd), n_features = length(cols),
      R2_real = out$real$R2, R2_perm = out$perm$R2,
      R2b_real = out$real$R2_bounded, R2b_perm = out$perm$R2_bounded,
      z_real = out$real$max_abs_z, z_perm = out$perm$max_abs_z,
      f10_real = out$real$frac_z_gt10, f10_perm = out$perm$frac_z_gt10,
      lam_real = out$real$lambda, lam_perm = out$perm$lambda
    )
  }))
}

res <- rbindlist(lapply(seq_len(nrow(pick)), function(i) {
  id <- pick$ID[i]
  r <- tryCatch(one_event(id), error = function(e) {
    cat(sprintf("  event %d ERROR: %s\n", id, conditionMessage(e)))
    NULL
  })
  if (!is.null(r)) cat(sprintf("  %d done (%s)\n", id, pick$`Event Type`[i]))
  r
}), fill = TRUE)

# Persist the raw fits BEFORE any post-processing -- 12 events x 2 spaces x 2 fits is
# minutes of compute and a downstream bug must not discard it (it did, first run).
# Per-transcript_filter, like everything from 02-2 on: event `ID` is a row position in that
# filter's event_gr, so a single shared path is overwritten by whichever filter ran last and
# reads back as silently wrong events -- the bug that forced aggregated_dt to become per-filter.
.raw_out <- Sys.getenv(
  "PERMUTED_RAW_OUT",
  sprintf("processed_data/permuted_groups_control_raw_%s.csv", tf)
)
dir.create(dirname(.raw_out), recursive = TRUE, showWarnings = FALSE)
fwrite(res, .raw_out)

# NB `i.` cannot prefix a backticked/spaced column name in a data.table update join --
# `i.`Event Type`` is a PARSE error. Alias in the join table's j. Documented in CLAUDE.md;
# this is its FOURTH occurrence in this repo.
res[pick[, .(ID, k_et = `Event Type`)], on = "ID", event_type := i.k_et]
res[, dR2b := R2b_perm - R2b_real]

cat("\n=== per event (bounded R2; raw is blow-up prone, 50.1% of rows |R2|>1) ===\n")
print(res[, .(ID, event_type, feature_set, n_samples,
  R2b_real = round(R2b_real, 4), R2b_perm = round(R2b_perm, 4), dR2b = round(dR2b, 4),
  f10_real = round(f10_real, 3), f10_perm = round(f10_perm, 3))])

cat("\n=== summary by Event Type x feature space ===\n")
print(res[, .(n = .N,
  med_R2b_real = round(median(R2b_real), 4),
  med_R2b_perm = round(median(R2b_perm), 4),
  med_gain = round(median(dR2b), 4),
  pct_perm_better = round(100 * mean(dR2b > 0), 1)),
  by = .(event_type, feature_set)][order(event_type, feature_set)])

cat("\n=== gain stratified by out-of-support severity of the REAL split ===\n")
res[, oos := cut(f10_real, c(-Inf, 0.01, 0.1, 0.5, Inf),
  labels = c("<=1%", "1-10%", "10-50%", ">50%"))]
print(res[, .(n = .N, med_gain = round(median(dR2b), 4),
  med_z_real = round(median(z_real), 1), med_z_perm = round(median(z_perm), 1)),
  by = .(feature_set, oos)][order(feature_set, oos)])

# Snakemake owns the path (PERMUTED_OUT is {output}); the sprintf default keeps the script
# runnable standalone. Same contract as 09f-floor-sample.R's FLOOR_OUT.
f <- Sys.getenv(
  "PERMUTED_OUT",
  sprintf("processed_data/permuted_groups_control_%s.csv", tf)
)
dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
fwrite(res, f)
cat(sprintf("\nwritten: %s\n", f))
cat("\nNOTE: a positive med_gain means permuted folds outperform real ones, i.e. group\n")
cat("structure is exploitable. It does NOT establish that the structure is cell-type\n")
cat("biology, and part of it is interpolation-vs-extrapolation -- see the oos table.\n")
