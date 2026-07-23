## ---------------------------------------------------------------------------
## Tier-1 fit-free ridge screen — one event per invocation
##
## "Is this event's PSI predictable from its local epigenome, beyond
## expression/protocol?" — answered with a CLOSED-FORM group-leave-one-ontology-
## out ridge partial-R²/CCC (no tuning loop, no model refit), tested against a
## feature-rotation null over many matched control events (each null = one cheap
## matrix computation, so hundreds are affordable → exact empirical p → FDR in
## the aggregator).
##
## This REPLACES the significance-testing role the Tier-2 elastic-net (09zz) used
## to do with 10 rotations. 09zz now runs on Tier-1 hits only, real PSI only.
##
## CLI (mirrors 09-1-ml-local-array.sh dispatch of 09zz):
##   Rscript 09s-ridge-screen.R <cfg_rds> <event_id> [cores]
## Reuses the SAME cfg + session rds that 09-1 writes for the event models.
##
## Output: processed_data/event_models/<tf>/screen/<ID>_screen.csv.gz (one row:
##   real R²/CCC + null mean/sd + p_emp + effect), plus a sidecar
##   <ID>_screen_null.csv.gz (the full per-control null R²/CCC vector) for the
##   09-2 per-event null-distribution histograms — near-free (values already
##   computed, only their summary is kept in the main row).
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  cfg_path <- args[1]
  this_id <- as.integer(args[2])
  cores <- if (length(args) >= 3L) as.integer(args[3L]) else 1L

  cfg <- readRDS(cfg_path)
  setwd(cfg$project_dir)
  data.table::setDTthreads(cores)

  source("09-ml-shared.R")

  # --- config -------------------------------------------------------------
  feature_table_dir <- cfg$feature_table_dir
  event_dir <- cfg$event_dir
  grouping_col <- cfg$grouping_col
  nfolds <- cfg$nfolds
  seed_base <- if (!is.null(cfg$seed)) cfg$seed else getOption("EpiATLAS_AS_SEED", 42L)
  # Tier-1 rotation count: cheap, so default to a few hundred (p-floor ≈ 1/(R+1)).
  n_rotations <- if (!is.null(cfg$screen_rotations)) {
    cfg$screen_rotations
  } else {
    getOption("EpiATLAS_AS_SCREEN_ROTATIONS", 200L)
  }

  screen_dir <- file.path(event_dir, "screen")
  dir.create(screen_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(screen_dir, paste0(this_id, "_screen.csv.gz"))
  if (file.exists(out_file)) {
    message("Already computed: ", this_id)
    quit(save = "no", status = 0L)
  }

  # --- session (event metadata + PSI matrix for control matching) --------
  sess <- readRDS(cfg$session_rds)
  psi_table <- sess$psi_table
  event_dt <- sess$event_dt
  rm(sess)
  gc()

  # NA-safe single-row writer (so every event appears in aggregation).
  write_row <- function(row) {
    data.table::fwrite(data.table::as.data.table(row), out_file)
  }
  na_row <- function(n_samples = NA_integer_, note = NA_character_) {
    list(
      ID = this_id, n_samples = n_samples, n_features = NA_integer_,
      R_used = 0L, lambda = NA_real_, screen_df = NA_real_,
      screen_R2 = NA_real_, screen_CCC = NA_real_,
      null_R2_mean = NA_real_, null_R2_sd = NA_real_, null_CCC_mean = NA_real_,
      p_emp = NA_real_, effect = NA_real_, note = note
    )
  }

  # --- load this event's feature table -----------------------------------
  feature_data <- data.table::fread(
    file.path(feature_table_dir, paste0("feature_table_", this_id, ".csv.gz"))
  )
  feature_data <- feature_data[!is.na(PSI)]
  if (nrow(feature_data) < 6L) {
    write_row(na_row(nrow(feature_data), "too_few_samples"))
    quit(save = "no", status = 0L)
  }

  # --- CV supergroups (leave-one-ontology-group-out) ---------------------
  groups <- resolve_supergroup_folds(feature_data[[grouping_col]], nfolds)
  if (is.null(groups)) {
    write_row(na_row(nrow(feature_data), "collapsed_to_one_group"))
    quit(save = "no", status = 0L)
  }

  # --- partition columns: X (epigenetic) vs Z (confounds) ----------------
  parts <- screen_partition_columns(names(feature_data), grouping_col)
  this_uuids <- feature_data[["uuid"]]
  confound_df <- as.data.frame(
    feature_data[, parts$confound_cols, with = FALSE]
  )
  prep <- prep_event(feature_data[["PSI"]], confound_df)

  X_real <- as.matrix(feature_data[, parts$x_cols, with = FALSE])
  real <- ridge_screen_stat(prep, X_real, groups)
  if (!is.finite(real$R2)) {
    write_row(na_row(nrow(feature_data), "real_stat_na"))
    quit(save = "no", status = 0L)
  }

  # --- matched controls (same selection as 09zz's feature-rotation) ------
  subset_psi_matrix <- psi_table[this_uuids, , drop = FALSE]
  other_ids <- as.integer(colnames(subset_psi_matrix)[
    colSums(is.na(subset_psi_matrix)) == 0 &
      apply(subset_psi_matrix, 2L, sd, na.rm = TRUE) > 0
  ])
  this_event <- event_dt[ID == this_id]
  other_ids <- other_ids[
    other_ids != this_id &
      other_ids %in% event_dt[`Event Type` == this_event$`Event Type`, ID] &
      other_ids %in% event_dt[seqnames != this_event$seqnames, ID] &
      other_ids %in% event_dt[Variability == this_event$Variability, ID] &
      other_ids %in% event_dt[transcript_filter == this_event$transcript_filter, ID]
  ]
  rm(subset_psi_matrix, psi_table)
  gc()

  # --- feature-rotation null: swap in each control's epigenetic matrix,
  #     keep THIS event's PSI + confounds (prep) + groups fixed --------------
  load_control_X <- function(cid) {
    cdt <- tryCatch(
      data.table::fread(
        file.path(feature_table_dir, paste0("feature_table_", cid, ".csv.gz"))
      ),
      error = function(e) NULL
    )
    if (is.null(cdt)) return(NULL)
    cdt <- cdt[uuid %in% this_uuids]
    cdt <- cdt[match(this_uuids, uuid)]
    if (anyNA(cdt$uuid) || !all(cdt$uuid == this_uuids)) return(NULL)
    cparts <- screen_partition_columns(names(cdt), grouping_col)
    as.matrix(cdt[, cparts$x_cols, with = FALSE])
  }

  n_use <- min(n_rotations, length(other_ids))
  null_R2 <- numeric(0)
  null_CCC <- numeric(0)
  if (n_use > 0L) {
    set.seed(seed_base + this_id)
    sampled <- sample(other_ids, n_use)
    null_stats <- pbmcapply::pbmclapply(
      sampled,
      function(cid) {
        Xc <- load_control_X(cid)
        if (is.null(Xc)) return(c(NA_real_, NA_real_))
        s <- ridge_screen_stat(prep, Xc, groups)
        c(s$R2, s$CCC)
      },
      mc.cores = cores
    )
    null_mat <- do.call(rbind, null_stats)
    null_R2 <- null_mat[is.finite(null_mat[, 1L]), 1L]
    null_CCC <- null_mat[is.finite(null_mat[, 2L]), 2L]

    # Persist the full per-control null vector (sidecar) for the 09-2 per-event
    # null-distribution histograms. These values are already computed above —
    # the main screen row only keeps their mean/sd — so this is a near-free
    # write, not extra compute. One row per sampled control (order matches
    # `sampled`); `control_id` records which locus was borrowed. NAs kept so the
    # count reflects attempted rotations; 09-2 filters to finite for the p-value
    # geometry. Filename is `_screen_null.csv.gz` (NOT `_screen.csv.gz`) so
    # 09s-aggregate.R's `_screen\.csv\.gz$` glob does not pick it up.
    data.table::fwrite(
      data.table::data.table(
        ID = this_id,
        control_id = sampled,
        null_R2 = null_mat[, 1L],
        null_CCC = null_mat[, 2L]
      ),
      file.path(screen_dir, paste0(this_id, "_screen_null.csv.gz"))
    )
  }

  R_used <- length(null_R2)
  p_emp <- if (R_used > 0L) {
    (1 + sum(null_R2 >= real$R2)) / (1 + R_used)
  } else {
    NA_real_
  }
  effect <- if (R_used > 0L) real$R2 - mean(null_R2) else NA_real_

  write_row(list(
    ID = this_id,
    n_samples = nrow(feature_data),
    n_features = real$n_features,
    R_used = R_used,
    lambda = real$lambda,
    screen_df = real$df,
    screen_R2 = real$R2,
    screen_CCC = real$CCC,
    null_R2_mean = if (R_used > 0L) mean(null_R2) else NA_real_,
    null_R2_sd = if (R_used > 0L) sd(null_R2) else NA_real_,
    null_CCC_mean = if (length(null_CCC) > 0L) mean(null_CCC) else NA_real_,
    p_emp = p_emp,
    effect = effect,
    note = NA_character_
  ))
  message("Done screen: ", this_id, " R2=", round(real$R2, 4),
          " p=", signif(p_emp, 3), " (R_used=", R_used, ")")
}
