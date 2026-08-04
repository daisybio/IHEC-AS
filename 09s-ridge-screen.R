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

  # cfg_path <- "processed_data/event_glmnet_cfg_biotype_filtered.rds"
  # this_id <- 20541
  # cores <- 10L

  cfg <- readRDS(cfg_path)
  setwd(cfg$project_dir)
  data.table::setDTthreads(cores)

  source("09-ml-shared.R")

  # --- config -------------------------------------------------------------
  feature_table_dir <- cfg$feature_table_dir
  event_dir <- cfg$event_dir
  grouping_col <- cfg$grouping_col
  nfolds <- cfg$nfolds
  seed_base <- if (!is.null(cfg$seed)) {
    cfg$seed
  } else {
    getOption("EpiATLAS_AS_SEED", 42L)
  }
  # Tier-1 rotation count: cheap, so default to a few hundred (p-floor ≈ 1/(R+1)).
  # PER EVENT TYPE as of 2026-08-04 — resolved further down, once this event's
  # Event Type is known (see resolve_screen_rotations in 09-ml-shared.R for why a
  # single scalar cannot serve both RI and SE).
  rotation_spec <- if (!is.null(cfg$screen_rotations)) {
    cfg$screen_rotations
  } else {
    getOption("EpiATLAS_AS_SCREEN_ROTATIONS", 200L)
  }

  screen_dir <- file.path(event_dir, "screen")
  dir.create(screen_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(screen_dir, paste0(this_id, "_screen.csv.gz"))
  null_file <- file.path(screen_dir, paste0(this_id, "_screen_null.csv.gz"))

  # Screen is now run PER FEATURE SET (long/short/local) — the omnibus single
  # ridge over all epigenetic features under-powers spatially-local signal in
  # p>>n (a strong exon-local signal diluted by hundreds of far chromHMM windows).
  # `long`=all epigenetic; `short`=drop far chromHMM; `local`=drop all chromHMM.
  # Screen feature sets: prefer the screen-specific cfg field (lets the screen run
  # short/local-only in long-reuse mode without touching Tier-2's feature_sets).
  feature_sets <- if (!is.null(cfg$screen_feature_sets)) {
    cfg$screen_feature_sets
  } else if (!is.null(cfg$feature_sets)) {
    cfg$feature_sets
  } else {
    c("long", "short", "local")
  }

  # Resume gate is PER FEATURE SET, not per file: a pre-existing output from the
  # earlier omnibus run has no `feature_set` column at all — those rows ARE the
  # `long` set, so they're retained (tagged) and only the missing sets recomputed.
  # Same for the null sidecar, so 09-2's Plot 4 keeps its long facets.
  # Rows are only reusable if they were produced by the CURRENT statistic. The
  # stamp (SCREEN_STAT_VERSION, defined next to prep_event/ridge_screen_stat in
  # 09-ml-shared.R) is compared here; anything stamped with an older version -- or
  # unstamped, which means version 1, the pre-2026-07-29 leaky statistic -- is
  # DISCARDED and recomputed. Without this check a correctness fix is invisible to
  # the reuse logic: Snakemake re-runs every per-event job because the script
  # changed, but each exits "Already computed" on file existence and keeps the
  # stale numbers.
  read_existing <- function(path) {
    if (!file.exists(path)) {
      return(NULL)
    }
    dt <- tryCatch(data.table::fread(path), error = function(e) NULL)
    if (is.null(dt) || !nrow(dt)) {
      return(NULL)
    }
    ver <- if ("stat_version" %in% names(dt)) {
      suppressWarnings(as.integer(dt$stat_version))
    } else {
      rep(1L, nrow(dt))
    }
    stale <- is.na(ver) | ver != SCREEN_STAT_VERSION
    if (all(stale)) {
      message(
        "Discarding ", basename(path), ": stat_version ",
        paste(unique(ifelse(is.na(ver), "NA", ver)), collapse = "/"),
        " != current ", SCREEN_STAT_VERSION, " -- recomputing"
      )
      return(NULL)
    }
    if (any(stale)) {
      message(
        "Dropping ", sum(stale), " stale row(s) from ", basename(path),
        " (stat_version != ", SCREEN_STAT_VERSION, ")"
      )
      dt <- dt[!stale]
    }
    if (!"feature_set" %in% names(dt)) {
      dt[, feature_set := "long"]
    }
    dt[]
  }
  existing_rows <- read_existing(out_file)
  existing_null <- read_existing(null_file)
  if (!is.null(existing_rows)) {
    feature_sets <- base::setdiff(
      feature_sets,
      unique(as.character(existing_rows$feature_set))
    )
    if (length(feature_sets) == 0L) {
      message("Already computed (all feature sets): ", this_id)
      quit(save = "no", status = 0L)
    }
    message(
      "Resuming ",
      this_id,
      " — missing sets: ",
      paste(feature_sets, collapse = ",")
    )
  }

  # --- session (event metadata + PSI matrix for control matching) --------
  sess <- readRDS(cfg$session_rds)
  psi_table <- sess$psi_table
  event_dt <- sess$event_dt
  # chromHMM-vicinity objects — needed for the long/short/local partition
  # (mirrors 09zz's build_explanatory_vars): `short` drops far-vicinity chromHMM,
  # `local` drops all chromHMM.
  chromhmm_hits_smaller <- sess$chromhmm_hits_smaller
  keep_rows_manual <- sess$keep_rows_manual
  rm(sess)
  gc()

  this_et <- as.character(event_dt[ID == this_id, `Event Type`][1L])
  this_tf <- as.character(event_dt[ID == this_id, transcript_filter][1L])

  # Rotation count for THIS event type. RI needs its (near-)full ~818-control pool
  # to push the empirical-p floor under BH's bound; SE clears at 200 and would cost
  # ~10 days at RI's setting. Errors loudly on an unknown Event Type rather than
  # quietly using a value that decides FDR admissibility.
  n_rotations <- resolve_screen_rotations(rotation_spec, this_et)

  # Writer + NA filler emit ONE ROW PER FEATURE SET, carrying the
  # feature_set / Event Type / transcript_filter keys the aggregator groups FDR by
  # (`by = .(transcript_filter, Event Type, feature_set)`).
  # Retained rows from an earlier (omnibus/partial) run are re-emitted alongside
  # the newly computed sets. Legacy rows also predate the Event Type /
  # transcript_filter keys the aggregator groups FDR by, so backfill them here.
  # ATOMIC: in resume mode this OVERWRITES a file whose retained rows are the only
  # copy of an already-completed screen. Write to a temp sibling + rename (atomic on
  # the same filesystem) so a crash/kill mid-write can't truncate existing results.
  # NOTE: the temp name MUST keep the .csv.gz extension and compress= must be
  # explicit — fwrite picks its compression from the file extension, so a temp path
  # ending in ".tmp<pid>" silently writes PLAIN CSV that then gets renamed to
  # .csv.gz (fread tolerates it; zcat/gzfile do not).
  write_atomic <- function(dt, path) {
    tmp <- paste0(path, ".tmp", Sys.getpid(), ".csv.gz")
    data.table::fwrite(dt, tmp, compress = "gzip")
    if (!file.rename(tmp, path)) {
      unlink(tmp)
      stop("Failed to rename ", tmp, " -> ", path)
    }
  }
  write_rows <- function(dt) {
    if (!is.null(existing_rows)) {
      if (!"Event Type" %in% names(existing_rows)) {
        existing_rows[, `Event Type` := this_et]
      }
      if (!"transcript_filter" %in% names(existing_rows)) {
        existing_rows[, transcript_filter := this_tf]
      }
      dt <- data.table::rbindlist(list(existing_rows, dt), fill = TRUE)
    }
    write_atomic(dt, out_file)
  }
  na_rows <- function(n_samples = NA_integer_, note = NA_character_) {
    data.table::rbindlist(lapply(feature_sets, function(fs) {
      list(
        ID = this_id,
        feature_set = fs,
        `Event Type` = this_et,
        transcript_filter = this_tf,
        n_samples = n_samples,
        n_features = NA_integer_,
        R_used = 0L,
        n_rotations_requested = NA_integer_,
        n_eligible_controls = NA_integer_,
        lambda = NA_real_,
        screen_df = NA_real_,
        screen_R2 = NA_real_,
        screen_R2_bounded = NA_real_,
        screen_CCC = NA_real_,
        null_R2_mean = NA_real_,
        null_R2_sd = NA_real_,
        null_CCC_mean = NA_real_,
        p_emp = NA_real_,
        effect = NA_real_,
        effect_bounded = NA_real_,
        effect_ccc = NA_real_,
        max_abs_oof = NA_real_,
        max_abs_z = NA_real_,
        frac_z_gt10 = NA_real_,
        lambda_min = NA_real_,
        note = note,
        stat_version = SCREEN_STAT_VERSION
      )
    }))
  }

  # --- load this event's feature table -----------------------------------
  feature_data <- data.table::fread(
    file.path(feature_table_dir, paste0("feature_table_", this_id, ".csv.gz"))
  )
  feature_data <- feature_data[!is.na(PSI)]
  if (nrow(feature_data) < 6L) {
    write_rows(na_rows(nrow(feature_data), "too_few_samples"))
    quit(save = "no", status = 0L)
  }

  # --- CV supergroups (leave-one-ontology-group-out) ---------------------
  groups <- resolve_supergroup_folds(feature_data[[grouping_col]], nfolds)
  if (is.null(groups)) {
    write_rows(na_rows(nrow(feature_data), "collapsed_to_one_group"))
    quit(save = "no", status = 0L)
  }

  # --- confounds (Z), shared across all feature sets --------------------
  parts <- screen_partition_columns(names(feature_data), grouping_col)
  this_uuids <- feature_data[["uuid"]]
  confound_df <- as.data.frame(
    feature_data[, parts$confound_cols, with = FALSE]
  )
  prep <- prep_event(feature_data[["PSI"]], confound_df, groups)

  # long/short/local column sets for a feature table, given its event id.
  # Mirrors 09zz's build_explanatory_vars: long = all epigenetic X; short = drop
  # the far-vicinity chromHMM windows (keep the near ones from chromhmm_hits_smaller);
  # local = drop ALL chromHMM (exon-proximal marks / DNAm / splice-site / RBP only).
  feature_set_columns <- function(x_cols, ev_id) {
    chromhmm <- x_cols[grepl("chromhmm", x_cols, fixed = TRUE)]
    smaller_ids <- chromhmm_hits_smaller[
      chromhmm_hits_smaller[, "queryHits"] == which(keep_rows_manual == ev_id),
      "subjectHits"
    ]
    far <- chromhmm[Reduce(
      `&`,
      lapply(sprintf("chromhmm_%d", smaller_ids), function(s) {
        !endsWith(chromhmm, s)
      }),
      rep(TRUE, length(chromhmm))
    )]
    list(
      long = x_cols,
      short = base::setdiff(x_cols, far),
      local = base::setdiff(x_cols, chromhmm)
    )
  }
  ridge_R2 <- function(cols, dt) {
    if (length(cols) == 0L) {
      return(list(
        R2 = NA_real_,
        CCC = NA_real_,
        lambda = NA_real_,
        df = NA_real_,
        n_features = 0L
      ))
    }
    ridge_screen_stat(prep, as.matrix(dt[, cols, with = FALSE]), groups)
  }
  this_fs_cols <- feature_set_columns(parts$x_cols, this_id)

  # --- matched controls (same selection as 09zz's feature-rotation) ------
  # PSI-INDEPENDENT as of 2026-08-04 (SCREEN_STAT_VERSION 4). The candidate universe
  # is every modelable event -- psi_table's columns -- with NO condition on the
  # control's own PSI. Eligibility is now the matching criteria alone (below).
  #
  # WHY. control_fs_R2 uses only the control's X, scored against THIS event's
  # PSI/confounds/folds (`prep` is fixed), and PSI/IJC/SJC are blocked out of x_cols
  # entirely -- so the control's PSI is never read. Two clauses used to gate on it
  # anyway, both computed on `psi_table[this_uuids, ]`:
  #
  #   colSums(is.na(...)) == 0   the control had to have observed PSI on EVERY focal
  #                              sample. Pure artifact of v1 feature tables being
  #                              built only over rows where the control's own PSI was
  #                              observed: eligibility tracked the control's RNA
  #                              coverage, not its suitability. Measured cost: RI's
  #                              usable pool ~92 against a matched-criteria ceiling
  #                              of 791-887, which is the difference between clearing
  #                              FDR and not. FEATURE_TABLE_VERSION 2 builds every
  #                              table over the full cohort, so it is now moot.
  #
  #   sd(...) > 0                the control's PSI had to VARY across the focal
  #                              samples. Also PSI-derived, and redundant: 09-1
  #                              already admits only events with global sd(PSI) > 0
  #                              (its ids_to_build filter), and psi_table is built
  #                              from that same set, so every candidate column is
  #                              non-constant by construction. Keeping it would ALSO
  #                              be a live bug once the coverage clause goes: sd of a
  #                              column with <2 non-NA values is NA, and `TRUE & NA`
  #                              is NA, so NA would enter the logical index and then
  #                              other_ids. The coverage clause is what masked that
  #                              (`FALSE & NA` is FALSE).
  #
  # Exchangeability is preserved exactly, not traded away: every control is still
  # scored on the focal event's exact `this_uuids`, so `ss_tot` is identical across
  # controls. (That caveat belongs to the DIFFERENT, weaker ">=90% overlap"
  # relaxation, where controls would be scored on differing sample subsets.)
  #
  # To restore the old behaviour, re-intersect other_ids with the two clauses above;
  # note that doing so returns RI to a ~92-control pool and to failing FDR.
  other_ids <- as.integer(colnames(psi_table))
  this_event <- event_dt[ID == this_id]
  other_ids <- other_ids[
    other_ids != this_id &
      other_ids %in% event_dt[`Event Type` == this_event$`Event Type`, ID] &
      other_ids %in% event_dt[seqnames != this_event$seqnames, ID] &
      other_ids %in% event_dt[Variability == this_event$Variability, ID] &
      other_ids %in%
        event_dt[transcript_filter == this_event$transcript_filter, ID]
  ]
  stopifnot(!anyNA(other_ids))
  rm(psi_table)
  gc()

  # --- feature-rotation null, PER FEATURE SET. Each matched control contributes
  #     one R² per feature set, computed from the CONTROL's OWN long/short/local
  #     columns (same as 09zz), on THIS event's PSI/confounds/groups (prep fixed). --
  # Returns BOTH R2 and CCC per feature set, as `R2_<fs>` / `CCC_<fs>`. CCC used to be
  # skipped here on the grounds that a second per-control statistic would double the null
  # cost -- that was wrong: ridge_screen_stat computes CCC from the same fitted model, in
  # the same call, so carrying it out is free. Having a real null CCC is what makes a
  # CCC-based effect size possible, and CCC is the bounded, stratum-stable quantity
  # (0 sign flips across |R2| strata on real data, vs 1 for every R2-derived arm).
  null_names <- c(paste0("R2_", feature_sets), paste0("CCC_", feature_sets))
  control_fs_R2 <- function(cid) {
    na <- setNames(rep(NA_real_, length(null_names)), null_names)
    cdt <- tryCatch(
      data.table::fread(
        file.path(feature_table_dir, paste0("feature_table_", cid, ".csv.gz"))
      ),
      error = function(e) NULL
    )
    if (is.null(cdt)) {
      return(na)
    }
    cdt <- cdt[uuid %in% this_uuids]
    cdt <- cdt[match(this_uuids, uuid)]
    if (anyNA(cdt$uuid) || !all(cdt$uuid == this_uuids)) {
      return(na)
    }
    cparts <- screen_partition_columns(names(cdt), grouping_col)
    cfs <- feature_set_columns(cparts$x_cols, cid)
    res <- lapply(feature_sets, function(fs) ridge_R2(cfs[[fs]], cdt))
    setNames(
      c(
        vapply(res, function(z) z$R2, numeric(1)),
        vapply(res, function(z) z$CCC, numeric(1))
      ),
      null_names
    )
  }

  n_use <- min(n_rotations, length(other_ids))
  if (n_use > 0L) {
    set.seed(seed_base + this_id)
    sampled <- sample(other_ids, n_use)
    null_list <- pbmcapply::pbmclapply(sampled, control_fs_R2, mc.cores = cores)
    null_mat <- do.call(rbind, null_list) # n_use × feature_sets (named cols)
  } else {
    sampled <- integer(0)
    null_mat <- matrix(
      numeric(0),
      nrow = 0L,
      ncol = length(null_names),
      dimnames = list(NULL, null_names)
    )
  }

  # --- per feature set: real stat + empirical p + effect ----------------
  rows <- list()
  null_side <- list()
  for (fs in feature_sets) {
    real <- ridge_R2(this_fs_cols[[fs]], feature_data)
    nR2 <- null_mat[, paste0("R2_", fs)]
    nR2 <- nR2[is.finite(nR2)]
    nCCC <- null_mat[, paste0("CCC_", fs)]
    nCCC <- nCCC[is.finite(nCCC)]
    R_used <- length(nR2)
    ok <- is.finite(real$R2) && R_used > 0L
    rows[[fs]] <- list(
      ID = this_id,
      feature_set = fs,
      `Event Type` = this_et,
      transcript_filter = this_tf,
      n_samples = nrow(feature_data),
      n_features = real$n_features,
      # R_used counts controls that actually returned a FINITE statistic. The two
      # columns after it are what make a shortfall diagnosable instead of just
      # visible: n_eligible_controls is the matched pool before the cap,
      # n_rotations_requested is the cap. R_used << eligible means controls are
      # failing inside control_fs_R2; eligible < requested means the pool, not the
      # cap, is the binding constraint (the whole question for RI).
      R_used = R_used,
      n_rotations_requested = as.integer(n_rotations),
      n_eligible_controls = length(other_ids),
      lambda = real$lambda,
      screen_df = real$df,
      screen_R2 = real$R2,
      # bounded companion, R2/(1+|R2|); strictly increasing, so it cannot change any
      # ranking or p. For tables and effect sizes only.
      screen_R2_bounded = real$R2_bounded,
      screen_CCC = real$CCC,
      null_R2_mean = if (R_used > 0L) mean(nR2) else NA_real_,
      null_R2_sd = if (R_used > 0L) stats::sd(nR2) else NA_real_,
      # Real value now (was hard-coded NA): control_fs_R2 returns CCC alongside R2 at no
      # extra cost, since both come out of the same fitted model.
      null_CCC_mean = if (length(nCCC)) mean(nCCC) else NA_real_,
      # p_emp is computed from the RAW R2 and is unchanged. It is rank-based, so the
      # unbounded magnitude is irrelevant to it -- and any strictly monotone rescaling
      # (e.g. screen_R2_bounded) would give a numerically identical p.
      p_emp = if (ok) (1 + sum(nR2 >= real$R2)) / (1 + R_used) else NA_real_,
      # Retained for continuity, but mean-based on an unbounded quantity, so it is
      # outlier-dominated and must NOT be quoted as an effect size (measured: 61%
      # positive for purely numerical reasons, median flipping to -222 in the worst
      # |R2| stratum). The `effect > 0` clause has been dropped from the hit rule; the
      # one-sided p_emp already encodes direction.
      effect = if (ok) real$R2 - mean(nR2) else NA_real_,
      # Usable effect sizes: both are differences of BOUNDED quantities, so a single
      # extreme null cannot dominate the mean.
      effect_bounded = if (ok) {
        real$R2_bounded - mean(.squash_r2(nR2))
      } else {
        NA_real_
      },
      effect_ccc = if (is.finite(real$CCC) && length(nCCC)) {
        real$CCC - mean(nCCC)
      } else {
        NA_real_
      },
      # B3 out-of-support diagnostics: the cross-cell-type extrapolation is the finding,
      # so it is recorded per event rather than smoothed away.
      max_abs_oof = real$max_abs_oof,
      max_abs_z = real$max_abs_z,
      frac_z_gt10 = real$frac_z_gt10,
      lambda_min = real$lambda_min,
      # prep_event's per-fold notes (rank deficiency, a confound constant or all-NA in a
      # fold's train rows, a factor level present only in held-out rows). Previously an
      # event voided this way carried note = NA and was indistinguishable from a clean one.
      note = if (length(prep$notes)) {
        paste(prep$notes, collapse = ";")
      } else {
        NA_character_
      },
      # Provenance stamp: which version of the statistic produced this row. The
      # resume gate above discards anything not matching SCREEN_STAT_VERSION.
      stat_version = SCREEN_STAT_VERSION
    )
    # Per-event null sidecar (feature_set-tagged) for 09-2's null histograms AND
    # for the aggregator's floor-free sensitivity p -- stamped for the same reason.
    if (nrow(null_mat) > 0L) {
      null_side[[fs]] <- data.table::data.table(
        ID = this_id,
        feature_set = fs,
        control_id = sampled,
        null_R2 = null_mat[, paste0("R2_", fs)],
        null_CCC = null_mat[, paste0("CCC_", fs)],
        stat_version = SCREEN_STAT_VERSION
      )
    }
  }
  write_rows(data.table::rbindlist(rows))
  if (length(null_side) > 0L) {
    null_out <- data.table::rbindlist(null_side)
    if (!is.null(existing_null)) {
      null_out <- data.table::rbindlist(
        list(existing_null, null_out),
        fill = TRUE
      )
    }
    write_atomic(null_out, null_file)
  }
  message(
    "Done screen: ",
    this_id,
    " (",
    this_et,
    ") — ",
    paste(
      vapply(
        feature_sets,
        function(fs) {
          sprintf(
            "%s R2=%.3f p=%.3g",
            fs,
            rows[[fs]]$screen_R2,
            rows[[fs]]$p_emp
          )
        },
        character(1)
      ),
      collapse = " | "
    )
  )
}
