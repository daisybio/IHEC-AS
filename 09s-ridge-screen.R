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

  # --- DETECTION-FLOOR INJECTION (opt-in; unset => ZERO behaviour change) -----
  # Positive control for the screen: replace this event's PSI with a synthetic response
  # built from its OWN epigenetic features at a controlled effect size, keep confounds,
  # folds and control pool byte-identical, and run the real worker. Answers "how small an
  # effect could this screen have found?", which the 0-hit result cannot be interpreted
  # without. See revision/EXECUTION-PLAN.md item 1.
  #
  # EPIATLAS_AS_SCREEN_INJECT="rho=0.05;s=10;mode=rankmap;fs=local;seed=1"
  #   rho  target in-sample signal variance fraction. The two sanity anchors:
  #        rho=1 is pure signal and MUST be detected -- if it is not, the harness is
  #        broken, not the data.
  #        rho=0 is pure noise rank-mapped onto the event's own PSI values, i.e. a
  #        PERMUTATION of PSI across samples. It is NOT the real screen rerun (the real
  #        PSI-sample pairing is destroyed), so do not expect identical numbers -- expect
  #        a calibrated null: p_emp ~ Uniform and no hits. The "identical to production"
  #        check is simply leaving the env var unset, which this block is gated on.
  #   s    number of non-zero coefficients in the sparse random direction (default 10)
  #   mode Three response shapes, and the differences between them decompose WHY PSI is a
  #        hard response for a linear screen:
  #        "gaussian" leaves the latent signal untransformed -- unbounded and exactly linear
  #          in the features, so it is the ceiling (measured: R2 = 1.000 at rho=1).
  #        "logit" maps it through plogis(signal * logit_scale) -- bounded in (0,1), smooth,
  #          NO boundary mass. Isolates the cost of a monotone bounding transform per se.
  #          No epsilon is needed here because the response is GENERATED on the logit scale
  #          rather than transformed from observed PSI, so it never reaches 0 or 1 exactly
  #          (unlike logit(observed PSI), where the epsilon choice moves the event ranking --
  #          see CLAUDE.md).
  #        "rankmap" (default) maps it onto this event's OWN observed PSI values -- preserves
  #          the real marginal exactly, boundary mass included, still monotone (measured:
  #          R2 = 0.518 at rho=1).
  #        gaussian - logit  = cost of bounding.  logit - rankmap = cost of PSI's specific,
  #        boundary-heavy marginal. Together they say how much detectability the response
  #        distribution itself costs.
  #   logit_scale  multiplier inside plogis for mode=logit (default 2; larger = wider spread)
  #   fs   feature space the signal is drawn from: "local" (default) or "long"
  #   rotations  override the per-Event-Type rotation count for THIS run only. The floor
  #        experiment scores on p_pooled_z, which needs the null's MOMENTS (mean, sd), not the
  #        exceedance depth RI's 818 exists to buy -- so 200 for both event types is adequate
  #        there and ~4x cheaper for RI. Without this, cfg$screen_rotations is always set and
  #        the getOption fallback never fires.
  #   outdir  output directory, used VERBATIM when given, so Snakemake can OWN the path
  #        (rule floor_inject_one) without reimplementing the derived-name scheme in Python.
  #        Omitted => derived from event_dir plus a parameter-encoded subdirectory, as before.
  #        NB with outdir given, uniqueness is the CALLER's responsibility -- Snakemake's
  #        output paths are unique by construction; a human passing one outdir for two rho
  #        values would collide.
  #        REFUSED if it resolves to <event_dir>/screen: without the derived-path guarantee
  #        an injection run could otherwise clobber the real screen and be picked up by the
  #        resume gate as genuine. The parameter-derived subdirectory is always appended
  #        UNDER outdir, so two rho values can never collide inside one outdir either.
  .inject <- local({
    spec <- Sys.getenv("EPIATLAS_AS_SCREEN_INJECT", "")
    if (!nzchar(spec)) {
      return(NULL)
    }
    kv <- strsplit(strsplit(spec, ";", fixed = TRUE)[[1]], "=", fixed = TRUE)
    p <- stats::setNames(vapply(kv, function(x) x[2L], character(1)),
      vapply(kv, function(x) x[1L], character(1)))
    g <- function(k, default) if (k %in% names(p) && nzchar(p[[k]])) p[[k]] else default
    if (!"rho" %in% names(p)) stop("EPIATLAS_AS_SCREEN_INJECT needs rho=")
    rho <- as.numeric(p[["rho"]])
    if (!is.finite(rho) || rho < 0 || rho > 1) stop("inject rho must be in [0, 1]")
    mode <- g("mode", "rankmap")
    if (!mode %in% c("rankmap", "gaussian", "logit")) {
      stop("inject mode must be rankmap|gaussian|logit")
    }
    fs <- g("fs", "local")
    if (!fs %in% c("local", "long")) stop("inject fs must be local|long")
    list(rho = rho, s = as.integer(g("s", "10")), mode = mode, fs = fs,
      seed = as.integer(g("seed", "1")),
      logit_scale = as.numeric(g("logit_scale", "2")),
      outdir = g("outdir", ""),
      rotations = {
        r <- g("rotations", "")
        if (nzchar(r)) as.integer(r) else NA_integer_
      })
  })

  # An injection run MUST NOT write into the production screen directory -- out_file is
  # screen/<id>_screen.csv.gz and would silently clobber the real result (and be picked up
  # by the resume gate below as if it were real). Redirect before read_existing() runs.
  screen_dir <- if (is.null(.inject)) {
    file.path(event_dir, "screen")
  } else {
    # outdir given => VERBATIM (Snakemake owns the path). Absent => event_dir + derived name.
    if (nzchar(.inject$outdir)) {
      .inject$outdir
    } else {
    file.path(event_dir, sprintf(
      "screen_inject_rho%s_s%d_%s_%s%s",
      sub("\\.", "p", format(.inject$rho, trim = TRUE)),
      .inject$s, .inject$mode, .inject$fs,
      if (identical(.inject$mode, "logit")) {
        paste0("_ls", sub("\\.", "p", format(.inject$logit_scale, trim = TRUE)))
      } else {
        ""
      }
    ))
    }
  }
  # HARD GUARD, applied to the FINAL path whichever branch produced it.
  # normalizePath(mustWork = FALSE) so a not-yet-created dir still compares, and so a
  # relative path like ./a/../a/screen is caught as well as the literal one.
  if (!is.null(.inject) && identical(
    normalizePath(screen_dir, mustWork = FALSE),
    normalizePath(file.path(event_dir, "screen"), mustWork = FALSE)
  )) {
    stop("refusing to write injected results into the real screen directory")
  }
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
  # to push the empirical-p floor under BH's bound. SE stays at 200 on COST grounds,
  # not because 200 suffices (see 09-1's cfg comment): SE needs k >= 1,611 at R=200,
  # and escalation is deferred to Stage 2, which lifts only floor-tied events. Raising
  # every SE event instead would cost ~11 days against 1.5. Errors loudly on an unknown
  # Event Type rather than quietly using a value that decides FDR admissibility.
  n_rotations <- if (!is.null(.inject) && !is.na(.inject$rotations)) {
    .inject$rotations
  } else {
    resolve_screen_rotations(rotation_spec, this_et)
  }

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
    # provenance for injection runs -- without this an injected result is
    # indistinguishable from a real one once it leaves this directory
    if (!is.null(.inject)) {
      dt[, `:=`(
        inject_rho = .inject$rho, inject_s = .inject$s,
        inject_mode = .inject$mode, inject_fs = .inject$fs,
        inject_seed = .inject$seed, inject_logit_scale = .inject$logit_scale
      )]
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
  # resolve_supergroup_folds() guarantees >= 2 groups but says nothing about their SIZES.
  # If one supergroup holds all but one sample, that fold's TRAIN partition is a single
  # row: no standardisation, no confound fit, no kernel -- statistically meaningless, and
  # it used to crash (stats::sd() of length 1 is NA, and .build_fold_z's
  # `if (sd <= 1e-12)` on NA is a hard error). That killed event 28871 and, with
  # restart-times: 0, took the whole 34k-job run down with it. Such an event is
  # unscreenable, so record it as such and exit 0 rather than aborting the run.
  # NB this changes no number for any event that already succeeded (those all have
  # train >= 2 in every fold), so it needs no SCREEN_STAT_VERSION bump.
  .min_train <- length(groups) - max(table(groups))
  if (.min_train < 2L) {
    write_rows(na_rows(nrow(feature_data), sprintf(
      "fold_train_too_small:%d", .min_train
    )))
    quit(save = "no", status = 0L)
  }

  # --- confounds (Z), shared across all feature sets --------------------
  parts <- screen_partition_columns(names(feature_data), grouping_col)

  # Replace PSI with the synthetic response BEFORE prep_event, so every downstream step --
  # confound projection, folds, kernel, control rotations, p_emp -- sees it exactly as it
  # would see real PSI. Placed after `parts` because it needs x_cols, and before `prep`
  # because prep is built from PSI.
  if (!is.null(.inject)) {
    .xa <- parts$x_cols
    .chm <- .xa[grepl("chromhmm", .xa, fixed = TRUE)]
    .cols <- if (.inject$fs == "long") .xa else base::setdiff(.xa, .chm)
    .X <- as.matrix(feature_data[, .cols, with = FALSE])
    storage.mode(.X) <- "double"
    # same standardisation the screen itself applies, so the injected direction lives on
    # the scale the ridge will actually see
    .sd0 <- apply(.X, 2L, stats::sd, na.rm = TRUE)
    .X <- .X[, is.finite(.sd0) & .sd0 > 1e-8, drop = FALSE]
    if (ncol(.X) < 1L) stop("inject: no usable feature columns for event ", this_id)
    .mu <- colMeans(.X, na.rm = TRUE)
    .mu[!is.finite(.mu)] <- 0
    .X <- sweep(.X, 2L, .mu, "-")
    .X[!is.finite(.X)] <- 0
    .sdv <- sqrt(colSums(.X^2) / max(1L, nrow(.X) - 1L))
    .sdv[!is.finite(.sdv) | .sdv <= 1e-8] <- 1
    .X <- sweep(.X, 2L, .sdv, "/")

    set.seed(seed_base + this_id + .inject$seed +
      as.integer(round(1e6 * .inject$rho)))
    .s <- min(.inject$s, ncol(.X))
    .beta <- numeric(ncol(.X))
    .beta[sample.int(ncol(.X), .s)] <- sample(c(-1, 1), .s, replace = TRUE)
    .sig <- as.numeric(.X %*% .beta)
    .sdsig <- stats::sd(.sig)
    if (!is.finite(.sdsig) || .sdsig <= 1e-12) {
      stop("inject: degenerate signal direction for event ", this_id)
    }
    .sig <- (.sig - mean(.sig)) / .sdsig
    .ylat <- sqrt(.inject$rho) * .sig +
      sqrt(1 - .inject$rho) * stats::rnorm(length(.sig))
    feature_data[, PSI := switch(.inject$mode,
      gaussian = .ylat,
      # bounded in (0,1), smooth, no boundary mass -- isolates the cost of bounding alone
      logit = stats::plogis(.ylat * .inject$logit_scale),
      # rank-map onto this event's OWN observed PSI values: exact marginal, monotone
      rankmap = sort(feature_data[["PSI"]])[rank(.ylat, ties.method = "first")]
    )]
    rm(.X)
  }

  this_uuids <- feature_data[["uuid"]]
  confound_df <- as.data.frame(
    feature_data[, parts$confound_cols, with = FALSE]
  )
  prep <- prep_event(feature_data[["PSI"]], confound_df, groups)
  # Confound-free baseline prep (Z = intercept only), for two diagnostics added
  # 2026-08-20: r2_confounds (how much do confounds alone explain, vs the
  # intercept-only null?) and r2_x_alone (how much does X alone explain, with
  # NO confound residualisation at all?) -- both real-event-only, never
  # computed for the null rotation (cost would scale with R_used for no
  # diagnostic benefit; the null's job is the beyond-confounds comparison
  # already in screen_R2). Cheap: prep_event() itself is a few qr() calls, and
  # ridge_screen_stat(prep_noconf, ...) below is 3 extra kernel fits (one per
  # feature set) against an event that already pays for 1 + R_used x 3 fits.
  # See revision/file-changes/09-confound-projection-instability.md for why
  # these numbers matter (found confounds-alone underperforms the null on the
  # majority of real events -- an OLS-instability artifact, not real
  # anti-signal, but worth surfacing per-event rather than only in a scratch
  # script). Bounded via the same R2/(1+|R2|) squash as screen_R2_bounded
  # (v3 decision, see project_tier1_v3_scoring_decision.md) -- same
  # non-injectivity-avoidance reasoning applies to these two new numbers.
  prep_noconf <- prep_event(feature_data[["PSI"]], NULL, groups)

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
  # X alone, no confound residualisation at all -- same shape as ridge_R2 but
  # closes over prep_noconf instead of prep. Real-event-only (see comment at
  # prep_noconf's definition above).
  ridge_R2_noconf <- function(cols, dt) {
    if (length(cols) == 0L) {
      return(list(
        R2 = NA_real_,
        CCC = NA_real_,
        lambda = NA_real_,
        df = NA_real_,
        n_features = 0L
      ))
    }
    ridge_screen_stat(prep_noconf, as.matrix(dt[, cols, with = FALSE]), groups)
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
    # ONE data.table thread per fork. The parent set setDTthreads(cores) for its own
    # reads, but this function runs inside pbmclapply's mc.cores forks, so inheriting
    # that gives cores x cores threads on a `cores`-CPU allocation (4 x 4 = 16 on the
    # analysis rule). 09-1's Phase-1 worker already does this; the screen did not.
    # Measured: BLAS thread count makes NO difference to the dominant `long` ridge
    # (0.53 s at 1 thread vs 0.55 s at 8), so serialising here costs nothing and only
    # removes contention. BLAS/OMP are capped to 1 in the Snakefile rule for the same
    # reason -- OpenBLAS otherwise sizes itself to the NODE's core count (80), not the
    # cgroup's, inside every fork.
    data.table::setDTthreads(1L)
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
  # r2_confounds: how much do confounds (Z) alone explain, vs the
  # intercept-only null? Depends only on y/Z/groups, not on X or feature_set --
  # computed once outside the loop, same value stored on every feature_set's
  # row (cheap, matches the existing convention of repeating event-level
  # columns like ID/Event Type across feature_set rows).
  r2_confounds <- if (is.finite(prep_noconf$ss_tot) && prep_noconf$ss_tot > 0) {
    1 - prep$ss_tot / prep_noconf$ss_tot
  } else {
    NA_real_
  }
  r2_confounds_bounded <- .squash_r2(r2_confounds)
  for (fs in feature_sets) {
    real <- ridge_R2(this_fs_cols[[fs]], feature_data)
    real_noconf <- ridge_R2_noconf(this_fs_cols[[fs]], feature_data)
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
      # Added 2026-08-20 (v5): how much do confounds alone, and X alone with NO
      # confound residualisation, explain vs the intercept-only null?
      # r2_confounds is the same value across every feature_set row for this
      # event (see its computation above the loop). Both bounded via the same
      # R2/(1+|R2|) squash as screen_R2_bounded, same reasoning (strictly
      # increasing, so p_emp-style rank comparisons on these would be
      # unaffected -- though neither of these two feeds p_emp/q; they are
      # diagnostic-only, sitting alongside screen_R2/screen_R2_bounded, not
      # replacing them). See revision/file-changes/09-confound-projection-instability.md.
      r2_confounds = r2_confounds,
      r2_confounds_bounded = r2_confounds_bounded,
      r2_x_alone = real_noconf$R2,
      r2_x_alone_bounded = real_noconf$R2_bounded,
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
