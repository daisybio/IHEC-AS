# =============================================================================
# Session data and config
# =============================================================================
load("processed_data/aggregating.rda", verbose = TRUE)

aggregated_dt <- fread(
  "processed_data/aggregated_dt_filtered.csv.gz",
  stringsAsFactors = TRUE
)

response <- "PSI"
grouping_col <- "ontology" # CV fold grouping column
nfolds <- 5L # number of CV folds = ontology supergroups
nrotations <- 10L # negative-control rotations per event
feature_sets <- c("long", "short", "local") # explanatory-set names and dir suffixes

source("07-ml-event-glmnet-tidymodels.R")

# =============================================================================
# Output directory
# =============================================================================
event_dir <- file.path("processed_data", "event_models")
# setNames(
#   file.path("processed_data", paste0("event_models_", feature_sets)),
#   feature_sets
# )

already_computed_ids <- as.integer(gsub(
  "(_robust)*.rds$",
  "",
  basename(list.files(event_dir, pattern = ".rds"))
))
# lapply(event_dirs, function(d) {
#   dir.create(d, showWarnings = FALSE)
#   as.integer(gsub(
#     "(_robust)*.rds$",
#     "",
#     basename(list.files(d, pattern = ".rds"))
#   ))
# })

# All three dirs must agree on which events are complete
# stopifnot(identical(
#   already_computed_ids[[1]],
#   Reduce(base::intersect, already_computed_ids)
# ))
# already_computed_ids <- Reduce(base::intersect, already_computed_ids)

# =============================================================================
# Load shared data needed for feature table construction
# =============================================================================
file_table <- fread("processed_data/file_table.csv.gz")
setkey(file_table, "local_file")

# WGBS chromHMM coverage — cached as .fst for fast re-read
wgbs_chromhmm_file <- file.path(sample_dt_dir, "WGBS_chromhmm.fst")
if (file.exists(wgbs_chromhmm_file)) {
  wgbs_chromhmm <- fst::read_fst(wgbs_chromhmm_file, as.data.table = TRUE)
} else {
  col_idx <- c(2L, 3L, 4L, 5L)
  wgbs_chromhmm <- fread(
    cmd = sprintf(
      "zcat %s | egrep 'chromhmm_'",
      file.path(sample_dt_dir, "WGBS_agg.csv.gz")
    ),
    select = col_idx,
    stringsAsFactors = TRUE
  )
  col_names <- strsplit(
    system(
      sprintf(
        "zcat %s | head -n 1",
        file.path(sample_dt_dir, "WGBS_agg.csv.gz")
      ),
      intern = TRUE
    ),
    ","
  )[[1]]
  setnames(wgbs_chromhmm, col_names[col_idx])
  wgbs_chromhmm <- melt(
    wgbs_chromhmm,
    id.vars = c("IHEC", "name"),
    value.name = "score",
    variable.name = "experiment_type",
    measure.vars = c("score", "n")
  )
  wgbs_chromhmm[experiment_type == "score", experiment_type := "DNAm"]
  wgbs_chromhmm[experiment_type == "n", experiment_type := "CpGs"]
  wgbs_chromhmm[, IHEC := sub("\\.[0-9]+$", "", IHEC)]
  fst::write_fst(wgbs_chromhmm, wgbs_chromhmm_file)
}
setkey(wgbs_chromhmm, name, IHEC)

# =============================================================================
# Determine events to process
# =============================================================================
# Events in aggregated_dt that have sufficient samples and non-zero PSI variance.
# Sanity check: events in keep_rows_manual but not aggregated_dt should have no
# sample-level data (they were filtered out upstream for a legitimate reason).
stopifnot(all(is.na(event_dt[
  ID %in% setdiff(keep_rows_manual, aggregated_dt[, unique(ID)]),
  ..sample_cols
])))

ids_to_build <- aggregated_dt[, unique(ID)]
ids_to_build <- ids_to_build[
  ids_to_build %in%
    aggregated_dt[ID %in% ids_to_build, .N, by = ID][N >= minimum_events, ID]
]
ids_to_build <- ids_to_build[
  ids_to_build %in%
    aggregated_dt[ID %in% ids_to_build, .(sd = sd(get(response))), by = ID][
      sd > 0,
      ID
    ]
]
event_dt[ID %in% ids_to_build, table(transcript_filter, `Event Type`)]

# PSI matrix (events × samples) used inside each SLURM job to select
# rotation-control events that match on type, chromosome, variability, filter.
psi_table <- dcast(
  aggregated_dt[ID %in% base::intersect(keep_rows_manual, ids_to_build)],
  uuid ~ ID,
  value.var = response
)

ids_to_build <- setdiff(ids_to_build, already_computed_ids)

# =============================================================================
# Build chip_matrix — rows = chromHMM regions, cols = ChIP-Seq files (cached)
# =============================================================================
# Pre-slicing the full matrix once avoids re-reading ~2,262 files per event.
feature_table_dir <- file.path("processed_data", "event_feature_tables")
dir.create(feature_table_dir, showWarnings = FALSE)

chip_files <- file_table[assay_type == "ChIP-Seq", unique(local_file)]
active_chrom_ids <- sort(unique(to(chromhmm_hits[
  from(chromhmm_hits) %in% which(keep_rows_manual %in% ids_to_build)
])))

chip_matrix_cache <- file.path("processed_data", "chip_matrix.rds")
if (file.exists(chip_matrix_cache)) {
  chip_matrix <- readRDS(chip_matrix_cache)
} else {
  all_cols <- pbmcapply::pbmclapply(seq_along(chip_files), function(i) {
    dt <- fread(
      cmd = paste0("zcat ", chip_files[i], " | grep -P '^chromhmm_'"),
      select = c(1L, 6L),
      col.names = c("name", "mean")
    )
    dt[, id_int := as.integer(sub("chromhmm_", "", name, fixed = TRUE))]
    row_idx <- match(dt$id_int, active_chrom_ids)
    col_vec <- rep(NA_real_, length(active_chrom_ids))
    col_vec[row_idx[!is.na(row_idx)]] <- dt$mean[!is.na(row_idx)]
    col_vec
  })
  chip_matrix <- do.call(cbind, all_cols)
  rm(all_cols)
  dimnames(chip_matrix) <- list(
    sprintf("chromhmm_%d", active_chrom_ids),
    chip_files
  )
  saveRDS(chip_matrix, chip_matrix_cache, compress = FALSE)
}

# chromHMM hits within a narrower window (vicinity/10) — defines the
# "short" feature sets as the subset of regions within this tighter radius.
chromhmm_hits_smaller <- findOverlaps(
  event_gr[keep_rows_manual],
  activeChromHMM,
  maxgap = vicinity / 10,
  ignore.strand = TRUE
)

# =============================================================================
# Phase 1: Build feature tables (local, I/O-bound)
# =============================================================================
# Idempotent — already-built tables are skipped.  Completes before Phase 2.
# Requires all large shared objects above; Phase 2 workers need none of them.
pbmcapply::pbmclapply(ids_to_build, function(id) {
  data.table::setDTthreads(1L)
  feature_table_file <- file.path(
    feature_table_dir,
    paste0("feature_table_", id, ".csv.gz")
  )
  if (file.exists(feature_table_file)) {
    return(invisible(NULL))
  }

  feature_data <- aggregated_dt[ID == id]
  chromhmm_ids <- to(chromhmm_hits[
    from(chromhmm_hits) == which(keep_rows_manual == id)
  ])
  event_files <- file_table[
    assay_type == "ChIP-Seq" &
      epirr_id_without_version %in% feature_data[, unique(IHEC)],
    local_file
  ]

  # Slice pre-built chip_matrix for this event's regions and samples
  mat_sub <- chip_matrix[
    sprintf("chromhmm_%d", chromhmm_ids),
    event_files,
    drop = FALSE
  ]
  chromhmm_data <- melt(
    as.data.table(mat_sub, keep.rownames = "name"),
    id.vars = "name",
    variable.name = "local_file",
    value.name = "score"
  )
  chromhmm_data[,
    IHEC := as.factor(file_table[local_file, epirr_id_without_version])
  ]
  chromhmm_data[,
    experiment_type := as.factor(file_table[local_file, experiment_type])
  ]
  chromhmm_data[, local_file := NULL]

  chromhmm_features <- rbindlist(
    list(
      chromhmm_data,
      wgbs_chromhmm[
        CJ(
          name = sprintf("chromhmm_%d", chromhmm_ids),
          IHEC = feature_data[, unique(IHEC)],
          unique = TRUE
        ),
        nomatch = NULL
      ]
    ),
    use.names = TRUE
  )

  # Pivot to wide: one row per IHEC sample, one col per (mark, region) pair
  chromhmm_features_wide <- dcast(
    chromhmm_features,
    formula = IHEC ~ experiment_type + name,
    value.var = "score",
    sep = ";"
  )
  non_all_na_cols <- chromhmm_features_wide[, sapply(.SD, function(col) {
    !all(is.na(col))
  })]
  chromhmm_features_wide <- chromhmm_features_wide[, ..non_all_na_cols]

  cols_to_add <- names(non_all_na_cols)[non_all_na_cols != "IHEC"]
  feature_data[
    chromhmm_features_wide,
    on = .(IHEC),
    (cols_to_add) := mget(cols_to_add)
  ]
  fwrite(feature_data, feature_table_file)
  invisible(NULL)
})

# =============================================================================
# Phase 2: ML runs — one SLURM job per event, minimal memory footprint
# =============================================================================
# Each worker loads only:
#   • feature_table_{id}.csv.gz  — the event's pre-built wide feature table
#   • session_rds                — psi_table + event_dt + hits objects (~300 MB)
# aggregated_dt, chip_matrix, wgbs_chromhmm, file_table are NOT needed.

session_rds <- file.path("processed_data", "session_09_1_ml_local.rds")
saveRDS(
  list(
    psi_table = psi_table,
    event_dt = event_dt,
    chromhmm_hits_smaller = chromhmm_hits_smaller,
    keep_rows_manual = keep_rows_manual
  ),
  session_rds
)

.slurm_cfg <- list(
  project_dir = normalizePath("."),
  session_rds = normalizePath(session_rds),
  feature_table_dir = normalizePath(feature_table_dir),
  event_dir = normalizePath(event_dir), # setNames(normalizePath(event_dirs), names(event_dirs)),
  feature_sets = feature_sets,
  response = response,
  grouping_col = grouping_col,
  nfolds = nfolds,
  nrotations = nrotations
)

library(future.batchtools)
library(future.apply)
plan(
  batchtools_slurm,
  template = normalizePath("scripts/batchtools.slurm.tmpl"),
  resources = list(
    partition = "shared-cpu",
    memory = "6G",
    ncpus = 1L,
    walltime = "2:00:00",
    conda_env = "ihec-as"
  )
)

## in dev make it sequential for easier debugging; in prod, switch to batchtools_slurm
# plan(multisession, workers = 10L)

event_models <- future_lapply(
  ids_to_build[1:10],
  function(id) {
    cfg <- .slurm_cfg
    setwd(cfg$project_dir)
    source(file.path(cfg$project_dir, "07-ml-event-glmnet-tidymodels.R"))

    sess <- readRDS(cfg$session_rds)
    psi_table <- sess$psi_table
    event_dt <- sess$event_dt
    chromhmm_hits_smaller <- sess$chromhmm_hits_smaller
    keep_rows_manual <- sess$keep_rows_manual
    rm(sess)

    feature_table_dir <- cfg$feature_table_dir
    feature_sets <- cfg$feature_sets
    event_dir <- cfg$event_dir
    response <- cfg$response
    grouping_col <- cfg$grouping_col
    nfolds <- cfg$nfolds
    nrotations <- cfg$nrotations

    tryCatch(
      {
        feature_data <- data.table::fread(file.path(
          feature_table_dir,
          paste0("feature_table_", id, ".csv.gz")
        ))

        # All columns that are not metadata / response are candidate predictors
        explanatory <- names(feature_data)[
          !names(feature_data) %in%
            c(
              "IHEC",
              "ID",
              "Event Type",
              "Variability",
              "gene_id",
              "uuid",
              "transcript_filter",
              "project",
              grouping_col,
              response
            )
        ]
        chromhmm_explanatory <- explanatory[grepl(
          "chromhmm",
          explanatory,
          fixed = TRUE
        )]

        # "local" set: chromHMM regions within the narrow window (vicinity/10)
        # "short" set: removes narrow-window regions, keeping only wider ones
        # "long"  set: all chromHMM regions including distal (full vicinity)
        smaller_chromhmm_ids <- S4Vectors::to(chromhmm_hits_smaller[
          S4Vectors::from(chromhmm_hits_smaller) ==
            which(keep_rows_manual == id)
        ])
        old_chromhmm_explanatory <- chromhmm_explanatory[Reduce(
          `&`,
          lapply(sprintf("chromhmm_%d", smaller_chromhmm_ids), function(suff) {
            !endsWith(chromhmm_explanatory, suff)
          })
        )]

        explanatory_vars <- setNames(
          list(
            explanatory,
            explanatory[!explanatory %in% old_chromhmm_explanatory],
            explanatory[!explanatory %in% chromhmm_explanatory]
          ),
          feature_sets
        )

        # Rotation controls: events on different chromosomes that match this
        # event on type, variability, and transcript filter
        subset_psi_matrix <- as.matrix(
          psi_table[feature_data[, uuid]],
          rownames = "uuid"
        )
        other_ids <- as.integer(colnames(subset_psi_matrix)[
          colSums(is.na(subset_psi_matrix)) == 0 &
            apply(subset_psi_matrix, 2, sd, na.rm = TRUE) > 0
        ])
        this_event <- event_dt[ID == id]
        other_ids <- other_ids[
          other_ids != id &
            other_ids %in%
              event_dt[`Event Type` == this_event$`Event Type`, ID] &
            other_ids %in% event_dt[seqnames != this_event$seqnames, ID] &
            other_ids %in% event_dt[Variability == this_event$Variability, ID] &
            other_ids %in%
              event_dt[transcript_filter == this_event$transcript_filter, ID]
        ]

        this_rotations <- min(nrotations, length(other_ids))
        if (this_rotations < nrotations) {
          warning(sprintf(
            "Not enough rotation controls for %d, using %d",
            id,
            this_rotations
          ))
        }
        stopifnot(
          rownames(subset_psi_matrix) == feature_data[, as.character(uuid)]
        )
        set.seed(id)
        rotated_psis <- subset_psi_matrix[, as.character(sample(
          other_ids,
          this_rotations
        ))]

        wflow_res <- run_event_glmnet(
          this_feature_data = cbind(feature_data, rotated_psis),
          explanatory_vars,
          response,
          rotated_psis,
          grouping_col,
          nfolds,
          seed = id
        )

        saveRDS(
          wflow_res,
          file.path(event_dir, paste0(id, ".rds")),
          compress = TRUE
        )
        # Save every workflow result to its feature-set directory.
        # Filename convention (all under event_dirs[[set_name]]/):
        #   {id}.rds            — primary PSI, non-logit
        #   logit_{id}.rds      — primary PSI, logit scale
        #   rotated_{rot}_{id}.rds — negative-control (rotated PSI response)
        # for (set_name in feature_sets) {
        #   suffix <- paste0("_", set_name, "_glmnet")
        #   set_wfs <- wflow_res$workflow_results[
        #     endsWith(names(wflow_res$workflow_results), suffix)
        #   ]
        #   for (wf_name in names(set_wfs)) {
        #     fname <- if (wf_name == paste0(response, suffix)) {
        #       paste0(id, ".rds")
        #     } else if (wf_name == paste0("logit_", response, suffix)) {
        #       paste0("logit_", id, ".rds")
        #     } else {
        #       rot_id <- sub(suffix, "", wf_name, fixed = TRUE)
        #       paste0("rotated_", rot_id, "_", id, ".rds")
        #     }
        #     saveRDS(
        #       set_wfs[[wf_name]],
        #       file.path(event_dirs[[set_name]], fname)
        #     )
        #   }
        # }
        invisible(NULL)
      },
      error = function(e) e$message
    )
  },
  future.globals = ".slurm_cfg",
  future.seed = TRUE
)

# =============================================================================
# Post-hoc error audit
# =============================================================================
names(event_models) <- ids_to_build
failed <- event_models[!sapply(event_models, is.null)]
message(sprintf("%d / %d events failed", length(failed), length(ids_to_build)))
if (length(failed) > 0L) {
  table(unlist(failed))
}

# sd = 0: PSI constant across samples — expected to fail, verify that is the cause
sd_check <- aggregated_dt[
  ID %in% names(failed)[failed == "not_all_na[response] is not TRUE"],
  .(sd = sd(get(response))),
  by = ID
]
stopifnot(sd_check[, all(sd == 0, na.rm = TRUE)])

# Too few ontology groups for k-fold CV — expected for events with rare tissues
group_check <- aggregated_dt[
  ID %in% names(failed)[startsWith(unlist(failed), "`k` should be less than")],
  .(unique_ontology = uniqueN(get(grouping_col))),
  by = ID
]
stopifnot(group_check[, all(unique_ontology < nfolds)])

# Constant y within folds — PSI has no within-fold variance
constant_check <- aggregated_dt[
  ID %in%
    names(failed)[
      failed == "y is constant; gaussian glmnet fails at standardization step"
    ],
  .(uniqueResponses = unique(get(response))),
  by = .(ID, get(grouping_col))
]
stopifnot(constant_check[,
  .(groupsWithVar = sum(uniqueResponses > 1)),
  by = ID
][, all(groupsWithVar < nfolds)])
