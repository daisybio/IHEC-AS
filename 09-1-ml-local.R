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
already_computed_ids <- as.integer(
  sub(
    "^event_summary_(.+)\\.csv\\.gz$",
    "\\1",
    basename(list.files(event_dir, pattern = "^event_summary_.*\\.csv\\.gz$"))
  )
)
# already_computed_ids <- as.integer(gsub(
#   "(_robust)*.rds$",
#   "",
#   basename(list.files(event_dir, pattern = ".rds"))
# ))
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
# first only make biotype_filtered ids
ids_to_build <- ids_to_build[
  ids_to_build %in% event_dt[transcript_filter == "biotype_filtered", ID]
]
event_dt[ID %in% ids_to_build, table(transcript_filter, `Event Type`)]


# PSI matrix (events × samples) used inside each SLURM job to select
# rotation-control events that match on type, chromosome, variability, filter.
psi_table <- as.matrix(
  dcast(
    aggregated_dt[ID %in% base::intersect(keep_rows_manual, ids_to_build)],
    uuid ~ ID,
    value.var = response
  ),
  rownames = "uuid"
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
    event_dt = event_dt[, .(
      ID,
      `Event Type`,
      seqnames,
      Variability,
      transcript_filter
    )],
    chromhmm_hits_smaller = as.matrix(chromhmm_hits_smaller),
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

dir.create(event_dir, showWarnings = FALSE)
dir.create("event_glmnet_logs", showWarnings = FALSE)

ids_file <- normalizePath(
  file.path("processed_data", "event_glmnet_ids.txt"),
  mustWork = FALSE
)
cfg_file <- normalizePath(
  file.path("processed_data", "event_glmnet_cfg.rds"),
  mustWork = FALSE
)
writeLines(as.character(ids_to_build), ids_file)
saveRDS(.slurm_cfg, cfg_file)

n <- length(ids_to_build)
system(sprintf(
  'sbatch --array=0-%d%%10 "%s" "%s" "%s" "%s"',
  10, #n - 1,
  normalizePath("09-1-ml-local-array.sh"),
  ids_file,
  cfg_file,
  normalizePath(".")
))
message(sprintf("Submitted %d array jobs", n))
