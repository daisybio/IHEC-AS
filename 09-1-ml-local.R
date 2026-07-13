# =============================================================================
# Session data and config
# =============================================================================
# Per-filter (matches 02-2 .. 05): one transcript_filter per invocation, chosen
# by TRANSCRIPT_FILTER (Snakemake `event_models` rule sets it; standalone falls
# back to the primary filter). `aggregating.rda` is GONE (removed by the
# per-filter refactor); every object it used to supply is now read explicitly
# from the per-filter artifacts 02-3/03/05 write.
tf <- Sys.getenv(
  "TRANSCRIPT_FILTER",
  getOption("EpiATLAS_AS_PRIMARY_FILTER", "biotype_filtered")
)

aggregated_dt <- fread(
  sprintf("processed_data/aggregated_dt_filtered_%s.csv.gz", tf),
  stringsAsFactors = TRUE
)
# event-level metadata (ID, Event Type, seqnames, Variability, transcript_filter,
# is_autosome, cluster_representative_annotated, ...) — replaces the old
# aggregating.rda wide event_dt for everything EXCEPT per-sample PSI columns
# (those now live in psi_long_dt, long format, read below).
event_annotations_dt <- fread(
  sprintf("processed_data/event_annotations_dt_%s.csv.gz", tf)
)
psi_long_dt <- fread(sprintf("processed_data/psi_long_dt_%s.csv.gz", tf))
keep_rows_manual <- readRDS(
  sprintf("processed_data/keep_rows_manual_%s.rds", tf)
)
sample_cols <- readRDS(sprintf("processed_data/sample_cols_%s.rds", tf))
# chromHMM-vicinity objects (09-1's own feature engineering, separate from the
# pooled aggregated_dt/splicing_ml route) — see 03's saveRDS comment for why
# these 3 must come from the SAME 03 run (internal index consistency: event_gr
# row position == ID; chromhmm_hits' from() indexes into
# event_gr[keep_rows_manual] positions).
event_gr <- readRDS(sprintf("processed_data/event_gr_%s.rds", tf))
activeChromHMM <- readRDS(sprintf("processed_data/activeChromHMM_%s.rds", tf))
chromhmm_hits <- readRDS(sprintf("processed_data/chromhmm_hits_%s.rds", tf))

response <- "PSI"
grouping_col <- "ontology" # CV fold grouping column
nfolds <- 5L # number of CV folds = ontology supergroups
nrotations <- 10L # negative-control rotations per event
feature_sets <- c("long", "short", "local") # explanatory-set names and dir suffixes

source("09zz-ml-event-glmnet-tidymodels.R")

# =============================================================================
# Output directory
# =============================================================================
# Per-filter: matches the Snakemake event_models rule's declared
# processed_data/event_models/{transcript_filter}/.done sentinel.
event_dir <- file.path("processed_data", "event_models", tf)
dir.create(event_dir, recursive = TRUE, showWarnings = FALSE)
# setNames(
#   file.path("processed_data", paste0("event_models_", feature_sets)),
#   feature_sets
# )
already_computed_ids <- as.integer(
  sub(
    "^(\\d+)_all_metrics\\.csv\\.gz$",
    "\\1",
    basename(list.files(event_dir, pattern = "^\\d+_all_metrics\\.csv\\.gz$"))
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
# file_table.csv.gz (cohort-wide, unchanged) has NO local_file column itself —
# construct it the same way 05-create-aggregated-dt.Rmd does, pointing at the
# per-filter ChIP aggregation output (whose window set — aggregateOver_{tf}.bed
# — already includes activeChromHMM[chromhmm_in_vicinity], so these .tab.gz
# files carry the chromhmm_* rows this script needs, same files 05 uses).
file_table <- fread("processed_data/file_table.csv.gz")
chip_agg_dir <- file.path(sample_dt_dir, sprintf("chip_agg_%s", tf))
file_table[,
  local_file := file.path(chip_agg_dir, paste0(basename(file_path), ".tab.gz"))
]
setkey(file_table, "local_file")

# WGBS chromHMM coverage — cached as .fst for fast re-read (per-filter)
wgbs_chromhmm_file <- file.path(sample_dt_dir, sprintf("WGBS_chromhmm_%s.fst", tf))
if (file.exists(wgbs_chromhmm_file)) {
  wgbs_chromhmm <- fst::read_fst(wgbs_chromhmm_file, as.data.table = TRUE)
} else {
  # WGBS_agg_{tf}.csv.gz columns (fixed order): ID, ihec, score, n, name.
  # col 1 (ID) dropped, cols 2:5 kept, renamed to IHEC (capital — required by
  # the CJ(name=..., IHEC=...) join in the Phase 1 feature-table chunk below;
  # the on-disk header uses lowercase "ihec", unlike the old bare WGBS_agg.csv.gz
  # this replaced). egrep on a header-bearing stream drops the header line (no
  # "chromhmm_" text in it), so select-by-name isn't possible — select by
  # position + hardcode names instead of the old separate header-fetch dance.
  wgbs_chromhmm <- fread(
    cmd = sprintf(
      "zcat %s | egrep 'chromhmm_'",
      file.path(sample_dt_dir, sprintf("WGBS_agg_%s.csv.gz", tf))
    ),
    select = c(2L, 3L, 4L, 5L),
    col.names = c("IHEC", "score", "n", "name"),
    stringsAsFactors = TRUE
  )
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
# Sanity check: events in keep_rows_manual but not aggregated_dt should have
# INSUFFICIENT sample-level data (< 2 non-NA PSI samples), not necessarily ZERO
# — 03's keep_rows_manual (is_autosome & cluster_representative_annotated) does
# not itself require >=2 non-NA samples, but 05 additionally drops events with
# unset Variability (== NA sd_psi == <2 non-NA PSI samples after the 02-3
# read-coverage/VST masking; see 05-create-aggregated-dt.Rmd.md's NA-Variability
# note). A dropped event can have exactly 1 non-NA sample (sd undefined, not
# all-NA), so the old all-NA stopifnot no longer holds — relaxed to <2 here,
# checked against psi_long_dt (same 02-3 masking as aggregated_dt's source).
.dropped_ids <- setdiff(keep_rows_manual, aggregated_dt[, unique(ID)])
.n_non_na <- psi_long_dt[ID %in% .dropped_ids & !is.na(psi), .N, by = ID]
stopifnot(all(.n_non_na$N < 2L))

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
# NB: no "first only make biotype_filtered ids" subset needed anymore —
# aggregated_dt is already scoped to exactly `tf` (per-filter architecture),
# unlike the old combined-filters event_dt this line used to subset from.
event_annotations_dt[ID %in% ids_to_build, table(transcript_filter, `Event Type`)]


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
feature_table_dir <- file.path(
  "processed_data",
  sprintf("event_feature_tables_%s", tf)
)
dir.create(feature_table_dir, showWarnings = FALSE)

chip_files <- file_table[assay_type == "ChIP-Seq", unique(local_file)]
active_chrom_ids <- sort(unique(to(chromhmm_hits[
  from(chromhmm_hits) %in% which(keep_rows_manual %in% ids_to_build)
])))

chip_matrix_cache <- file.path(
  "processed_data",
  sprintf("chip_matrix_%s.rds", tf)
)
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

session_rds <- file.path(
  "processed_data",
  sprintf("session_09_1_ml_local_%s.rds", tf)
)
saveRDS(
  list(
    psi_table = psi_table,
    # 07 reads sess$event_dt with exactly these 5 columns (ID/Event
    # Type/seqnames/Variability/transcript_filter) — event_annotations_dt has
    # them all under the same names, so 07 itself needs no change here.
    event_dt = event_annotations_dt[, .(
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
  file.path("processed_data", sprintf("event_glmnet_ids_%s.txt", tf)),
  mustWork = FALSE
)
cfg_file <- normalizePath(
  file.path("processed_data", sprintf("event_glmnet_cfg_%s.rds", tf)),
  mustWork = FALSE
)
writeLines(as.character(ids_to_build), ids_file)
saveRDS(.slurm_cfg, cfg_file)

# =============================================================================
# Phase 2 dispatch: local test OR SLURM array
# =============================================================================
# Set n_local_test > 0 to run a small subset locally instead of sbatch.
#   n_outer_cores — how many events run in parallel
#   n_inner_cores — workflows parallelised within each event
#                   (passed as 3rd CLI arg; SLURM always uses 1)
n_local_test <- 0L # set > 0 to bypass sbatch
n_outer_cores <- 5L
n_inner_cores <- 11L

if (n_local_test > 0L) {
  test_ids <- head(ids_to_build, n_local_test)
  message(sprintf(
    "Local test: %d events, %d outer x %d inner cores",
    length(test_ids),
    n_outer_cores,
    n_inner_cores
  ))
  dir.create("event_glmnet_logs", showWarnings = FALSE)
  pbmcapply::pbmclapply(
    test_ids,
    function(id) {
      system2(
        "Rscript",
        args = c(
          shQuote(normalizePath("09zz-ml-event-glmnet-tidymodels.R")),
          shQuote(cfg_file),
          as.character(id),
          as.character(n_inner_cores)
        ),
        wait = TRUE,
        stdout = file.path("event_glmnet_logs", sprintf("local_%d.log", id)),
        stderr = file.path("event_glmnet_logs", sprintf("local_%d.err", id))
      )
    },
    mc.cores = n_outer_cores
  )
  message("Local test runs complete.")
} else {
  n <- length(ids_to_build)
  system(sprintf(
    'sbatch --array=0-%d%%20 "%s" "%s" "%s" "%s"',
    1000, #n - 1L,
    normalizePath("09-1-ml-local-array.sh"),
    ids_file,
    cfg_file,
    normalizePath(".")
  ))
  message(sprintf("Submitted %d array jobs", n))
}
