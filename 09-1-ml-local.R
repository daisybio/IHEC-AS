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
# RBP Step 4 (§3.4): wide per-RBP expression (one rbp_<name> col per POSTAR3 RBP)
# + the per-event nearby-RBP map. Used to attach each event's binding-site-local
# RBP columns to its feature table (Phase 1 below). (The 3-col rbp_score_* / rbp_n
# aggregate already rides in via aggregated_dt from 05 §4.12f.)
rbp_wide <- fread(sprintf("processed_data/rbp_wide_expression_%s.csv.gz", tf))
rbp_per_event <- readRDS(sprintf("processed_data/rbp_per_event_%s.rds", tf))

response <- "PSI"
grouping_col <- "ontology" # CV fold grouping column
nfolds <- 5L # number of CV folds = ontology supergroups
# Tier-2 elastic-net no longer runs its own rotation null — significance is now
# owned by the Tier-1 fit-free ridge screen (09s-ridge-screen.R). 09zz fits the
# real PSI only, so nrotations = 0 here (belt-and-suspenders alongside 09zz's
# stripped rotation loop). Screen rotation count is separate (cfg$screen_rotations).
nrotations <- 0L
feature_sets <- c("long", "short", "local") # explanatory-set names and dir suffixes

source("09zz-ml-event-glmnet-tidymodels.R")
# FEATURE_TABLE_VERSION / feature_table_dir_for() live next to SCREEN_STAT_VERSION,
# so the build stamp and the statistic stamp are maintained side by side. Safe to
# source here: 09-ml-shared.R is pure definitions with no top-level side effects.
source("09-ml-shared.R")
tidymodels::tidymodels_prefer()

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
wgbs_chromhmm_file <- file.path(
  sample_dt_dir,
  sprintf("WGBS_chromhmm_%s.fst", tf)
)
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
.dropped_ids <- base::setdiff(keep_rows_manual, aggregated_dt[, unique(ID)])
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
event_annotations_dt[
  ID %in% ids_to_build,
  table(transcript_filter, `Event Type`)
]


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

# Full modelable id set (min-samples + non-constant-PSI filtered) — the Tier-1
# ridge screen runs on ALL of these; written to disk below for 09s-dispatch.R.
# (Kept before the already-computed setdiff so the screen covers every event.)
all_modelable_ids <- ids_to_build

ids_to_build <- base::setdiff(ids_to_build, already_computed_ids)

# =============================================================================
# Build chip_matrix — rows = chromHMM regions, cols = ChIP-Seq files (cached)
# =============================================================================
# Pre-slicing the full matrix once avoids re-reading ~2,262 files per event.
# Version-suffixed (FEATURE_TABLE_VERSION, 09-ml-shared.R): the Phase-1 loop skips
# any event whose table already exists, so writing v2 tables into the v1 directory
# would silently skip all ~34k of them and exit 0. A fresh directory also keeps the
# v1 tables readable for side-by-side verification instead of deleting 72 GB.
feature_table_dir <- feature_table_dir_for(tf)
dir.create(feature_table_dir, recursive = TRUE, showWarnings = FALSE)
message(sprintf(
  "Feature tables: v%d -> %s", FEATURE_TABLE_VERSION, feature_table_dir
))

# =============================================================================
# PSI-independent row set for the feature tables (FEATURE_TABLE_VERSION 2)
# =============================================================================
# WHY. Tier-1's matched-control null takes a CONTROL event's epigenetic features and
# scores them against the FOCAL event's PSI/confounds/folds (`prep` is fixed; see
# 09s-ridge-screen.R's control_fs_R2). The control's own PSI is never read -- it is
# blocked out of x_cols entirely. But a control could only ever be used if it had a
# row for every focal sample, and at v1 a feature table was built from
# `aggregated_dt_filtered[ID == id]`, i.e. ONLY the samples where that event's own
# PSI was observed (mean 263 of 415 cohort-wide). So control eligibility was gated on
# a quantity the statistic never uses -- an artifact of the row set, not a design
# choice. It bit RI hardest: RI's usable pool measured ~92 against a matched-criteria
# ceiling of ~818, which is the difference between clearing FDR and not.
#
# FIX. Build every table over the FULL cohort sample set. PSI/IJC/SJC stay NA where
# unobserved (that is the point -- they are blocked from x_cols and the focal fit
# re-applies `!is.na(PSI)` at load in 09s-ridge-screen.R, so the FOCAL statistic is
# untouched); every epigenetic and expression feature is filled for all samples.
#
# HOW, and why it is safe. Observed rows are taken VERBATIM from aggregated_dt --
# not reassembled -- and only the previously-missing rows are built from the
# PSI-independent sources. A join bug in the assembly can therefore degrade a null,
# but it structurally cannot perturb the focal statistic, which reads only the
# verbatim rows. Sources are chosen for the same reason: every per-uuid and
# per-event-constant column is lifted from aggregated_dt itself, so its values and
# its factor levels agree with the observed rows by construction rather than by a
# join that has to be checked.
all_uuids <- sample_cols # 415 for biotype_filtered (03's PSI-matrix sample columns)

# --- classify every aggregated_dt column by the key it is constant over ---------
# classify_feature_columns() (09-ml-shared.R) errors on any column it cannot place,
# so a new column added by 05 stops the build instead of being silently left NA on
# every added row.
feature_cols <- classify_feature_columns(names(aggregated_dt))
message(sprintf(
  paste0(
    "Feature-table column classes: %d key / %d per-(ID,IHEC) / %d per-uuid / ",
    "%d per-(uuid,gene) / %d per-(ID,uuid) / %d per-event"
  ),
  length(feature_cols$key), length(feature_cols$per_epigenome),
  length(feature_cols$per_uuid), length(feature_cols$per_gene),
  length(feature_cols$per_event_sample), length(feature_cols$per_event)
))

# --- per-uuid covariates, taken from aggregated_dt itself -----------------------
# Verified 1:1 over all 415 uuids, so this is the complete cohort covariate table
# without touching file_table/metadata/qc_flag_covariates again -- and therefore
# without any chance of disagreeing with the observed rows on a value or a factor
# level. anyDuplicated() would catch a uuid whose covariates were not in fact
# constant (e.g. a qc_flag_count join gone wrong upstream).
sample_covariates <- unique(
  aggregated_dt[, c("uuid", "IHEC", feature_cols$per_uuid), with = FALSE]
)
stopifnot(!anyDuplicated(sample_covariates$uuid))
.no_cov <- base::setdiff(all_uuids, as.character(sample_covariates$uuid))
if (length(.no_cov)) {
  # A uuid with no observed PSI anywhere in this filter has no covariate row and no
  # epigenetic data to attach, so it cannot contribute a usable row to any event.
  message(sprintf(
    "%d of %d cohort uuids have no row in aggregated_dt - excluded from the full row set",
    length(.no_cov), length(all_uuids)
  ))
  all_uuids <- base::intersect(all_uuids, as.character(sample_covariates$uuid))
}
setkey(sample_covariates, uuid)


chip_files <- file_table[assay_type == "ChIP-Seq", unique(local_file)]
# all_modelable_ids, matching the Phase-1 loop below: chip_matrix's rows are the only
# source for each event's chromHMM slice, so a region set narrowed to ids_to_build
# would make that slice fail for any event Tier-2 had already fitted.
active_chrom_ids <- sort(unique(to(chromhmm_hits[
  from(chromhmm_hits) %in% which(keep_rows_manual %in% all_modelable_ids)
])))

chip_matrix_cache <- file.path(
  "processed_data",
  sprintf("chip_matrix_%s.rds", tf)
)
# STALENESS GUARD. This cache is 13.8 GB and was previously reused on file existence
# alone, with no comparison against its inputs -- and chromHMM is ~2,370 of the ~2,460
# control x_cols, so silently reusing a cache older than the ChIP tabs it was built
# from would mean nearly every control feature is stale. The copy on disk at
# implementation time was dated 2026-07-21 against 2026-07-17 data, i.e. already in the
# regime where a production rerun of the ChIP aggregation invalidates it.
#
# Rebuilt (not just warned about) when any input .tab.gz is newer, or when the region
# set it was built for no longer covers what this run needs -- the latter matters now
# that active_chrom_ids follows all_modelable_ids.
.chip_cache_ok <- FALSE
if (file.exists(chip_matrix_cache)) {
  .cache_mtime <- file.mtime(chip_matrix_cache)
  .newest_input <- suppressWarnings(max(file.mtime(chip_files), na.rm = TRUE))
  if (is.finite(.newest_input) && .newest_input > .cache_mtime) {
    message(sprintf(
      "chip_matrix cache (%s) is OLDER than its newest ChIP tab (%s) - rebuilding",
      format(.cache_mtime, "%Y-%m-%d %H:%M"),
      format(.newest_input, "%Y-%m-%d %H:%M")
    ))
  } else {
    chip_matrix <- readRDS(chip_matrix_cache)
    .want_rows <- sprintf("chromhmm_%d", active_chrom_ids)
    .missing_rows <- length(base::setdiff(.want_rows, rownames(chip_matrix)))
    .missing_cols <- length(base::setdiff(chip_files, colnames(chip_matrix)))
    if (.missing_rows > 0L || .missing_cols > 0L) {
      message(sprintf(
        "chip_matrix cache lacks %d needed region(s) and %d file(s) - rebuilding",
        .missing_rows, .missing_cols
      ))
      rm(chip_matrix)
      gc()
    } else {
      .chip_cache_ok <- TRUE
      message("chip_matrix cache reused (newer than all ChIP tabs, covers all regions)")
    }
  }
}
if (.chip_cache_ok) {
  invisible(NULL)
} else {
  # BATCHED FILL, not pbmclapply-over-everything-then-cbind.
  #
  # The old shape was `all_cols <- pbmclapply(seq_along(chip_files), ...)` followed by
  # `do.call(cbind, all_cols)`. The final matrix is 707,997 x 2,430 doubles = 12.8 GB
  # (exactly the size of the cache on disk), so that shape needs 12.8 GB for the list
  # PLUS 12.8 GB for the cbind result -- and much worse in practice, because
  # mclapply/pbmclapply return results through serialization pipes: each of the 16
  # workers holds its own ~860 MB chunk and a serialized copy of it while the parent
  # accumulates the full 12.8 GB. Measured: OOM-killed at 62.3 GB under a 64000 limit
  # (job 6433872) and again at 93.6 GB under 96000 (job 6433873) -- it grows to fill
  # whatever it is given, so raising mem_mb is not the fix.
  #
  # Filling a preallocated matrix one batch at a time bounds the peak at
  # (final matrix) + (one batch), i.e. ~12.8 GB + ~1 GB instead of 26 GB+ of
  # transients. Batches still use all mc.cores, so throughput is unchanged.
  chip_matrix <- matrix(
    NA_real_,
    nrow = length(active_chrom_ids),
    ncol = length(chip_files),
    dimnames = list(sprintf("chromhmm_%d", active_chrom_ids), chip_files)
  )
  .batch <- 200L
  .starts <- seq.int(1L, length(chip_files), by = .batch)
  message(sprintf(
    "Building chip_matrix: %d regions x %d files (%.1f GB), %d batches of <=%d",
    length(active_chrom_ids), length(chip_files),
    length(active_chrom_ids) * length(chip_files) * 8 / 1024^3,
    length(.starts), .batch
  ))
  for (.bi in seq_along(.starts)) {
    .idx <- seq.int(
      .starts[.bi], min(.starts[.bi] + .batch - 1L, length(chip_files))
    )
    .cols <- pbmcapply::pbmclapply(.idx, function(i) {
      # one thread per worker: 16 forks x data.table's own 16 threads would otherwise
      # multiply fread's buffers (the event loop below sets this for the same reason)
      data.table::setDTthreads(1L)
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
    # Diagnose a dead fork HERE, per batch. mclapply returns NULL for a worker the
    # OOM killer took and warns "scheduled cores ... did not deliver results"; the
    # old code only found out at cbind, where it surfaced as "length of 'dimnames'
    # [2] not equal to array extent" -- saying nothing about the real cause.
    .bad <- !vapply(
      .cols,
      function(z) is.numeric(z) && length(z) == length(active_chrom_ids),
      logical(1)
    )
    if (any(.bad)) {
      stop(sprintf(
        paste0(
          "chip_matrix batch %d/%d: %d of %d tabs returned no usable column -- ",
          "almost certainly the OOM killer taking forked workers. Check MaxRSS ",
          "(sacct -j <id>) against this rule's mem_mb; the matrix alone is %.1f GB."
        ),
        .bi, length(.starts), sum(.bad), length(.idx),
        length(active_chrom_ids) * length(chip_files) * 8 / 1024^3
      ))
    }
    for (.k in seq_along(.idx)) chip_matrix[, .idx[.k]] <- .cols[[.k]]
    rm(.cols)
    gc()
  }
  stopifnot(
    nrow(chip_matrix) == length(active_chrom_ids),
    ncol(chip_matrix) == length(chip_files)
  )
  saveRDS(chip_matrix, chip_matrix_cache, compress = FALSE)
  message("chip_matrix built and cached")
}

# =============================================================================
# PSI-independent sources for the full-cohort row set (LOADED LATE, ON PURPOSE)
# =============================================================================
# These three reads sit AFTER the chip_matrix build, not before it. The build forks
# mc.cores workers over 2,430 ChIP tabs and already peaks hard: `all_cols` is
# 2,430 x ~709k doubles (~13.8 GB) and `do.call(cbind, ...)` transiently doubles it.
# R's GC dirties pages inside the forks, so copy-on-write does NOT keep a large
# parent object free -- every extra GB held here is multiplied by the fork width.
# Loading the grid (~5 GB) + rbp_score + gene_expr before that step is what OOM-killed
# a real run at 62.3 GB against a 62.5 GB limit (job 6433872, 2026-08-04), with the
# symptom surfacing as a cryptic `dimnames` length error after the workers died.
# Nothing between here and the Phase-1 loop needs them, so keep them last.

# --- the unfiltered (IHEC x ID) grid: event-proximal features for ALL samples ---
# Written by 05 immediately before the PSI merge, so it is PSI-independent by
# construction. Its DNAm columns are M-value-transformed there with the identical
# formula 05 applies to the filtered table -- WITHOUT that, control and focal DNAm
# would sit on different scales inside one ridge fit, silently, since raw beta
# (0-100) and M-values (-9.97..9.97) are both plausible numbers.
#
# Path is PER FILTER. The grid's *content* is filter-specific (event IDs are row
# positions in that filter's event_gr, so the same integer denotes a different event
# under a different filter), so a single shared aggregated_dt.csv.gz would be
# overwritten by whichever 05 ran last and read back here as silently wrong features
# under matching ID values. Nothing else in the repo reads this file.
grid_file <- sprintf("processed_data/aggregated_dt_%s.csv.gz", tf)
if (!file.exists(grid_file)) {
  stop(
    "Missing ", grid_file, " -- the PSI-independent grid 05 writes before the PSI ",
    "merge. Re-run 05 for this transcript_filter (a pre-2026-08-04 05 wrote it to ",
    "the unsuffixed processed_data/aggregated_dt.csv.gz, in RAW beta; that file is ",
    "NOT usable here)."
  )
}
# stringsAsFactors = FALSE deliberately: only the join key needs to be factor-aligned
# with aggregated_dt (below), and every other character column is reconciled by
# align_types() at assembly time. Reading as factors instead would create columns
# whose integer codes index a DIFFERENT level set than aggregated_dt's -- which joins
# and rbindlist can silently mismatch.
grid <- fread(grid_file, stringsAsFactors = FALSE)
grid <- grid[ID %in% all_modelable_ids]
stopifnot(all(c("ID", "IHEC", feature_cols$per_epigenome) %in% names(grid)))
.grid_missing_ids <- base::setdiff(all_modelable_ids, grid[, unique(ID)])
if (length(.grid_missing_ids)) {
  stop(
    length(.grid_missing_ids), " events to build are absent from ", grid_file,
    " (e.g. ", paste(head(.grid_missing_ids, 3L), collapse = ", "),
    ") -- the grid is stale or was written for a different transcript_filter."
  )
}
# Guard the scale explicitly rather than trusting the 05 comment: M-values are
# negative for any window below 50% methylation, raw beta never is.
.grid_dnam <- grep("^DNAm;", names(grid), value = TRUE)
if (length(.grid_dnam)) {
  .dnam_min <- suppressWarnings(min(vapply(
    .grid_dnam, function(cn) min(grid[[cn]], na.rm = TRUE), numeric(1)
  )))
  if (is.finite(.dnam_min) && .dnam_min >= 0) {
    stop(
      "DNAm in ", grid_file, " has no negative values (min ", .dnam_min,
      ") -- it is RAW beta, not M-values. Re-run 05: mixing this with the ",
      "M-value-transformed filtered table puts control and focal DNAm on ",
      "different scales inside the same ridge fit."
    )
  }
}
# Align the join key with aggregated_dt's IHEC factor so the (ID, IHEC) joins below
# compare like with like. An IHEC in the grid but absent from aggregated_dt's levels
# has no sample in the modelled cohort, so it can never be joined to and dropping it
# to NA is correct -- but report the count rather than assume it is zero.
.grid_ihec_chr <- as.character(grid$IHEC)
grid[, IHEC := factor(.grid_ihec_chr, levels = levels(aggregated_dt$IHEC))]
.grid_ihec_unmatched <- sum(is.na(grid$IHEC) & !is.na(.grid_ihec_chr))
if (.grid_ihec_unmatched > 0L) {
  message(sprintf(
    "%d grid rows carry an IHEC not present in aggregated_dt (no modelled sample) - not joinable",
    .grid_ihec_unmatched
  ))
}
rm(.grid_ihec_chr)
setkey(grid, ID, IHEC)

# --- per-(ID, uuid) and per-(uuid, gene_id) sources ----------------------------
# The only two families that cannot be recycled from aggregated_dt: they vary with
# BOTH the event and the sample, so the added rows have no observed row to copy.
rbp_score_dt <- fread(sprintf("processed_data/rbp_score_dt_%s.csv.gz", tf))
rbp_score_dt <- rbp_score_dt[
  uuid %in% all_uuids & ID %in% all_modelable_ids,
  c("ID", "uuid", intersect(feature_cols$per_event_sample, names(rbp_score_dt))),
  with = FALSE
]
setkey(rbp_score_dt, ID, uuid)

gene_expr_dt <- fread(
  sprintf("processed_data/gene_expression_normalised_%s.csv.gz", tf)
)
gene_expr_dt <- gene_expr_dt[
  uuid %in% all_uuids,
  c("gene_id", "uuid", feature_cols$per_gene),
  with = FALSE
]
# STRIP THE ENSEMBL VERSION SUFFIX. gene_expression_normalised carries versioned ids
# ("ENSG00000000419.12") while aggregated_dt/event_annotations_dt carry bare ones
# ("ENSG00000000419"), so keying on the raw value matches NOTHING and every added row
# silently gets NA expression. Same mismatch that bit 05's housekeeping-gene check.
gene_expr_dt[, gene_id := sub("\\.\\d+$", "", gene_id)]
if (anyDuplicated(gene_expr_dt, by = c("gene_id", "uuid"))) {
  # two versions of one gene collapsing onto the same key would make the match()
  # below pick an arbitrary one
  stop(
    "gene_expression_normalised_", tf,
    " has duplicate (gene_id, uuid) rows after stripping version suffixes"
  )
}
setkey(gene_expr_dt, gene_id, uuid)
.ge_hit <- uniqueN(gene_expr_dt$gene_id[
  gene_expr_dt$gene_id %in% as.character(unique(aggregated_dt$gene_id))
])
message(sprintf(
  "gene_expression_normalised: %d of %d aggregated_dt genes matched after version strip",
  .ge_hit, uniqueN(aggregated_dt$gene_id)
))
if (.ge_hit == 0L) {
  stop("No gene_id overlap with aggregated_dt -- expression join would be all-NA")
}

# IJC/SJC are deliberately left NA on the added rows. They are blocked from x_cols
# (they are what PSI is computed from), and a row with no PSI has no junction-count
# reading to report either, so NA is the honest value rather than a lookup worth a
# 76 MB join.

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
#
# Loops all_modelable_ids, NOT ids_to_build. ids_to_build has already had
# already_computed_ids (events with a Tier-2 _all_metrics.csv.gz) subtracted, but that
# subtraction is about skipping Tier-2 FITS, not table builds: the Tier-1 screen runs
# on every modelable event, and every event is also a candidate matched CONTROL, so a
# missing table costs a control silently (control_fs_R2 returns NA on an unreadable
# file). This was masked while the v1 tables existed for all events; with a
# FEATURE_TABLE_VERSION bump writing to an empty directory it would leave a real hole
# for exactly the events Tier-2 had already fitted.
.build_res <- pbmcapply::pbmclapply(all_modelable_ids, function(id) {
  data.table::setDTthreads(1L)
  feature_table_file <- file.path(
    feature_table_dir,
    paste0("feature_table_", id, ".csv.gz")
  )
  if (file.exists(feature_table_file)) {
    return(invisible(NULL))
  }

  # Observed rows VERBATIM, then the full-cohort row set built around them. The
  # verbatim half is what 09s-ridge-screen.R's focal fit sees (it re-applies
  # `!is.na(PSI)` at load), so the focal statistic is bit-identical to v1 regardless
  # of anything the assembly does; the added rows exist only so this event can serve
  # as a matched CONTROL for events whose samples it does not itself cover.
  feature_data <- build_full_event_rows(
    aggregated_dt[ID == id], id, all_uuids, sample_covariates, feature_cols,
    grid, rbp_score_dt, gene_expr_dt
  )
  chromhmm_ids <- to(chromhmm_hits[
    from(chromhmm_hits) == which(keep_rows_manual == id)
  ])
  # Now spans every cohort sample, not just the ones with observed PSI, because
  # feature_data does -- so the chromHMM slice and the WGBS CJ below widen with it
  # automatically and need no separate change.
  event_files <- file_table[
    assay_type == "ChIP-Seq" &
      epirr_id_without_version %in% as.character(feature_data[, unique(IHEC)]),
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

  # RBP Step 4 (§3.4): attach the wide per-RBP expression columns for the RBPs
  # with a binding site near THIS event (rbp_per_event). Folded in with the other
  # predictors (long/short/local); events with no nearby RBP add nothing.
  nearby <- rbp_per_event[ID == id, rbp]
  nearby_rbps <- if (length(nearby) == 1L) unlist(nearby[[1]]) else character(0)
  # base::intersect — `conflicted` is active (09zz attaches dplyr + GenomicRanges,
  # both export intersect); bare intersect() errors inside the worker.
  rbp_cols <- base::intersect(paste0("rbp_", nearby_rbps), names(rbp_wide))
  if (length(rbp_cols) > 0L) {
    feature_data[
      rbp_wide[, c("uuid", "transcript_filter", rbp_cols), with = FALSE],
      on = .(uuid, transcript_filter),
      (rbp_cols) := mget(paste0("i.", rbp_cols))
    ]
  }

  # ATOMIC: this job can hit its 1440-min wall mid-write, and a truncated .csv.gz left
  # at the final path would be SKIPPED by the file.exists() gate above on the resubmit
  # -- a silently corrupt table that fread may still partially accept. Write to a temp
  # sibling and rename (atomic on one filesystem). The temp name MUST keep the .csv.gz
  # extension and pass compress= explicitly: fwrite picks compression from the
  # extension, so a ".tmp<pid>" suffix writes PLAIN CSV that then gets renamed to
  # .csv.gz (fread tolerates it, zcat/gzfile do not) -- same trap already documented
  # for 09s-ridge-screen.R's write_atomic.
  tmp_file <- file.path(
    dirname(feature_table_file),
    sprintf(".tmp%d_%s", Sys.getpid(), basename(feature_table_file))
  )
  fwrite(feature_data, tmp_file, compress = "gzip")
  if (!file.rename(tmp_file, feature_table_file)) {
    unlink(tmp_file)
    stop("Failed to rename ", tmp_file, " -> ", feature_table_file)
  }
  invisible(NULL)
})

# --- Phase 1 completeness gate -------------------------------------------------
# pbmclapply/mclapply do NOT abort on a worker error: the failing element comes back
# as a try-error (or NULL if the worker was OOM-killed) and the loop reports success.
# Previously the result was neither assigned nor inspected, so any per-event failure
# left that table missing and this script still exited 0 -- and the screen would then
# silently lose the event (a missing control table makes control_fs_R2 return NA, and a
# missing focal table makes its own job fail much later, ~34k jobs in). That is the
# exact "silently successful" mode the version stamps exist to prevent, so check it.
# FILE EXISTENCE is the authority here, not the try-error count. With
# mc.preschedule = TRUE (the default) mclapply hands each core a contiguous chunk and a
# single failure poisons the RETURN VALUES of that entire chunk -- measured: asking for
# errors at elements 3 and 5 of 6 over 2 cores reports 1, 3 and 5. At 34k events over 16
# cores one real failure would claim ~2,100 "errored" events whose tables were in fact
# written fine. So the error list is only used to surface a diagnostic message; whether
# the build is complete is decided by which files are actually on disk.
.expected <- file.path(
  feature_table_dir, sprintf("feature_table_%d.csv.gz", all_modelable_ids)
)
.missing <- !file.exists(.expected)
.build_failed <- vapply(
  .build_res,
  function(z) inherits(z, "try-error") || inherits(z, "condition"),
  logical(1)
)
if (any(.build_failed)) {
  .first <- which(.build_failed)[1L]
  message(sprintf(
    paste0(
      "Worker error(s) reported (count %d is inflated by mc.preschedule chunk ",
      "poisoning -- treat the missing-file count below as authoritative). First: %s"
    ),
    sum(.build_failed),
    conditionMessage(attr(.build_res[[.first]], "condition"))
  ))
}
if (any(.missing)) {
  stop(sprintf(
    "Feature-table build INCOMPLETE: %d of %d tables missing (e.g. event %s). %s",
    sum(.missing), length(.expected),
    paste(head(all_modelable_ids[.missing], 3L), collapse = ", "),
    "Re-run to resume -- the build is idempotent, existing tables are skipped."
  ))
}
message(sprintf(
  "Feature tables complete: %d/%d present in %s",
  sum(!.missing), length(.expected), feature_table_dir
))

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
  nrotations = nrotations,
  # Tier-1 ridge-screen rotation count (feature-rotation controls per event);
  # cheap → a few hundred gives an empirical-p floor ≈ 1/(R+1). Consumed by
  # 09s-ridge-screen.R, not by the elastic-net.
  #
  # PER EVENT TYPE, and it has to be. 09s-ridge-screen.R caps usage at
  # min(screen_rotations, eligible_controls), so the PSI-independent control pool
  # (FEATURE_TABLE_VERSION 2) is only an upper bound -- at a flat 200 the whole
  # change buys nothing. BH admits a p at rank k only if p <= k*q/m, and the
  # empirical p is floored at 1/(R+1):
  #   RI  m ~= 1,776 per (tf, Event Type, feature_set) family, matched pool measured
  #       791-887. At R=818 the floor is 1.22e-3, needing k >= 22 against a real
  #       (v2, partial) k of 25/20/11 for local/short/long -- so `local` plausibly
  #       clears and `long` does not. At R=200 the floor is 4.98e-3, needing k >= 88.
  #       RI cannot clear at 200 for arithmetic reasons alone.
  #   SE  m ~= 32,370, pool ~15,000, already inside its bound at 200. Raising SE to
  #       match RI would cost ~10 days of screen wall-clock for no inferential gain,
  #       which is why this is not a single global number.
  # An Event Type with no entry here is a hard error in resolve_screen_rotations()
  # rather than a silent 200 -- see 09-ml-shared.R.
  screen_rotations = getOption(
    "EpiATLAS_AS_SCREEN_ROTATIONS", c(RI = 818L, SE = 200L)
  ),
  # Feature sets the Tier-1 screen runs per event (long/short/local) — SEPARATE
  # from `feature_sets` above (which 09zz/Tier-2 always needs all three of). The
  # screen splits so a spatially-local signal isn't diluted by far chromHMM in
  # one omnibus ridge; FDR is then computed per (tf, Event Type, feature_set).
  # Default = all three, and it is ALSO the right setting for reusing a prior
  # omnibus screen: 09s-ridge-screen.R resumes PER FEATURE SET, tagging an
  # untagged legacy per-event file as `long` (identical quantity — one ridge over
  # all epigenetic X) and recomputing only the missing sets. Narrow this only to
  # deliberately skip a set entirely.
  screen_feature_sets = getOption(
    "EpiATLAS_AS_SCREEN_FEATURE_SETS", c("long", "short", "local")
  ),
  # §4.9: thread the global seed to the array workers (09zz uses seed_base +
  # this_id per event; falls back to the option if this field is absent).
  seed = getOption("EpiATLAS_AS_SEED", 42L)
)

dir.create(event_dir, showWarnings = FALSE)
dir.create("event_glmnet_logs", showWarnings = FALSE)

ids_file <- normalizePath(
  file.path("processed_data", sprintf("event_glmnet_ids_%s.txt", tf)),
  mustWork = FALSE
)
all_ids_file <- normalizePath(
  file.path("processed_data", sprintf("event_glmnet_all_ids_%s.txt", tf)),
  mustWork = FALSE
)
cfg_file <- normalizePath(
  file.path("processed_data", sprintf("event_glmnet_cfg_%s.rds", tf)),
  mustWork = FALSE
)
# All modelable events → the Tier-1 ridge screen (09s-dispatch.R fires an array
# over these; the screen has its own per-event idempotency).
writeLines(as.character(all_modelable_ids), all_ids_file)
saveRDS(.slurm_cfg, cfg_file)

# Snakemake `build_feature_tables` runs this with EPIATLAS_AS_ML_PHASE=build to
# stop here (feature tables + session + cfg + all-events id list produced), so
# the build stays a distinct DAG step ahead of the screen. `event_models` runs
# it with no phase set → falls through to the hits-gated elastic-net dispatch.
if (identical(Sys.getenv("EPIATLAS_AS_ML_PHASE"), "build")) {
  message(
    "Build-only phase — feature tables + session + cfg + all_ids written; ",
    "skipping elastic-net dispatch."
  )
  quit(save = "no", status = 0L)
}

# =============================================================================
# Phase 2 dispatch — Tier-2 elastic-net, HITS ONLY
# =============================================================================
# The EN interpreter runs only on events the Tier-1 ridge screen flagged
# (q < threshold, effect > 0), listed in tier1_hits_<tf>.txt by 09s-aggregate.R.
# Until that file exists (screen not run yet) this invocation is BUILD-ONLY:
# feature tables + session + cfg + all-events id list are produced and no EN is
# dispatched — so `rule build_feature_tables` and `rule event_models` can both
# just run this script (build is idempotent; EN fires once hits exist).
hits_file <- file.path(event_dir, sprintf("tier1_hits_%s.txt", tf))
if (!file.exists(hits_file)) {
  message(
    "Build complete. No Tier-1 hits file yet (", hits_file, ") — run the ridge ",
    "screen + aggregate (09s-dispatch.R) first, then re-run this to dispatch ",
    "the elastic-net on hits."
  )
  quit(save = "no", status = 0L)
}
hit_ids <- as.integer(readLines(hits_file))
ids_for_en <- base::intersect(ids_to_build, hit_ids)
writeLines(as.character(ids_for_en), ids_file)
message(sprintf(
  "Tier-1 hits: %d total; %d not-yet-computed → dispatching EN on %d events",
  length(hit_ids), length(ids_for_en), length(ids_for_en)
))
if (length(ids_for_en) == 0L) {
  message("Nothing to dispatch (all hits already computed).")
  quit(save = "no", status = 0L)
}

# Set n_local_test > 0 to run a small subset locally instead of sbatch.
#   n_outer_cores — how many events run in parallel
#   n_inner_cores — workflows parallelised within each event
#                   (passed as 3rd CLI arg; SLURM always uses 1)
n_local_test <- 0L # set > 0 to bypass sbatch
n_outer_cores <- 5L
n_inner_cores <- 11L

if (n_local_test > 0L) {
  test_ids <- head(ids_for_en, n_local_test)
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
  n <- length(ids_for_en)
  system(sprintf(
    'sbatch --array=0-%d%%20 "%s" "%s" "%s" "%s"',
    n - 1L,
    normalizePath("09-1-ml-local-array.sh"),
    ids_file,
    cfg_file,
    normalizePath(".")
  ))
  message(sprintf("Submitted %d array jobs", n))
}
