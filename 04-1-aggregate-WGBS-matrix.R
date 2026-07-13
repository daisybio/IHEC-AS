# this script aggregates the methylation matrices to the regions of interest
set.seed(getOption("EpiATLAS_AS_SEED", 42L))
# Per-filter invocation (PLAN §4.12e per-filter incremental outputs): aggregate
# over this transcript_filter's event regions and write a per-filter output.
# Snakemake passes TRANSCRIPT_FILTER; interactive runs fall back to primary.
this_transcript_filter <- Sys.getenv(
  "TRANSCRIPT_FILTER",
  unset = getOption("EpiATLAS_AS_PRIMARY_FILTER", "biotype_filtered")
)

# Same standalone append_sanity() pattern as 02-2/02-3/03 — append rows to the
# shared qc/sanity_summary.csv so 04-1's β checks land in the central ledger.
append_sanity <- function(check_name, status, value, threshold, notes = "") {
  dir.create("qc", showWarnings = FALSE)
  sanity_file <- "qc/sanity_summary.csv"
  new_row <- data.table::data.table(check_name, status, value, threshold, notes)
  if (file.exists(sanity_file)) {
    # upsert by check_name: drop any prior row for this check before writing
    # the new one, so sanity_summary.csv doesn't accumulate stale/duplicate
    # rows across reruns (append=TRUE never truncated -> the file only grew).
    existing <- data.table::fread(sanity_file)
    existing <- existing[check_name != new_row$check_name]
    data.table::fwrite(rbind(existing, new_row, fill = TRUE), sanity_file)
  } else {
    data.table::fwrite(new_row, sanity_file)
  }
}
# first locate the wgbs matrices
chr_files <- list.files(
  wgbs_matrices_data_dir,
  pattern = ".meth10.csv.gz$",
  full.names = TRUE
)
chr_seqnames <- tstrsplit(
  basename(list.files(
    wgbs_matrices_data_dir,
    pattern = ".meth10.csv.gz$",
    full.names = TRUE
  )),
  split = ".",
  fixed = TRUE,
  keep = 1
)[[1]]
# load the image with the aggregated regions (per transcript_filter)
aggregateOver <- rtracklayer::import(
  sprintf("processed_data/aggregateOver_%s.bed", this_transcript_filter)
)

# function to aggregate the matrices.
# EFFICIENCY (2026-07-08): melt ONLY the CpGs that fall inside an aggregateOver
# window, not the whole chromosome. The old version melted the full nCpG x
# nSamples matrix (chr1 ~ tens of GB) and then index-subset the in-window cells
# — that full melt was the OOM (jobs 6289903/6289920) and the runtime cost.
# Event windows cover a small fraction of each chromosome, so subsetting the
# CpG rows before the melt is ~3x faster and cuts peak memory by ~an order of
# magnitude. Output is byte-identical to the old code (validated on chr22).
aggregate_matrix <- function(file, seqname) {
  chr_matrix <- fread(file)
  sample_cols <- setdiff(names(chr_matrix), "l")
  # CpG position -> overlapping aggregateOver windows
  chr_gr <- chr_matrix[, GRanges(
    seqnames = seqname,
    IRanges(start = l, end = l)
  )]
  agg_hits <- findOverlaps(aggregateOver, chr_gr, ignore.strand = TRUE)
  if (length(agg_hits) == 0L) {
    return(NULL)
  }
  # (window ID, CpG row) pairs — a CpG may fall in several (overlapping) windows
  hit_dt <- data.table(ID = from(agg_hits), row = to(agg_hits))
  # melt only the in-window CpG rows (bounded by in-window CpGs x samples,
  # not the whole chromosome)
  keep_rows <- unique(hit_dt$row)
  sub <- chr_matrix[keep_rows]
  sub[, row := keep_rows]
  long <- melt(
    sub,
    id.vars = "row",
    measure.vars = sample_cols,
    variable.name = "ihec",
    value.name = "score"
  )
  long <- long[score != -1]
  # expand each CpG's per-sample scores to every window it overlaps, then
  # aggregate to (window, sample): mean beta + covered-CpG count
  agg_dt <- merge(long, hit_dt, by = "row", allow.cartesian = TRUE)[
    ,
    .(score = mean(score), n = .N),
    by = .(ID, ihec)
  ]
  agg_dt[, name := as.factor(aggregateOver$name[ID])]
  return(agg_dt)
}
# DT threads must be 1 in parent before fork (OpenMP + fork = crash)
data.table::setDTthreads(1L)
# bind all chromosomes.
# MEMORY: cap fork count. Each fork loads one chromosome's full methylation
# matrix (all ~600 samples x that chr's CpGs) + melts it; at mc.cores=8 the
# large chromosomes (chr1/chr2) coincide and OOM-kill the step at 64G (job
# 6289903, MaxRSS 65G). Fewer forks = lower simultaneous peak; the per-chr
# read+overlap is I/O/CPU-bound so the throughput hit is small. Override via
# EpiATLAS_AS_WGBS_CORES. Per-fork peak is LARGE on big chromosomes: chr1 has
# ~1.15M in-window CpGs x 645 samples, so the melt + merge + copies peak ~48G in
# ONE fork (measured: 2 forks hit 97G at the 96G cap). So default serial
# (mc.cores=1) -> ~48G peak, safe under 96G. Raise only with proportionally more
# mem (2 forks need ~128-160G). A chunked merge would cap this — see Open Items.
result <- pbmcapply::pbmcmapply(
  aggregate_matrix,
  file = chr_files,
  seqname = chr_seqnames,
  SIMPLIFY = FALSE,
  mc.cores = getOption("EpiATLAS_AS_WGBS_CORES", 1L)
)
stopifnot(length(Reduce(intersect, lapply(result, function(x) x$name))) == 0)
# create dir if not exists
dir.create(sample_dt_dir, showWarnings = FALSE)
dir.create(file.path(sample_dt_dir, "error_logs"), showWarnings = FALSE)
dir.create(file.path(sample_dt_dir, "logs"), showWarnings = FALSE)
# bind all chromosomes into single table
full_result <- rbindlist(result)

# β distribution sanity checks (§4.13c) — meth10 score is % methylation (0-100), not fraction
stopifnot(all(full_result$score >= 0 & full_result$score <= 100))

# §4.13c item 4 — NA/computability. A (ID,ihec) row exists in full_result iff
# that window had >=1 covered CpG (n>=1; score is their mean, never NA by
# construction). So "CpGs>0 => β computable" reduces to: no NA score among the
# covered rows. (Events with CpGs=0 have no row -> become NA at the 05 join,
# which is expected and checked there.)
stopifnot(!anyNA(full_result$score))
stopifnot(all(full_result$n >= 1L))
append_sanity(
  sprintf("dnam_beta_na_covered_%s", this_transcript_filter),
  "PASS", 0, 0,
  "NA β among CpG-covered windows (must be 0; uncovered windows have no row)"
)

full_result[, feature_set := sub("_[0-9]+$", "", name)]

# median CpG count >= 5 per window type
for (fs in unique(full_result$feature_set)) {
  med_n <- full_result[feature_set == fs, median(n)]
  if (med_n < 5) {
    warning(sprintf("median CpGs = %.1f for feature_set '%s' (< 5)", med_n, fs))
  }
  append_sanity(
    sprintf("dnam_cpg_median_%s__%s", fs, this_transcript_filter),
    if (med_n >= 5) "PASS" else "WARN", med_n, 5,
    "median CpG count per window (>= 5 target)"
  )
}

# coverage gap: fraction of aggregateOver regions with no WGBS data
n_total <- length(aggregateOver)
n_covered <- full_result[, uniqueN(ID)]
gap_frac <- 1 - n_covered / n_total
if (gap_frac > 0.1) {
  warning(sprintf("%.1f%% of windows have no WGBS coverage", gap_frac * 100))
}
append_sanity(
  sprintf("dnam_coverage_gap_%s", this_transcript_filter),
  if (gap_frac <= 0.1) "PASS" else "WARN", gap_frac, 0.1,
  "fraction of aggregateOver windows with no WGBS data"
)

# β histograms per window type
dir.create("qc/distributions", recursive = TRUE, showWarnings = FALSE)
pdf(
  sprintf("qc/distributions/dnam_beta_%s.pdf", this_transcript_filter),
  width = 9, height = 4
)
for (fs in unique(full_result$feature_set)) {
  hist(
    full_result[feature_set == fs, score],
    main = paste("beta distribution -", fs),
    xlab = "mean beta",
    breaks = 50,
    col = "steelblue",
    border = "white"
  )
}
# §4.13c item 5 — M-value preview. score is % methylation (0-100); the M-value
# formula assumes beta in [0,1], so divide by 100 first. M = log2((b+e)/(1-b+e))
# with e=0.001. Plotted alongside beta for shape comparison; NOT written to the
# output (schema unchanged) — a preview for the 07/09 M-value feature choice.
for (fs in unique(full_result$feature_set)) {
  b <- full_result[feature_set == fs, score] / 100
  mval <- log2((b + 0.001) / (1 - b + 0.001))
  hist(
    mval,
    main = paste("M-value preview -", fs),
    xlab = "M = log2((beta+e)/(1-beta+e)), e=0.001",
    breaks = 50,
    col = "darkorange",
    border = "white"
  )
}
dev.off()

full_result[, feature_set := NULL]

# write the result to file (per transcript_filter)
fwrite(
  full_result,
  file.path(
    sample_dt_dir,
    sprintf("WGBS_agg_%s.csv.gz", this_transcript_filter)
  )
)
