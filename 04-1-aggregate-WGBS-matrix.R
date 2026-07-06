# this script aggregates the methylation matrices to the regions of interest
set.seed(getOption("EpiATLAS_AS_SEED", 42L))
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
# load the image with the aggregated regions
aggregateOver <- rtracklayer::import("processed_data/aggregateOver.bed")

# function to aggregate the matrices
aggregate_matrix <- function(file, seqname) {
  chr_matrix <- fread(file)
  # transform to GRanges
  chr_gr <- chr_matrix[, GRanges(
    seqnames = seqname,
    IRanges(start = l, end = l)
  )]
  # find overlaps
  agg_hits <- findOverlaps(aggregateOver, chr_gr, ignore.strand = TRUE)
  # melt matrix
  melt_mat <- melt(
    chr_matrix,
    id.vars = 'l',
    variable.name = 'ihec',
    value.name = 'score'
  )
  stopifnot(nrow(melt_mat) / melt_mat[, nlevels(ihec)] == length(chr_gr))
  # need to repeat ID rows times ihec entries
  melt_agg_from_hits <- rep(from(agg_hits), melt_mat[, nlevels(ihec)])
  # now need to multiply the target id with the corresponding entry
  melt_agg_to_hits <- unlist(lapply(
    seq.int(melt_mat[, nlevels(ihec)]),
    function(ihec_i) {
      to(agg_hits) + (rep(nrow(chr_matrix), length(agg_hits)) * (ihec_i - 1))
    }
  ))
  # aggregate by ID and IHEC
  agg_dt <- data.table(
    ID = melt_agg_from_hits,
    ihec = melt_mat[melt_agg_to_hits, ihec],
    score = melt_mat[melt_agg_to_hits, score]
  )[score != -1, .(score = mean(score), n = length(score)), by = .(ID, ihec)]

  agg_dt[, name := as.factor(aggregateOver$name[ID])]
  return(agg_dt)
}
# DT threads must be 1 in parent before fork (OpenMP + fork = crash)
data.table::setDTthreads(1L)
# bind all chromosomes
result <- pbmcapply::pbmcmapply(
  aggregate_matrix,
  file = chr_files,
  seqname = chr_seqnames,
  SIMPLIFY = FALSE
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

full_result[, feature_set := sub("_[0-9]+$", "", name)]

# median CpG count >= 5 per window type
for (fs in unique(full_result$feature_set)) {
  med_n <- full_result[feature_set == fs, median(n)]
  if (med_n < 5) {
    warning(sprintf("median CpGs = %.1f for feature_set '%s' (< 5)", med_n, fs))
  }
}

# coverage gap: fraction of aggregateOver regions with no WGBS data
n_total <- length(aggregateOver)
n_covered <- full_result[, uniqueN(ID)]
gap_frac <- 1 - n_covered / n_total
if (gap_frac > 0.1) {
  warning(sprintf("%.1f%% of windows have no WGBS coverage", gap_frac * 100))
}

# β histograms per window type
dir.create("qc/distributions", recursive = TRUE, showWarnings = FALSE)
pdf("qc/distributions/dnam_beta.pdf", width = 9, height = 4)
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
dev.off()

full_result[, feature_set := NULL]

# write the result to file
fwrite(full_result, file.path(sample_dt_dir, 'WGBS_agg.csv.gz'))
