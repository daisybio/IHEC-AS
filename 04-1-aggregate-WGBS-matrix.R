# this script aggregates the methylation matrices to the regions of interest
# first locate the wgbs matrices
chr_files <- list.files(wgbs_matrices_data_dir, pattern = ".meth10.csv.gz$", full.names = TRUE)
chr_seqnames <- tstrsplit(basename(list.files(wgbs_matrices_data_dir, pattern = ".meth10.csv.gz$", full.names = TRUE)), split =".", fixed=TRUE, keep = 1)[[1]]
# load the image with the aggregated regions
load("processed_data/aggregating.rda")

# function to aggregate the matrices
aggregate_matrix <- function(file, seqname){
  # read the matrix
  chr_matrix <- fread(file)
  # transform to GRanges
  chr_gr <- chr_matrix[, GRanges(seqnames = seqname, IRanges(start=l, end=l))]
  # find overlaps
  agg_hits <- findOverlaps(aggregateOver, chr_gr, ignore.strand=TRUE)
  # melt matrix
  melt_mat <- melt(chr_matrix, id.vars = 'l', variable.name = 'ihec', value.name = 'score')
  stopifnot(nrow(melt_mat)/melt_mat[, nlevels(ihec)] == length(chr_gr))
  # need to repeat ID rows times ihec entries
  melt_agg_from_hits <- rep(from(agg_hits), melt_mat[, nlevels(ihec)])
  # now need to multiply the target id with the corresponding entry
  melt_agg_to_hits <- unlist(lapply(seq.int(melt_mat[, nlevels(ihec)]), function(ihec_i) to(agg_hits) + (rep(nrow(chr_matrix), length(agg_hits)) * (ihec_i - 1))))
  # aggregate by ID and IHEC
  agg_dt <- data.table(ID = melt_agg_from_hits,
                       ihec = melt_mat[melt_agg_to_hits, ihec],
                       score = melt_mat[melt_agg_to_hits, score])[score != -1, .(score=mean(score)),
                                                                    by = .(ID, ihec)]
  agg_dt[, name:=as.factor(aggregateOver$name[ID])]
  gc()
  return(agg_dt)
}
# bind all chromosomes
result <- pbmcapply::pbmcmapply(aggregate_matrix, file=chr_files, seqname=chr_seqnames, SIMPLIFY=FALSE)
stopifnot(length(Reduce(intersect, lapply(result, function(x)x$name))) == 0)
# create dir if not exists
dir.create(sample_dt_dir, showWarnings = FALSE)
dir.create(file.path(sample_dt_dir, "error_logs"), showWarnings = FALSE)
dir.create(file.path(sample_dt_dir, "logs"), showWarnings = FALSE)
# write the result to file
fwrite(rbindlist(result), file.path(sample_dt_dir, 'WGBS_agg.csv.gz'))
