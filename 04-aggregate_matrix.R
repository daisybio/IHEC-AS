chr_files <- list.files(wgbs_matrices_data_dir, pattern = ".meth10.csv.gz$", full.names = TRUE)
chr_seqnames <- tstrsplit(basename(list.files(wgbs_matrices_data_dir, pattern = ".meth10.csv.gz$", full.names = TRUE)), split =".", fixed=TRUE, keep = 1)[[1]]
load("aggregating.rda")

aggregate_matrix <- function(file, seqname){
  chr_matrix <- fread(file)
  chr_gr <- chr_matrix[, GRanges(seqnames = seqname, IRanges(start=l, end=l))]
  agg_hits <- findOverlaps(aggregateOver, chr_gr, ignore.strand=TRUE)
  # melt matrix
  melt_mat <- melt(chr_matrix, id.vars = 'l', variable.name = 'ihec', value.name = 'score')
  stopifnot(nrow(melt_mat)/melt_mat[, nlevels(ihec)] == length(chr_gr))
  # need to repeat ID rows times ihec entries
  melt_agg_from_hits <- rep(from(agg_hits), melt_mat[, nlevels(ihec)])
  # now need to multiply the target id with the corresponding entry
  melt_agg_to_hits <- unlist(lapply(seq.int(melt_mat[, nlevels(ihec)]), function(ihec_i) to(agg_hits) + (rep(nrow(chr_matrix), length(agg_hits)) * ihec_i)))
  # aggregate by ID and IHEC
  agg_dt <- data.table(ID = melt_agg_from_hits,
                       ihec = melt_mat[melt_agg_to_hits, ihec],
                       score = melt_mat[melt_agg_to_hits, score])[score != -1, .(score=mean(score)),
                                                                    by = .(ID, ihec)]
  agg_dt[, name:=as.factor(aggregateOver$name[ID])]
  gc()
  return(agg_dt)
}
#chr1 <- aggregate_matrix(chr_files[1], chr_seqnames[1])
# library(clustermq)
result <- pbmcapply::pbmcmapply(aggregate_matrix, file=chr_files, seqname=chr_seqnames, SIMPLIFY=FALSE)
stopifnot(length(Reduce(intersect, lapply(result, function(x)x$name))) == 0)
fwrite(rbindlist(result), 'sample_dts/WGBS_agg.csv.gz')