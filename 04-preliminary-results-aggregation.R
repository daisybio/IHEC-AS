# setwd('~/hiwi/IHEC-AS/')
load('aggregating.rda')
source('04-aggregating_helper.R')
dir.create(sample_dt_dir, recursive=TRUE)
already_computed_files <- list.files(sample_dt_dir, pattern = 'merge_dt.csv.gz')
if (length(already_computed_files) == 0) {
  samples_to_consider <- sample_cols
} else {
  samples_to_consider <- sample_cols[!sample_cols %in% tstrsplit(already_computed_files, '-', fixed = TRUE , keep = 1)[[1]]]
}
aggregated_data <- aggregate_multiple_samples(samples_to_consider, event_dt, event_gr, upstream_gr, downstream_gr, make_cCREs = TRUE)
# names(aggregated_data) <- sample_cols
# saveRDS(aggregated_data, file.path(psi_input_dir, 'aggregated_reference.rds'))