setwd('~/hiwi/IHEC-AS/')
load('aggregating.rda')
source('04-aggregating_helper.R')
aggregated_data <- aggregate_multiple_samples(sample_cols, event_dt, event_gr, upstream_gr, downstream_gr)
# names(aggregated_data) <- sample_cols
saveRDS(aggregated_data, file.path(psi_input_dir, 'aggregated_reference_no_cCRE.rds'))