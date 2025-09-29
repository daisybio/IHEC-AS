# load image and table
load('processed_data/aggregating.rda')
aggregated_dt <- fread('processed_data/aggregated_dt_filtered.csv.gz', stringsAsFactors=TRUE)
# remove manual_ids, i.e., filter IDs in keep_rows
aggregated_dt <- aggregated_dt[ID %in% keep_rows]
# check if there are still manual_ids in the dt
stopifnot(!any(manual_ids %in% aggregated_dt$ID))
aggregated_dt[, `Event Type`:=factor(`Event Type`, levels=to_analyze)]

minimum_events_per_group <- 10
n_permutations <- 1000
agg_long_tmp <- fread("processed_data/agg_long.csv.gz")

# perm_pcor_long <- rbindlist(
pbmclapply(1:n_permutations, function(i) {
  set.seed(i)
  # shuffle the annotation column
  perm_annotation_mapping <- metadata[, .(IHEC = epirr_id_without_version, 
                                          perm_annotation = sample(harmonized_sample_ontology_term_high_order_fig1))]
  agg_long_tmp[perm_annotation_mapping, on=.NATURAL, annotation := perm_annotation]
  rand_cor_dt <- agg_long_tmp[, {if (.N >= minimum_events_per_group & sd(value) != 0 & sd(PSI) != 0){
    if (identical(value, gene_expression)) { # gene expression
      res <- cor.test(PSI, value, method = 'spearman')
    } else {
      res <- ppcor::pcor.test(PSI, value, gene_expression, method = 'spearman')
    }
    list(corr=res$estimate, p_val=res$p.value, n = .N)
  }}, by = .(feature, annotation, `Event Type`, ID, gene_id, Variability)]
  fwrite(rand_cor_dt, file.path("processed_data", "random_correlations", sptrinf(paste0("rand_cor_dt%0", nchar(n_permutations), "d.csv.gz"), i)))
}, mc.cores = 20L)
# , idcol = "i")
# fwrite(perm_pcor_long, "processed_data/perm_pcor_long.csv.gz")