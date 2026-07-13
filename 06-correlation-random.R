# ⚠️ LEGACY / NOT WIRED (not a Snakefile rule). Three blockers before this runs:
#   (1) it reads processed_data/agg_long.csv.gz, which NOTHING in the current
#       pipeline produces (06-correlation.Rmd does not write it) — provenance
#       must be resolved / the long table regenerated first;
#   (2) the fwrite below had a typo `sptrinf` (fixed to sprintf) — so it never
#       completed a run. Migrated the load/path/column bits here for per-filter
#       consistency, but this is UNTESTED end-to-end.
#   (3) RETRACTED 2026-07-13: this script's grouping
#       (by=.(feature, annotation, `Event Type`, ID, gene_id, Variability), i.e.
#       per-event, ontology-stratified) does NOT match 06-correlation.Rmd's
#       actual current aggregation_columns (c("Event Type", "transcript_filter",
#       "uuid", "Variability") -- cross-event, WITHIN one sample, no ontology/ID
#       term at all). Shuffling ontology here would not even perturb that
#       statistic. Do not revive as-written -- see
#       revision/file-changes/06-correlation.Rmd.md's retraction section for
#       what a correct null would need to permute instead.
# Per-filter migration (aggregating.rda is gone): read keep_rows + metadata
# explicitly; to_analyze comes from .Rprofile.
this_transcript_filter <- Sys.getenv(
  "TRANSCRIPT_FILTER",
  getOption("EpiATLAS_AS_PRIMARY_FILTER", "biotype_filtered")
)
keep_rows_manual <- readRDS(sprintf("processed_data/keep_rows_manual_%s.rds", this_transcript_filter))
metadata <- fread(sample_metadata_file)
aggregated_dt <- fread(sprintf('processed_data/aggregated_dt_filtered_%s.csv.gz', this_transcript_filter), stringsAsFactors=TRUE)
# keep_rows_manual IS the modelled set (autosomal + cluster_representative);
# the old manual_ids exclusion is already folded into it.
aggregated_dt <- aggregated_dt[ID %in% keep_rows_manual]
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
    if (identical(value, gene_expression_vst)) { # gene expression (§4.12b: vst)
      res <- cor.test(PSI, value, method = 'spearman')
    } else {
      res <- ppcor::pcor.test(PSI, value, gene_expression_vst, method = 'spearman')
    }
    list(corr=res$estimate, p_val=res$p.value, n = .N)
  }}, by = .(feature, annotation, `Event Type`, ID, gene_id, Variability)]
  fwrite(rand_cor_dt, file.path("processed_data", "random_correlations", sprintf(paste0("rand_cor_dt%0", nchar(n_permutations), "d.csv.gz"), i)))
}, mc.cores = 20L)
# , idcol = "i")
# fwrite(perm_pcor_long, "processed_data/perm_pcor_long.csv.gz")