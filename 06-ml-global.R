source('06-ml.R')

for (this_event in aggregated_dt[, levels(`Event Type`)]){
  print(this_event)
  feature_data <- aggregated_dt[this_event == `Event Type` & ID %in% keep_rows, -'Event Type']
  # feature_data[, cCRE_near:=ID %in% ids_with_cCREs]
  #names(feature_data)[names(feature_data) != response & keep_cols(feature_data, aggregation) & !grepl('tile', names(feature_data), fixed = TRUE) & !endsWith(names(feature_data), 'percentage') & !endsWith(names(feature_data), 'count') & !endsWith(names(feature_data), 'presence')]
  
  # cCRE_rows <- feature_data[, ID %in% ids_with_cCREs]
  # no_cCRE_rows <- feature_data[, !ID %in% ids_with_cCREs]
  
  all_var <- feature_data[, rep(TRUE, .N)]
  high_var <- feature_data[, Variability == 'High Variability']
  low_var <- feature_data[, Variability == 'Low Variability']
  var_list <- list(all_var=all_var, low_var=low_var, high_var=high_var)
  
  feature_data[, ID := NULL]
  feature_data[, IHEC := NULL]
  feature_data[, Variability := NULL]
  feature_data[, gene_id := NULL]
  
  explanatory_all <- names(feature_data)[!names(feature_data) %in% c(response, grouping_cols)]
  
  # make multiple models: one if an enhancer is adjacent and one without enhancers
  #1 all data
  
  #2 without enhancers
  
  #3 enhancer without enhancer variables
  
  #4 enhancer with variables
  
  # non_enhancer_cols <- setdiff(explanatory_all, enhancer_cols)
  # cCRE_type_cols <- sapply(cCRE_cols, function(cCRE_col_i) c(non_enhancer_cols, cCRE_col_i))
  # cCRE_type_explanatory <- unlist(replicate(length(family), cCRE_type_cols, simplify = FALSE), recursive = FALSE)
  # cCRE_type_explanatory <- c() #cCRE_type_explanatory[order(names(cCRE_type_explanatory))] #TODO add all cCRE Types here
  explanatory <-
    c(
      replicate(length(family), list(all = explanatory_all)),
      replicate(length(family), list(noEpigenetic = explanatory_all[Reduce(`&`, lapply(c(histone_marks, 'DNAm', 'summed_enhancer', 'max_promoter'), function(prefix)
        ! startsWith(explanatory_all, prefix)))])),
      # replicate(length(family), list(nocCRE = non_enhancer_cols)),
      replicate(length(family), list(onlyEpigenetic = explanatory_all[Reduce(`|`, lapply(c(histone_marks, 'DNAm', 'summed_enhancer', 'max_promoter'), function(prefix)
        startsWith(explanatory_all, prefix)))]))
      # replicate(length(family), list(all = explanatory_all)),
      # cCRE_type_explanatory
    )#, replicate(2, explanatory, simplify = FALSE))
  
  family <- c('binomial')
  alpha <- list(c(1))
  # data_rows <- list(all=TRUE)#c(replicate(length(family)*2, list(all=TRUE)), replicate(length(family), list(nocCRE=no_cCRE_rows)), replicate(length(family) * length(cCRE_type_explanatory), list(cCRE=cCRE_rows)))#replicate(4, cCRE_rows, simplify = FALSE))
  filter_rows <- var_list #unlist(sapply(var_list, function(var) sapply(data_rows, function(data) data & var, simplify = FALSE), simplify = FALSE), recursive = FALSE)
  # run_glmnet(data = feature_data, explanatory = explanatory$all, response = response, filter_rows = filter_rows$high_var.nocCRE, family = family, parallel = parallel, nfolds = nfolds, alpha = alpha[[1]], grouping_cols = grouping_cols)
  
  combination_ids <- expand.grid(filter_rows=seq_along(filter_rows), explanatory=seq_along(explanatory))
  
  result_names <- paste(this_event, names(filter_rows)[combination_ids$filter_rows], names(explanatory)[combination_ids$explanatory], family, sep='_')
  res_list <- mapply(run_glmnet, SIMPLIFY = FALSE, USE.NAMES = TRUE, model_name=result_names, save_model=TRUE, filter_rows=filter_rows[combination_ids$filter_rows], explanatory=explanatory[combination_ids$explanatory], family=family, response=response, data=list(feature_data), parallel=parallel, nfolds=nfolds, alpha=alpha, grouping_cols=list(grouping_cols))
  # test_id <- which(result_names == 'RI_high_var_noEpigenetic_binomial')
  # run_glmnet(
  #   model_name = result_names[test_id],
  #   save_model = TRUE,
  #   filter_rows = filter_rows[[combination_ids$filter_rows[test_id]]],
  #   explanatory = explanatory[[combination_ids$explanatory[test_id]]],
  #   family = family,
  #   response = response,
  #   data = feature_data,
  #   parallel = parallel,
  #   nfolds = nfolds,
  #   alpha = alpha[[1]],
  #   grouping_cols = grouping_cols
  # )
  # res <-
  #   slurmR::Slurm_Map(
  #     f = run_glmnet,
  #     model_name = result_names,
  #     data = list(feature_data),
  #     explanatory = explanatory[combination_ids$explanatory],
  #     filter_rows = filter_rows[combination_ids$filter_rows],
  #     family = family,
  #     response = response,
  #     nfolds = nfolds,
  #     parallel = parallel,
  #     save_model = TRUE,
  #     alpha = alpha,
  #     grouping_cols = list(grouping_cols),
  #     njobs = length(result_names),
  #     mc.cores = ncores,
  #     plan = "submit",
  #     sbatch_opt = list(
  #       time = "01:00:00",
  #       mem = "50G",
  #       c = "40",
  #       p = "exbio-cpu"
  #     ),
  #     export = ls()
  #   )
  # clustermq::Q(run_glmnet, model_name=result_names, save_model=TRUE, filter_rows=filter_rows, explanatory=explanatory, family=family, response=response, data=list(feature_data), parallel=TRUE, nfolds=nfolds, alpha=alpha, grouping_cols=list(grouping_cols), n_jobs = 3, log_worker=TRUE, memory="60G", template=list(cores=40))
  # names(res_list) <- 
  
  # saveRDS(res_list, paste0(this_event, '_models.rds'))
}