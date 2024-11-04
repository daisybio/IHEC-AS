source('07-ml-helper.R')

# create directory to store genome-wide models
global_model_path <- file.path('processed_data', 'global_models')
dir.create(global_model_path, showWarnings = FALSE)

for (this_event in to_analyze){
  print(this_event)
  # set feature data
  feature_data <- aggregated_dt[this_event == `Event Type` & ID %in% keep_rows, 
                                -'Event Type']
  
  # create filters for variability
  all_var <- feature_data[, rep(TRUE, .N)]
  high_var <- feature_data[, Variability == 'High Variability']
  low_var <- feature_data[, Variability == 'Low Variability']
  filter_rows <- list(all_var=all_var, low_var=low_var, high_var=high_var)
  
  # remove columns that are not needed
  feature_data[, ID := NULL]
  feature_data[, IHEC := NULL]
  feature_data[, Variability := NULL]
  feature_data[, gene_id := NULL]
  
  # get all explanatory columns
  explanatory_all <- names(feature_data)[!names(feature_data) %in%
                                           c(response, grouping_cols)]
  
  # set some fixed values
  family <- c('binomial')
  alpha <- list(c(1))
  
  epigenetic_columns <- c(histone_marks, 'DNAm', 'summed_enhancer', 
                          'max_promoter')
  # make the list of explanatory columns for the different models
  explanatory <-
    c(
      replicate(length(family), 
                list(all = explanatory_all)),
      replicate(length(family), 
                list(noEpigenetic = 
                     explanatory_all[Reduce(`&`, lapply(epigenetic_columns,
                    function(prefix) !startsWith(explanatory_all, prefix)))])),
      replicate(length(family), 
                list(onlyEpigenetic = 
                     explanatory_all[Reduce(`|`, lapply(epigenetic_columns, 
                     function(prefix) startsWith(explanatory_all, prefix)))])))
  # make all combinations of filter_rows and explanatory columns
  combination_ids <- expand.grid(filter_rows=seq_along(filter_rows), 
                                 explanatory=seq_along(explanatory))
  # set the names for the result model
  result_names <- paste(this_event, 
                        names(filter_rows)[combination_ids$filter_rows], 
                        names(explanatory)[combination_ids$explanatory], 
                        family, 
                        sep='_')
  result_names <- file.path(global_model_path, result_names)
  # run the models
  res_list <-
    mapply(
      run_ML,
      SIMPLIFY = FALSE,
      USE.NAMES = TRUE,
      model_name = result_names,
      save_model = TRUE,
      filter_rows = filter_rows[combination_ids$filter_rows], 
      explanatory = explanatory[combination_ids$explanatory], 
      family = family,
      response = response,
      data = list(feature_data),
      parallel = parallel,
      nfolds = nfolds,
      alpha = alpha,
      grouping_cols = list(grouping_cols)
    )
}
