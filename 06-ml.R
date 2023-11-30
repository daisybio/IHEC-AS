# setwd('~/hiwi/IHEC-AS/')
# source('.Rprofile')

library(glmnet)
library(ranger)
library(caret)
library(RSNNS)
library(MLmetrics)

run_glmnet <- function(model_name, data, explanatory, response, filter_rows = TRUE, scale_explanatory = TRUE, scale_response = FALSE, family = 'gaussian', type.measure = ifelse(family == 'gaussian', 'mse', 'Balanced Accuracy'), nfolds=5, parallel=FALSE, alpha = 1, seed = 1234, grouping_cols = c("seqnames", "harmonized_sample_ontology_intermediate"), save_model=FALSE) {
  # set seed
  set.seed(seed)
  # make copy so we don't change original data by reference
  this_data <- copy(subset(data, filter_rows, c(explanatory, response, grouping_cols)))
  # remove cols that are all NA, maybe the mark is missing or no overlap was there
  not_all_na <- this_data[, sapply(.SD, function(x) !(all(is.na(x)) | sd(x, na.rm = TRUE) == 0)), .SDcols = c(explanatory, response)]
  stopifnot(not_all_na[response])
  explanatory <- explanatory[explanatory %in% names(not_all_na)[not_all_na]]
  not_all_na <- c(not_all_na, setNames(rep(TRUE, length(grouping_cols)), grouping_cols))
  this_data <- this_data[, ..not_all_na]
  
  if (family == 'binomial') {
    # split in train and test
    perc_test <- .1
    
    test_by_group <- sapply(grouping_cols, function(group) {
      if (group == 'seqnames')
        return(c('chr1'))
      if (group == 'harmonized_sample_ontology_intermediate')
        return(c('T cell', 'lymphocyte of B lineage', 'B cell derived cell line'))
      # return(this_data[, sort(table(factor(get(group)))/.N)])
      group_counts <-
        sample(this_data[, table(factor(get(group)))])
      # this_data[get(group) %in% names(group_counts)[cumsum(group_counts) <= (perc_test*this_data[, .N])], which = TRUE]
      names(group_counts)[seq.int(min(which(
        cumsum(group_counts) >= (perc_test * this_data[, .N])
      )))]
    }, simplify = FALSE)
    sapply(grouping_cols, function(group) message(paste0(group, ": removing ", paste(test_by_group[[group]], collapse = ', '), " for testing. This accounts for ", this_data[, sum(get(group) %in% test_by_group[[group]])/.N])))
    test_ids_by_group <- sapply(names(test_by_group), function(group) this_data[get(group) %in% test_by_group[[group]], which=TRUE])
    test_ids <- Reduce(union, test_ids_by_group)
    test_data <- copy(this_data[test_ids])
    this_data <- this_data[-test_ids]
    
    # split in included and excluded
    dichotomization_thresholds <- c(1/3, 2/3) # this_data[, quantile(get(response), c(1/3, 2/3))]
    message(paste('dichotomization thresholds:', dichotomization_thresholds[1], dichotomization_thresholds[2]))
    
    old_response <- response
    response <- 'binary'
    for (dt in list(this_data, test_data)) {
      dt[get(old_response) <= dichotomization_thresholds[1], (response):=control_class]
      dt[get(old_response) >= dichotomization_thresholds[2], (response):=case_class]
      dt[, (response) := ordered(get(response), levels=class_levels)]
    }
    this_data <- na.omit(this_data, cols=response)
    test_data <- na.omit(test_data, cols=response)
    
    real_test_ids <-
      Reduce(intersect, lapply(names(test_by_group), function(grouping_col)
        test_data[get(grouping_col) %in% test_by_group[[grouping_col]], which = TRUE]))
    message(sprintf("Number of train samples: %d, Number of test samples: %d, Percentage %.2f", nrow(this_data), nrow(test_data), nrow(test_data)/(nrow(this_data) + nrow(test_data))))
    message(sprintf("Number of train samples: %d, Number of real test samples: %d, Percentage %.2f", nrow(this_data), nrow(test_data[real_test_ids]), nrow(test_data[real_test_ids])/(nrow(this_data) + nrow(test_data[real_test_ids]))))
    message(paste('train classes', paste(this_data[, 100 * table(get(response))/.N], collapse = ' ')))
    message(paste('test classes', paste(test_data[, 100 * table(get(response))/.N], collapse = ' ')))
    message(paste('real test classes', paste(test_data[real_test_ids, 100 * table(get(response))/.N], collapse = ' ')))
    # ggplot(this_data, aes(x=seqnames, fill=binary)) + geom_bar(position = 'fill') + facet_wrap(~ harmonized_sample_ontology_intermediate)
  }
  
  # now replace NA by mean value and 
  for (j in which(names(this_data) %in% explanatory)) {
    this_mean <- mean(this_data[[j]], na.rm = TRUE)
    set(this_data, which(is.na(this_data[[j]])), j, this_mean)
    if (family == 'binomial')
      set(test_data, which(is.na(test_data[[j]])), j, this_mean)
  }
  means_sds <- melt(this_data[, c(sapply(explanatory, function(x) c(mean(get(x), na.rm = TRUE), sd(get(x), na.rm=TRUE)), simplify = FALSE), list(name=c('mean', 'sd')))], id.vars = 'name')
  
  # scale all explanatory vars to mean = 0 and sd = 1
  if (scale_explanatory ) {
    this_data[, (explanatory) := lapply(explanatory, function(x) scale(get(x)))]
    if (family == 'binomial')
      test_data <- test_data[, (explanatory):=lapply(explanatory, function(x) (get(x) - means_sds[name == 'mean' & variable == x, value])/means_sds[name == 'sd' & variable == x, value])]
  }
  if (scale_response & family == 'gaussian') this_data[, (response) := scale(get(response))]
  # make cv
  #TODO: find way to split into foldid somthing like this: sapply(grouping_cols, function(group) caret::groupKFold(this_data[, factor(get(group))], nfolds))
  cvfit <- list()
  
  if (is.null(grouping_cols)) {
    folds <- createFolds(this_data[, get(response)], k = nfolds, returnTrain = TRUE)
    # Initialize a vector to store test fold assignment
    foldid <- integer(this_data[, .N])
    
    # Assign each observation to the fold where it's in the test set
    for (i in seq_along(folds)) {
      test_indices <- setdiff(seq.int(this_data[, .N]), folds[[i]])
      foldid[test_indices] <- i
    }
  } else {
    folds <- list()
    i <- 0
    max_tries <- 10
    while(length(folds) < nfolds){
      if (i >= max_tries) stop(paste("Stopping... tried to find a grouped fold many times:", i))
      i <- i + 1
      folds <- caret::groupKFold(this_data[, get(grouping_cols[1])], nfolds)
      
    }
    # Initialize a vector to store test fold assignment
    foldid <- integer(this_data[, .N])
    
    # Assign each observation to the fold where it's in the test set
    for (i in seq_along(folds)) {
      test_indices <- setdiff(seq.int(this_data[, .N]), folds[[i]])
      foldid[test_indices] <- i
    }
  }
  # foldid <- sample(seq.int(nfolds), size = nrow(this_data), replace = TRUE)
  # lapply(folds, function(ids) this_data[-ids, table(get(response), droplevels(seqnames))])
  
  for (random in c(FALSE, TRUE)){
    if (random) this_data[, (response):=sample(get(response))]
    
    if(family == 'gaussian') {
      if(random) browser()
      this_cvfit <- lapply(alpha, function(a) cv.glmnet(as.matrix(this_data[, ..explanatory]), this_data[, get(response)], family = family, type.measure = type.measure, foldid = foldid, parallel = parallel, keep = FALSE, alpha = a))
      names(this_cvfit) <- paste(random, alpha, sep = "::")
      coefs <- coef(this_cvfit[[paste(random, alpha, sep = "::")]], s = my_lambda) #model_glmnet$bestTune$lambda)
      nonzero_explanatory <- rownames(coefs)[as.matrix(coefs != 0 & rownames(coefs) != '(Intercept)')]
      if(length(nonzero_explanatory) > 0){
        this_cvfit[[paste(random, "OLS", sep = "::")]] <- glm(formula = formula(paste(response, '~', paste0('`', nonzero_explanatory, '`', collapse = ' + '))), data = this_data, family = family, model = FALSE, x = FALSE)
        this_cvfit[[paste(random, "OLS", sep = "::")]]$data <- NULL
      # this_cvfit[[paste(random, "FULLOLS", sep = "::")]] <- glm(formula = formula(paste(response, '~', paste0('`', explanatory, '`', collapse = ' + '))), data = this_data, family = family, model = FALSE, y = FALSE)
      }
    } else if(family == 'binomial') {
      myControl <- caret::trainControl(
        method = "cv",
        index = folds,
        search = 'grid',
        verboseIter = TRUE,
        classProbs = TRUE,
        allowParallel = parallel,
        returnData = FALSE,
        summaryFunction = function(data, lev = NULL, model = NULL) {
          cm <- caret::confusionMatrix(data$pred, data$obs)
          out <- c(cm$overall, cm$byClass, twoClassSummary(data, lev, model), prSummary(data, lev, model), MCC=ModelMetrics::mcc(ifelse(data$obs == case_class, 1, 0), data[[case_class]], 0.5))
          return(out)
        },
        selectionFunction = "oneSE"
      )
      
      this_cvfit <- list()
      
      
      
      
      # this_cvfit <- lapply(alpha, function(a) cv.glmnet(as.matrix(this_data[, ..explanatory]), this_data[, get(response)], family = family, type.measure = type.measure, foldid = foldid, parallel = parallel, keep = FALSE, alpha = a))
      # names(this_cvfit) <- paste(random, alpha, sep = "::")
      # cvfit_train_results <- data.table(obs = this_data[, get(response)],
      #                             pred = ordered(predict(this_cvfit[[1]], newx=as.matrix(this_data[, ..explanatory]), s=my_lambda, type="class")[, my_lambda], levels = class_levels))
      # cvfit_train_results[[case_class]] <- predict(this_cvfit[[1]], newx=as.matrix(this_data[, ..explanatory]), s=my_lambda, type="response")[, my_lambda]
      # cvfit_train_results[[control_class]] <- 1 - cvfit_train_results[[case_class]]
      # twoClassSummary(cvfit_train_results, lev = class_levels)
      # caret::confusionMatrix(data=cvfit_train_results[, pred], reference=this_data[, get(response)], positive=case_class, mode='everything')
      # 
      # cvfit_test_results <- data.table(obs = test_data[, get(response)],
      #                             pred = ordered(predict(this_cvfit[[1]], newx=as.matrix(test_data[, ..explanatory]), s=my_lambda, type="class")[, my_lambda], levels = class_levels))
      # cvfit_test_results[[case_class]] <- predict(this_cvfit[[1]], newx=as.matrix(test_data[, ..explanatory]), s=my_lambda, type="response")[, my_lambda]
      # cvfit_test_results[[control_class]] <- 1 - cvfit_test_results[[case_class]]
      # twoClassSummary(cvfit_test_results, lev = class_levels)
      # caret::confusionMatrix(data=cvfit_test_results[, pred], reference=test_data[, get(response)], positive=case_class, mode='everything')
      model_glmnet <-
        caret::train(
          x = this_data[, ..explanatory],
          y = this_data[, get(response)],
          method = 'glmnet',
          metric = type.measure,
          trControl = myControl,
          tuneGrid = expand.grid(alpha = alpha,
                                 lambda = 10^-(rev(seq(1,4,length.out = 21)))))
      this_cvfit[[paste(random, alpha, sep = "::")]] <- model_glmnet
      
      my_s <-
        max(model_glmnet$results[model_glmnet$results[[type.measure]] >= 
                                   (max(model_glmnet$results[, type.measure], na.rm=TRUE) - 
                                      model_glmnet$results[which.max(model_glmnet$results[, type.measure]), paste0(type.measure, 'SD')]/sqrt(nfolds)),
                                 'lambda'], na.rm = TRUE)

      stopifnot(model_glmnet$bestTune$lambda == my_s)
      # cvfit_train_results <- data.table(obs = this_data[, get(response)],
      #                             pred = ordered(predict(model_glmnet$finalModel, newx=as.matrix(this_data[, ..explanatory]), s=my_s, type="class")[, 1], levels = class_levels))
      # cvfit_train_results[[case_class]] <- predict(model_glmnet$finalModel, newx=as.matrix(this_data[, ..explanatory]), s=my_s, type="response")[, 1]
      # cvfit_train_results[[control_class]] <- 1 - cvfit_train_results[[case_class]]
      # twoClassSummary(cvfit_train_results, lev = class_levels)
      # caret::confusionMatrix(data=cvfit_train_results[, pred], reference=this_data[, get(response)], positive=case_class, mode='everything')
      # 
      # cvfit_test_results <- data.table(obs = test_data[seqnames == 'chr1', get(response)],
      #                             pred = ordered(predict(model_glmnet$finalModel, newx=as.matrix(test_data[seqnames == 'chr1', ..explanatory]), s=my_s, type="class")[, 1], levels = class_levels))
      # cvfit_test_results[[case_class]] <- predict(model_glmnet$finalModel, newx=as.matrix(test_data[seqnames == 'chr1', ..explanatory]), s=my_s, type="response")[, 1]
      # cvfit_test_results[[control_class]] <- 1 - cvfit_test_results[[case_class]]
      # twoClassSummary(cvfit_test_results, lev = class_levels)
      # caret::confusionMatrix(data=cvfit_test_results[, pred], reference=test_data[seqnames == 'chr1', get(response)], positive=case_class, mode='everything')
      
      coefs <- coef(model_glmnet$finalModel, s = my_s) #model_glmnet$bestTune$lambda)
      nonzero_explanatory <- rownames(coefs)[as.matrix(coefs != 0 & rownames(coefs) != '(Intercept)')]
      if(length(nonzero_explanatory) > 0) {
        this_cvfit[[paste(random, "OLS", sep = "::")]] <- glm(formula = formula(paste(response, '~', paste0('`', nonzero_explanatory, '`', collapse = ' + '))), data = this_data, family = family, model = FALSE, y = FALSE, x = FALSE)
        this_cvfit[[paste(random, "OLS", sep = "::")]]$data <- NULL
      }
      
      model_ranger <-
        caret::train(
          x = this_data[, ..explanatory],
          y = this_data[, get(response)],
          method = 'ranger',
          metric = type.measure,
          # class.weights = this_data[, table(get(response))/.N],
          #sample.
          # num.threads = ifelse(parallel, max(floor(ncores/nfolds), 1L), 1L),
          num.trees = 200,
          trControl = myControl,
          tuneGrid = expand.grid(
            mtry = unique(floor(seq.int(1, sqrt(length(explanatory)), length.out = 5))), #,
            splitrule = "gini",
            min.node.size = 10
          ),
          importance = "impurity"
        )
      this_cvfit[[paste(random, "RF", sep = "::")]] <- model_ranger
      
      # train_results <- cbind(data.table(obs = this_data[, get(response)], data.table(model_ranger$finalModel$predictions)))
      # train_results[, pred := ordered(ifelse(model_ranger$finalModel$predictions[, control_class] > model_ranger$finalModel$predictions[, case_class], control_class, case_class))]
      # twoClassSummary(train_results, lev = class_levels)
      # caret::confusionMatrix(data=train_results[, pred], reference=this_data[, get(response)], positive=case_class, mode='everything')
      # 
      # rf_predictions <- predict(model_ranger$finalModel, data = test_data[seqnames == 'chr1', ..explanatory])
      # test_results <- cbind(data.table(obs = test_data[seqnames == 'chr1', get(response)], data.table(rf_predictions$predictions)))
      # test_results[, pred := ordered(ifelse(rf_predictions$predictions[, control_class] > rf_predictions$predictions[, case_class], control_class, case_class))]
      # twoClassSummary(test_results, lev = class_levels)
      # caret::confusionMatrix(data=test_results[, pred], reference=test_data[seqnames == 'chr1', get(response)], positive=case_class, mode='everything')
      
      # numTrees <- pbmclapply(1:200, function(i){
      #   rf_predictions <- predict(model_ranger$finalModel, data = test_data[seqnames == 'chr1', ..explanatory], num.trees = i)
      #   test_results <- cbind(data.table(obs = test_data[seqnames == 'chr1', get(response)], data.table(rf_predictions$predictions)))
      #   test_results[, pred := ordered(ifelse(rf_predictions$predictions[, control_class] > rf_predictions$predictions[, case_class], control_class, case_class))]
      #   twoClassSummary(test_results, lev = class_levels)
      #   caret::confusionMatrix(data=test_results[, pred], reference=test_data[seqnames == 'chr1', get(response)], positive=case_class, mode='everything')
      # })
      
      model_mlp <- caret::train(
        x = this_data[, ..explanatory],
        y = this_data[, get(response)],
        method = 'mlp',
        metric = type.measure,
        trControl = myControl,
        # tuneLength = floor(length(explanatory) / 3)
        tuneGrid = expand.grid(size = unique(floor(seq.int(1, length(explanatory)*2/3, length.out = 5))))
        #   layer1 = ((1:3) * 2) - 1,
        #   layer2 = ((1:3) * 2) - 1,
        #   layer3 = ((1:3) * 2) - 1
        # )
      )
      this_cvfit[[paste(random, "MLP", sep = "::")]] <- model_mlp
      
      # preds <- predict(model_mlp$finalModel)
      # colnames(preds) <- class_levels
      # train_results <- cbind(data.table(obs = this_data[, get(response)], data.table(preds)))
      # train_results[, pred := ordered(ifelse(get(control_class) > get(case_class), control_class, case_class))]
      # twoClassSummary(train_results, lev = class_levels)
      # caret::confusionMatrix(data=train_results[, pred], reference=this_data[, get(response)], positive=case_class, mode='everything')
      # 
      # mlp_predictions <- predict(model_mlp$finalModel, newdata = test_data[seqnames == 'chr1', ..explanatory])
      # colnames(mlp_predictions) <- class_levels
      # test_results <- cbind(data.table(obs = test_data[seqnames == 'chr1', get(response)], mlp_predictions))
      # test_results[, pred := ordered(ifelse(get(control_class) > get(case_class), control_class, case_class))]
      # twoClassSummary(test_results, lev = class_levels)
      # caret::confusionMatrix(data=test_results[, pred], reference=test_data[seqnames == 'chr1', get(response)], positive=case_class, mode='everything')
      
    } else stop('unknown family')
    cvfit <- c(cvfit, this_cvfit)
  }
  res_list <- list(cvfit=cvfit, nsamples=nrow(this_data))
  
  if (family == 'binomial') {
    if (scale_explanatory) {
      res_list <- c(res_list, list(test_data=test_data, test_by_group=test_by_group, explanatory=explanatory, dichotomization_thresholds=dichotomization_thresholds, means_sds=means_sds))
    }
    res_list <- c(res_list, list(test_data=test_data,test_by_group=test_by_group, explanatory=explanatory, dichotomization_thresholds=dichotomization_thresholds))
  }
  if (save_model)
    saveRDS(object = res_list, file = paste0(model_name, '.rds'))
  return(res_list)
}

aggregated_dt <- fread('aggregated_dt_filtered.csv.gz', stringsAsFactors = TRUE)
response <- 'PSI'
nfolds <- 10
parallel <- TRUE

load('aggregating.rda', verbose=TRUE)
# cCREs <- data.table::fread('data/GRCh38-cCREs.bed')
# names(cCREs) <- c('seqnames', 'start', 'end', 'some_id', 'accession', 'cCRE_type')
# cCRE_gr <- cCREs[, GRanges(seqnames = seqnames, IRanges(start = start, end = end))]
# cCRE_hits <- findOverlaps(event_gr, cCRE_gr, maxgap = vicinity, ignore.strand=TRUE)

# ids_with_cCREs <- aggregated_dt[, unique(ID)[unique(ID) %in% from(cCRE_hits)]]

# if (parallel) {
#   require(doMC)
#   registerDoMC(nfolds)
# }

# enhancer_cols <- names(aggregated_dt)[grepl('(CTCF)|(DNase)|(ELS)|(PLS)', names(aggregated_dt))]
# 
# cCRE_regions <- c('CTCF-only', 'dELS', 'pELS', 'PLS', 'DNase-H3K4me3')
# cCRE_cols <- sapply(cCRE_regions, function(cCRE_type) names(aggregated_dt)[grepl(cCRE_type, names(aggregated_dt))])

aggregated_dt[event_dt, on=c(ID="ID"), seqnames:= seqnames]
sample_metadata <- data.table::fread(sample_metadata_file, stringsAsFactors = TRUE)
aggregated_dt[sample_metadata, on=c(IHEC='epirr_id_without_version'), harmonized_sample_ontology_intermediate:=harmonized_sample_ontology_intermediate]
grouping_cols <- c('seqnames','harmonized_sample_ontology_intermediate')
# 
# for (this_event in aggregated_dt[, levels(`Event Type`)]){
#   print(this_event)
#   feature_data <- aggregated_dt[this_event == `Event Type`, -'Event Type']
#   # feature_data[, cCRE_near:=ID %in% ids_with_cCREs]
#   #names(feature_data)[names(feature_data) != response & keep_cols(feature_data, aggregation) & !grepl('tile', names(feature_data), fixed = TRUE) & !endsWith(names(feature_data), 'percentage') & !endsWith(names(feature_data), 'count') & !endsWith(names(feature_data), 'presence')]
#   
#   # cCRE_rows <- feature_data[, ID %in% ids_with_cCREs]
#   # no_cCRE_rows <- feature_data[, !ID %in% ids_with_cCREs]
#   
#   all_var <- feature_data[, ID %in% var_events[[this_event]]$`0`]
#   # high_var <- feature_data[, ID %in% var_events[[this_event]]$`0.5`]
#   # low_var <- feature_data[, ID %in% setdiff(var_events[[this_event]]$`0`, var_events[[this_event]]$`0.5`)]
#   var_list <- list(all_var=all_var)#, low_var=low_var, high_var=high_var)
#   
#   feature_data[, ID := NULL]
#   feature_data[, IHEC := NULL]
#   
#   explanatory_all <- names(feature_data)[!names(feature_data) %in% c(response, grouping_cols)]
#   
#   # make multiple models: one if an enhancer is adjacent and one without enhancers
#   #1 all data
#   
#   #2 without enhancers
#   
#   #3 enhancer without enhancer variables
#   
#   #4 enhancer with variables
#   
#   # non_enhancer_cols <- setdiff(explanatory_all, enhancer_cols)
#   # cCRE_type_cols <- sapply(cCRE_cols, function(cCRE_col_i) c(non_enhancer_cols, cCRE_col_i))
#   # cCRE_type_explanatory <- unlist(replicate(length(family), cCRE_type_cols, simplify = FALSE), recursive = FALSE)
#   # cCRE_type_explanatory <- c() #cCRE_type_explanatory[order(names(cCRE_type_explanatory))] #TODO add all cCRE Types here
#   explanatory <-
#     c(
#       replicate(length(family), list(all = explanatory_all)),
#       replicate(length(family), list(noEpigenetic = explanatory_all[Reduce(`&`, lapply(c(histone_marks, 'DNAm'), function(prefix) !startsWith(explanatory_all, prefix)))])),
#       # replicate(length(family), list(nocCRE = non_enhancer_cols)),
#       replicate(length(family), list(onlyEpigenetic = explanatory_all[Reduce(`|`, lapply(c(histone_marks, 'DNAm'), function(prefix) startsWith(explanatory_all, prefix)))]))
#       # replicate(length(family), list(all = explanatory_all)),
#       # cCRE_type_explanatory
#     )#, replicate(2, explanatory, simplify = FALSE))
#   
#   family <- c('binomial')
#   alpha <- list(c(1))
#   data_rows <- list(all=TRUE)#c(replicate(length(family)*2, list(all=TRUE)), replicate(length(family), list(nocCRE=no_cCRE_rows)), replicate(length(family) * length(cCRE_type_explanatory), list(cCRE=cCRE_rows)))#replicate(4, cCRE_rows, simplify = FALSE))
#   filter_rows <- unlist(sapply(var_list, function(var) sapply(data_rows, function(data) data & var, simplify = FALSE), simplify = FALSE), recursive = FALSE)
#   # run_glmnet(data = feature_data, explanatory = explanatory$all, response = response, filter_rows = filter_rows$high_var.nocCRE, family = family, parallel = parallel, nfolds = nfolds, alpha = alpha[[1]], grouping_cols = grouping_cols)
#   result_names <- paste(this_event, names(filter_rows), names(explanatory), family, sep='_')
#   res_list <- mapply(run_glmnet, SIMPLIFY = FALSE, USE.NAMES = TRUE, model_name=result_names, save_model=TRUE, filter_rows=filter_rows, explanatory=explanatory, family=family, response=response, data=list(feature_data), parallel=parallel, nfolds=nfolds, alpha=alpha, grouping_cols=list(grouping_cols))
#   
#   # names(res_list) <- 
#   
#   # saveRDS(res_list, paste0(this_event, '_models.rds'))
# }
# 
# # load the aggregated data with the cCRE-specific data as well
# 
# # event_specific_ids <- aggregated_dt[ID %in% ids_with_cCREs, #if(.N > 300) 
# #                                     .(N=.N), by=ID]
# glmnet.control(itrace = 0)
# # cCRE_files <- list.files(sample_dt_dir, '-cCRE_dt.csv.gz', full.names = TRUE)
# # cCRE_files <- cCRE_files[!startsWith(basename(cCRE_files), 'old')]
# # file2IHEC <- tstrsplit(basename(cCRE_files), '-', fixed=TRUE, keep = 1)[[1]]
# # names(file2IHEC) <- cCRE_files
# # value_cols <- c('wgbs;mean', 'H3K4me1;max', 'H3K27me3;max', 'H3K36me3;max', 'H3K27ac;max', 'H3K4me3;max', 'H3K9me3;max')
# # all_cCREs <- rbindlist(sapply(cCRE_files, fread, select = c('ID', value_cols), simplify = FALSE), idcol = 'IHEC', fill = TRUE)
# # all_cCREs[, IHEC:=as.factor(file2IHEC[IHEC])]
# 
# file_table <- fread('file_table.csv.gz')
# file_table[, filename:=file.path('sample_dts', basename(filename))]
# file_table[, filename:=paste0(filename, '.tab.gz')]
# file_table <- file_table[experiment_type != "DNAm"]
# file_table[, IHEC:=gsub(pattern = "\\..*", "", epirr_id)]
# setkey(file_table, "filename")
# 
# wgbs_chromhmm <- fread('sample_dts/WGBS_agg.csv.gz', skip = 'chromhmm', stringsAsFactors = TRUE)
# wgbs_chromhmm[, `V1`:=NULL]
# names(wgbs_chromhmm) <- c("IHEC", "score", "name")
# wgbs_chromhmm[, experiment_type:=factor('DNAm')]
# wgbs_chromhmm[, IHEC:=as.factor(gsub(pattern = "\\..*", "", IHEC))]
# 
# 
# # ids_to_check <- event_specific_ids[, ID]
# event_models <- pbmcapply::pbmclapply(keep_rows, function(id) {
#   tryCatch({
#     feature_data <- aggregated_dt[ID == id] #, -..enhancer_cols]
#     chromhmm_ids <- to(chromhmm_hits[from(chromhmm_hits) == id])
#     # cCREs_ids <- to(cCRE_hits)[from(cCRE_hits) == id]
#     # cCREs_this <- all_cCREs[ID %in% cCREs_ids & IHEC %in% feature_data[, IHEC]]
#     chromhmm_data <- rbindlist(
#       pbapply::pbsapply(
#         file_table[IHEC %in% feature_data[, unique(IHEC)], filename],
#         function(file) {
#           dt <-
#             fread(cmd = paste0("zcat ", file, " | egrep '", paste(paste0("chromhmm_", chromhmm_ids, "\\s"), collapse = "|"), "'"))
#           names(dt) <-
#             c('name', 'size', 'covered', 'sum', 'mean0', 'mean', 'min', 'max')
#           dt <- dt[, c('name', 'mean')]
#           dt
#         },
#         USE.NAMES = TRUE,
#         simplify = FALSE
#       ),
#       idcol = 'filename'
#     )
# 
#     chromhmm_data[, IHEC:=as.factor(file_table[filename, IHEC])]
#     chromhmm_data[, experiment_type:=as.factor(file_table[filename, experiment_type])]
#     chromhmm_data[, filename:=NULL]
#     setnames(chromhmm_data, old = 'mean', new= 'score')
#     
#     chromhmm_features <- rbindlist(list(chromhmm_data, wgbs_chromhmm[IHEC %in% feature_data[, IHEC] & name %in% paste0("chromhmm_", chromhmm_ids)]), use.names = TRUE)
#     chromhmm_features_wide <- dcast(chromhmm_features, formula = IHEC ~ experiment_type + name, value.var = "score", sep = ';')
#     
#     non_all_na_cols <- chromhmm_features_wide[, sapply(.SD, function(column) !all(is.na(column)))]
#     chromhmm_features_wide <- chromhmm_features_wide[, ..non_all_na_cols]
#     # for (j in which(!names(chromhmm_features_wide) == 'IHEC' | startsWith(names(chromhmm_features_wide), "DNAm")))
#     #   set(cCREs_this,NULL,j,log2(min(cCREs_this[[j]][cCREs_this[[j]] != 0], na.rm = TRUE) + cCREs_this[[j]]))
#     # cCREs_wide <- dcast(cCREs_this, formula = IHEC ~ ID, value.var = value_cols, fill = NA)
#     cols_to_add <- names(non_all_na_cols)[non_all_na_cols != 'IHEC']
#       #names(cCREs_wide)[Reduce(`|`, lapply(value_cols, function(prefix) startsWith(names(cCREs_wide), prefix)))]
#     feature_data[chromhmm_features_wide, on=.(IHEC), (cols_to_add):=mget(cols_to_add)]
#     # files_to_read <- list.files(file.path(psi_input_dir, 'sample_dts'), paste0('(', paste(feature_data[, IHEC], collapse = '|'), ')-cCRE_dt.csv.gz'), full.names = TRUE)
#     # cCRE_data <- rbindlist(sapply(files_to_read, function(file) fread(cmd = paste0('zgrep -e ', '\"^\\(ID\\|', paste(cCREs_ids, collapse = '\\|'), '\\),\" ', file)), simplify = FALSE), idcol = 'IHEC')
#     # cCRE_data <- rbindlist(sapply(files_to_read, function(file) fread(file)[ID %in% cCREs_ids], simplify = FALSE), idcol = 'IHEC')
#     
#     cvfit <- run_glmnet(data = feature_data, 
#                         explanatory = names(feature_data)[!names(feature_data) %in% c('IHEC', 'ID', 'Event Type', 'seqnames', 'harmonized_sample_ontology_intermediate', response)], 
#                         response = response,
#                         # filter_rows = aggregated_dt[, ID == id],
#                         nfolds = 5,
#                         grouping_cols = NULL,
#                         parallel = FALSE
#                         # alpha = c(.5, 1)
#     )
#     
#     # feature_data[, pred:=predict(cvfit$cvfit$`FALSE::OLS`)]
#     # ggplot(feature_data, aes(x=pred, y=PSI, color=harmonized_sample_ontology_intermediate)) + geom_abline(slope = 1, intercept = 0) + geom_point() + theme_bw() + theme(legend.position = 'bottom')
#   
#     cvfit}, error = function(e){e}
#   )
# })
# names(event_models) <- ids_to_check
# saveRDS(event_models, 'event_models.rds')
