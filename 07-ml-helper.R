library(caret)
library(glmnet)
library(ranger)
library(ModelMetrics)

# helper file for training the ML models
run_ML <-
  function(model_name,
           data,
           explanatory,
           response,
           test_set = TRUE,
           filter_rows = TRUE,
           scale_explanatory = TRUE,
           scale_response = FALSE,
           family = 'gaussian',
           type.measure = ifelse(family == 'gaussian', 'mse', 'BalancedAccuracy'),
           nfolds = 5,
           parallel = FALSE,
           alpha = 1,
           seed = 1234,
           grouping_cols = c("seqnames", "harmonized_sample_ontology_term_high_order_fig1"),
           save_model = FALSE,
           verbose = TRUE,
           robust_features = FALSE) {
    if (parallel) {
      require(doMC)
      registerDoMC(min(ncores, nfolds))
    }
    if (verbose) print(model_name)
    # set seed
    set.seed(seed)
    # make copy so we don't change original data by reference
    this_data <-
      copy(subset(data, filter_rows, c(explanatory, response, grouping_cols)))
    # remove cols that are all NA, maybe the mark is missing or no overlap was there
    not_all_na <-
      this_data[, sapply(.SD, function(x)
        ! (all(is.na(x)) |
             sd(x, na.rm = TRUE) == 0)), .SDcols = c(explanatory, response)]
    stopifnot(not_all_na[response])
    # filter explanatory variables
    explanatory <-
      explanatory[explanatory %in% names(not_all_na)[not_all_na]]
    # keep grouping cols 
    not_all_na <-
      c(not_all_na, setNames(rep(TRUE, length(grouping_cols)), grouping_cols))
    this_data <- this_data[, ..not_all_na]
    
    if (test_set) {
    # split in train and test
    perc_test <- .1
    
    test_by_group <- sapply(grouping_cols, function(group) {
      # fixed test sets for seqnames and cell type
      if (group == 'seqnames')
        return(c('chr1'))
      if (group == 'harmonized_sample_ontology_term_high_order_fig1')
        return(c('T lymphocyte', 'B lymphocyte'))
      # otherwise by percentage
      group_counts <-
        sample(this_data[, table(factor(get(group)))])
      names(group_counts)[seq.int(min(which(
        cumsum(group_counts) >= (perc_test * this_data[, .N])
      )))]
    }, simplify = FALSE)
    # print percentages per group
    if (verbose)
      sapply(grouping_cols, function(group)
        message(
          paste0(
            group,
            ": removing ",
            paste(test_by_group[[group]], collapse = ', '),
            " for testing. This accounts for ",
            this_data[, sum(get(group) %in% test_by_group[[group]]) / .N]
          )
        ))
    # get the test_ids per group and overall
    test_ids_by_group <-
      sapply(names(test_by_group), function(group)
        this_data[get(group) %in% test_by_group[[group]], which = TRUE])
    test_ids <- Reduce(union, test_ids_by_group)
    # create test dataset
    test_data <- copy(this_data[test_ids])
    if (length(test_ids) > 0)
      this_data <- this_data[-test_ids]
    # for binary predictions dichotomize response
    if (family == 'binomial') {
      # split in included and excluded
      dichotomization_thresholds <-
        c(1 / 3, 2 / 3) 
      if (verbose)
        message(
          paste(
            'dichotomization thresholds:',
            dichotomization_thresholds[1],
            dichotomization_thresholds[2]
          )
        )
      
      # swap response to new classes
      old_response <- response
      response <- 'binary'
      for (dt in list(this_data, test_data)) {
        dt[get(old_response) <= dichotomization_thresholds[1], (response) := control_class]
        dt[get(old_response) >= dichotomization_thresholds[2], (response) :=
             case_class]
        dt[, (response) := ordered(get(response), levels = class_levels)]
      }
      # drop events without class
      this_data <- na.omit(this_data, cols = response)
      test_data <- na.omit(test_data, cols = response)
      
      # get the test ids that were excluded from both groups
      real_test_ids <-
        Reduce(intersect, lapply(names(test_by_group), function(grouping_col)
          test_data[get(grouping_col) %in% test_by_group[[grouping_col]], which = TRUE]))
      if (verbose) {
        message(
          sprintf(
            "Number of train samples: %d, Number of test samples: %d, Percentage %.2f",
            nrow(this_data),
            nrow(test_data),
            nrow(test_data) / (nrow(this_data) + nrow(test_data))
          )
        )
        message(
          sprintf(
            "Number of train samples: %d, Number of real test samples: %d, Percentage %.2f",
            nrow(this_data),
            nrow(test_data[real_test_ids]),
            nrow(test_data[real_test_ids]) / (nrow(this_data) + nrow(test_data[real_test_ids]))
          )
        )
        message(paste('train classes', paste(this_data[, 100 * table(get(response)) /
                                                         .N], collapse = ' ')))
        message(paste('test classes', paste(test_data[, 100 * table(get(response)) /
                                                        .N], collapse = ' ')))
        message(paste('real test classes', paste(test_data[real_test_ids, 100 * table(get(response)) /
                                                             .N], collapse = ' ')))
      }
    } else {
      # remove cols that are all NA, maybe the mark is missing or no overlap was there
      # gotta do this again, because we dropped some samples
      not_all_na <-
        this_data[, sapply(.SD, function(x)
          ! (all(is.na(x)) |
               sd(x, na.rm = TRUE) == 0)), .SDcols = c(explanatory, response)]
      not_all_na[is.na(not_all_na)] <- FALSE
      stopifnot(not_all_na[response])
      explanatory <-
        explanatory[explanatory %in% names(not_all_na)[not_all_na]]
      not_all_na <-
        c(not_all_na, setNames(rep(TRUE, length(grouping_cols)), grouping_cols))
      this_data <- this_data[, ..not_all_na]
      test_data <- test_data[, ..not_all_na]
    }
    } else {
      test_by_group <- NULL
    }
    # now replace NA by mean value and scale explanatory vars
    for (j in which(names(this_data) %in% explanatory)) {
      this_mean <- mean(this_data[[j]], na.rm = TRUE)
      set(this_data, which(is.na(this_data[[j]])), j, this_mean)
      if (test_set) set(test_data, which(is.na(test_data[[j]])), j, this_mean)
    }
    # get means and sds for scaling
    means_sds <-
      melt(this_data[, c(sapply(explanatory, function(x)
        c(
          mean(get(x), na.rm = TRUE), sd(get(x), na.rm = TRUE)
        ), simplify = FALSE),
        list(name = c('mean', 'sd')))], id.vars = 'name')
    
    # scale all explanatory vars to mean = 0 and sd = 1, use training mean and sd for test set as well
    if (scale_explanatory) {
      this_data[, (explanatory) := lapply(explanatory, function(x)
        scale(get(x)))]
      if (test_set)
      test_data <-
        test_data[, (explanatory) := lapply(explanatory, function(x)
          (get(x) - means_sds[name == 'mean' &
                                variable == x, value]) / means_sds[name == 'sd' &
                                                                     variable == x, value])]
    }
    if (scale_response &
        family == 'gaussian')
      this_data[, (response) := scale(get(response))]
    
    cvfit <- list()
    if (verbose) message('finding folds')
    # here we create folds
    # if no grouping is given, we use createFolds
    if (is.null(grouping_cols)) {
      folds <-
        createFolds(this_data[, get(response)], k = nfolds, returnTrain = TRUE)
      # Initialize a vector to store test fold assignment
      foldid <- integer(this_data[, .N])
      
      # Assign each observation to the fold where it's in the test set
      for (i in seq_along(folds)) {
        test_indices <- setdiff(seq.int(this_data[, .N]), folds[[i]])
        foldid[test_indices] <- i
      }
    } else {
      # if  grouping is given, we use groupKFold and try a few times
      folds <- list()
      i <- 0
      max_tries <- 50
      while (length(folds) < nfolds) {
        if (i >= max_tries)
          stop(paste("Stopping... tried to find a grouped fold many times:", i))
        i <- i + 1
        folds <-
          caret::groupKFold(this_data[, get(grouping_cols[1])], nfolds)
      }
      # Initialize a vector to store test fold assignment
      foldid <- integer(this_data[, .N])
      
      # Assign each observation to the fold where it's in the test set
      for (i in seq_along(folds)) {
        test_indices <- setdiff(seq.int(this_data[, .N]), folds[[i]])
        foldid[test_indices] <- i
      }
    }
    
    # save the fits to a list
    fits <- list()
    
    # now do the training, we do this for the real and the randomized response vector
    for (random in c(FALSE, TRUE)) {
      if (verbose) message(sprintf("random: %s", random))
      if (random)
        this_data[, (response) := sample(get(response))] # randomize response
      
      # depending on the response we use different models
      # gaussian means event specific while binomial are the genome-wide models
      if (family == 'gaussian') {
          if (verbose) message(sprintf("glmnet"))
          # set penalty factors, such that gene expression in not penalized and always selected by the models
          penalty.factor <- rep(1, length(explanatory))
          penalty.factor[explanatory == 'gene_expression'] <- 0
          
          # fit the original models
          orig_glmnet <-
            cv.glmnet(
              as.matrix(this_data[, ..explanatory]),
              this_data[, get(response)],
              family = 'gaussian',
              type.measure = type.measure,
              foldid = foldid,
              penalty.factor = penalty.factor,
              parallel = parallel,
              keep = FALSE,
              alpha = alpha
            )
          # save the fit
          fits[[paste(random, alpha, sep = "::")]] <- orig_glmnet
          
          # get the coefficients for the best fit
          coefs <-
            coef(fits[[paste(random, alpha, sep = "::")]], s = my_lambda)
          # get the nonzero features, except Intercept and gene_expression
          nonzero_explanatory <-
            rownames(coefs)[as.matrix(coefs != 0 &
                                        !rownames(coefs) %in% c('(Intercept)', 'gene_expression'))]
          
          # if there are any nonzero features, we fit the models with only these features
          if (length(nonzero_explanatory) > 0) {
            # first per fold, to get the selected features per fold
            fold_features <- lapply(folds, function(ids) {
              net <- glmnet(as.matrix(this_data[ids, ..explanatory]),
                              this_data[ids, get(response)],
                              family = 'gaussian', #(link = link),
                              type.measure = type.measure,
                              penalty.factor = penalty.factor,
                              keep = FALSE,
                              alpha = alpha,
                              lambda = orig_glmnet$lambda)
              single_coefs <- 
                coef(net, s = orig_glmnet$lambda.1se)[
                  coef(net, s = orig_glmnet$lambda.1se)[, 1] != 0, ][-1]
              names(single_coefs)
              })
            fold_features$full <- nonzero_explanatory
            # save the features for later
            fits[[paste(random, "FEATURES", sep = "::")]] <- fold_features
            # filter the variables for the ones in at least 2 folds and the full one
            if (robust_features) {
              features_per_fold <- rbindlist(lapply(fold_features, function(x) list(feature = x)), idcol = 'fold')
              robust_feature_dt <- features_per_fold[, if('full' %in% fold & uniqueN(fold) >= 3) .(robust = uniqueN(fold)), by=.(feature)]
              nonzero_explanatory <- nonzero_explanatory[nonzero_explanatory %in% robust_feature_dt[, feature]]
              model_name <- paste(model_name, 'robust', sep = '_')
            }
            # remove gene expression, because we add it manually
            nonzero_explanatory <- 
              nonzero_explanatory[nonzero_explanatory != 'gene_expression']
            
            # fit the OLS models
            # fit the base model with only the intercept first
            fits[[paste(random, "OLS_BASE", sep = "::")]] <-
              glm(
                formula = formula(paste(
                  response,
                  '~',
                  1
                  )),
                data = this_data,
                family = 'gaussian', #(link = link),
                model = FALSE,
                x = FALSE
              )
            fits[[paste(random, "OLS_BASE", sep = "::")]]$data <-
              NULL
            # fit one with only epigenetic features
            fits[[paste(random, "OLS", sep = "::")]] <-
              glm(
                formula = formula(paste(
                  response,
                  '~',
                  ifelse(length(nonzero_explanatory) == 0, 
                         1, 
                         paste0('`', nonzero_explanatory, '`', collapse = ' + ')
                ))),
                data = this_data,
                family = 'gaussian', #(link = link),
                model = FALSE,
                x = FALSE
              )
            fits[[paste(random, "OLS", sep = "::")]]$data <-
              NULL
            # fit one with epigenetic features and gene expression
            fits[[paste(random, "OLS_GE", sep = "::")]] <-
              glm(
                formula = formula(paste(
                  response,
                  '~',
                  paste0('`', unique(c(nonzero_explanatory, 'gene_expression')),
                         '`', collapse = ' + ')
                )),
                data = this_data,
                family = "gaussian",#(link = link),
                model = FALSE,
                x = FALSE
              )
            fits[[paste(random, "OLS_GE", sep = "::")]]$data <-
              NULL
            # fit one with logit link
            fits[[paste(random, "LOGIT", sep = "::")]] <-
              glm(
                formula = formula(paste(
                  response,
                  '~',
          ifelse(length(nonzero_explanatory) == 0,
                 1, 
                 paste0('`', nonzero_explanatory, '`', collapse = ' + ')
                ))),
                data = this_data,
                family = 'quasibinomial', #(link = 'logit'),
                # weights = this_data[, 2^gene_expression],
                model = FALSE,
                x = FALSE
              )
            fits[[paste(random, "LOGIT", sep = "::")]]$data <-
              NULL
          }
      } else if (family == 'binomial') {
        
        # set up the control for caret
        myControl <- caret::trainControl(
          method = "cv",
          index = folds,
          search = 'grid',
          verboseIter = verbose,
          classProbs = TRUE,
          allowParallel = parallel,
          returnData = FALSE,
          summaryFunction = function(data,
                                     lev = NULL,
                                     model = NULL) {
            cm <- caret::confusionMatrix(data$pred, data$obs)
            # add multiple metrics
            out <-
              c(
                cm$overall,
                cm$byClass,
                twoClassSummary(data, lev, model),
                prSummary(data, lev, model),
                MCC = ModelMetrics::mcc(ifelse(data$obs == case_class, 1, 0), 
                                        data[[case_class]], 0.5)
              )
            names(out) <- gsub(' ', '', names(out), fixed = TRUE)
            return(out)
          },
          selectionFunction = "oneSE"
        )
        # set up the case weights for the imbalanced classes
        case_weights <- 
          this_data[, as.numeric(1 / (table(get(response))
                                      [get(response)] / length(get(response))))]
        if (verbose) message(sprintf("glmnet"))
        # fit the glmnet models
        model_glmnet <-
          caret::train(
            x = this_data[, ..explanatory],
            y = this_data[, get(response)],
            method = 'glmnet',
            metric = type.measure,
            trControl = myControl,
            weights = case_weights,
            tuneGrid = expand.grid(alpha = alpha,
                                   lambda = 10 ^ -(rev(
                                     seq(1, 4, length.out = 21)
                                   )))
          )
        fits[[paste(random, alpha, sep = "::")]] <- model_glmnet
        
        # double check the best lambda is the same as the one we calculated
        my_s <- max(model_glmnet$results[model_glmnet$results[[type.measure]] >= 
              (max(model_glmnet$results[, type.measure], na.rm = TRUE) -
                 model_glmnet$results[which.max(
                   model_glmnet$results[, type.measure]), 
                   paste0(type.measure, 'SD')] / sqrt(nfolds)), 'lambda'], 
              na.rm = TRUE)
        stopifnot(model_glmnet$bestTune$lambda == my_s)
        
        # get the coefficients
        coefs <-
          coef(model_glmnet$finalModel, s = my_s)
        # get the nonzero explanatory variables
        nonzero_explanatory <-
          rownames(coefs)[as.matrix(coefs != 0 &
                                      rownames(coefs) != '(Intercept)')]
        if (length(nonzero_explanatory) > 0) {
          # fit the OLS model
          fits[[paste(random, "OLS", sep = "::")]] <-
            glm(
              formula = formula(paste(
                response,
                '~',
                paste0('`', nonzero_explanatory, '`', collapse = ' + ')
              )),
              data = this_data,
              family = family,
              model = FALSE,
              y = FALSE,
              x = FALSE
            )
          fits[[paste(random, "OLS", sep = "::")]]$data <- NULL
        }
        if (verbose) message(sprintf("ranger"))
        # fit the random forest model
        model_ranger <-
          caret::train(
            x = this_data[, ..explanatory],
            y = this_data[, get(response)],
            method = 'ranger',
            metric = type.measure,
            sample.fraction = rep(this_data[, min(table(get(response)))/.N], 2),
            num.threads = ifelse(parallel, max(floor(ncores/nfolds), 1L), 1L),
            num.trees = 200,
            trControl = myControl,
            tuneGrid = expand.grid(
              mtry =  caret::var_seq(
                length(explanatory),
                classification = TRUE,
                len = 3
              ),
              splitrule = "gini",
              min.node.size = 10
            ),
            importance = "impurity"
          )
        fits[[paste(random, "RF", sep = "::")]] <- model_ranger
        
      } else
        stop('unknown family')
    }
    if (!test_set | 
        (family == 'gaussian' & !any(endsWith(names(fits), 'OLS')))) {
      test_data <- NULL
    }
    # return the results
    res_list <-
      list(
        cvfit = fits, # the models
        nsamples = nrow(this_data), # the number of samples
        test_data = test_data, # the test data
        test_by_group = test_by_group, # the test ids by group
        explanatory = explanatory, # the explanatory variables
        folds = folds, # the folds
        foldid = foldid # the fold ids
      )
    if (family == 'binomial') {
      # add the thresholds for the dichotomization
      res_list <-
        c(res_list,
          list(dichotomization_thresholds = dichotomization_thresholds))
    }
    if (scale_explanatory) {
      # add the means and sds for the explanatory variables
      res_list <- c(res_list, list(means_sds = means_sds))
    }
    if (save_model) { # save the model if requested
      if (verbose) message(sprintf("writing model ..."))
      saveRDS(object = res_list, file = paste0(model_name, '.rds'))
    }
    return(res_list)
  }

# set some default values
aggregated_dt <-
  fread('processed_data/aggregated_dt_filtered.csv.gz', stringsAsFactors = TRUE)
response <- 'PSI'
nfolds <- 10
parallel <- TRUE

# load the data
load('processed_data/aggregating.rda', verbose = TRUE)

# add the groups to the data
aggregated_dt[event_dt, on =.(ID), seqnames := seqnames] # chromosomes
aggregated_dt[metadata, # add the biological groups
              on = c(IHEC = 'epirr_id_without_version'),
              harmonized_sample_ontology_term_high_order_fig1 :=
                harmonized_sample_ontology_term_high_order_fig1]
grouping_cols <-
  c('seqnames', 'harmonized_sample_ontology_term_high_order_fig1')