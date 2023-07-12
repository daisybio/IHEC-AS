# setwd('~/hiwi/IHEC-AS/')
source('.Rprofile')

library(glmnet)

run_glmnet <- function(data, explanatory, response, filter_rows = TRUE, scale_explanatory = TRUE, scale_response = FALSE, family = 'gaussian', type.measure = ifelse(family == 'gaussian', 'mse', 'auc'), nfolds=10, parallel=FALSE, name=NULL, alpha = 1, seed = 1234) {
  # make copy so we don't change original data by reference
  this_data <- copy(subset(data, filter_rows, c(explanatory, response)))
  # remove cols that are all NA, maybe the mark is missing or no overlap was there
  not_all_na <- this_data[, sapply(.SD, function(x) !all(is.na(x)) & sd(x, na.rm = TRUE) != 0)]
  stopifnot(not_all_na[response])
  explanatory <- explanatory[explanatory %in% names(not_all_na)[not_all_na]]
  this_data <- this_data[, ..not_all_na]
  
  # get only rows with full data
  # do I want this?
  # this_data <- na.omit(this_data)
  if (family == 'binomial') {
    # split in train and test
    nrows <- this_data[, .N]
    set.seed(seed)
    test_ids <- sample(seq(nrows), .2*nrows)
    test_data <- this_data[test_ids]
    this_data <- this_data[-test_ids]
    
    # split in included and excluded
    dichotomization_thresholds <- this_data[, quantile(get(response), c(1/3, 2/3))]
    message(paste('dichotomization thresholds:', dichotomization_thresholds[1], dichotomization_thresholds[2]))
    
    old_response <- response
    response <- 'binary'
    for (dt in list(this_data, test_data)) {
      dt[get(old_response) <= dichotomization_thresholds[1], (response):=factor('0')]
      dt[get(old_response) >= dichotomization_thresholds[2], (response):=factor('1')]
    }
    this_data <- na.omit(this_data, cols=response)
    test_data <- na.omit(test_data, cols=response)
    
    message(paste('train classes', paste(this_data[, 100 * table(get(response))/.N], collapse = ' ')))
    message(paste('test classes', paste(test_data[, 100 * table(get(response))/.N], collapse = ' ')))
    
    # now replace NA by mean value and 
    for (j in which(names(this_data) %in% explanatory)) {
      this_mean <- mean(this_data[[j]], na.rm = TRUE)
      set(this_data, which(is.na(this_data[[j]])), j, this_mean)
      set(test_data, which(is.na(test_data[[j]])), j, this_mean)
    }
    means_sds <- melt(this_data[, c(sapply(explanatory, function(x) c(mean(get(x), na.rm = TRUE), sd(get(x), na.rm=TRUE)), simplify = FALSE), list(name=c('mean', 'sd')))], id.vars = 'name')
  }
  
  # scale all explanatory vars to mean = 0 and sd = 1
  if (scale_explanatory) {
    this_data[, (explanatory) := lapply(explanatory, function(x) scale(get(x)))]
  }
  if (scale_response & family == 'gaussian') this_data[, (response) := scale(get(response))]
  # make cv
  set.seed(seed)
  foldid <- sample(seq.int(nfolds), size = nrow(this_data), replace = TRUE)
  cvfit <- lapply(alpha, function(a) cv.glmnet(as.matrix(this_data[, ..explanatory]), this_data[, get(response)], family = family, type.measure = type.measure, foldid = foldid, parallel = parallel, keep = FALSE, alpha = a))
  names(cvfit) <- alpha
  browser()
  nonzero_explanatory <- rownames(coef(cvfit[[max(alpha)]], s = my_lambda))[as.matrix(coef(cvfit[[max(alpha)]], s = my_lambda)) != 0 & rownames(coef(cvfit[[max(alpha)]], s = my_lambda)) != '(Intercept)']
  cvfit$OLS <- glm(formula = formula(paste(response, '~', paste0('`', nonzero_explanatory, '`', collapse = ' + '))), data = this_data, family = family, model = FALSE, y = FALSE)
  res_list <- list(cvfit=cvfit, nsamples=nrow(this_data))
  if (family == 'binomial') {
    if (scale_explanatory) {
      test_data <- test_data[, (explanatory):=lapply(explanatory, function(x) (get(x) - means_sds[name == 'mean' & variable == x, value])/means_sds[name == 'sd' & variable == x, value])]
      return(c(res_list, list(test_data=test_data, explanatory=explanatory, dichotomization_thresholds=dichotomization_thresholds, 
    return(c(res_list, list(test_data=test_data, explanatory=explanatory, dichotomization_thresholds=dichotomization_thresholds, means_sds=means_sds))))))
    }
    return(c(res_list, list(test_data=test_data, explanatory=explanatory, dichotomization_thresholds=dichotomization_thresholds)))
  }
  return(res_list)
}

aggregated_dt <- fread('aggregated_dt_filtered.csv.gz', stringsAsFactors = TRUE)
response <- 'PSI'
nfolds <- 10
parallel <- TRUE

load('aggregating.rda', verbose=TRUE)
cCREs <- data.table::fread('data/GRCh38-cCREs.bed')
names(cCREs) <- c('seqnames', 'start', 'end', 'some_id', 'accession', 'cCRE_type')
cCRE_gr <- cCREs[, GRanges(seqnames = seqnames, IRanges(start = start, end = end))]
cCRE_hits <- findOverlaps(event_gr, cCRE_gr, maxgap = 5000, ignore.strand=TRUE)

ids_with_cCREs <- aggregated_dt[, unique(ID)[unique(ID) %in% from(cCRE_hits)]]

if (parallel) {
  require(doParallel)
  registerDoParallel(min(nfolds, options()$mc.cores))
}

enhancer_cols <- names(aggregated_dt)[grepl('(CTCF)|(DNase)|(ELS)|(PLS)', names(aggregated_dt))]

cCRE_regions <- c('CTCF-only', 'dELS', 'pELS', 'PLS', 'DNase-H3K4me3')
cCRE_cols <- sapply(cCRE_regions, function(cCRE_type) names(aggregated_dt)[grepl(cCRE_type, names(aggregated_dt))])

for (this_event in aggregated_dt[, levels(`Event Type`)]){
  print(this_event)
  feature_data <- aggregated_dt[this_event == `Event Type`, -c('IHEC', 'Event Type')]
  # feature_data[, cCRE_near:=ID %in% ids_with_cCREs]
  #names(feature_data)[names(feature_data) != response & keep_cols(feature_data, aggregation) & !grepl('tile', names(feature_data), fixed = TRUE) & !endsWith(names(feature_data), 'percentage') & !endsWith(names(feature_data), 'count') & !endsWith(names(feature_data), 'presence')]
  
  cCRE_rows <- feature_data[, ID %in% ids_with_cCREs]
  no_cCRE_rows <- feature_data[, !ID %in% ids_with_cCREs]
  
  all_var <- feature_data[, ID %in% var_events[[this_event]]$`0`]
  high_var <- feature_data[, ID %in% var_events[[this_event]]$`0.5`]
  low_var <- feature_data[, ID %in% setdiff(var_events[[this_event]]$`0`, var_events[[this_event]]$`0.5`)]
  var_list <- list(all_var=all_var, low_var=low_var, high_var=high_var)
  
  feature_data[, ID := NULL]
  
  
  explanatory_all <- names(feature_data)[names(feature_data) != response]
  
  # make multiple models: one if an enhancer is adjacent and one without enhancers
  #1 all data
  
  #2 without enhancers
  
  #3 enhancer without enhancer variables
  
  #4 enhancer with variables
  
  non_enhancer_cols <- setdiff(explanatory_all, enhancer_cols)
  cCRE_type_cols <- sapply(cCRE_cols, function(cCRE_col_i) c(non_enhancer_cols, cCRE_col_i))
  cCRE_type_explanatory <- unlist(replicate(length(family), cCRE_type_cols, simplify = FALSE), recursive = FALSE)
  cCRE_type_explanatory <- cCRE_type_explanatory[order(names(cCRE_type_explanatory))]
  explanatory <- c(replicate(length(family), list(nocCRE=non_enhancer_cols)), replicate(length(family)*2, list(all=explanatory_all)), cCRE_type_explanatory)#, replicate(2, explanatory, simplify = FALSE))
  
  family <- c('binomial')
  alpha <- list(c(1))
  data_rows <- c(replicate(length(family)*2, list(all=TRUE)), replicate(length(family), list(nocCRE=no_cCRE_rows)), replicate(length(family) * length(cCRE_type_explanatory), list(cCRE=cCRE_rows)))#replicate(4, cCRE_rows, simplify = FALSE))
  filter_rows <- unlist(sapply(var_list, function(var) sapply(data_rows, function(data) data & var, simplify = FALSE), simplify = FALSE), recursive = FALSE)
  
  res_list <- mapply(run_glmnet, SIMPLIFY = FALSE, USE.NAMES = TRUE, filter_rows=filter_rows, explanatory=explanatory, family=family, response=response, data=list(feature_data), parallel=parallel, nfolds=nfolds, alpha=alpha)
  
  names(res_list) <- paste(names(filter_rows), names(explanatory), family, sep='_')
  saveRDS(res_list, paste0(this_event, '_models.rds'))
}

foreach::registerDoSEQ()

# load the aggregated data with the cCRE-specific data as well

event_specific_ids <- aggregated_dt[ID %in% ids_with_cCREs, #if(.N > 300) 
                                    .(N=.N), by=ID]
glmnet.control(itrace = 0)
cCRE_files <- list.files(sample_dt_dir, '-cCRE_dt.csv.gz', full.names = TRUE)
cCRE_files <- cCRE_files[!startsWith(basename(cCRE_files), 'old')]
file2IHEC <- tstrsplit(basename(cCRE_files), '-', fixed=TRUE, keep = 1)[[1]]
names(file2IHEC) <- cCRE_files
value_cols <- c('wgbs;mean', 'H3K4me1;max', 'H3K27me3;max', 'H3K36me3;max', 'H3K27ac;max', 'H3K4me3;max', 'H3K9me3;max')
all_cCREs <- rbindlist(sapply(cCRE_files, fread, select = c('ID', value_cols), simplify = FALSE), idcol = 'IHEC', fill = TRUE)
all_cCREs[, IHEC:=as.factor(file2IHEC[IHEC])]


ids_to_check <- event_specific_ids[, ID]
event_models <- pbmcapply::pbmclapply(ids_to_check, function(id) {
  tryCatch({
    feature_data <- aggregated_dt[ID == id, -..enhancer_cols]
    cCREs_ids <- to(cCRE_hits)[from(cCRE_hits) == id]
    cCREs_this <- all_cCREs[ID %in% cCREs_ids & IHEC %in% feature_data[, IHEC]]
    non_all_na_cols <- cCREs_this[, sapply(.SD, function(column) !all(is.na(column)))]
    cCREs_this <- cCREs_this[, ..non_all_na_cols]
    for (j in which(!names(cCREs_this) %in% c('IHEC', 'ID')))
      set(cCREs_this,NULL,j,log2(min(cCREs_this[[j]][cCREs_this[[j]] != 0], na.rm = TRUE) + cCREs_this[[j]]))
    cCREs_wide <- dcast(cCREs_this, formula = IHEC ~ ID, value.var = value_cols, fill = NA)
    cols_to_add <- names(cCREs_wide)[Reduce(`|`, lapply(value_cols, function(prefix) startsWith(names(cCREs_wide), prefix)))]
    feature_data[cCREs_wide, on=.(IHEC), (cols_to_add):=mget(cols_to_add)]
    # files_to_read <- list.files(file.path(psi_input_dir, 'sample_dts'), paste0('(', paste(feature_data[, IHEC], collapse = '|'), ')-cCRE_dt.csv.gz'), full.names = TRUE)
    # cCRE_data <- rbindlist(sapply(files_to_read, function(file) fread(cmd = paste0('zgrep -e ', '\"^\\(ID\\|', paste(cCREs_ids, collapse = '\\|'), '\\),\" ', file)), simplify = FALSE), idcol = 'IHEC')
    # cCRE_data <- rbindlist(sapply(files_to_read, function(file) fread(file)[ID %in% cCREs_ids], simplify = FALSE), idcol = 'IHEC')
    cvfit <- run_glmnet(data = feature_data, 
                        explanatory = names(feature_data)[!names(feature_data) %in% c('IHEC', 'ID', 'Event Type', response)], 
                        response = response,
                        # filter_rows = aggregated_dt[, ID == id],
                        nfolds = 5,
                        alpha = c(.5, 1))
    cvfit}, error = function(e){e}
  )
})
names(event_models) <- ids_to_check
saveRDS(event_models, 'event_models.rds')
