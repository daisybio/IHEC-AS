setwd('~/hiwi/IHEC-AS/')
source('.Rprofile')

library(glmnet)
glmnet.control(itrace = 1)

run_glmnet <- function(data, explanatory, response, filter_rows = TRUE, scale_explanatory = TRUE, scale_response = TRUE, family = 'gaussian', type.measure = ifelse(family == 'gaussian', 'mse', 'auc'), nfolds=10, parallel=FALSE, name=NULL, alpha = 1) {
  # make copy so we don't change original data by reference
  this_data <- copy(subset(data, filter_rows, c(explanatory, response)))
  
  # remove cols that are all NA, maybe the mark is missing or no overlap was there
  not_all_na <- this_data[, sapply(.SD, function(x) !all(is.na(x)))]
  explanatory <- explanatory[explanatory %in% names(not_all_na)[not_all_na]]
  this_data <- this_data[, ..not_all_na]
  
  # get only rows with full data
  # do I want this?
  this_data <- na.omit(this_data)
  if (family == 'binomial') {
    dichotomization_thresholds <- this_data[, quantile(get(response), seq(0, 1, 1/3))]
    message(paste('dichotomization thresholds:', dichotomization_thresholds[2], dichotomization_thresholds[3]))
    this_data[get(response) < dichotomization_thresholds[2], binary:='0']
    this_data[get(response) > dichotomization_thresholds[3], binary:='1']
    this_data <- na.omit(this_data)
    response <- 'binary'
    this_data[, binary:=as.factor(binary)]
    # Create test set here
    nrows <- this_data[, .N]
    test_ids <- sample(seq(nrows), .2*nrows)
    test_data <- this_data[test_ids]
    this_data <- this_data[-test_ids]
    message(paste('train classes', paste(this_data[, 100 * table(binary)/(nrows - length(test_ids))], collapse = ' ')))
    message(paste('test classes', paste(test_data[, 100 * table(binary)/length(test_ids)], collapse = ' ')))
    means_sds <- melt(this_data[, c(sapply(explanatory, function(x) c(mean(get(x)), sd(get(x))), simplify = FALSE), list(name=c('mean', 'sd')))], id.vars = 'name')
  }
  
  # scale all explanatory vars to mean = 0 and sd = 1
  if (scale_explanatory) this_data[, (explanatory) := lapply(explanatory, function(x) scale(get(x)))]
  if (scale_response & family == 'gaussian') this_data[, (response) := scale(get(response))]
  # make cv
  foldid <- sample(seq.int(nfolds), size = nrow(this_data), replace = TRUE)
  cvfit <- lapply(alpha, function(a) cv.glmnet(as.matrix(this_data[, ..explanatory]), this_data[, get(response)], family = family, type.measure = type.measure, foldid = foldid, parallel = parallel, keep = TRUE, alpha = a))
  names(cvfit) <- alpha
  # nonzero_explanatory <- rownames(coef(cvfit[[max(alpha)]], s = 'lambda.1se'))[as.matrix(coef(cvfit$`1`, s = 'lambda.1se')) != 0 & rownames(coef(cvfit$`1`, s = 'lambda.1se')) != '(Intercept)']
  # if (length(nonzero_explanatory) == 0) 
  nonzero_explanatory <- rownames(coef(cvfit[[max(alpha)]], s = 'lambda.1se'))[as.matrix(coef(cvfit$`1`, s = 'lambda.min')) != 0 & rownames(coef(cvfit$`1`, s = 'lambda.min')) != '(Intercept)']
  cvfit$OLS <- glm(formula = formula(paste(response, '~', paste0('`', nonzero_explanatory, '`', collapse = ' + '))), data = this_data, family = family)
  res_list <- list(cvfit=cvfit, nsamples=nrow(this_data))
  if (family == 'binomial') {
    if (scale_explanatory) {
      test_data <- test_data[, (explanatory):=lapply(explanatory, function(x) (get(x) - means_sds[name == 'mean' & variable == x, value])/means_sds[name == 'sd' & variable == x, value])]
    }
    return(c(res_list, list(test_data=test_data, explanatory=explanatory)))
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

#
for (this_event in aggregated_dt[, levels(event_name)]){
  print(this_event)
  feature_data <- aggregated_dt[this_event == event_name, -c('IHEC', 'event_name')]

    #names(feature_data)[names(feature_data) != response & keep_cols(feature_data, aggregation) & !grepl('tile', names(feature_data), fixed = TRUE) & !endsWith(names(feature_data), 'percentage') & !endsWith(names(feature_data), 'count') & !endsWith(names(feature_data), 'presence')]

  cCRE_rows <- feature_data[, ID %in% ids_with_cCREs]
  no_cCRE_rows <- feature_data[, !ID %in% ids_with_cCREs]
  feature_data[, ID := NULL]

  explanatory <- names(feature_data)[names(feature_data) != response]

  # make multiple models: one if an enhancer is adjacent and one without enhancers
  #1 all data

  #2 without enhancers

  #3 enhancer without enhancer variables

  #4 enhancer with variables
  family <- c('gaussian', 'binomial')
  explanatory <- c(replicate(6, intersect(explanatory, names(feature_data)[!names(feature_data) %in% enhancer_cols]), simplify = FALSE))#, replicate(2, explanatory, simplify = FALSE))
  filter_rows <- c(replicate(2, TRUE, simplify = FALSE), replicate(2, no_cCRE_rows, simplify = FALSE), replicate(2, cCRE_rows, simplify = FALSE))#replicate(4, cCRE_rows, simplify = FALSE))
  alpha <- list(c(.5, 1))
  res_list <- mapply(run_glmnet, USE.NAMES = TRUE, filter_rows=filter_rows, explanatory=explanatory, family=family, response=response, data=list(feature_data), parallel=parallel, nfolds=nfolds, alpha=alpha)
  names(res_list) <- paste(c(rep('all', 2), rep('noEnhancer', 2), rep('enhancer', 4)), c(rep('noEnhancer', 6), rep('enhancer', 2)), c('psi', 'binary'), sep = '_')[1:length(explanatory)]
  saveRDS(res_list, paste0(this_event, '_models.rds'))
}

foreach::registerDoSEQ()

# load the aggregated data with the cCRE-specific data as well

event_specific_ids <- aggregated_dt[ID %in% ids_with_cCREs, #if(.N > 300) 
  .(N=.N), by=ID]
glmnet.control(itrace = 0)
cCRE_files <- list.files(file.path(psi_input_dir, 'sample_dts'), '-cCRE_dt.csv.gz', full.names = TRUE)
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
    cCREs_wide <- dcast(cCREs_this, formula = IHEC ~ ID, value.var = value_cols)
    cols_to_add <- names(cCREs_wide)[Reduce(`|`, lapply(value_cols, function(prefix) startsWith(names(cCREs_wide), prefix)))]
    feature_data[cCREs_wide, on=.(IHEC), (cols_to_add):=mget(cols_to_add)]
    # files_to_read <- list.files(file.path(psi_input_dir, 'sample_dts'), paste0('(', paste(feature_data[, IHEC], collapse = '|'), ')-cCRE_dt.csv.gz'), full.names = TRUE)
    # cCRE_data <- rbindlist(sapply(files_to_read, function(file) fread(cmd = paste0('zgrep -e ', '\"^\\(ID\\|', paste(cCREs_ids, collapse = '\\|'), '\\),\" ', file)), simplify = FALSE), idcol = 'IHEC')
    # cCRE_data <- rbindlist(sapply(files_to_read, function(file) fread(file)[ID %in% cCREs_ids], simplify = FALSE), idcol = 'IHEC')
    cvfit <- run_glmnet(data = feature_data, 
                        explanatory = names(feature_data)[!names(feature_data) %in% c(c('IHEC', 'ID', 'event_name'), response)], 
                        response = response,
                        # filter_rows = aggregated_dt[, ID == id],
                        nfolds = 5,
                        alpha = c(.5, 1))
    cvfit}, error = function(e){e}
  )
})
names(event_models) <- ids_to_check
saveRDS(event_models, 'event_models.rds')