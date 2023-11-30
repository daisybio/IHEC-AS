source('06-ml.R')

file_table <- fread('file_table.csv.gz')
file_table[, filename:=file.path('sample_dts', basename(filename))]
file_table[, filename:=paste0(filename, '.tab.gz')]
file_table <- file_table[experiment_type != "DNAm"]
file_table[, IHEC:=gsub(pattern = "\\..*", "", epirr_id)]
setkey(file_table, "filename")

wgbs_chromhmm <- fread('sample_dts/WGBS_agg.csv.gz', skip = 'chromhmm', stringsAsFactors = TRUE)
wgbs_chromhmm[, `V1`:=NULL]
names(wgbs_chromhmm) <- c("IHEC", "score", "name")
wgbs_chromhmm[, experiment_type:=factor('DNAm')]
wgbs_chromhmm[, IHEC:=as.factor(gsub(pattern = "\\..*", "", IHEC))]


# ids_to_check <- event_specific_ids[, ID]
already_computed <- sapply(list.files('event_models'), basename)
already_computed_ids <- as.integer(sub(pattern = '.rds', replacement = '', already_computed, fixed = TRUE))
event_models <- #pbmcapply::pbmc
lapply(keep_rows[!keep_rows %in% already_computed_ids], function(id) {
  tryCatch({
    feature_data <- aggregated_dt[ID == id] #, -..enhancer_cols]
    chromhmm_ids <- to(chromhmm_hits[from(chromhmm_hits) == id])
    # cCREs_ids <- to(cCRE_hits)[from(cCRE_hits) == id]
    # cCREs_this <- all_cCREs[ID %in% cCREs_ids & IHEC %in% feature_data[, IHEC]]
    if (feature_data[, uniqueN(PSI) > 1]) browser()
    else return(feature_data[, PSI])
    chromhmm_data <- rbindlist(
      sapply(
        file_table[IHEC %in% feature_data[, unique(IHEC)], filename],
        function(file) {
          dt <-
            fread(cmd = paste0("zcat ", file, " | egrep '", paste(paste0("chromhmm_", chromhmm_ids, "\\s"), collapse = "|"), "'"))
          names(dt) <-
            c('name', 'size', 'covered', 'sum', 'mean0', 'mean', 'min', 'max')
          dt <- dt[, c('name', 'mean')]
          dt
        },
        USE.NAMES = TRUE,
        simplify = FALSE
      ),
      idcol = 'filename'
    )

    chromhmm_data[, IHEC:=as.factor(file_table[filename, IHEC])]
    chromhmm_data[, experiment_type:=as.factor(file_table[filename, experiment_type])]
    chromhmm_data[, filename:=NULL]
    setnames(chromhmm_data, old = 'mean', new= 'score')
    
    chromhmm_features <- rbindlist(list(chromhmm_data, wgbs_chromhmm[IHEC %in% feature_data[, IHEC] & name %in% paste0("chromhmm_", chromhmm_ids)]), use.names = TRUE)
    chromhmm_features_wide <- dcast(chromhmm_features, formula = IHEC ~ experiment_type + name, value.var = "score", sep = ';')
    
    non_all_na_cols <- chromhmm_features_wide[, sapply(.SD, function(column) !all(is.na(column)))]
    chromhmm_features_wide <- chromhmm_features_wide[, ..non_all_na_cols]
    # for (j in which(!names(chromhmm_features_wide) == 'IHEC' | startsWith(names(chromhmm_features_wide), "DNAm")))
    #   set(cCREs_this,NULL,j,log2(min(cCREs_this[[j]][cCREs_this[[j]] != 0], na.rm = TRUE) + cCREs_this[[j]]))
    # cCREs_wide <- dcast(cCREs_this, formula = IHEC ~ ID, value.var = value_cols, fill = NA)
    cols_to_add <- names(non_all_na_cols)[non_all_na_cols != 'IHEC']
      #names(cCREs_wide)[Reduce(`|`, lapply(value_cols, function(prefix) startsWith(names(cCREs_wide), prefix)))]
    feature_data[chromhmm_features_wide, on=.(IHEC), (cols_to_add):=mget(cols_to_add)]
    # files_to_read <- list.files(file.path(psi_input_dir, 'sample_dts'), paste0('(', paste(feature_data[, IHEC], collapse = '|'), ')-cCRE_dt.csv.gz'), full.names = TRUE)
    # cCRE_data <- rbindlist(sapply(files_to_read, function(file) fread(cmd = paste0('zgrep -e ', '\"^\\(ID\\|', paste(cCREs_ids, collapse = '\\|'), '\\),\" ', file)), simplify = FALSE), idcol = 'IHEC')
    # cCRE_data <- rbindlist(sapply(files_to_read, function(file) fread(file)[ID %in% cCREs_ids], simplify = FALSE), idcol = 'IHEC')
    
    cvfit <- run_glmnet(model_name = file.path('event_models', id),
                        data = feature_data, 
                        explanatory = names(feature_data)[!names(feature_data) %in% c('IHEC', 'ID', 'Event Type', 'seqnames', 'harmonized_sample_ontology_intermediate', response)], 
                        response = response,
                        # filter_rows = aggregated_dt[, ID == id],
                        nfolds = 5,
                        grouping_cols = NULL,
                        parallel = FALSE,
                        save_model=TRUE
                        # alpha = c(.5, 1)
    )
    
    # feature_data[, pred:=predict(cvfit$cvfit$`FALSE::OLS`)]
    # ggplot(feature_data, aes(x=pred, y=PSI, color=harmonized_sample_ontology_intermediate)) + geom_abline(slope = 1, intercept = 0) + geom_point() + theme_bw() + theme(legend.position = 'bottom')
  
    cvfit}, error = function(e){e}
  )
})
