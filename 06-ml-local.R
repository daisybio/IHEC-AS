source("06-ml.R")

# file_table <- fread("file_table.csv.gz")
# file_table[, filename := file.path("sample_dts", basename(filename))]
# file_table[, filename := paste0(filename, ".tab.gz")]
# file_table <- file_table[experiment_type != "DNAm"]
# file_table[, IHEC := gsub(pattern = "\\..*", "", epirr_id)]
# setkey(file_table, "filename")
# 
# wgbs_chromhmm <- fread("sample_dts/WGBS_agg.csv.gz", skip = "chromhmm", stringsAsFactors = TRUE)
# wgbs_chromhmm[, `V1` := NULL]
# names(wgbs_chromhmm) <- c("IHEC", "score", "name")
# wgbs_chromhmm[, experiment_type := factor("DNAm")]
# wgbs_chromhmm[, IHEC := as.factor(gsub(pattern = "\\..*", "", IHEC))]

nfolds_events <- 5

psi_table <- dcast(aggregated_dt, IHEC ~ ID, value.var = response)

already_computed <- sapply(list.files("event_models", pattern = ".rds"), basename)
already_computed_ids <- as.integer(gsub(pattern = "(_robust)*.rds$", replacement = "", already_computed))

data.table::setDTthreads(10)
# why not use keep_rows instead of aggregated_dt?
stopifnot(all(is.na(event_dt[ID %in% setdiff(keep_rows, aggregated_dt[, unique(ID)]), ..sample_cols])))
event_models <- # pbapply::pblapply(
  pbmcapply::pbmclapply(#mc.cores = 10,
                        X = aggregated_dt[, unique(ID)]#[!aggregated_dt[, unique(ID)] %in% already_computed_ids]
    , function(id) {
      tryCatch(
        {
          feature_table_file <- file.path("event_models", paste0("feature_table_", id, ".csv.gz"))
          if (!file.exists(feature_table_file)) {
            feature_data <- aggregated_dt[ID == id] # , -..enhancer_cols]
            chromhmm_ids <- to(chromhmm_hits[from(chromhmm_hits) == id])
            # cCREs_ids <- to(cCRE_hits)[from(cCRE_hits) == id]
            # cCREs_this <- all_cCREs[ID %in% cCREs_ids & IHEC %in% feature_data[, IHEC]]
            chromhmm_data <- rbindlist(
              pbmcmapply(mc.cores = 10,
                file_table[IHEC %in% feature_data[, unique(IHEC)], filename],
                  FUN=function(file) {
                  dt <-
                    fread(cmd = paste0("zcat ", file, " | egrep '", paste(paste0("chromhmm_", chromhmm_ids, "\\s"), collapse = "|"), "'"))
                  names(dt) <-
                    c("name", "size", "covered", "sum", "mean0", "mean", "min", "max")
                  dt <- dt[, c("name", "mean")]
                  dt
                },
                SIMPLIFY = FALSE
              ),
              idcol = "filename"
            )

            chromhmm_data[, IHEC := as.factor(file_table[filename, IHEC])]
            chromhmm_data[, experiment_type := as.factor(file_table[filename, experiment_type])]
            chromhmm_data[, filename := NULL]
            setnames(chromhmm_data, old = "mean", new = "score")

            chromhmm_features <- rbindlist(list(chromhmm_data, wgbs_chromhmm[IHEC %in% feature_data[, IHEC] & name %in% paste0("chromhmm_", chromhmm_ids)]), use.names = TRUE)
            chromhmm_features_wide <- dcast(chromhmm_features, formula = IHEC ~ experiment_type + name, value.var = "score", sep = ";")

            non_all_na_cols <- chromhmm_features_wide[, sapply(.SD, function(column) !all(is.na(column)))]
            chromhmm_features_wide <- chromhmm_features_wide[, ..non_all_na_cols]
            # for (j in which(!names(chromhmm_features_wide) == 'IHEC' | startsWith(names(chromhmm_features_wide), "DNAm")))
            #   set(cCREs_this,NULL,j,log2(min(cCREs_this[[j]][cCREs_this[[j]] != 0], na.rm = TRUE) + cCREs_this[[j]]))
            # cCREs_wide <- dcast(cCREs_this, formula = IHEC ~ ID, value.var = value_cols, fill = NA)
            cols_to_add <- names(non_all_na_cols)[non_all_na_cols != "IHEC"]
            # names(cCREs_wide)[Reduce(`|`, lapply(value_cols, function(prefix) startsWith(names(cCREs_wide), prefix)))]
            feature_data[chromhmm_features_wide, on = .(IHEC), (cols_to_add) := mget(cols_to_add)]
            # files_to_read <- list.files(file.path(psi_input_dir, 'sample_dts'), paste0('(', paste(feature_data[, IHEC], collapse = '|'), ')-cCRE_dt.csv.gz'), full.names = TRUE)
            # cCRE_data <- rbindlist(sapply(files_to_read, function(file) fread(cmd = paste0('zgrep -e ', '\"^\\(ID\\|', paste(cCREs_ids, collapse = '\\|'), '\\),\" ', file)), simplify = FALSE), idcol = 'IHEC')
            # cCRE_data <- rbindlist(sapply(files_to_read, function(file) fread(file)[ID %in% cCREs_ids], simplify = FALSE), idcol = 'IHEC')
            fwrite(feature_data, feature_table_file)
          } else {
            feature_data <- fread(feature_table_file)
          }
          
          explanatory <- names(feature_data)[!names(feature_data) %in% c("IHEC", "ID", "Event Type", grouping_cols, response)]
          
          # Apply the function to each column specified in explanatory and create a new data.table

          # cor_results <- feature_data[, rbindlist(lapply(.SD, function(x) {
          #   rmv <- is.na(x) | is.na(PSI) | is.na(gene_expression)
          #   if (all(rmv) || sd(x[!rmv]) == 0 || sd(PSI[!rmv]) == 0) {
          #     return(list(method = NA, cor = NA, pcor = NA, pcor_p = NA))
          #   } # Return NA if all values are missing
          #   # do the following for pearson and spearman correlation
          #   cor_methods <- c("pearson", "spearman")
          #   names(cor_methods) <- cor_methods
          #   rbindlist(lapply(cor_methods, function(method) {
          #     cor <- cor(x[!rmv], PSI[!rmv], use = "na.or.complete", method = method)
          #     par_cor <- ppcor::pcor.test(x = x[!rmv], y = PSI[!rmv], z = gene_expression[!rmv], method = method)
          #     list(
          #       cor = cor,
          #       pcor = par_cor$estimate,
          #       pcor_p = par_cor$p.value
          #     )
          #   }), idcol = "method")
          # }), idcol = "feature"), .SDcols = explanatory]
          # 
          # fwrite(cor_results, file.path("event_models", paste0("cor_table_", id, ".csv.gz")))
          # ggplot(cor_results[!is.na(method)], aes(x = cor, y = pcor, label = feature, color = p.adjust(pcor_p, "fdr") < .05)) +
          #   geom_abline(intercept = 0, slope = 1, color = "grey") +
          #   geom_point() +
          #   geom_label_repel() +
          #   facet_wrap(~ method) +
          #   theme_bw() + theme(legend.position = "bottom")

          # cvfit <- run_glmnet(
          #   model_name = file.path("event_models", id),
          #   data = feature_data,
          #   explanatory = explanatory, #explanatory[!explanatory %in% old_chromhmm_explanatory], 
          #   response = response,
          #   test_set = FALSE,
          #   # filter_rows = aggregated_dt[, ID == id],
          #   nfolds = nfolds_events,
          #   grouping_cols = c("harmonized_sample_ontology_term_high_order_fig1"),
          #   parallel = FALSE,
          #   save_model = TRUE,
          #   v = FALSE,
          #   robust_features = TRUE
          #   # alpha = c(.5, 1)
          # )
          
          # use smaller window:
          chromhmm_explanatory  <- explanatory[grepl('chromhmm', explanatory, fixed = TRUE)]
          
          chromhmm_hits_smaller <- findOverlaps(event_gr, activeChromHMM, maxgap = vicinity/10, ignore.strand=TRUE)
          smaller_chromhmm_ids <- to(chromhmm_hits_smaller[from(chromhmm_hits_smaller) == id])
          
          # now filter chromhmm_explanatory for the smaller chromhmm ids, i.e., only chromhmm_explanatory columns that end with chromhmm_{smaller_chromhmm_ids}
          old_chromhmm_explanatory <- chromhmm_explanatory[Reduce(`&`, lapply(paste0('chromhmm_', smaller_chromhmm_ids), function(suff) !endsWith(chromhmm_explanatory, suff)))]
          
          # cvfit_smaller <- run_glmnet(
          #   model_name = file.path("event_models_smaller", id),
          #   data = feature_data,
          #   explanatory = explanatory[!explanatory %in% old_chromhmm_explanatory], 
          #   response = response,
          #   test_set = FALSE,
          #   # filter_rows = aggregated_dt[, ID == id],
          #   nfolds = nfolds_events,
          #   grouping_cols = c("harmonized_sample_ontology_term_high_order_fig1"),
          #   parallel = FALSE,
          #   save_model = TRUE,
          #   v = FALSE,
          #   robust_features = TRUE
          #   # alpha = c(.5, 1)
          # )
          
          cvfit_local <- run_glmnet(
            model_name = file.path("event_models_local", id),
            data = feature_data,
            explanatory = explanatory[!explanatory %in% chromhmm_explanatory], 
            response = response,
            test_set = FALSE,
            # filter_rows = aggregated_dt[, ID == id],
            nfolds = nfolds_events,
            grouping_cols = c("harmonized_sample_ontology_term_high_order_fig1"),
            parallel = FALSE,
            save_model = TRUE,
            v = FALSE,
            robust_features = TRUE
            # alpha = c(.5, 1)
          )
          
          subset_psi_table <- psi_table[IHEC %in% feature_data[, IHEC]]
          other_ids <- as.integer(names(subset_psi_table[, .SD, .SDcols = colSums(is.na(subset_psi_table)) == 0])[-1])
          other_ids <- other_ids[other_ids != id & # not this id
                                   other_ids %in% event_dt[`Event Type` == event_dt[ID == id, `Event Type`], ID] & # same event type
                                   other_ids %in% event_dt[seqnames != event_dt[ID == id, seqnames], ID]] # different chromosome
          stopifnot(all(subset_psi_table[, as.character(IHEC)] == feature_data[, as.character(IHEC)]))
          set.seed(id)
          set(feature_data, j = response, value = subset_psi_table[, as.character(sample(other_ids, 1)), with = FALSE])

          # cvfit_rotated <- run_glmnet(
          #   model_name = file.path("event_models_rotated", id),
          #   data = feature_data,
          #   explanatory = explanatory,
          #   response = response,
          #   test_set = FALSE,
          #   # filter_rows = aggregated_dt[, ID == id],
          #   nfolds = nfolds_events,
          #   grouping_cols = c("harmonized_sample_ontology_term_high_order_fig1"),
          #   parallel = FALSE,
          #   save_model = TRUE,
          #   v = FALSE,
          #   robust_features = TRUE
          #   # alpha = c(.5, 1)
          # )
          # 
          # cvfit_smaller_rotated <- run_glmnet(
          #   model_name = file.path("event_models_smaller_rotated", id),
          #   data = feature_data,
          #   explanatory = explanatory[!explanatory %in% old_chromhmm_explanatory],
          #   response = response,
          #   test_set = FALSE,
          #   # filter_rows = aggregated_dt[, ID == id],
          #   nfolds = nfolds_events,
          #   grouping_cols = c("harmonized_sample_ontology_term_high_order_fig1"),
          #   parallel = FALSE,
          #   save_model = TRUE,
          #   v = FALSE,
          #   robust_features = TRUE
          #   # alpha = c(.5, 1)
          # )
          
          cvfit_local_rotated <- run_glmnet(
            model_name = file.path("event_models_local_rotated", id),
            data = feature_data,
            explanatory = explanatory[!explanatory %in% chromhmm_explanatory], 
            response = response,
            test_set = FALSE,
            # filter_rows = aggregated_dt[, ID == id],
            nfolds = nfolds_events,
            grouping_cols = c("harmonized_sample_ontology_term_high_order_fig1"),
            parallel = FALSE,
            save_model = TRUE,
            v = FALSE,
            robust_features = TRUE
            # alpha = c(.5, 1)
          )
          # feature_data[, pred:=predict(cvfit$cvfit$`FALSE::OLS`)]
          # ggplot(feature_data, aes(x=pred, y=PSI, b=harmonized_sample_ontology_intermediate)) + geom_abline(slope = 1, intercept = 0) + geom_point() + theme_bw() + theme(legend.position = 'bottom')

          NULL
        },
        error = function(e) {
          e$message
        }
      )
    }
  )

# check for the models that did not work
# names(event_models) <- aggregated_dt[, unique(ID)][!aggregated_dt[, unique(ID)] %in% already_computed_ids]
# unique(unlist(event_models))

# sd_problems <- aggregated_dt[ID %in% names(event_models)[event_models == "not_all_na[response] is not TRUE"], .(sd=sd(get(response))), by=ID]
# sd_problems[, table(sd, useNA = "ifany")]
# stopifnot(sd_problems[, all(sd == 0, na.rm = TRUE)])
# 
# group_problems <- aggregated_dt[ID %in% names(event_models)[startsWith(unlist(event_models), "`k` should be less than")], .(unique_ontology=uniqueN(harmonized_sample_ontology_intermediate)), by=ID]
# group_problems[, table(unique_ontology)]
# stopifnot(group_problems[, all(unique_ontology < nfolds_events)])
# 
# constant_y <- aggregated_dt[ID %in% names(event_models)[event_models == "y is constant; gaussian glmnet fails at standardization step"], uniqueN(get(response)), by = ID]
# aggregated_dt[ID == 18425, .(unique(get(response)), uniqueN(get(response))), by=.(harmonized_sample_ontology_intermediate)]
# aggregated_dt[ID == 15100, .(unique(get(response)), uniqueN(get(response))), by=.(harmonized_sample_ontology_intermediate)]
# aggregated_dt[ID == 13433, .(unique(get(response)), uniqueN(get(response))), by=.(harmonized_sample_ontology_intermediate)]

