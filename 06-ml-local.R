source("06-ml.R")

# # only needed when I need to build new feauture_tables
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
nrotations <- 10

event_dirs <- c('event_models', 'event_models_smaller', 'event_models_local')
already_computed_ids <- lapply(event_dirs, function(event_dir) {
  # first create event_dir, without showing warnings
  dir.create(event_dir, showWarnings = FALSE)
  already_computed <- sapply(list.files(event_dir, pattern = ".rds"), basename)
  already_computed_ids <- as.integer(gsub(pattern = "(_robust)*.rds$", replacement = "", already_computed))
  already_computed_ids
})
# make sure that the intersection of all those is the same as the first one
stopifnot(identical(already_computed_ids[[1]], Reduce(intersect, already_computed_ids)))


# data.table::setDTthreads(10)

# why not use keep_rows_manual instead of aggregated_dt?
stopifnot(all(is.na(event_dt[ID %in% setdiff(keep_rows_manual, aggregated_dt[, unique(ID)]), ..sample_cols])))

# set up ids to build models from: 
ids_to_build <- aggregated_dt[, unique(ID)]
# remove ids with less that minimum_events in aggregated_dt
ids_to_build <- ids_to_build[ids_to_build %in% aggregated_dt[ID %in% ids_to_build, .N, by=ID][N >= minimum_events, ID]]
# remove ids where the responses standard deviation is not zero
ids_to_build <- ids_to_build[ids_to_build %in% aggregated_dt[ID %in% ids_to_build, .(sd = sd(get(response))), by=ID][sd > 0, ID]]
event_dt[ID %in% ids_to_build, table(`Event Type`)]

psi_table <- dcast(aggregated_dt[ID %in% intersect(keep_rows, ids_to_build)], IHEC ~ ID, value.var = response)

# recompute only the models that have not been computed yet
ids_to_build <- setdiff(ids_to_build, already_computed_ids[[1]])

event_models <- #pbapply::pblapply(
  pbmcapply::pbmclapply(#mc.cores = 10,
                        X = ids_to_build
    , function(id) {
      tryCatch(
        {
          feature_table_file <- file.path("event_feature_tables", paste0("feature_table_", id, ".csv.gz"))
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
          
          explanatory <- names(feature_data)[!names(feature_data) %in% c("IHEC", "ID", "Event Type", "Variability", "gene_id", grouping_cols, response)]
          ### distal elements of multiple distance:
          chromhmm_explanatory  <- explanatory[grepl('chromhmm', explanatory, fixed = TRUE)]
          # use smaller window:
          chromhmm_hits_smaller <- findOverlaps(event_gr, activeChromHMM, maxgap = vicinity/10, ignore.strand=TRUE)
          smaller_chromhmm_ids <- to(chromhmm_hits_smaller[from(chromhmm_hits_smaller) == id])
          
          # now filter chromhmm_explanatory for the smaller chromhmm ids, i.e., only chromhmm_explanatory columns that end with chromhmm_{smaller_chromhmm_ids}
          old_chromhmm_explanatory <- chromhmm_explanatory[Reduce(`&`, lapply(paste0('chromhmm_', smaller_chromhmm_ids), function(suff) !endsWith(chromhmm_explanatory, suff)))]
          
          ### Rotated responses:
          # same IHEC IDs available for all events, also remove IHEC column
          subset_psi_matrix <- as.matrix(psi_table[IHEC %in% feature_data[, IHEC]], rownames = "IHEC")
          # all IDs with no NAs and with sd > 0
          other_ids <- as.integer(colnames(subset_psi_matrix)[colSums(is.na(subset_psi_matrix)) == 0 &
                                          apply(subset_psi_matrix, 2, sd, na.rm = TRUE) > 0])
          other_ids <- other_ids[other_ids != id & # not this id
                                   other_ids %in% event_dt[`Event Type` == event_dt[ID == id, `Event Type`], ID] & # same event type
                                   other_ids %in% event_dt[seqnames != event_dt[ID == id, seqnames], ID] & # different chromosome
                                   other_ids %in% event_dt[Variability == event_dt[ID == id, Variability], ID] # same variability
                                 ]
          this_rotations <- nrotations
          if (length(other_ids) < nrotations) {
            this_rotations <- length(other_ids)
            warning(paste0('Not enough other IDs for ', id, ', using ', this_rotations, ' rotations'))
          }
          # stopifnot(length(other_ids) >= nrotations)
          stopifnot(all(rownames(subset_psi_matrix) == feature_data[, as.character(IHEC)]))
          set.seed(id)
          rotated_psis <- subset_psi_matrix[, as.character(sample(other_ids, this_rotations))]
          orig_psi <- feature_data[, get(response)]
          
          explanatory_vars <- list(explanatory,
                                   explanatory[!explanatory %in% old_chromhmm_explanatory], 
                                   explanatory[!explanatory %in% chromhmm_explanatory])
          
          ### actually compute models
          mapply(event_dirs, explanatory_vars, FUN = function(event_dir, this_explanatory) {
            
            if (!dir.exists(event_dir))
              dir.create(event_dir)
            filtered_dir <- paste0(event_dir, '_filtered')
            if (!dir.exists(filtered_dir))
              dir.create(filtered_dir)
            rotated_dir <- paste0(event_dir, '_rotated')
            if (!dir.exists(rotated_dir))
              dir.create(rotated_dir)
            for (rotation_i in seq.int(this_rotations))
              if (!dir.exists(file.path(rotated_dir, rotation_i)))
                dir.create(file.path(rotated_dir, rotation_i))
            
            trans_features <- c()
            for (rotation_i in seq.int(this_rotations)) {
              
              set(feature_data, j = response, value = rotated_psis[, rotation_i])
              new_id <- colnames(rotated_psis)[rotation_i]
              tryCatch({
              cvfit_rotated <- run_glmnet(
                model_name = file.path(rotated_dir, rotation_i, paste(id, new_id, sep = 'to')),
                data = feature_data,
                explanatory = this_explanatory,
                response = response,
                test_set = FALSE,
                nfolds = nfolds_events,
                grouping_cols = c("harmonized_sample_ontology_term_high_order_fig1"),
                parallel = FALSE,
                save_model = TRUE,
                verbose = FALSE,
                robust_features = TRUE
              )
              if (!is.null(cvfit_rotated$cvfit$`FALSE::FEATURES`)) {
                trans_features <- union(trans_features, cvfit_rotated$cvfit$`FALSE::FEATURES`$full)
              }
              },
              error = function(e) {
                print(paste0('Error in ', id, ' rotation ', rotation_i, ': ', e$message))
              })
              
            }
            
            set(feature_data, j = response, value = orig_psi)
            cvfit <- run_glmnet(
              model_name = file.path(event_dir, id),
              data = feature_data,
              explanatory = this_explanatory,
              response = response,
              test_set = FALSE,
              nfolds = nfolds_events,
              grouping_cols = c("harmonized_sample_ontology_term_high_order_fig1"),
              parallel = FALSE,
              save_model = TRUE,
              verbose = FALSE,
              robust_features = TRUE
            )

            cvfit_filtered <- run_glmnet(
              model_name = file.path(filtered_dir, id),
              data = feature_data,
              explanatory = this_explanatory[!this_explanatory %in% trans_features[trans_features != 'gene_expression']],
              response = response,
              # lambda_list = list('FALSE' = cvfit$cvfit$`FALSE::1`$lambda, 'TRUE' = cvfit$cvfit$`TRUE::1`$lambda),
              test_set = FALSE,
              nfolds = nfolds_events,
              grouping_cols = c("harmonized_sample_ontology_term_high_order_fig1"),
              parallel = FALSE,
              save_model = TRUE,
              verbose = FALSE,
              robust_features = TRUE
            )
            
          })
          
          
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
names(event_models) <- ids_to_build
table(unlist(event_models))

sd_problems <- aggregated_dt[ID %in% names(event_models)[event_models == "not_all_na[response] is not TRUE"], .(sd=sd(get(response))), by=ID]
sd_problems[, table(sd, useNA = "ifany")]
stopifnot(sd_problems[, all(sd == 0, na.rm = TRUE)])

group_problems <- aggregated_dt[ID %in% names(event_models)[startsWith(unlist(event_models), "`k` should be less than")], .(unique_ontology=uniqueN(harmonized_sample_ontology_term_high_order_fig1)), by=ID]
group_problems[, table(unique_ontology)]
stopifnot(group_problems[, all(unique_ontology < nfolds_events)])

constant_y <- aggregated_dt[ID %in% names(event_models)[event_models == "y is constant; gaussian glmnet fails at standardization step"],
                            .(uniqueResponses=unique(get(response))),
                            by = .(ID, harmonized_sample_ontology_term_high_order_fig1)]
constant_y[, .(groupsWithVar=sum(uniqueResponses > 1)), by = ID][, table(groupsWithVar)]
# make sure that per ID in constant_y there are not 5 groups that have uniqueRespnse > 1
stopifnot(all(constant_y[, .(groupsWithVar=sum(uniqueResponses > 1)), by = ID][, groupsWithVar < 5]))

