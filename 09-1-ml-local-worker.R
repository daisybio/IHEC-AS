#!/usr/bin/env Rscript
# Called by 09-1-ml-local-array.sh: <cfg_rds> <event_id>
args <- commandArgs(trailingOnly = TRUE)
cfg_path <- args[1]
id <- as.integer(args[2])

## dev:
# cfg_path <- "processed_data/event_glmnet_cfg.rds"
# id <- 29728

cfg <- readRDS(cfg_path)
setwd(cfg$project_dir)
source(file.path(cfg$project_dir, "07-ml-event-glmnet-tidymodels.R"))

sess <- readRDS(cfg$session_rds)
psi_table <- sess$psi_table
event_dt <- sess$event_dt
chromhmm_hits_smaller <- sess$chromhmm_hits_smaller
keep_rows_manual <- sess$keep_rows_manual
rm(sess)
gc()

feature_table_dir <- cfg$feature_table_dir
feature_sets <- cfg$feature_sets
event_dir <- cfg$event_dir
response <- cfg$response
grouping_col <- cfg$grouping_col
nfolds <- cfg$nfolds
nrotations <- cfg$nrotations

tryCatch(
  {
    feature_data <- data.table::fread(file.path(
      feature_table_dir,
      paste0("feature_table_", id, ".csv.gz")
    ))

    explanatory <- names(feature_data)[
      !names(feature_data) %in%
        c(
          "IHEC",
          "ID",
          "Event Type",
          "Variability",
          "gene_id",
          "uuid",
          "transcript_filter",
          "project",
          "harmonized_sample_ontology_term_high_order_fig1",
          grouping_col,
          response
        )
    ]
    chromhmm_explanatory <- explanatory[grepl(
      "chromhmm",
      explanatory,
      fixed = TRUE
    )]

    smaller_chromhmm_ids <- chromhmm_hits_smaller[
      chromhmm_hits_smaller[, "queryHits"] == which(keep_rows_manual == id),
      "subjectHits"
    ]
    old_chromhmm_explanatory <- chromhmm_explanatory[Reduce(
      `&`,
      lapply(sprintf("chromhmm_%d", smaller_chromhmm_ids), function(suff) {
        !endsWith(chromhmm_explanatory, suff)
      })
    )]

    explanatory_vars <- setNames(
      list(
        explanatory,
        explanatory[!explanatory %in% old_chromhmm_explanatory],
        explanatory[!explanatory %in% chromhmm_explanatory]
      ),
      feature_sets
    )

    subset_psi_matrix <- psi_table[feature_data[, uuid], , drop = FALSE]
    other_ids <- as.integer(colnames(subset_psi_matrix)[
      colSums(is.na(subset_psi_matrix)) == 0 &
        apply(subset_psi_matrix, 2, sd, na.rm = TRUE) > 0
    ])
    this_event <- event_dt[ID == id]
    other_ids <- other_ids[
      other_ids != id &
        other_ids %in%
          event_dt[`Event Type` == this_event$`Event Type`, ID] &
        other_ids %in% event_dt[seqnames != this_event$seqnames, ID] &
        other_ids %in% event_dt[Variability == this_event$Variability, ID] &
        other_ids %in%
          event_dt[transcript_filter == this_event$transcript_filter, ID]
    ]

    this_rotations <- min(nrotations, length(other_ids))
    if (this_rotations < nrotations) {
      warning(sprintf(
        "Not enough rotation controls for %d, using %d",
        id,
        this_rotations
      ))
    }
    stopifnot(
      rownames(subset_psi_matrix) == feature_data[, as.character(uuid)]
    )
    set.seed(id)
    print(glue::glue(
      "Running event {id} with {this_rotations} rotations (out of {nrotations})..."
    ))
    rotated_psis <- subset_psi_matrix[, as.character(sample(
      other_ids,
      this_rotations
    ))]

    rm(psi_table, subset_psi_matrix, event_dt, chromhmm_hits_smaller, keep_rows_manual)
    gc()

    parallel <- 1L

    ## dev:
    # this_feature_data = cbind(feature_data, rotated_psis)
    # seed = id
    # parallel = 40L

    wflow_res <- run_event_glmnet(
      this_feature_data = cbind(feature_data, rotated_psis),
      explanatory_vars,
      response,
      rotated_psis,
      grouping_col,
      nfolds,
      seed = id,
      parallel = parallel
    )

    # Save all per-event result tables as CSVs; no RDS written
    wrs <- wflow_res$workflow_results
    bind_wrs <- function(fn) data.table::rbindlist(lapply(names(wrs), fn), fill = TRUE)

    # 1. all_metrics — per-fold CV metrics across hyperparameter grid
    am_dt <- bind_wrs(function(wid) {
      m <- wrs[[wid]]$all_metrics
      if (is.null(m)) return(NULL)
      cbind(wflow_id = wid, data.table::as.data.table(m))
    })
    am_dt[, ID := id]
    data.table::fwrite(am_dt, file.path(event_dir, paste0("all_metrics_", id, ".csv.gz")))

    # 2. best_params — selected hyperparameters + CV summary metrics
    bp_dt <- bind_wrs(function(wid) {
      cbind(wflow_id = wid, data.table::as.data.table(wrs[[wid]]$best_params))
    })
    bp_dt[, ID := id]
    data.table::fwrite(bp_dt, file.path(event_dir, paste0("best_params_", id, ".csv.gz")))

    # 3. fit_metrics — in-sample metrics on full training data
    fm_dt <- bind_wrs(function(wid) {
      m <- wrs[[wid]]$fit_metrics
      if (is.null(m)) return(NULL)
      cbind(wflow_id = wid, data.table::as.data.table(m))
    })
    fm_dt[, ID := id]
    data.table::fwrite(fm_dt, file.path(event_dir, paste0("fit_metrics_", id, ".csv.gz")))

    # 4. event_summary — model type + robust feature count per workflow
    es_dt <- bind_wrs(function(wid) {
      wr <- wrs[[wid]]
      rf <- wr$robust_features
      data.table::data.table(
        wflow_id = wid,
        model_type = wr$model_type,
        n_robust = length(rf),
        has_model = length(rf) > 0L,
        expression_bias = data.table::fifelse(
          "gene_expression" %in% rf, "Expression Bias", "Epigenetic Only"
        )
      )
    })
    es_dt[, ID := id]
    data.table::fwrite(es_dt, file.path(event_dir, paste0("event_summary_", id, ".csv.gz")))

    # 5. fold_features — selected features per CV fold per workflow
    ff_dt <- bind_wrs(function(wid) {
      ff <- wrs[[wid]]$fold_features
      if (is.null(ff) || length(ff) == 0L) return(NULL)
      cbind(
        wflow_id = wid,
        data.table::rbindlist(
          lapply(names(ff), function(fold) data.table::data.table(fold = fold, feature = ff[[fold]])),
          fill = TRUE
        )
      )
    })
    ff_dt[, ID := id]
    data.table::fwrite(ff_dt, file.path(event_dir, paste0("fold_features_", id, ".csv.gz")))

    # 6. nonzero_coefs — robust non-zero coefficients per workflow
    nc_dt <- bind_wrs(function(wid) {
      nc <- wrs[[wid]]$nonzero_coefs
      if (is.null(nc) || nrow(nc) == 0L) return(NULL)
      data.table::data.table(wflow_id = wid, feature = rownames(nc), coef = nc[, 1L])
    })
    nc_dt[, ID := id]
    data.table::fwrite(nc_dt, file.path(event_dir, paste0("nonzero_coefs_", id, ".csv.gz")))

    message("Done: ", id)
  },
  error = function(e) {
    message("ERROR for id ", id, ": ", e$message)
  }
)
