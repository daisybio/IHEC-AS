## ---------------------------------------------------------------------------
## Event-specific regularised regression via tidymodels
##
## Fits glmnet (elastic net) models for a single splicing event using
## biology-guided grouped CV folds.  One workflow per (feature set × response)
## combination; response = primary PSI or a rotated PSI negative control.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidymodels)
  library(glmnet) # required by parsnip's glmnet engine
})

tidymodels::tidymodels_prefer()
conflicted::conflicts_prefer(base::setdiff)


# ---------------------------------------------------------------------------
# Ontology hierarchy — biology-guided CV supergroup assignment
# ---------------------------------------------------------------------------
# Paths are coarse-to-fine across 4 levels:
#   [super] > [germ_layer/lineage] > [subdomain] > [leaf]
#
# Three super-groups reflect biologically distinct axes:
#   "blood"         — circulating immune / haematopoietic cells
#   "somatic"       — differentiated solid-tissue cells (ecto/endo/mesoderm)
#   "developmental" — uncommitted / de-differentiated states (stem, cancer)
#
# Pairwise distances: inter-super = 1.0; intra-super inter-germ = 0.75;
#                     intra-germ inter-subdomain = 0.5; within subdomain = 0.25.
.ontology_bio_paths <- list(
  # Blood — circulating immune / haematopoietic
  "b lymphocyte" = c("blood", "immune", "lymphoid", "b_lymphocyte"),
  "t lymphocyte" = c("blood", "immune", "lymphoid", "t_lymphocyte"),
  "natural killer cell" = c("blood", "immune", "lymphoid", "natural_killer"),
  "myeloid cell" = c("blood", "immune", "myeloid", "myeloid_cell"),
  "monocyte" = c("blood", "immune", "myeloid", "monocyte"),
  "macrophage" = c("blood", "immune", "myeloid", "macrophage"),
  "neutrophil" = c("blood", "immune", "myeloid", "neutrophil"),
  "eosinophil" = c("blood", "immune", "myeloid", "eosinophil"),
  "dendritic cell" = c("blood", "immune", "myeloid", "dendritic"),
  "mononuclear cell" = c("blood", "immune", "myeloid", "mononuclear"),
  "hematopoietic cell" = c("blood", "immune", "hematopoietic", "hematopoietic"),
  "erythroid lineage cell" = c("blood", "immune", "hematopoietic", "erythroid"),
  "peripheral blood" = c(
    "blood",
    "immune",
    "hematopoietic",
    "peripheral_blood"
  ),
  # Somatic — ectodermal (CNS + neural-crest)
  "brain" = c("somatic", "ectodermal", "nervous_system", "brain"),
  "nervous system" = c(
    "somatic",
    "ectodermal",
    "nervous_system",
    "nervous_system"
  ),
  "neural" = c("somatic", "ectodermal", "nervous_system", "neural"),
  "melanocyte" = c("somatic", "ectodermal", "neural_crest", "melanocyte"),
  # Somatic — endodermal (gut, liver/pancreas, lung)
  "digestive system" = c(
    "somatic",
    "endodermal",
    "digestive",
    "digestive_system"
  ),
  "colon" = c("somatic", "endodermal", "digestive", "colon"),
  "mucosa" = c("somatic", "endodermal", "digestive", "mucosa"),
  "epithelial" = c("somatic", "endodermal", "digestive", "epithelial"),
  "liver" = c("somatic", "endodermal", "hepatopancreatic", "liver"),
  "pancreas" = c("somatic", "endodermal", "hepatopancreatic", "pancreas"),
  "endoderm-derived structure" = c(
    "somatic",
    "endodermal",
    "general",
    "endodermal_structure"
  ),
  "lung" = c("somatic", "endodermal", "respiratory", "lung"),
  # Somatic — mesodermal (connective tissue, muscle, kidney)
  "connective tissue cell" = c(
    "somatic",
    "mesodermal",
    "connective",
    "connective_tissue"
  ),
  "mesoderm-derived structure" = c(
    "somatic",
    "mesodermal",
    "connective",
    "mesodermal_structure"
  ),
  "muscle" = c("somatic", "mesodermal", "muscle", "muscle"),
  "kidney" = c("somatic", "mesodermal", "renal", "kidney"),
  # Developmental — stem / progenitor / embryonic
  "stem cell" = c("developmental", "stem_progenitor", "stem", "stem_cell"),
  "embryonic cell (metazoa)" = c(
    "developmental",
    "stem_progenitor",
    "embryonic",
    "embryonic_cell"
  ),
  "extraembryonic cell" = c(
    "developmental",
    "stem_progenitor",
    "extraembryonic",
    "extraembryonic_cell"
  ),
  "placenta" = c(
    "developmental",
    "stem_progenitor",
    "extraembryonic",
    "placenta"
  ),
  # Developmental — transformed (shares stem_progenitor subdomain so the
  # developmental internal distance = 0.5 < somatic germ-layer distance = 0.75,
  # preventing cancer from splitting off as a lone tiny fold before somatic separates)
  "cancer cell line" = c(
    "developmental",
    "stem_progenitor",
    "transformed",
    "cancer_cell_line"
  )
)


.path_distance <- function(a, b) {
  if (identical(a, b)) {
    return(0.0)
  }
  n <- min(length(a), length(b))
  lcp <- if (n == 0L) 0L else sum(cumprod(a[seq_len(n)] == b[seq_len(n)]) == 1L)
  raw <- (length(a) - lcp) + (length(b) - lcp)
  raw / max(1L, length(a) + length(b))
}


# Names a cluster by the longest common path prefix of its members.
.cluster_name_from_paths <- function(paths, cluster_id) {
  if (length(paths) == 0L) {
    return(paste0("hc_cluster_", cluster_id))
  }
  prefix <- paths[[1L]]
  for (p in paths[-1L]) {
    n <- min(length(prefix), length(p))
    lcp <- if (n == 0L) {
      0L
    } else {
      sum(cumprod(prefix[seq_len(n)] == p[seq_len(n)]) == 1L)
    }
    prefix <- prefix[seq_len(lcp)]
    if (length(prefix) == 0L) break
  }
  if (length(prefix) > 0L) {
    paste0(
      "hc_",
      paste(prefix[seq_len(min(2L, length(prefix)))], collapse = "_")
    )
  } else {
    paste0("hc_cluster_", cluster_id)
  }
}


# Maps a character vector of ontology labels to biology-guided supergroup names.
# Uses hierarchical clustering (average linkage) on path distances to form
# exactly n_groups groups.  Unknown labels are assigned to "other".
.map_ontology_to_supergroups <- function(ontology_vec, n_groups = 5L) {
  labels <- sort(names(.ontology_bio_paths))
  n <- length(labels)
  stopifnot(n_groups >= 3L, n_groups <= n)

  dist_mat <- matrix(0.0, n, n, dimnames = list(labels, labels))
  for (i in seq_len(n - 1L)) {
    for (j in seq.int(i + 1L, n)) {
      d <- .path_distance(
        .ontology_bio_paths[[labels[i]]],
        .ontology_bio_paths[[labels[j]]]
      )
      dist_mat[i, j] <- d
      dist_mat[j, i] <- d
    }
  }

  cluster_ids <- cutree(
    hclust(as.dist(dist_mat), method = "average"),
    k = n_groups
  )

  # Name each cluster by its longest common path prefix
  members <- split(labels, cluster_ids)
  used_names <- character(0L)
  id_to_name <- vapply(
    names(members),
    function(cid) {
      paths <- lapply(.ontology_bio_paths[members[[cid]]], identity)
      base_name <- .cluster_name_from_paths(paths, cluster_id = cid)
      if (base_name %in% used_names) {
        base_name <- paste0(base_name, "_", sum(used_names == base_name) + 1L)
      }
      used_names <<- c(used_names, base_name)
      base_name
    },
    character(1L)
  )

  norm_ont <- tolower(trimws(as.character(ontology_vec)))
  mapped <- cluster_ids[norm_ont]
  ifelse(is.na(mapped), "other", id_to_name[as.character(mapped)])
}


# ---------------------------------------------------------------------------
# Metrics  (module-level: defined once, reused per call)
# ---------------------------------------------------------------------------
# All responses are PSI-scale ([0, 1]); standard yardstick metrics suffice.
# Model selection uses rsq_trad (traditional R²: 1 − SS_res/SS_tot).
all_metrics <- yardstick::metric_set(
  yardstick::rmse,
  yardstick::rsq_trad,
  yardstick::ccc
)


# ---------------------------------------------------------------------------
# Main function
# ---------------------------------------------------------------------------

#' Fit glmnet for a single splicing event
#'
#' One workflow per (feature set × response): primary PSI + rotated PSI
#' negative controls.  Tuned via biology-guided grouped CV; returns robust
#' non-zero coefficients per workflow.
#'
#' @param this_feature_data  data.table: predictors, response, grouping col.
#' @param explanatory_vars   Named list of predictor-name vectors (one per
#'   feature set, e.g. \code{list(long = ..., short = ..., local = ...)}).
#' @param response      Primary PSI response column name.
#' @param rotated_psis  Matrix/data.table of rotated-PSI negative-control
#'   columns appended to \code{this_feature_data}.
#' @param grouping_col  Ontology column for biology-guided grouped CV.
#' @param nfolds        Number of CV folds (ontology supergroups).
#' @param seed          Base random seed; per-workflow offset added inside loop.
#' @param parallel      \code{mc.cores} for \code{pbmclapply}.
#' @param cv_folds_path Path for the fold-assignment summary CSV (gzipped).
#'   If \code{NULL}, nothing is written.
#'
#' @return Named list (one element per workflow) where each element contains:
#'   \code{wflow_id}, \code{best_params}, \code{all_metrics}, \code{fit_metrics},
#'   \code{model_type}, \code{nonzero_coefs}, \code{fold_features},
#'   \code{robust_features}.
run_event_glmnet <- function(
  this_feature_data,
  explanatory_vars,
  response,
  rotated_psis = NULL,
  grouping_col = "ontology",
  nfolds = 5L,
  seed = 1234L,
  parallel = 1L,
  cv_folds_path = NULL
) {
  # -- Ontology-guided grouped CV folds ---------------------------------------
  set.seed(seed)
  n_unique_ont <- length(unique(tolower(trimws(
    as.character(this_feature_data[[grouping_col]])
  ))))
  actual_nfolds <- max(3L, min(nfolds, n_unique_ont))
  if (actual_nfolds != nfolds) {
    warning(sprintf(
      "CV folds reduced %d → %d (only %d unique ontology groups in data)",
      nfolds, actual_nfolds, n_unique_ont
    ))
  }
  this_feature_data[,
    hc_group := .map_ontology_to_supergroups(
      this_feature_data[[grouping_col]],
      n_groups = actual_nfolds
    )
  ]
  n_unique_hc <- data.table::uniqueN(this_feature_data[["hc_group"]])
  if (n_unique_hc < actual_nfolds) {
    warning(sprintf(
      "CV folds further reduced %d → %d (mapping collapsed to %d supergroups)",
      actual_nfolds, n_unique_hc, n_unique_hc
    ))
    actual_nfolds <- n_unique_hc
  }
  cv_folds <- rsample::group_vfold_cv(
    this_feature_data,
    group = "hc_group",
    v = actual_nfolds,
    balance = "observations"
  )

  if (!is.null(cv_folds_path)) {
    cv_folds_dt <- data.table::rbindlist(mapply(
      cv_folds$splits,
      cv_folds$id,
      FUN = function(split, id) {
        split |>
          rsample::assessment() |>
          dplyr::select(dplyr::all_of(grouping_col), hc_group) |>
          dplyr::mutate(resample = id) |>
          dplyr::count(dplyr::across(dplyr::everything()))
      },
      SIMPLIFY = FALSE
    ))
    data.table::fwrite(cv_folds_dt, cv_folds_path)
  }

  # -- Recipes ----------------------------------------------------------------
  # One base recipe per response (primary PSI + rotated PSI columns), then one
  # variant per explanatory set.
  responses_to_build <- c(response, colnames(rotated_psis))

  # Base: all columns are "misc" (excluded) except the nominated outcome.
  base_recipes <- sapply(
    responses_to_build,
    function(this_response) {
      recipe(
        this_feature_data,
        roles = rep("misc", length(this_feature_data))
      ) |>
        update_role(all_of(this_response), new_role = "outcome")
    },
    simplify = FALSE
  )

  # For each (response, explanatory set) pair: promote predictors, impute,
  # drop zero-variance columns, dummy-encode nominals, and add
  # protocol × gene_expression_vst interaction terms. (§4.12b: 05 no longer
  # emits a bare `gene_expression` column.)
  explanatory_recipe_list <- unlist(
    sapply(
      explanatory_vars,
      function(var) {
        sapply(
          base_recipes,
          function(base_recipe) {
            base_recipe |>
              update_role(all_of(var), new_role = "predictor") |>
              step_rm(has_role("misc")) |>
              step_impute_mean(all_numeric_predictors()) |>
              step_dummy(all_nominal_predictors(), one_hot = TRUE) |>
              step_interact(
                terms = ~ starts_with("protocol"):gene_expression_vst
              ) |>
              step_zv(all_predictors())
          },
          simplify = FALSE
        )
      },
      simplify = FALSE
    ),
    recursive = FALSE
  )

  # unlist() joins nested list names with "." — normalise to "_" so that
  # workflow_set produces predictable IDs like "long_PSI_glmnet".
  names(explanatory_recipe_list) <- gsub(
    ".",
    "_",
    names(explanatory_recipe_list),
    fixed = TRUE
  )

  # -- Model specifications ---------------------------------------------------
  # glmnet: mixture and penalty both tuned.  The submodel trick lets tune() fit
  # ONE path per (mixture, fold) and evaluate all 50 lambda values via
  # multi_predict() — ~6 fits per fold yield 300 evaluations.
  glmnet_spec <- linear_reg(penalty = tune(), mixture = tune()) |>
    set_engine("glmnet")

  # -- Workflow set & grids ---------------------------------------------------
  wflow_set <- workflow_set(
    preproc = explanatory_recipe_list,
    models = list(glmnet = glmnet_spec)
  )

  # glmnet: structured crossing grid — 6 alpha × 50 lambda values.
  # Alpha values skew toward lasso (mixture ≈ 1) for sparsity and feature
  # identification.  Intermediate values (0.5–0.9) prevent arbitrary selection
  # among correlated chromHMM features where pure lasso is unstable.
  glmnet_grid <- tidyr::crossing(
    penalty = 10^seq(-5, 0, length.out = 50),
    mixture = c(0.5, 0.7, 0.9)
  ) |>
    arrange(mixture, penalty) # ensure lambda path is contiguous for each alpha

  glmnet_params <- extract_parameter_set_dials(glmnet_spec)

  for (wid in wflow_set$wflow_id) {
    wflow_set <- wflow_set |>
      option_add(param_info = glmnet_params, id = wid) |>
      option_add(grid = glmnet_grid, id = wid)
  }

  wflow_to_tune <- wflow_set
  rm(wflow_set, explanatory_recipe_list, base_recipes)

  n_wflows <- nrow(wflow_to_tune)

  # -- Feature extraction helpers ---------------------------------------------

  # Non-zero glmnet coefficients at the selected penalty, excluding the
  # intercept and gene_expression_vst (used as an interaction base term).
  .extract_glmnet_coefs <- function(best_params, final_fit) {
    penalty_val <- best_params$penalty[[1L]]
    engine <- workflows::extract_fit_engine(final_fit)
    coef_mat <- as.matrix(coef(engine, s = penalty_val))
    mask <- coef_mat[, 1L] != 0 &
      !rownames(coef_mat) %in% c("(Intercept)", "gene_expression_vst")
    coef_mat[mask, 1L, drop = FALSE]
  }

  # Robustness filter: a feature must be non-zero in the full-data fit AND
  # appear in >= 2 CV folds (uniqueN(fold) >= 3 counts the "full" label too).
  # Replicates the filter used in 07-ml-helper.R.
  .robust_glmnet_coefs <- function(
    best_params,
    final_fit,
    final_wf,
    cv_folds,
    this_feature_data
  ) {
    full_coefs <- .extract_glmnet_coefs(best_params, final_fit)
    full_features <- rownames(full_coefs)

    if (length(full_features) == 0L) {
      return(list(
        nonzero_coefs = full_coefs,
        fold_features = list(),
        robust_features = character(0L)
      ))
    }

    penalty_val <- best_params$penalty[[1L]]
    fold_features <- lapply(cv_folds$splits, function(split) {
      fold_fit <- parsnip::fit(final_wf, data = rsample::analysis(split))
      fold_eng <- workflows::extract_fit_engine(fold_fit)
      fold_cm <- as.matrix(coef(fold_eng, s = penalty_val))
      rownames(fold_cm)[
        fold_cm[, 1L] != 0 &
          !rownames(fold_cm) %in% c("(Intercept)")
      ]
    })
    # Name each fold element so rbindlist's idcol produces distinct IDs.
    # lapply() returns an unnamed list; without explicit names the idcol column
    # gets empty strings for all numeric folds, making uniqueN(fold) = 2
    # regardless of how many folds a feature appears in.
    names(fold_features) <- seq_along(fold_features)
    fold_features[["0"]] <- full_features

    feat_dt <- data.table::rbindlist(
      lapply(fold_features, function(x) list(feature = x)),
      idcol = "fold"
    )

    robust_dt <- feat_dt[,
      if ("0" %in% fold && data.table::uniqueN(fold) >= 3L) {
        list(n_folds = data.table::uniqueN(fold))
      },
      by = "feature"
    ]

    robust_feats <- full_features[full_features %in% robust_dt$feature]
    list(
      nonzero_coefs = full_coefs[robust_feats, , drop = FALSE],
      fold_features = fold_features,
      robust_features = robust_feats
    )
  }

  # -- Tune, finalize, fit, and extract features in a single pass -------------
  workflow_results <- pbmcapply::pbmclapply(
    seq_len(n_wflows),
    function(i) {
      wid <- wflow_to_tune$wflow_id[[i]]
      tryCatch(
        {
          print(glue::glue(
            "Tuning workflow {i} / {n_wflows} ({wid})..."
          ))
          set.seed(seed + i) # per-workflow offset for reproducibility
          wf <- workflowsets::extract_workflow(wflow_to_tune, id = wid)
          this_grid <- wflow_to_tune$option[[i]]$grid

          tuned <- tune::tune_grid(
            wf,
            resamples = cv_folds,
            grid = this_grid,
            metrics = all_metrics,
            control = tune::control_grid(
              verbose = FALSE,
              allow_par = FALSE # prevent nested forking inside outer pbmclapply
            )
          )

          all_metrics_dt <- data.table::as.data.table(
            tuned |> collect_metrics(summarize = FALSE)
          )

          # Finalize, fit, and extract features for one lambda selection criterion.
          .process_selection <- function(best_params) {
            if (nrow(best_params) == 0L) {
              return(NULL)
            }
            best_cv_summary <- tuned |>
              collect_metrics(summarize = TRUE) |>
              dplyr::filter(.config == best_params$.config) |>
              dplyr::select(.metric, mean, std_err, n) |>
              tidyr::pivot_wider(
                names_from = .metric,
                values_from = c(mean, std_err, n),
                names_glue = "cv_{.value}_{.metric}"
              )
            best_params <- dplyr::bind_cols(best_params, best_cv_summary)
            final_wf <- tune::finalize_workflow(wf, best_params)
            final_fit <- parsnip::fit(final_wf, data = this_feature_data)

            wf_response <- workflows::extract_preprocessor(final_wf)$var_info |>
              dplyr::filter(role == "outcome") |>
              dplyr::pull(variable)

            fit_preds <- predict(final_fit, new_data = this_feature_data)$.pred
            fit_metrics <- all_metrics(
              data = tibble(
                truth = this_feature_data[[wf_response]],
                estimate = fit_preds
              ),
              truth = truth,
              estimate = estimate
            )

            feat <- .robust_glmnet_coefs(
              best_params,
              final_fit,
              final_wf,
              cv_folds,
              this_feature_data
            )
            rm(final_fit, final_wf)
            c(
              list(
                best_params = best_params,
                fit_metrics = fit_metrics,
                model_type = "glmnet"
              ),
              feat
            )
          }

          sel_min <- .process_selection(select_best(tuned, metric = "rsq_trad"))
          sel_1se <- .process_selection(select_by_one_std_err(
            tuned,
            metric = "rsq_trad",
            desc(penalty)
          ))
          rm(tuned)

          list(
            wflow_id = wid,
            all_metrics = all_metrics_dt,
            selections = list(min = sel_min, `1se` = sel_1se)
          )
        },
        error = function(e) {
          call_str <- if (!is.null(conditionCall(e))) {
            paste0(" [", deparse(conditionCall(e))[[1L]], "]")
          } else {
            ""
          }
          message(sprintf(
            "Workflow %s error: %s%s",
            wid,
            conditionMessage(e),
            call_str
          ))
          structure(conditionMessage(e), class = "try-error")
        }
      )
    },
    mc.cores = parallel,
    ignore.interactive = TRUE
  )

  names(workflow_results) <- wflow_to_tune$wflow_id
  workflow_results
}


# ---------------------------------------------------------------------------
# CLI entry point — invoked by 09-1-ml-local-array.sh as:
#   Rscript 09zz-ml-event-glmnet-tidymodels.R <cfg_rds> <event_id>
# ---------------------------------------------------------------------------
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  cfg_path <- args[1]
  this_id <- as.integer(args[2])

  ## dev:
  # cfg_path <- "processed_data/event_glmnet_cfg.rds"
  # this_id <- 29728

  cfg <- readRDS(cfg_path)
  setwd(cfg$project_dir)

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
        paste0("feature_table_", this_id, ".csv.gz")
      ))

      explanatory <- names(feature_data)[
        !names(feature_data) %in%
          c(
            "IHEC",
            "ID",
            "Event Type",
            "Variability",
            "seqnames",
            "gene_id",
            "uuid",
            "transcript_filter",
            "project",
            "harmonized_sample_ontology_term_high_order_fig1",
            # §4.12b: 05 emits BOTH gene_expression_getmm + gene_expression_vst;
            # event-specific models use vst (single-gene, cross-sample —
            # matches 06/09-1 routing), getmm is the pooled splicing_ml copy —
            # exclude it here so it isn't also fit as a second expression
            # predictor alongside vst.
            "gene_expression_getmm",
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
        chromhmm_hits_smaller[, "queryHits"] ==
          which(keep_rows_manual == this_id),
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
      this_event <- event_dt[ID == this_id]
      other_ids <- other_ids[
        other_ids != this_id &
          other_ids %in% event_dt[`Event Type` == this_event$`Event Type`, ID] &
          other_ids %in% event_dt[seqnames != this_event$seqnames, ID] &
          other_ids %in% event_dt[Variability == this_event$Variability, ID] &
          other_ids %in%
            event_dt[transcript_filter == this_event$transcript_filter, ID]
      ]

      this_rotations <- min(nrotations, length(other_ids))
      if (this_rotations < nrotations) {
        warning(sprintf(
          "Not enough rotation controls for %d, using %d",
          this_id,
          this_rotations
        ))
      }
      stopifnot(
        rownames(subset_psi_matrix) == feature_data[, as.character(uuid)]
      )
      set.seed(this_id)
      print(glue::glue(
        "Running event {this_id} with {this_rotations} rotations (out of {nrotations})..."
      ))
      rotated_psis <- subset_psi_matrix[, as.character(sample(
        other_ids,
        this_rotations
      ))]

      rm(
        psi_table,
        subset_psi_matrix,
        event_dt,
        chromhmm_hits_smaller,
        keep_rows_manual
      )
      gc()

      parallel <- if (length(args) >= 3L) as.integer(args[3L]) else 1L

      wrs <- run_event_glmnet(
        this_feature_data = cbind(feature_data, rotated_psis),
        explanatory_vars,
        response,
        rotated_psis,
        grouping_col,
        nfolds,
        seed = this_id,
        parallel = parallel,
        cv_folds_path = file.path(
          event_dir,
          paste0(this_id, "_cv_folds.csv.gz")
        )
      )
      rm(feature_data, rotated_psis, explanatory_vars)
      gc()

      # Log any pbmclapply worker errors before writing CSVs.
      worker_errors <- Filter(
        inherits_try_error <- function(x) inherits(x, "try-error"),
        wrs
      )
      if (length(worker_errors) > 0L) {
        message(sprintf(
          "WARNING: %d/%d workflows errored for event %d: %s",
          length(worker_errors),
          length(wrs),
          this_id,
          paste(names(worker_errors), collapse = ", ")
        ))
      }

      # Write one file per table; rows tagged with lambda_selection ("min"/"1se").
      # Each per-workflow-selection fn call is guarded: a failure logs to stderr
      # and returns NULL rather than aborting the whole write pass.
      write_wrs_sel <- function(name, fn) {
        dt <- data.table::rbindlist(
          lapply(names(wrs), function(wid) {
            wr <- wrs[[wid]]
            if (inherits(wr, "try-error")) {
              return(NULL)
            }
            data.table::rbindlist(
              lapply(names(wr$selections), function(sel) {
                wr_sel <- wr$selections[[sel]]
                if (is.null(wr_sel)) {
                  return(NULL)
                }
                tryCatch(
                  {
                    res <- fn(wid, sel, wr_sel)
                    if (!is.null(res)) {
                      res[, lambda_selection := sel][]
                    } else {
                      NULL
                    }
                  },
                  error = function(e) {
                    message(sprintf(
                      "write_wrs_sel[%s/%s/%s] error: %s",
                      name,
                      wid,
                      sel,
                      conditionMessage(e)
                    ))
                    NULL
                  }
                )
              }),
              fill = TRUE
            )
          }),
          fill = TRUE
        )
        if (is.null(dt) || nrow(dt) == 0L) {
          dt <- data.table::data.table()
        }
        dt[, ID := this_id]
        data.table::fwrite(
          dt,
          file.path(event_dir, paste0(this_id, "_", name, ".csv.gz"))
        )
      }

      # all_metrics: shared across selections — write once (completion sentinel)
      write_wrs_once <- function(name, fn) {
        dt <- data.table::rbindlist(
          lapply(names(wrs), function(wid) {
            tryCatch(fn(wid), error = function(e) {
              message(sprintf(
                "write_wrs_once[%s/%s] error: %s",
                name,
                wid,
                conditionMessage(e)
              ))
              NULL
            })
          }),
          fill = TRUE
        )
        if (is.null(dt) || nrow(dt) == 0L) {
          dt <- data.table::data.table()
        }
        dt[, ID := this_id]
        data.table::fwrite(
          dt,
          file.path(event_dir, paste0(this_id, "_", name, ".csv.gz"))
        )
      }

      write_wrs_sel("best_params", function(wid, sel, wr_sel) {
        bp <- wr_sel$best_params
        if (is.null(bp)) {
          return(NULL)
        }
        data.table::as.data.table(bp)[, wflow_id := wid][]
      })

      write_wrs_sel("fit_metrics", function(wid, sel, wr_sel) {
        m <- wr_sel$fit_metrics
        if (is.null(m)) {
          return(NULL)
        }
        data.table::as.data.table(m)[, wflow_id := wid][]
      })

      write_wrs_sel("event_summary", function(wid, sel, wr_sel) {
        rf <- wr_sel$robust_features
        data.table::data.table(
          wflow_id = wid,
          model_type = wr_sel$model_type,
          n_robust = length(rf),
          has_model = length(rf) > 0L,
          expression_bias = data.table::fifelse(
            length(rf) == 0L,
            NA_character_,
            data.table::fifelse(
              any(grepl("gene_expression", rf, fixed = TRUE)),
              "Expression Bias",
              "Epigenetic Only"
            )
          )
        )
      })

      write_wrs_sel("fold_features", function(wid, sel, wr_sel) {
        ff <- wr_sel$fold_features
        if (is.null(ff) || length(ff) == 0L) {
          return(NULL)
        }
        data.table::rbindlist(
          lapply(names(ff), function(fold) {
            data.table::data.table(
              wflow_id = wid,
              fold = fold,
              feature = ff[[fold]]
            )
          }),
          fill = TRUE
        )
      })

      write_wrs_sel("nonzero_coefs", function(wid, sel, wr_sel) {
        nc <- wr_sel$nonzero_coefs
        if (is.null(nc) || nrow(nc) == 0L) {
          return(NULL)
        }
        data.table::data.table(
          wflow_id = wid,
          feature = rownames(nc),
          coef = nc[, 1L]
        )
      })

      # all_metrics written last — its existence is the completion sentinel
      write_wrs_once("all_metrics", function(wid) {
        wr <- wrs[[wid]]
        if (inherits(wr, "try-error")) {
          return(NULL)
        }
        m <- wr$all_metrics
        if (is.null(m)) {
          return(NULL)
        }
        data.table::as.data.table(m)[, wflow_id := wid][]
      })

      rm(wrs)
      gc()
      message("Done: ", this_id)
    },
    error = function(e) {
      call_str <- if (!is.null(conditionCall(e))) {
        paste0(" [", deparse(conditionCall(e))[[1L]], "]")
      } else {
        ""
      }
      message("ERROR for id ", this_id, ": ", conditionMessage(e), call_str)
    }
  )
}
