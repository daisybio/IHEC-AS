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
# base::intersect over dplyr::/GenomicRanges:: (09-1 RBP-col join + control matching
# use plain vector intersect). Set here since 09-1 source()s this file before Phase 1.
conflicted::conflicts_prefer(base::intersect)


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
# Uses hierarchical clustering (average linkage) on path distances to form up
# to n_groups groups. Unknown labels are assigned to "other".
#
# Clusters only the labels PRESENT in `ontology_vec` (not the full fixed
# vocabulary in .ontology_bio_paths) -- fixed 2026-07-22. Clustering the full
# 34-label universe every time is data-independent: this dataset's ontology
# labels are majority blood/immune, and at the configured distances all 13
# blood/immune leaves sit close enough together (subdomain gaps of 0.5) that
# any k in [2,33] merges them into ONE cluster before the 21 non-blood labels
# (whose germ-layer gaps are wider, 0.75) finish splitting -- so a blood-only
# event always collapsed to 1 supergroup regardless of n_groups, which then
# crashed rsample::group_vfold_cv(v=1) ("v must be >= 2"). Clustering only the
# present labels lets blood's own lymphoid/myeloid/hematopoietic subdomains
# split against each other instead of being swamped by irrelevant somatic/
# developmental leaves that this event's samples don't even contain.
.map_ontology_to_supergroups <- function(ontology_vec, n_groups = 5L) {
  norm_ont <- tolower(trimws(as.character(ontology_vec)))
  labels <- sort(unique(norm_ont[norm_ont %in% names(.ontology_bio_paths)]))
  n <- length(labels)

  if (n == 0L) {
    return(rep("other", length(norm_ont)))
  }

  effective_groups <- max(1L, min(n_groups, n))

  if (n == 1L || effective_groups == 1L) {
    # Nothing to split: either only one known label is present, or n_groups
    # collapsed to 1 -- every known label falls in a single named group.
    id_to_name <- setNames(
      .cluster_name_from_paths(.ontology_bio_paths[labels], cluster_id = 1L),
      "1"
    )
    mapped <- ifelse(norm_ont %in% labels, "1", NA_character_)
    return(ifelse(is.na(mapped), "other", id_to_name[mapped]))
  }

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
    k = effective_groups
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

  # -- Nominal-predictor level fixing -----------------------------------------
  # Cast every nominal (character) predictor to a factor whose levels come from
  # the FULL event data, ONCE, before any CV split. R factor subsetting keeps the
  # parent's full `levels`, so a per-fold training split that happens to contain
  # only one observed value (e.g. all `observed` for a mark's `H3K*_source` flag,
  # or a single `protocol`) still carries every level. Without this, `step_dummy()`
  # aborts at bake with "Only one factor level in <col>" on that fold — which,
  # under the outer pbmclapply's silent try(), sank every workflow for such an
  # event with no error in the logs (found via the 2026-07-22 event_models smoke).
  # Mirrors splicing_ml's "OHE categories from the full dataset, not per-fold"
  # rule. NA stays NA (→ all-zero dummy, unchanged behaviour). Restricted to
  # actual predictors so id/meta character cols (uuid, gene_id, …) — which are
  # role="misc" and step_rm'd before step_dummy anyway — aren't needlessly cast.
  predictor_cols <- base::intersect(
    unique(unlist(explanatory_vars, use.names = FALSE)),
    names(this_feature_data)
  )
  nominal_predictors <- predictor_cols[vapply(
    predictor_cols,
    function(cc) is.character(this_feature_data[[cc]]),
    logical(1)
  )]
  if (length(nominal_predictors) > 0L) {
    this_feature_data[,
      (nominal_predictors) := lapply(.SD, factor),
      .SDcols = nominal_predictors
    ]
  }

  # Drop GLOBALLY-constant nominal predictors (single level across the whole
  # event) — the cast above only rescues per-fold-single-but-globally-multi
  # columns; a column that is constant across every sample (e.g. an event whose
  # `H3K4me3_source` is `observed` for all samples) has just one level, still
  # trips `step_dummy`, and carries zero information anyway. Remove it from the
  # explanatory sets so it is never promoted to predictor (role stays "misc" →
  # step_rm'd before step_dummy). Note whether `protocol` survives: if protocol
  # is itself constant it gets dropped here, and its interaction term below must
  # be omitted (a constant-protocol × gene_expression_vst interaction is just a
  # rescaled gene_expression_vst, which is already a main-effect predictor —
  # no information lost).
  constant_nominals <- nominal_predictors[vapply(
    nominal_predictors,
    function(cc) nlevels(this_feature_data[[cc]]) < 2L,
    logical(1)
  )]
  if (length(constant_nominals) > 0L) {
    message(sprintf(
      "Event %s: dropping %d globally-constant nominal predictor(s): %s",
      this_feature_data[["ID"]][[1L]],
      length(constant_nominals),
      paste(constant_nominals, collapse = ", ")
    ))
    explanatory_vars <- lapply(
      explanatory_vars,
      function(v) base::setdiff(v, constant_nominals)
    )
  }
  protocol_varies <- "protocol" %in%
    unique(unlist(explanatory_vars, use.names = FALSE))

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
  # dummy-encode nominals, add the protocol × gene_expression_vst interaction
  # (only when protocol actually varies — see protocol_varies above), and drop
  # zero-variance columns. (§4.12b: 05 no longer emits a bare `gene_expression`
  # column.)
  explanatory_recipe_list <- unlist(
    sapply(
      explanatory_vars,
      function(var) {
        sapply(
          base_recipes,
          function(base_recipe) {
            rec <- base_recipe |>
              update_role(all_of(var), new_role = "predictor") |>
              step_rm(has_role("misc")) |>
              step_impute_mean(all_numeric_predictors()) |>
              step_dummy(all_nominal_predictors(), one_hot = TRUE)
            # Skip the interaction entirely for constant-protocol events — the
            # protocol column was dropped above, so starts_with("protocol")
            # would select nothing, and the interaction would be redundant with
            # the gene_expression_vst main effect regardless.
            if (protocol_varies) {
              rec <- rec |>
                step_interact(
                  terms = ~ starts_with("protocol"):gene_expression_vst
                )
            }
            rec |> step_zv(all_predictors())
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
# CLI entry point — invoked (per event) as:
#   Rscript 09zz-ml-event-glmnet-tidymodels.R <cfg_rds> <event_id> [cores]
# ---------------------------------------------------------------------------
# Guard on THIS file being the script Rscript actually invoked — NOT merely
# `!interactive()`. 09-1-ml-local.R source()s this file for run_event_glmnet(),
# and source() also runs under Rscript (!interactive() == TRUE there too); a bare
# !interactive() block would then execute with 09-1's own (cfg-less) args →
# `readRDS(NA)` crash during the build phase. `--file=` names the invoked script.
.is_09zz_script <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  length(f) == 1L && basename(f) == "09zz-ml-event-glmnet-tidymodels.R"
})
if (!interactive() && .is_09zz_script) {
  args <- commandArgs(trailingOnly = TRUE)
  cfg_path <- args[1]
  this_id <- as.integer(args[2])

  ## dev:
  # cfg_path <- "processed_data/event_glmnet_cfg.rds"
  # this_id <- 29728

  cfg <- readRDS(cfg_path)
  setwd(cfg$project_dir)

  sess <- readRDS(cfg$session_rds)
  # psi_table / event_dt were only used for rotation-control matching, now owned
  # by the Tier-1 screen — Tier-2 needs only the chromHMM-vicinity objects.
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
  # §4.9: base seed threaded from .slurm_cfg (fallback to the global option for
  # old cfgs). Per-event offset (this_id) added at use sites.
  seed_base <- if (!is.null(cfg$seed)) cfg$seed else getOption(
    "EpiATLAS_AS_SEED", 42L
  )

  # Build the three explanatory-variable sets (long/short/local) from a feature
  # table's columns. Factored out (was inline) so the FEATURE-rotation controls
  # can build their own sets from their own feature tables. `ev_id` selects the
  # per-event vicinity-narrowed chromHMM subset (for the "short" set); `resp` is
  # the outcome column to exclude.
  #   long  = all explanatory; short = drop vicinity-far chromHMM; local = drop
  #           all chromHMM.
  # Blocklist excludes id/meta/grouping/outcome, the pooled `gene_expression_getmm`
  # (event models use vst), AND `IJC`/`SJC`/`PSI` — the latter are ILLEGAL
  # features: PSI = IJC/(IJC+SJC), so using the junction counts as predictors is
  # trivial target leakage (same rule as splicing_ml's IJC/SJC/PSI drop).
  #
  # ---------------------------------------------------------------------------
  # AUDIT NOTE 2026-07-30 — KNOWN INCONSISTENCY WITH TIER 1. NOT FIXED (fixing it
  # would change results, and the current screen results were produced with this
  # behaviour). Documented here so nobody assumes the two tiers agree.
  #
  # This blocklist is LOOSER than Tier-1's `screen_partition_columns()`
  # (09-ml-shared.R). Columns Tier-1 deliberately excludes but that become
  # elastic-net PREDICTORS here:
  #   H3K*_source (6 cols)  — observed-vs-imputed PROVENANCE per mark
  #   qc_flag_count         — per-sample QC covariate
  #   width;* (3 cols)      — event geometry
  #   distance_TES, distance_gene_start, distance_gene_end
  #
  # Measured on 40 random real feature tables: `H3K*_source` VARIES within the
  # event in 100% of them (mean 5.75 of 6 columns), and `qc_flag_count` varies in
  # 100%. So these are live predictors for essentially every event, not a corner
  # case. The nominal-level fixing below only drops a `_source` column when it is
  # GLOBALLY constant across the event, which is almost never.
  #
  # Why it matters: `H3K*_source` encodes whether that mark's signal was imputed.
  # Imputation status tracks which epigenome a sample is (not all IHEC epigenomes
  # have all marks), so it is a direct cell-type/batch proxy — precisely the
  # fingerprint confound this project exists to expose — and it also flags that the
  # other features for that sample are model-derived. A "significant" chromatin
  # coefficient could be reading imputation status instead of biology.
  # `width;*` and the distance_* columns are constant within an event, so
  # `step_zv()` removes them; they are harmless, unlike the two above.
  #
  # Consequence for interpretation: when reading coef_table, check whether any
  # `H3K*_source_*` dummy or `qc_flag_count` carries weight before attributing a
  # model to chromatin. Currently dormant (Tier-2 only runs on Tier-1 hits, and the
  # last completed screen produced none), but live the moment there are hits.
  # ---------------------------------------------------------------------------
  build_explanatory_vars <- function(feat_dt, ev_id, resp) {
    explanatory <- names(feat_dt)[
      !names(feat_dt) %in%
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
          "gene_expression_getmm",
          "IJC",
          "SJC",
          "PSI",
          grouping_col,
          resp
        )
    ]
    chromhmm_explanatory <- explanatory[grepl(
      "chromhmm",
      explanatory,
      fixed = TRUE
    )]
    smaller_chromhmm_ids <- chromhmm_hits_smaller[
      chromhmm_hits_smaller[, "queryHits"] ==
        which(keep_rows_manual == ev_id),
      "subjectHits"
    ]
    # init = all-TRUE so an empty smaller-id set (possible for a control event)
    # doesn't error in Reduce and correctly marks every chromHMM col as "far".
    old_chromhmm_explanatory <- chromhmm_explanatory[Reduce(
      `&`,
      lapply(sprintf("chromhmm_%d", smaller_chromhmm_ids), function(suff) {
        !endsWith(chromhmm_explanatory, suff)
      }),
      rep(TRUE, length(chromhmm_explanatory))
    )]
    # AUDIT NOTE 2026-07-30 — this mapping is POSITIONAL and therefore fragile.
    # It assumes `feature_sets` is exactly c("long", "short", "local") IN THAT
    # ORDER (it is, defined at 09-1-ml-local.R:53), because the list below is
    # built as all / drop-far-chromHMM / drop-all-chromHMM. Reorder or subset
    # `feature_sets` and the labels silently swap — a `long` model would be
    # reported as `local`, with no error anywhere. NOT changed (a named-list
    # construction would be the fix, but it is a behaviour-neutral refactor only
    # if the vector really is in this order, and touching it now would put the
    # current results in question). If you ever make `feature_sets` configurable,
    # fix this first.
    setNames(
      list(
        explanatory,
        explanatory[!explanatory %in% old_chromhmm_explanatory],
        explanatory[!explanatory %in% chromhmm_explanatory]
      ),
      feature_sets
    )
  }

  tryCatch(
    {
      feature_data <- data.table::fread(file.path(
        feature_table_dir,
        paste0("feature_table_", this_id, ".csv.gz")
      ))

      explanatory_vars <- build_explanatory_vars(
        feature_data, this_id, response
      )

      # Event-level skip: require the full `nfolds` distinct ontology supergroups
      # (2026-07-22). An event whose present ontology labels cluster into fewer
      # than nfolds biology supergroups cannot support an nfolds-way grouped CV.
      # Rather than silently reduce the fold count and fit an underpowered model
      # — or crash on the degenerate columns such low-diversity events tend to
      # carry (globally-constant `H3K*_source`, all-NA-in-fold numeric marks) —
      # skip the whole event with a logged reason. Ontology is a per-sample
      # property, so this one check is authoritative for the event. A stub
      # `_all_metrics` is still written so 09-1's already_computed glob does not
      # re-dispatch it, plus a `_skipped` marker recording why.
      event_supergroups <- .map_ontology_to_supergroups(
        feature_data[[grouping_col]],
        n_groups = nfolds
      )
      n_surviving_groups <- length(unique(event_supergroups))
      if (n_surviving_groups < nfolds) {
        message(sprintf(
          paste0(
            "SKIP event %d: only %d ontology supergroup(s) survive (need %d) ",
            "-- too few for %d-fold biology-grouped CV"
          ),
          this_id, n_surviving_groups, nfolds, nfolds
        ))
        data.table::fwrite(
          data.table::data.table(
            ID = this_id,
            skipped = TRUE,
            reason = "insufficient_ontology_supergroups",
            n_supergroups = n_surviving_groups,
            required = nfolds
          ),
          file.path(event_dir, paste0(this_id, "_skipped.csv.gz"))
        )
        # sentinel: keeps 09-1 from re-dispatching (already_computed globs this)
        data.table::fwrite(
          data.table::data.table(ID = this_id),
          file.path(event_dir, paste0(this_id, "_all_metrics.csv.gz"))
        )
        quit(save = "no", status = 0)
      }

      parallel <- if (length(args) >= 3L) as.integer(args[3L]) else 1L
      set.seed(seed_base + this_id)
      print(glue::glue(
        "Fitting event {this_id} (real PSI only; Tier-2 interpretation)..."
      ))

      # Tier-2 interpretation: fit THIS event's real PSI on its own feature
      # matrix, across the long/short/local sets. Significance ("is this event
      # epigenetically predictable?") is decided upstream by the Tier-1 fit-free
      # ridge screen (09s-ridge-screen.R) — 09zz only runs on screen HITS, so
      # there is no rotation null here. Every wflow_id is "{fs}_PSI_glmnet"
      # (response_type == "real"); the parallelism is now across the 3 feature
      # sets inside run_event_glmnet (mc.cores = parallel).
      wrs <- run_event_glmnet(
        this_feature_data = feature_data,
        explanatory_vars = explanatory_vars,
        response = response,
        rotated_psis = NULL,
        grouping_col = grouping_col,
        nfolds = nfolds,
        seed = seed_base + this_id,
        parallel = parallel,
        cv_folds_path = file.path(event_dir, paste0(this_id, "_cv_folds.csv.gz"))
      )
      rm(feature_data, explanatory_vars, chromhmm_hits_smaller, keep_rows_manual)
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
