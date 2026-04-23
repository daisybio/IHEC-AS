## ---------------------------------------------------------------------------
## Event-specific regularised regression via tidymodels
##
## Fits penalised regression (glmnet) and MARS models for a single splicing
## event using biology-guided grouped CV folds.  Both PSI and logit-scale PSI
## response variants are fitted in parallel; the best workflow (highest CV CCC
## on the PSI scale) is returned along with robust non-zero lasso coefficients.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidymodels)
  library(glmnet) # required by parsnip's glmnet engine
  library(earth) # required by parsnip's MARS engine
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
# PSI / logit-scale metrics  (module-level: defined once, reused per call)
# ---------------------------------------------------------------------------
# Scale is auto-detected from truth: any value outside [0, 1] → logit scale.
#   _psi   variants: logit truth/estimate → plogis(); PSI → identity
#   _logit variants: PSI  truth/estimate → qlogis(.clamp()); logit → identity
# All workflows share one metric set; model selection uses rsq_trad_psi
# (traditional R²: 1 - SS_res/SS_tot) because it directly answers "what
# fraction of PSI variation do the epigenetic features explain?" and does not
# penalise calibration bias between PSI-recipe and logit-recipe workflows.
.clamp <- function(x, eps = 0.01) pmax(pmin(x, 1 - eps), eps)
.to_psi <- function(truth, x) {
  if (any(truth < 0 | truth > 1, na.rm = TRUE)) plogis(x) else x
}
.to_logit <- function(truth, x) {
  if (any(truth < 0 | truth > 1, na.rm = TRUE)) x else qlogis(.clamp(x))
}

rmse_psi_vec <- function(truth, estimate, na_rm = TRUE, ...) {
  yardstick::rmse_vec(
    .to_psi(truth, truth),
    .to_psi(truth, estimate),
    na_rm = na_rm
  )
}
rmse_psi <- function(data, ...) UseMethod("rmse_psi")
rmse_psi.data.frame <- function(data, truth, estimate, na_rm = TRUE, ...) {
  yardstick::numeric_metric_summarizer(
    "rmse_psi",
    rmse_psi_vec,
    data = data,
    truth = {{ truth }},
    estimate = {{ estimate }},
    na_rm = na_rm
  )
}
rmse_psi <- yardstick::new_numeric_metric(rmse_psi, direction = "minimize")

rsq_trad_psi_vec <- function(truth, estimate, na_rm = TRUE, ...) {
  yardstick::rsq_trad_vec(
    .to_psi(truth, truth),
    .to_psi(truth, estimate),
    na_rm = na_rm
  )
}
rsq_trad_psi <- function(data, ...) UseMethod("rsq_trad_psi")
rsq_trad_psi.data.frame <- function(data, truth, estimate, na_rm = TRUE, ...) {
  yardstick::numeric_metric_summarizer(
    "rsq_trad_psi",
    rsq_trad_psi_vec,
    data = data,
    truth = {{ truth }},
    estimate = {{ estimate }},
    na_rm = na_rm
  )
}
rsq_trad_psi <- yardstick::new_numeric_metric(
  rsq_trad_psi,
  direction = "maximize"
)

ccc_psi_vec <- function(truth, estimate, na_rm = TRUE, ...) {
  yardstick::ccc_vec(
    .to_psi(truth, truth),
    .to_psi(truth, estimate),
    na_rm = na_rm
  )
}
ccc_psi <- function(data, ...) UseMethod("ccc_psi")
ccc_psi.data.frame <- function(data, truth, estimate, na_rm = TRUE, ...) {
  yardstick::numeric_metric_summarizer(
    "ccc_psi",
    ccc_psi_vec,
    data = data,
    truth = {{ truth }},
    estimate = {{ estimate }},
    na_rm = na_rm
  )
}
ccc_psi <- yardstick::new_numeric_metric(ccc_psi, direction = "maximize")

rmse_logit_vec <- function(truth, estimate, na_rm = TRUE, ...) {
  yardstick::rmse_vec(
    .to_logit(truth, truth),
    .to_logit(truth, estimate),
    na_rm = na_rm
  )
}
rmse_logit <- function(data, ...) UseMethod("rmse_logit")
rmse_logit.data.frame <- function(data, truth, estimate, na_rm = TRUE, ...) {
  yardstick::numeric_metric_summarizer(
    "rmse_logit",
    rmse_logit_vec,
    data = data,
    truth = {{ truth }},
    estimate = {{ estimate }},
    na_rm = na_rm
  )
}
rmse_logit <- yardstick::new_numeric_metric(rmse_logit, direction = "minimize")

rsq_trad_logit_vec <- function(truth, estimate, na_rm = TRUE, ...) {
  yardstick::rsq_trad_vec(
    .to_logit(truth, truth),
    .to_logit(truth, estimate),
    na_rm = na_rm
  )
}
rsq_trad_logit <- function(data, ...) UseMethod("rsq_trad_logit")
rsq_trad_logit.data.frame <- function(
  data,
  truth,
  estimate,
  na_rm = TRUE,
  ...
) {
  yardstick::numeric_metric_summarizer(
    "rsq_trad_logit",
    rsq_trad_logit_vec,
    data = data,
    truth = {{ truth }},
    estimate = {{ estimate }},
    na_rm = na_rm
  )
}
rsq_trad_logit <- yardstick::new_numeric_metric(
  rsq_trad_logit,
  direction = "maximize"
)

ccc_logit_vec <- function(truth, estimate, na_rm = TRUE, ...) {
  yardstick::ccc_vec(
    .to_logit(truth, truth),
    .to_logit(truth, estimate),
    na_rm = na_rm
  )
}
ccc_logit <- function(data, ...) UseMethod("ccc_logit")
ccc_logit.data.frame <- function(data, truth, estimate, na_rm = TRUE, ...) {
  yardstick::numeric_metric_summarizer(
    "ccc_logit",
    ccc_logit_vec,
    data = data,
    truth = {{ truth }},
    estimate = {{ estimate }},
    na_rm = na_rm
  )
}
ccc_logit <- yardstick::new_numeric_metric(ccc_logit, direction = "maximize")

# rsq_trad_psi is the selection metric; _logit variants are diagnostic only.
all_metrics <- yardstick::metric_set(
  rmse_psi,
  rsq_trad_psi,
  ccc_psi,
  rmse_logit,
  rsq_trad_logit,
  ccc_logit
)


# ---------------------------------------------------------------------------
# Main function
# ---------------------------------------------------------------------------

#' Fit regularised regression for a single splicing event
#'
#' Builds a \code{workflow_set} of glmnet and MARS models over PSI and
#' logit-scale PSI responses, tunes each via biology-guided grouped CV, and
#' returns the best-performing workflow together with its robust features.
#'
#' @param this_feature_data  data.table: predictors, response, and grouping col.
#' @param explanatory_vars   Named list of predictor sets; each element is a
#'   character vector of column names and the name identifies the set
#'   (e.g. \code{list(local = c("h3k27ac", ...), chromhmm = c("s1", ...))}).
#' @param response      Name of the primary PSI response column.
#' @param rotated_psis  Optional data.table whose columns are additional response
#'   targets (e.g. rotated PSI components).  Recipes are built for each column
#'   but only the primary \code{response} workflows are tuned.
#' @param grouping_col  Ontology column used for biology-guided grouped CV.
#' @param nfolds        Number of CV folds (= number of ontology supergroups).
#' @param seed          Random seed; each workflow also gets a per-index offset
#'   for reproducibility.
#' @param verbose       Passed to \code{tune::control_grid(verbose = ...)}.
#'
#' @return A named list:
#'   \describe{
#'     \item{cv_results}{Workflow set with tuning results attached.}
#'     \item{workflow_results}{List of per-workflow result objects.}
#'     \item{best_workflow}{ID of the winning workflow.}
#'     \item{best_model_type}{\code{"glmnet"} or \code{"mars"}.}
#'     \item{final_fit}{Finalized fitted \code{workflows::workflow}.}
#'     \item{best_params}{Tibble with winning hyperparameters.}
#'     \item{cv_rsq}{Best-workflow CV R² on the PSI scale (traditional formula).}
#'     \item{nonzero_coefs}{Robust non-zero glmnet coefficients (NULL for MARS).}
#'     \item{fold_features}{Per-fold selected features (NULL for MARS).}
#'     \item{robust_features}{Feature names passing the robustness filter (NULL for MARS).}
#'     \item{importance}{MARS variable-importance matrix (NULL for glmnet).}
#'   }
run_event_glmnet <- function(
  this_feature_data,
  explanatory_vars,
  response,
  rotated_psis = NULL,
  grouping_col = "ontology",
  nfolds = 5L,
  seed = 1234L,
  parallel = FALSE
) {
  # -- Ontology-guided grouped CV folds ---------------------------------------
  set.seed(seed)
  this_feature_data[,
    hc_group := .map_ontology_to_supergroups(
      this_feature_data[[grouping_col]],
      n_groups = nfolds
    )
  ]
  cv_folds <- rsample::group_vfold_cv(
    this_feature_data,
    group = "hc_group",
    v = nfolds,
    balance = "observations"
  )

  # -- Recipes ----------------------------------------------------------------
  # Build one base recipe per response (primary + any rotated PSI columns), then
  # one variant per explanatory set, then duplicate with a logit outcome transform.

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
  # protocol × gene_expression interaction terms.
  explanatory_recipe_list <- unlist(
    sapply(
      base_recipes,
      function(base_recipe) {
        sapply(
          explanatory_vars,
          function(var) {
            base_recipe |>
              update_role(all_of(var), new_role = "predictor") |>
              step_rm(has_role("misc")) |>
              step_impute_mean(all_numeric_predictors()) |>
              step_zv(all_predictors()) |>
              step_dummy(all_nominal_predictors()) |>
              step_interact(terms = ~ starts_with("protocol"):gene_expression) # protocol × gene_expression
          },
          simplify = FALSE
        )
      },
      simplify = FALSE
    ),
    recursive = FALSE
  )

  # unlist() joins nested list names with "." — normalise to "_" so that
  # workflow_set produces predictable IDs like "PSI_long_glmnet".
  names(explanatory_recipe_list) <- gsub(
    ".",
    "_",
    names(explanatory_recipe_list),
    fixed = TRUE
  )

  # Logit-scale variants: step_logit(skip = TRUE) transforms the outcome during
  # prep but keeps predictions on the logit scale (back-transformed after tuning).
  logit_recipe_list <- lapply(explanatory_recipe_list, function(rec) {
    rec |> step_logit(all_outcomes(), offset = 0.01, skip = TRUE)
  })
  names(logit_recipe_list) <- paste0("logit_", names(explanatory_recipe_list))

  recipe_list <- c(explanatory_recipe_list, logit_recipe_list)

  # -- Model specifications ---------------------------------------------------
  # glmnet: mixture and penalty both tuned.  The submodel trick lets tune() fit
  # ONE path per (mixture, fold) and evaluate all 50 lambda values via
  # multi_predict() — ~6 fits per fold yield 300 evaluations.
  glmnet_spec <- linear_reg(penalty = tune(), mixture = tune()) |>
    set_engine("glmnet")

  # MARS: additive only (prod_degree = 1) with up to 50 basis-function terms.
  # Interactions (degree 2) are O(p²) in the forward pass and never outperformed
  # additive models in CV — cross-tissue folds prevent pairwise interactions from
  # generalising across tissue types.
  mars_spec <- mars(num_terms = tune(), prod_degree = 1L) |>
    set_engine("earth") |>
    set_mode("regression")

  mars_params <- mars_spec |>
    extract_parameter_set_dials() |>
    update(num_terms = num_terms(range = c(2L, 50L)))

  # -- Workflow set & grids ---------------------------------------------------
  wflow_set <- workflow_set(
    preproc = recipe_list,
    models = list(glmnet = glmnet_spec, mars = mars_spec)
  )

  # glmnet: structured crossing grid — 6 alpha × 50 lambda values.
  # Alpha values skew toward lasso (mixture ≈ 1) for sparsity and feature
  # identification.  Intermediate values (0.5–0.9) prevent arbitrary selection
  # among correlated chromHMM features where pure lasso is unstable.
  glmnet_grid <- tidyr::crossing(
    penalty = 10^seq(-5, 0, length.out = 50),
    mixture = c(0.1, 0.5, 0.7, 0.9, 0.95, 0.99)
  ) |>
    arrange(mixture, penalty) # ensure lambda path is contiguous for each alpha

  # MARS: 8 num_terms candidates (prod_degree fixed at 1 in the spec).
  mars_grid <- tibble::tibble(
    num_terms = c(2L, 5L, 10L, 15L, 20L, 30L, 40L, 50L)
  )

  params_list <- list(
    glmnet = extract_parameter_set_dials(glmnet_spec),
    mars = mars_params
  )

  for (wid in wflow_set$wflow_id) {
    model_key <- if (endsWith(wid, "glmnet")) "glmnet" else "mars"
    wflow_set <- wflow_set |>
      option_add(param_info = params_list[[model_key]], id = wid) |>
      option_add(
        grid = if (model_key == "glmnet") glmnet_grid else mars_grid,
        id = wid
      )
  }

  # -- Tune primary + rotated-response workflows (glmnet only) ----------------
  # all logit variants skipped; only glmnet built, no MARS
  wflow_to_tune <- wflow_set |>
    dplyr::filter(
      endsWith(wflow_id, "glmnet") &
        # (startsWith(wflow_id, "PSI_")) &
        (!startsWith(wflow_id, "logit_"))
    )

  n_wflows <- nrow(wflow_to_tune)

  # -- Feature extraction helpers ---------------------------------------------

  # Non-zero glmnet coefficients at the selected penalty, excluding the
  # intercept and gene_expression (used as an interaction base term).
  .extract_glmnet_coefs <- function(entry) {
    penalty_val <- entry$best_params$penalty[[1L]]
    engine <- workflows::extract_fit_engine(entry$final_fit)
    coef_mat <- as.matrix(coef(engine, s = penalty_val))
    mask <- coef_mat[, 1L] != 0 &
      !rownames(coef_mat) %in% c("(Intercept)", "gene_expression")
    coef_mat[mask, 1L, drop = FALSE]
  }

  # Robustness filter: a feature must be non-zero in the full-data fit AND
  # appear in >= 2 CV folds (uniqueN(fold) >= 3 counts the "full" label too).
  # Replicates the filter used in 07-ml-helper.R.
  .robust_glmnet_coefs <- function(entry, cv_folds, this_feature_data) {
    full_coefs <- .extract_glmnet_coefs(entry)
    full_features <- rownames(full_coefs)

    if (length(full_features) == 0L) {
      return(list(
        nonzero_coefs = full_coefs,
        fold_features = list(),
        robust_features = character(0L)
      ))
    }

    penalty_val <- entry$best_params$penalty[[1L]]
    fold_features <- lapply(cv_folds$splits, function(split) {
      fold_fit <- parsnip::fit(entry$final_wf, data = rsample::analysis(split))
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

  # MARS variable importance via evimp(); nsubsets = number of model subsets in
  # which each variable appeared (analogous to the fold-count robustness criterion).
  .extract_mars_importance <- function(entry) {
    tryCatch(
      earth::evimp(workflows::extract_fit_engine(entry$final_fit)),
      error = function(e) NULL
    )
  }

  # -- Tune, finalize, fit, and extract features in a single pass -------------
  # Sequential: only 6 workflows per event; nested forking inside the outer
  # pbmclapply (09-1-ml-local.R) wastes resources and conflicts with progress.
  workflow_results <- pbmcapply::pbmclapply(
    seq_len(n_wflows),
    function(i) {
      set.seed(seed + i) # per-workflow offset for reproducibility
      wid <- wflow_to_tune$wflow_id[[i]]
      wf <- workflowsets::extract_workflow(wflow_to_tune, id = wid)
      this_grid <- wflow_to_tune$option[[i]]$grid
      specs <- workflowsets::extract_spec_parsnip(wflow_to_tune, wid)

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

      # select_best for all model types: pick the hyperparameters with the highest
      # mean CV R² (rsq_trad_psi).  A null or near-null result is a meaningful
      # outcome — for the primary response it means no cross-tissue signal; for
      # rotated-response negative controls it is the expected result.
      # The 1-SE rule was previously applied to MARS only, but the asymmetry is
      # unjustified: both models have built-in regularisation (lasso penalty /
      # GCV pruning) that already controls complexity without a second heuristic.
      best_params <- select_best(tuned, metric = "rsq_trad_psi")
      final_wf <- tune::finalize_workflow(wf, best_params)
      final_fit <- parsnip::fit(final_wf, data = this_feature_data)

      # CV metric at the selected config (used for cross-workflow comparison)
      final_metrics <- collect_metrics(tuned) |>
        subset(.config == best_params$.config)
      cv_rsq <- subset(final_metrics, .metric == "rsq_trad_psi")$mean[[1L]]

      # In-sample metrics (diagnostic).  Logit-recipe predictions are on the logit
      # scale; back-transform to PSI so truth and estimate share the same scale
      # for the auto-detection logic inside all_metrics.
      fit_preds <- predict(final_fit, new_data = this_feature_data)$.pred
      if (grepl("logit", wid)) {
        fit_preds <- plogis(fit_preds)
      }
      fit_metrics <- all_metrics(
        data = tibble(
          truth = this_feature_data[[response]],
          estimate = fit_preds
        ),
        truth = truth,
        estimate = estimate
      )

      entry <- list(
        wflow_id = wid,
        specs = specs,
        best_params = best_params,
        final_wf = final_wf,
        final_fit = final_fit,
        cv_rsq = cv_rsq,
        tuned = tuned,
        final_metrics = final_metrics,
        fit_metrics = fit_metrics
      )

      if (inherits(specs, "linear_reg")) {
        feat <- .robust_glmnet_coefs(entry, cv_folds, this_feature_data)
        c(entry, list(model_type = "glmnet"), feat)
      } else if (inherits(specs, "mars")) {
        c(
          entry,
          list(
            model_type = "mars",
            importance = .extract_mars_importance(entry)
          )
        )
      } else {
        c(entry, list(model_type = "unknown"))
      }
    }
  )
  names(workflow_results) <- wflow_to_tune$wflow_id

  # Reconstruct grid_results workflowset (mirrors workflow_map output structure)
  # grid_results <- wflow_to_tune
  # grid_results[["result"]] <- lapply(workflow_results, `[[`, "tuned")
  workflow_results <- lapply(workflow_results, function(x) {
    x$tuned <- tuned |> collect_metrics(summarize = FALSE) # keep all metrics, not just the best config
    x
  })

  # -- Best workflow = highest CV R² among primary-response workflows only ----
  # Restrict to PSI / logit_PSI — rotated workflows are negative controls and
  # should not compete for "best" (their response is a different event's PSI).
  # primary_names <- names(workflow_results)[
  #   startsWith(names(workflow_results), response) |
  #     startsWith(names(workflow_results), paste0("logit_", response))
  # ]
  # best_idx <- which.max(vapply(
  #   workflow_results[primary_names],
  #   `[[`,
  #   numeric(1L),
  #   "cv_rsq"
  # ))
  # best_result <- workflow_results[primary_names][[best_idx]]

  # -- Return -----------------------------------------------------------------
  return(list(
    workflow_results = workflow_results,
    cv_folds = rbindlist(mapply(
      cv_folds$splits,
      cv_folds$id,
      FUN = function(split, id) {
        split |>
          assessment() |>
          select(all_of(grouping_col), hc_group) |>
          mutate(resample = id) |>
          count(across(everything()))
      },
      SIMPLIFY = FALSE
    ))
  ))
  # list(
  #   # cv_results = grid_results,
  #   workflow_results = workflow_results,
  #   best_workflow = best_result$wflow_id,
  #   best_model_type = best_result$model_type,
  #   final_fit = best_result$final_fit,
  #   best_params = best_result$best_params,
  #   cv_rsq = best_result$cv_rsq,
  #   # glmnet-specific (NULL for MARS)
  #   nonzero_coefs = best_result$nonzero_coefs,
  #   fold_features = best_result$fold_features,
  #   robust_features = best_result$robust_features,
  #   # MARS-specific (NULL for glmnet)
  #   importance = best_result$importance
  # )
}
