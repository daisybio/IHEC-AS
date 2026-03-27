# Shared feature preprocessing used by visualization and correlation scripts.

preprocess_mode_config <- function(mode = c("embedding", "correlation")) {
  mode <- match.arg(mode)
  if (mode == "embedding") {
    return(list(
      scale_features = TRUE,
      impute_after_scaling = TRUE,
      drop_sparse = FALSE,
      max_na_frac = 0.5,
      drop_constant = TRUE,
      clip_quantile = 0.99
    ))
  }

  list(
    scale_features = FALSE,
    impute_after_scaling = FALSE,
    drop_sparse = FALSE,
    max_na_frac = 0.5,
    drop_constant = TRUE,
    clip_quantile = 0.99
  )
}

build_feature_masks <- function(feature_names, histone_marks) {
  hist_pattern <- if (length(histone_marks)) {
    paste0("^(", paste(histone_marks, collapse = "|"), ");")
  } else {
    "^$"
  }

  list(
    hist_cols = grepl(hist_pattern, feature_names),
    log_cols = grepl("^width;|^distance_|^gene_expression$", feature_names)
  )
}

filter_feature_matrix <- function(
  feat_mat,
  drop_sparse = FALSE,
  max_na_frac = 0.5,
  drop_constant = TRUE
) {
  if (!is.matrix(feat_mat)) {
    feat_mat <- as.matrix(feat_mat)
  }

  original_cols <- colnames(feat_mat)

  if (drop_sparse && ncol(feat_mat) > 0L) {
    keep_sparse <- colMeans(is.na(feat_mat)) < max_na_frac
    feat_mat <- feat_mat[, keep_sparse, drop = FALSE]
  }

  if (drop_constant && ncol(feat_mat) > 0L) {
    col_var <- apply(feat_mat, 2, var, na.rm = TRUE)
    keep_var <- !is.na(col_var) & col_var > 0
    feat_mat <- feat_mat[, keep_var, drop = FALSE]
  }

  list(
    feat_mat = feat_mat,
    kept_cols = colnames(feat_mat),
    dropped_cols = setdiff(original_cols, colnames(feat_mat))
  )
}

clip_histone_features <- function(feat_mat, hist_cols, clip_quantile = 0.99) {
  if (ncol(feat_mat) == 0L || !any(hist_cols)) {
    return(feat_mat)
  }

  q_hi <- apply(
    feat_mat[, hist_cols, drop = FALSE],
    2,
    quantile,
    probs = clip_quantile,
    na.rm = TRUE
  )
  feat_mat[, hist_cols] <- t(pmin(t(feat_mat[, hist_cols, drop = FALSE]), q_hi))
  feat_mat
}

apply_log_transforms <- function(feat_mat, log_cols) {
  if (ncol(feat_mat) == 0L || !any(log_cols)) {
    return(feat_mat)
  }

  feat_mat[, log_cols] <- log1p(feat_mat[, log_cols])
  feat_mat
}

preprocess_feature_matrix <- function(
  feat_mat,
  histone_marks,
  mode = c("embedding", "correlation"),
  drop_sparse = NULL,
  max_na_frac = NULL,
  drop_constant = NULL,
  clip_quantile = NULL
) {
  mode <- match.arg(mode)
  cfg <- preprocess_mode_config(mode)

  if (!is.null(drop_sparse)) {
    cfg$drop_sparse <- drop_sparse
  }
  if (!is.null(max_na_frac)) {
    cfg$max_na_frac <- max_na_frac
  }
  if (!is.null(drop_constant)) {
    cfg$drop_constant <- drop_constant
  }
  if (!is.null(clip_quantile)) {
    cfg$clip_quantile <- clip_quantile
  }

  # Defensive coercion: mixed-type matrices (e.g. character/factor-like)
  # can reach quantile() in clipping and fail. Non-numeric values become NA.
  if (!is.numeric(feat_mat)) {
    original_colnames <- colnames(feat_mat)
    original_rownames <- rownames(feat_mat)
    feat_mat <- suppressWarnings(matrix(
      as.numeric(feat_mat),
      nrow = nrow(feat_mat),
      ncol = ncol(feat_mat),
      dimnames = list(original_rownames, original_colnames)
    ))
  }

  filtered <- filter_feature_matrix(
    feat_mat = feat_mat,
    drop_sparse = cfg$drop_sparse,
    max_na_frac = cfg$max_na_frac,
    drop_constant = cfg$drop_constant
  )
  feat_mat <- filtered$feat_mat

  if (ncol(feat_mat) == 0L) {
    return(list(
      feat_mat = feat_mat,
      kept_cols = character(),
      dropped_cols = filtered$dropped_cols,
      config = cfg
    ))
  }

  masks <- build_feature_masks(colnames(feat_mat), histone_marks)
  feat_mat <- clip_histone_features(
    feat_mat,
    masks$hist_cols,
    cfg$clip_quantile
  )
  feat_mat <- apply_log_transforms(feat_mat, masks$log_cols)

  if (cfg$scale_features) {
    feat_mat <- scale(feat_mat)
    if (cfg$impute_after_scaling) {
      feat_mat[is.na(feat_mat)] <- 0
    }
  }

  list(
    feat_mat = feat_mat,
    kept_cols = colnames(feat_mat),
    dropped_cols = filtered$dropped_cols,
    config = cfg
  )
}

preprocess_data_table_features <- function(
  dt,
  feature_cols,
  histone_marks,
  mode = c("embedding", "correlation"),
  drop_sparse = NULL,
  max_na_frac = NULL,
  drop_constant = NULL,
  clip_quantile = NULL
) {
  mode <- match.arg(mode)

  pre <- preprocess_feature_matrix(
    feat_mat = as.matrix(dt[, ..feature_cols]),
    histone_marks = histone_marks,
    mode = mode,
    drop_sparse = drop_sparse,
    max_na_frac = max_na_frac,
    drop_constant = drop_constant,
    clip_quantile = clip_quantile
  )

  out_dt <- data.table::copy(dt)

  if (length(pre$kept_cols) > 0L) {
    transformed <- data.table::as.data.table(pre$feat_mat)
    data.table::setnames(transformed, pre$kept_cols)
    out_dt[, (pre$kept_cols) := transformed]
  }

  if (length(pre$dropped_cols) > 0L) {
    out_dt[, (pre$dropped_cols) := NA_real_]
  }

  list(dt = out_dt, info = pre)
}
