## ---------------------------------------------------------------------------
## Shared helpers for the event-specific model pipeline
##
## Sourced by BOTH the Tier-1 ridge screen (09s-ridge-screen.R) and the Tier-2
## elastic-net interpreter (09zz-ml-event-glmnet-tidymodels.R). Pure function /
## data definitions only — no top-level side effects, safe to source repeatedly.
##
## Contents:
##   1. Ontology → biology-guided CV-supergroup mapping (was inline in 09zz).
##   2. screen_partition_columns() — split a feature table's columns into the
##      epigenetic design X, the unpenalised confounds Z, and everything excluded
##      (meta / geometry / QC / target-leakage).
##   3. resolve_supergroup_folds() — the actual_nfolds/collapse guard (was inline
##      in 09zz::run_event_glmnet) as a standalone helper.
##   4. Tier-1 closed-form statistic: prep_event() + ridge_screen_stat().
## ---------------------------------------------------------------------------

# ===========================================================================
# 1. Ontology hierarchy — biology-guided CV supergroup assignment
# ===========================================================================
# Paths are coarse-to-fine across 4 levels:
#   [super] > [germ_layer/lineage] > [subdomain] > [leaf]
# Three super-groups: "blood", "somatic", "developmental".
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
  "peripheral blood" = c("blood", "immune", "hematopoietic", "peripheral_blood"),
  # Somatic — ectodermal (CNS + neural-crest)
  "brain" = c("somatic", "ectodermal", "nervous_system", "brain"),
  "nervous system" = c("somatic", "ectodermal", "nervous_system", "nervous_system"),
  "neural" = c("somatic", "ectodermal", "nervous_system", "neural"),
  "melanocyte" = c("somatic", "ectodermal", "neural_crest", "melanocyte"),
  # Somatic — endodermal (gut, liver/pancreas, lung)
  "digestive system" = c("somatic", "endodermal", "digestive", "digestive_system"),
  "colon" = c("somatic", "endodermal", "digestive", "colon"),
  "mucosa" = c("somatic", "endodermal", "digestive", "mucosa"),
  "epithelial" = c("somatic", "endodermal", "digestive", "epithelial"),
  "liver" = c("somatic", "endodermal", "hepatopancreatic", "liver"),
  "pancreas" = c("somatic", "endodermal", "hepatopancreatic", "pancreas"),
  "endoderm-derived structure" = c("somatic", "endodermal", "general", "endodermal_structure"),
  "lung" = c("somatic", "endodermal", "respiratory", "lung"),
  # Somatic — mesodermal (connective tissue, muscle, kidney)
  "connective tissue cell" = c("somatic", "mesodermal", "connective", "connective_tissue"),
  "mesoderm-derived structure" = c("somatic", "mesodermal", "connective", "mesodermal_structure"),
  "muscle" = c("somatic", "mesodermal", "muscle", "muscle"),
  "kidney" = c("somatic", "mesodermal", "renal", "kidney"),
  # Developmental — stem / progenitor / embryonic
  "stem cell" = c("developmental", "stem_progenitor", "stem", "stem_cell"),
  "embryonic cell (metazoa)" = c("developmental", "stem_progenitor", "embryonic", "embryonic_cell"),
  "extraembryonic cell" = c("developmental", "stem_progenitor", "extraembryonic", "extraembryonic_cell"),
  "placenta" = c("developmental", "stem_progenitor", "extraembryonic", "placenta"),
  # Developmental — transformed (shares stem_progenitor subdomain so cancer's
  # internal distance = 0.5 < somatic germ-layer distance = 0.75)
  "cancer cell line" = c("developmental", "stem_progenitor", "transformed", "cancer_cell_line")
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
    paste0("hc_", paste(prefix[seq_len(min(2L, length(prefix)))], collapse = "_"))
  } else {
    paste0("hc_cluster_", cluster_id)
  }
}


# Maps a character vector of ontology labels to biology-guided supergroup names.
# Clusters only the labels PRESENT in `ontology_vec` (not the full 34-label
# vocabulary) — see the 2026-07-22 fix note in the git history: clustering the
# full universe collapsed blood-only events to 1 supergroup (v_fold crash).
.map_ontology_to_supergroups <- function(ontology_vec, n_groups = 5L) {
  norm_ont <- tolower(trimws(as.character(ontology_vec)))
  labels <- sort(unique(norm_ont[norm_ont %in% names(.ontology_bio_paths)]))
  n <- length(labels)

  if (n == 0L) {
    return(rep("other", length(norm_ont)))
  }

  effective_groups <- max(1L, min(n_groups, n))

  if (n == 1L || effective_groups == 1L) {
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


# Resolve the leave-one-group-out fold assignment for one event's samples,
# reproducing 09zz::run_event_glmnet's actual_nfolds/collapse guard.
# Returns a character vector of supergroup labels aligned to `ontology_vec`,
# or NULL if fewer than 2 usable groups result (screen cannot run → NA).
resolve_supergroup_folds <- function(ontology_vec, nfolds = 5L) {
  n_unique_ont <- length(unique(tolower(trimws(as.character(ontology_vec)))))
  actual_nfolds <- max(3L, min(nfolds, n_unique_ont))
  groups <- .map_ontology_to_supergroups(ontology_vec, n_groups = actual_nfolds)
  if (length(unique(groups)) < 2L) {
    return(NULL)
  }
  groups
}


# ===========================================================================
# 2. Feature-table column partition (screen)
# ===========================================================================
# Split a per-event feature table's columns into:
#   x_cols       — epigenetic design matrix to TEST (everything not excluded and
#                  not a confound)
#   confound_cols— unpenalised covariates for Z (protocol, gene_expression_vst,
#                  distance_TSS): partial these out so the screen R² is
#                  "epigenetic signal beyond confounds".
# Everything else is excluded: IDs/meta, event geometry (constant within an
# event → would drop under the variance filter anyway, excluded explicitly for
# clarity), QC/provenance, the pooled getmm expression, and the
# target-leakage columns IJC/SJC/PSI (PSI = IJC/(IJC+SJC)).
#
# Column selection is by an explicit blocklist (mirrors 09zz::build_explanatory_vars)
# plus pattern rules, so a genuinely new epigenetic column is auto-included, and
# `_source`/`width;`/geometry meta are auto-excluded.
screen_partition_columns <- function(col_names, grouping_col = "ontology") {
  confound_cols <- intersect(
    c("protocol", "gene_expression_vst", "distance_TSS"),
    col_names
  )

  exact_block <- c(
    "IHEC", "ID", "Event Type", "Variability", "seqnames", "gene_id",
    "uuid", "transcript_filter", "project",
    "harmonized_sample_ontology_term_high_order_fig1",
    "gene_expression_getmm", "gene_expression_vst",
    "protocol", "distance_TSS", "distance_TES",
    "distance_gene_start", "distance_gene_end",
    "qc_flag_count", "IJC", "SJC", "PSI",
    grouping_col
  )

  is_blocked <- col_names %in% exact_block |
    grepl("_source$", col_names) | # H3K*_source provenance
    grepl("^width;", col_names) # event-geometry widths (constant within event)

  x_cols <- col_names[!is_blocked]
  list(x_cols = x_cols, confound_cols = confound_cols)
}


# ===========================================================================
# 3. Tier-1 closed-form ridge screen statistic
# ===========================================================================
# Lin's concordance correlation coefficient (population / ÷n moments), matching
# yardstick::ccc(bias = TRUE).
.ccc <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 2L) {
    return(NA_real_)
  }
  mx <- mean(x)
  my <- mean(y)
  vx <- mean((x - mx)^2)
  vy <- mean((y - my)^2)
  cxy <- mean((x - mx) * (y - my))
  denom <- vx + vy + (mx - my)^2
  if (denom <= 0) {
    return(NA_real_)
  }
  2 * cxy / denom
}


# prep_event(): precompute the parts fixed across the real fit AND all null
# rotations (they share this event's PSI + confounds + samples). Only X changes
# per rotation, so this is computed ONCE per event.
#   y            — PSI vector (already NA-filtered, aligned to the sample rows)
#   confound_df  — data.frame of confound columns (protocol factor + numerics)
# Returns qrZ (QR of the confound design incl. intercept), ytil (PSI residualised
# on Z), ss_tot.
prep_event <- function(y, confound_df) {
  n <- length(y)
  Z <- matrix(1.0, nrow = n, ncol = 1L, dimnames = list(NULL, "(Intercept)"))
  if (!is.null(confound_df) && ncol(confound_df) > 0L) {
    for (nm in names(confound_df)) {
      v <- confound_df[[nm]]
      if (is.numeric(v)) {
        if (all(is.na(v))) next
        v[is.na(v)] <- mean(v, na.rm = TRUE)
        if (stats::sd(v) <= 1e-12) next # constant within event → absorbed by intercept
        Z <- cbind(Z, setNames(matrix(v, ncol = 1L), NULL))
        colnames(Z)[ncol(Z)] <- nm
      } else {
        f <- factor(ifelse(is.na(v), "NA", as.character(v)))
        if (nlevels(f) < 2L) next # constant → absorbed by intercept
        mm <- stats::model.matrix(~f)[, -1L, drop = FALSE]
        colnames(mm) <- paste0(nm, "_", levels(f)[-1L])
        Z <- cbind(Z, mm)
      }
    }
  }
  qrZ <- qr(Z)
  ytil <- qr.resid(qrZ, y)
  list(qrZ = qrZ, ytil = ytil, ss_tot = sum(ytil^2), n = n)
}


# ridge_screen_stat(): group-leave-one-out CV partial-R² and CCC of ridge
# predicting prep$ytil from X, in dual (kernel) form. λ chosen by GCV on a
# closed-form grid; the reported metric is the exact block leave-one-GROUP-out
# prediction at that λ. Returns list(R2, CCC, lambda, df, n_features) or an
# all-NA list if X has no usable columns / K is degenerate.
ridge_screen_stat <- function(prep, X, groups, grid = NULL) {
  na_out <- list(
    R2 = NA_real_, CCC = NA_real_, lambda = NA_real_,
    df = NA_real_, n_features = 0L
  )
  if (is.null(X) || ncol(X) == 0L || nrow(X) != prep$n) {
    return(na_out)
  }
  X <- as.matrix(X)
  storage.mode(X) <- "double"

  # center + mean-impute (NA → column mean → 0 after centering) + unit variance,
  # dropping constant / all-NA columns.
  mu <- colMeans(X, na.rm = TRUE)
  mu[!is.finite(mu)] <- 0
  X <- sweep(X, 2L, mu, "-")
  X[!is.finite(X)] <- 0
  sdv <- sqrt(colSums(X^2) / max(1L, (nrow(X) - 1L)))
  keep <- is.finite(sdv) & sdv > 1e-8
  if (!any(keep)) {
    return(na_out)
  }
  X <- sweep(X[, keep, drop = FALSE], 2L, sdv[keep], "/")
  p <- ncol(X)

  ytil <- prep$ytil
  MX <- qr.resid(prep$qrZ, X) # residualise features on Z
  K <- tcrossprod(MX) # n×n Gram
  eg <- eigen(K, symmetric = TRUE)
  d <- pmax(eg$values, 0)
  U <- eg$vectors
  n <- length(ytil)
  ystar <- crossprod(U, ytil) # Uᵀ ỹ

  if (is.null(grid)) {
    base <- sum(d) / n
    if (!is.finite(base) || base <= 0) base <- 1
    grid <- base * 10^seq(-3, 3, length.out = 25L)
  }
  gcv <- vapply(grid, function(lam) {
    ssr <- sum((lam / (d + lam) * ystar)^2)
    tr <- sum(d / (d + lam))
    denom <- (n - tr)^2
    if (denom <= 0) return(Inf)
    n * ssr / denom
  }, numeric(1))
  lam <- grid[which.min(gcv)]
  sfilt <- d / (d + lam)

  Sy <- U %*% (sfilt * ystar) # S ỹ  (smoother prediction)
  oof <- numeric(n)
  for (g in unique(groups)) {
    idx <- which(groups == g)
    Ub <- U[idx, , drop = FALSE]
    Sbb <- Ub %*% (sfilt * t(Ub)) # S[idx, idx], |B|×|B|
    resid <- tryCatch(
      solve(diag(length(idx)) - Sbb, ytil[idx] - Sy[idx]),
      error = function(e) {
        # (I − S_BB) singular (group perfectly fit) → fall back to a tiny ridge
        solve(diag(length(idx)) * (1 + 1e-8) - Sbb, ytil[idx] - Sy[idx])
      }
    )
    oof[idx] <- ytil[idx] - resid
  }

  list(
    R2 = 1 - sum((ytil - oof)^2) / prep$ss_tot,
    CCC = .ccc(ytil, oof),
    lambda = lam,
    df = sum(sfilt),
    n_features = p
  )
}
