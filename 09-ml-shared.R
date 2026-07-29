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


# SCREEN_STAT_VERSION — provenance stamp for the Tier-1 statistic below.
#
# BUMP THIS whenever prep_event() or ridge_screen_stat() changes in a way that
# alters the numbers they produce. It is written into every per-event screen
# output by 09s-ridge-screen.R and checked by that script's resume gate, so a
# changed statistic AUTO-INVALIDATES stale results instead of being silently
# reused. Without it, a correctness fix is invisible to the reuse logic: Snakemake
# re-runs all ~34k per-event jobs (their input script changed) but each one exits
# "Already computed" because its output file exists, and the run completes having
# recomputed nothing. That is not hypothetical — it is exactly what the
# 2026-07-29 leak fix would have hit.
#
# Rows with NO stamp are treated as version 1 (the pre-2026-07-29 statistic).
#
#   1 = original: confound residualisation, standardisation, kernel basis and GCV
#       lambda all fit GLOBALLY on all n rows, with an exact block-leave-one-group-
#       out formula applied on top -- a train/test LEAK, since every ingredient of
#       that "exact" formula depended on the held-out group. Results INVALID.
#   2 = 2026-07-29 fix: all of the above fit per fold on TRAIN rows only and
#       applied to held-out rows via the train-fitted transform/coefficients/lambda.
#
# A hash of the two function bodies was considered instead (no discipline needed)
# but rejected: it would also invalidate on a comment-only edit, and a false
# invalidation here costs a multi-hour 34k-job rerun.
SCREEN_STAT_VERSION <- 2L


# prep_event(): precompute the parts fixed across the real fit AND all null
# rotations (they share this event's PSI + confounds + samples + groups). Only
# X changes per rotation, so this is computed ONCE per event.
#   y            — PSI vector (already NA-filtered, aligned to the sample rows)
#   confound_df  — data.frame of confound columns (protocol factor + numerics)
#   groups       — leave-one-group-out CV fold label per row (same vector
#                  passed to ridge_screen_stat below)
# Confound (Z) adjustment is now fit PER FOLD, train rows only, and applied to
# that fold's test rows via the train-fitted coefficients — a held-out group's
# own PSI never contributes to the confound coefficient used to residualise
# its own target. (Fixed 2026-07-29: previously `qr(Z)`/`qr.resid` were fit
# once on ALL n rows, so a group's PSI leaked into its own adjustment.)
# Returns: Z (full confound design), group_idx (row indices per fold),
# fold_z (per-fold qr(Z_train), reused by ridge_screen_stat to residualise X
# the same way), ytil_train (per-fold in-sample residual of TRAIN y — the
# ridge target within that fold), ytil_oof (out-of-fold residual of every
# row's own held-out fold — the target the OOF prediction is scored against),
# ss_tot (built from ytil_oof, so R2's denominator matches its numerator).
prep_event <- function(y, confound_df, groups) {
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

  group_idx <- split(seq_len(n), groups)
  fold_z <- vector("list", length(group_idx))
  names(fold_z) <- names(group_idx)
  ytil_train <- vector("list", length(group_idx))
  names(ytil_train) <- names(group_idx)
  ytil_oof <- numeric(n)
  for (g in names(group_idx)) {
    test_idx <- group_idx[[g]]
    train_idx <- setdiff(seq_len(n), test_idx)
    qrZtr <- qr(Z[train_idx, , drop = FALSE])
    b_y <- qr.coef(qrZtr, y[train_idx])
    fold_z[[g]] <- qrZtr
    ytil_train[[g]] <- qr.resid(qrZtr, y[train_idx])
    ytil_oof[test_idx] <- y[test_idx] - Z[test_idx, , drop = FALSE] %*% b_y
  }

  list(
    Z = Z, n = n, group_idx = group_idx, fold_z = fold_z,
    ytil_train = ytil_train, ytil_oof = ytil_oof, ss_tot = sum(ytil_oof^2)
  )
}


# ridge_screen_stat(): group-leave-one-out CV partial-R² and CCC of ridge
# predicting PSI (residualised on confounds) from X, in dual (kernel) form.
# Returns list(R2, CCC, lambda, df, n_features) or an all-NA list if X has no
# usable columns.
#
# Rewritten 2026-07-29 to close every train/test leak in this function:
# standardisation (mu/sd), confound residualisation, the ridge kernel basis
# (K/eigendecomposition), and λ selection (GCV) are now ALL computed from
# TRAIN rows only, per fold, and applied to that fold's held-out rows purely
# by re-using the train-fitted transform/coefficients/λ — nothing about a
# held-out group's own X or y ever touches the model fit used to predict it.
# (The previous version built one global smoother from ALL n rows — via
# global standardisation, global confound-projection, and a single K/eigen
# basis — then used an exact block-CV formula on TOP of that shared
# smoother. That formula's exactness assumed the smoother itself didn't
# depend on the held-out group, which was false: every ingredient of it did.)
#
# Cost: this now eigendecomposes K_train once PER FOLD instead of once
# globally — G decompositions of a ~(n - n/G)-sized matrix instead of 1 of an
# n-sized one, a roughly G·((G-1)/G)^3 factor (≈1–3x for the G=3–5 folds used
# here), not the "G× worse" it looks like at first glance, because each
# fold's matrix is correspondingly smaller. No more block-inversion
# singularity edge case either (no (I − S_BB) solve) since prediction is a
# plain train→test kernel projection, not a same-smoother block update.
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

  # Global (X-only, never touches y) column filter: drop features with no
  # variance anywhere in this event. Unsupervised data cleaning, not a target
  # leak — just fixes the column SET so p doesn't ragged-vary fold to fold.
  sdv_all <- apply(X, 2L, function(col) stats::sd(col, na.rm = TRUE))
  keep <- is.finite(sdv_all) & sdv_all > 1e-8
  if (!any(keep)) {
    return(na_out)
  }
  X <- X[, keep, drop = FALSE]
  p <- ncol(X)

  oof <- numeric(prep$n)
  lambdas <- numeric(0)
  dfs <- numeric(0)
  for (g in names(prep$group_idx)) {
    test_idx <- prep$group_idx[[g]]
    train_idx <- setdiff(seq_len(prep$n), test_idx)
    Xtr <- X[train_idx, , drop = FALSE]
    Xte <- X[test_idx, , drop = FALSE]

    # center + mean-impute + unit variance from TRAIN rows only; the SAME
    # train-fitted mu/sd is then applied to the held-out test rows.
    mu <- colMeans(Xtr, na.rm = TRUE)
    mu[!is.finite(mu)] <- 0
    Xtr <- sweep(Xtr, 2L, mu, "-")
    Xte <- sweep(Xte, 2L, mu, "-")
    Xtr[!is.finite(Xtr)] <- 0
    Xte[!is.finite(Xte)] <- 0
    sdv <- sqrt(colSums(Xtr^2) / max(1L, (nrow(Xtr) - 1L)))
    # degenerate only within this fold's train subset (though not globally,
    # per the filter above) -> leave unscaled (still centered) rather than
    # dropping the column, so p stays fixed across folds.
    sdv[!is.finite(sdv) | sdv <= 1e-8] <- 1
    Xtr <- sweep(Xtr, 2L, sdv, "/")
    Xte <- sweep(Xte, 2L, sdv, "/")

    # residualise X on Z using this fold's TRAIN-fitted confound coefficients
    # (prep$fold_z[[g]] was fit on train rows only in prep_event) — applied to
    # test rows by explicit matrix projection, never refit on test.
    qrZtr <- prep$fold_z[[g]]
    b_X <- qr.coef(qrZtr, Xtr)
    MXtr <- qr.resid(qrZtr, Xtr)
    MXte <- Xte - prep$Z[test_idx, , drop = FALSE] %*% b_X

    ytil_tr <- prep$ytil_train[[g]]
    n_tr <- length(ytil_tr)

    Ktr <- tcrossprod(MXtr) # train × train Gram (train rows only)
    eg <- eigen(Ktr, symmetric = TRUE)
    d <- pmax(eg$values, 0)
    U <- eg$vectors
    ystar <- crossprod(U, ytil_tr)

    fold_grid <- grid
    if (is.null(fold_grid)) {
      base <- sum(d) / n_tr
      if (!is.finite(base) || base <= 0) base <- 1
      fold_grid <- base * 10^seq(-3, 3, length.out = 25L)
    }
    # standard trace-based GCV, but now legitimate: every quantity (d, ystar,
    # n_tr) comes from train rows only, so this is a train-internal criterion,
    # not a proxy contaminated by the held-out group.
    gcv <- vapply(fold_grid, function(lam) {
      ssr <- sum((lam / (d + lam) * ystar)^2)
      tr <- sum(d / (d + lam))
      denom <- (n_tr - tr)^2
      if (denom <= 0) return(Inf)
      n_tr * ssr / denom
    }, numeric(1))
    lam <- fold_grid[which.min(gcv)]
    sfilt <- d / (d + lam)

    alpha <- U %*% ((1 / (d + lam)) * ystar) # (K_train + λI)^-1 ỹ_train
    Ktest <- tcrossprod(MXte, MXtr) # test × train cross-Gram
    oof[test_idx] <- as.numeric(Ktest %*% alpha)

    lambdas <- c(lambdas, lam)
    dfs <- c(dfs, sum(sfilt))
  }

  list(
    R2 = 1 - sum((prep$ytil_oof - oof)^2) / prep$ss_tot,
    CCC = .ccc(prep$ytil_oof, oof),
    lambda = mean(lambdas),
    df = mean(dfs),
    n_features = p
  )
}
