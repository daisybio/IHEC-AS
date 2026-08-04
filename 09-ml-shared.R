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
# Column selection is by an explicit blocklist. NOTE it is deliberately STRICTER
# than 09zz::build_explanatory_vars, which blocks only the id/meta/target columns:
# this one additionally excludes `H3K*_source` (observed-vs-imputed provenance,
# which tracks which epigenome a sample is and is therefore a cell-type/batch
# proxy), `width;*` event geometry, distance_TES/gene_start/gene_end, and
# qc_flag_count. Tier-2 currently admits all of those as elastic-net predictors --
# an inconsistency, not a designed asymmetry (see the 09-event-models-redesign
# audit note). Do not "align" this by loosening it.
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
# 2b. Full-cohort feature-table assembly (FEATURE_TABLE_VERSION 2)
# ===========================================================================
# Used by 09-1-ml-local.R's Phase-1 build. Lives here, with every input passed
# explicitly rather than read from 09-1's globals, so a verification script can
# exercise the REAL assembly without also loading 09-1's 13.8 GB chip_matrix and
# ~9 M-row aggregated_dt.

# Classify a feature table's columns by the key each one is constant over. The
# assembly needs this because the added (previously-unobserved) rows have to source
# each column from wherever that column actually varies.
#
# NOTE what the partition check below does and does NOT catch. `per_event` is a
# catch-all setdiff, so a column 05 adds later ALWAYS lands there and the check cannot
# fail -- it only catches duplicates and a genuine set mismatch. The real guard against
# a misfiled column is the per-event constancy assertion in build_full_event_rows(),
# which tests the classification against the data. If you add a column that is not
# constant within an event, classify it explicitly here.
classify_feature_columns <- function(col_names) {
  key <- intersect(c("ID", "IHEC", "uuid"), col_names)
  # per (ID, IHEC) -- the event-proximal window features, i.e. exactly the unfiltered
  # grid's payload. `_source` is per (sample, mark) observed-vs-imputed provenance.
  per_epigenome <- grep(
    "^(H3K[^;]+|DNAm|CpGs);(3|5)(up|down)$|_source$", col_names, value = TRUE
  )
  # per uuid -- RNA-library-level covariates and whole-spliceosome expression
  per_uuid <- c(
    intersect(c("protocol", "ontology", "project", "qc_flag_count"), col_names),
    grep("^spliceosome_", col_names, value = TRUE)
  )
  # per (uuid, gene_id); gene_id is constant within an event, so per uuid per event
  per_gene <- intersect(
    c("gene_expression_vst", "gene_expression_getmm"), col_names
  )
  # per (ID, uuid) -- the only genuinely sample-by-event quantities
  per_event_sample <- intersect(
    c(
      "PSI", "IJC", "SJC",
      "rbp_score_sum", "rbp_score_mean", "rbp_score_max", "rbp_n"
    ),
    col_names
  )
  # everything else is constant within an event (geometry, annotation, splice-site /
  # Pangolin / GC scores) and is recycled from the event's observed rows
  per_event <- base::setdiff(
    col_names,
    c(key, per_epigenome, per_uuid, per_gene, per_event_sample)
  )
  out <- list(
    key = key, per_epigenome = per_epigenome, per_uuid = per_uuid,
    per_gene = per_gene, per_event_sample = per_event_sample,
    per_event = per_event
  )
  flat <- unlist(out, use.names = FALSE)
  if (anyDuplicated(flat) || !setequal(flat, col_names)) {
    stop(
      "classify_feature_columns(): groups are not an exact partition of the ",
      length(col_names), " columns (unclassified: ",
      paste(base::setdiff(col_names, flat), collapse = ", "),
      "; duplicated: ", paste(flat[duplicated(flat)], collapse = ", "), ")"
    )
  }
  out
}


# Coerce every shared column of `add` to `template`'s type. The added rows draw from
# four sources that each round-trip through fread independently, so a column can
# arrive as character where aggregated_dt holds a factor, or double where it holds an
# integer. rbindlist would paper over some of that and mis-code the rest (two factors
# with different level sets are the dangerous case: the integer codes mean different
# things). A value that does not survive the coercion -- a factor level aggregated_dt
# has never seen -- is an ERROR, not an NA.
align_types <- function(add, template) {
  for (nm in intersect(names(add), names(template))) {
    tmpl <- template[[nm]]
    cur <- add[[nm]]
    if (is.factor(tmpl)) {
      if (!is.factor(cur) || !identical(levels(cur), levels(tmpl))) {
        chr <- as.character(cur)
        new <- factor(chr, levels = levels(tmpl))
        lost <- which(is.na(new) & !is.na(chr))
        if (length(lost)) {
          stop(
            "Column '", nm, "' has ", length(lost),
            " value(s) outside the template's factor levels (e.g. '",
            chr[lost[1L]], "') -- source disagrees with aggregated_dt"
          )
        }
        data.table::set(add, j = nm, value = new)
      }
    } else if (is.integer(tmpl) && !is.integer(cur)) {
      if (is.numeric(cur) && !isTRUE(all.equal(cur, round(cur)))) {
        stop("Column '", nm, "' is integer in the template but non-integral here")
      }
      data.table::set(add, j = nm, value = as.integer(cur))
    } else if (is.numeric(tmpl) && !is.numeric(cur)) {
      data.table::set(add, j = nm, value = as.numeric(cur))
    }
  }
  add[]
}


# Expand ONE event's observed rows to the full cohort sample set.
#
#   obs        — aggregated_dt[ID == id], i.e. the rows where THIS event's PSI was
#                observed. Returned VERBATIM as the first rows of the result.
#   all_uuids  — character vector of every cohort uuid for this transcript_filter
#   sample_cov — one row per uuid: uuid, IHEC, and the per-uuid columns
#   cols       — classify_feature_columns(names(obs))
#   grid       — unfiltered (ID, IHEC) grid, keyed, IHEC factor-aligned to obs$IHEC
#   rbp_score  — (ID, uuid) RBP aggregates, keyed on ID first
#   gene_expr  — (gene_id, uuid) expression, keyed on gene_id first
#
# The observed rows are copied, never rebuilt, so a bug in the assembly can degrade a
# NULL rotation but structurally cannot perturb the focal statistic -- 09s-ridge-
# screen.R re-applies `!is.na(PSI)` at load, so the focal fit reads only these rows.
#
# Every lookup is a match() on character, not a data.table join: the sources' key
# columns are variously factor (aggregated_dt: stringsAsFactors = TRUE) and character
# (everything else), and an implicit factor-to-character join is exactly the coercion
# that fails quietly rather than loudly. At <=415 rows per side (after a keyed binary
# search for the slice) match() costs nothing.
build_full_event_rows <- function(obs, id, all_uuids, sample_cov, cols,
                                  grid, rbp_score, gene_expr) {
  add_uuids <- base::setdiff(all_uuids, as.character(obs$uuid))
  if (length(add_uuids) == 0L) {
    return(obs)
  }
  add <- data.table::copy(sample_cov[
    match(add_uuids, as.character(sample_cov$uuid)),
    c("uuid", "IHEC", cols$per_uuid),
    with = FALSE
  ])
  if (anyNA(add$uuid)) {
    stop("build_full_event_rows(): ", sum(is.na(add$uuid)),
         " uuid(s) absent from sample_cov for event ", id)
  }
  add[, ID := id]
  ihec_chr <- as.character(add$IHEC)
  uuid_chr <- as.character(add$uuid)

  # per-event constants: recycled from row 1 of the event's observed rows.
  #
  # VERIFY the constancy rather than trusting the classification. cols$per_event is a
  # catch-all setdiff in classify_feature_columns(), so a column 05 adds later lands
  # here by DEFAULT and no partition check can fail -- an earlier version of this file
  # claimed such a column would error, which was wrong. A genuinely per-sample column
  # misfiled here would be silently recycled, giving every added row the first observed
  # sample's value. Checking against this event's own rows costs ~30 columns x a few
  # hundred rows and turns that silent corruption into a loud failure. NA counts as a
  # level, so a column that is NA for some samples of one event also trips this -- which
  # is correct, since an event property cannot be present for only some of its samples.
  for (nm in cols$per_event) {
    v <- obs[[nm]]
    if (length(v) > 1L && data.table::uniqueN(v) > 1L) {
      stop(
        "Column '", nm, "' is classified per-event but varies WITHIN event ", id,
        " (", data.table::uniqueN(v), " distinct values over ", length(v),
        " observed rows). Classify it in classify_feature_columns() -- recycling it ",
        "would give every added row the first sample's value."
      )
    }
    data.table::set(add, j = nm, value = v[1L])
  }
  # PSI / IJC / SJC: unobserved on these rows -- which is the entire point. Indexing
  # by NA_integer_ yields a length-1 NA of the column's own type, preserving it.
  # IJC/SJC stay NA rather than being looked up: they are blocked from x_cols (PSI is
  # computed from them), and a sample with no PSI reading has no junction counts to
  # report either, so NA is the honest value.
  for (nm in cols$per_event_sample) {
    data.table::set(add, j = nm, value = obs[[nm]][NA_integer_])
  }
  # per (ID, uuid): RBP aggregates
  rbp_cols <- intersect(cols$per_event_sample, names(rbp_score))
  if (length(rbp_cols)) {
    rs <- rbp_score[.(id), nomatch = NULL]
    if (nrow(rs)) {
      i_rs <- match(uuid_chr, as.character(rs$uuid))
      for (nm in rbp_cols) data.table::set(add, j = nm, value = rs[[nm]][i_rs])
    }
  }
  # per (uuid, gene_id): expression for THIS event's gene. Columns are created
  # UNCONDITIONALLY (typed NA) before the fill -- if they were only created inside the
  # `nrow(ge)` branch, an event whose gene is absent from gene_expr would return a
  # table missing these columns entirely, and the caller's setcolorder/rbindlist would
  # then fail on that event rather than yielding NA for it.
  for (nm in cols$per_gene) {
    data.table::set(add, j = nm, value = obs[[nm]][NA_integer_])
  }
  if (length(cols$per_gene)) {
    # Strip the Ensembl version on the focal side too, not only when building
    # gene_expr. aggregated_dt's gene_id is bare only because 05 strips it (05:586)
    # AFTER joining expression on the versioned id (05:496) -- i.e. correctness here
    # depends on where that strip sits relative to the file write. Stripping both
    # sides makes the lookup work whichever format arrives, and is a no-op today.
    ge <- gene_expr[
      .(sub("\\.\\d+$", "", as.character(obs[["gene_id"]][1L]))),
      nomatch = NULL
    ]
    if (nrow(ge)) {
      i_ge <- match(uuid_chr, as.character(ge$uuid))
      for (nm in cols$per_gene) {
        data.table::set(add, j = nm, value = ge[[nm]][i_ge])
      }
    }
  }
  # per (ID, IHEC): event-proximal epigenetic features from the unfiltered grid
  gr <- grid[.(id), nomatch = NULL]
  if (!nrow(gr)) {
    stop("build_full_event_rows(): event ", id, " has no rows in the grid")
  }
  i_gr <- match(ihec_chr, as.character(gr$IHEC))
  for (nm in cols$per_epigenome) {
    data.table::set(add, j = nm, value = gr[[nm]][i_gr])
  }

  add <- align_types(add, obs)
  data.table::setcolorder(add, names(obs))
  # fill = FALSE: a column present in one table and not the other must fail loudly
  # here. (rbindlist(fill = TRUE) silently creating and filling a column is a
  # documented trap in this codebase.)
  data.table::rbindlist(list(obs, add), use.names = TRUE, fill = FALSE)
}


# ===========================================================================
# 3. Tier-1 closed-form ridge screen statistic
# ===========================================================================
# Lin's concordance correlation coefficient (population / ÷n moments), matching
# yardstick::ccc(bias = TRUE).
# Strictly increasing squash of an unbounded R2 onto (-1, 1). Used for REPORTING and for
# effect sizes only -- never for p_emp, which is computed from the raw value (and would
# be numerically identical either way, since the map is monotone).
.squash_r2 <- function(x) x / (1 + abs(x))


.ccc <- function(x, y) {
  # Fail together with R2, do NOT silently compute on the survivors. Previously this
  # opened with `ok <- is.finite(x) & is.finite(y); x <- x[ok]` so that a fold whose
  # confound design went rank-deficient (NA out-of-fold residuals) produced an NA
  # screen_R2 -- which has no `na.rm` -- alongside a FINITE screen_CCC computed on an
  # unknown subset of rows. A finite number computed on an unadvertised subset is worse
  # than NA, because it looks valid. Measured 1,137 rows / 379 events with exactly that
  # signature before the fix.
  if (length(x) != length(y)) {
    return(NA_real_)
  }
  if (anyNA(x) || anyNA(y) || !all(is.finite(x)) || !all(is.finite(y))) {
    return(NA_real_)
  }
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
#   3 = 2026-08-01: confound design built PER FOLD (.build_fold_z) -- NA imputation
#       uses the TRAIN mean only and constancy is decided on TRAIN rows only, which
#       closed a smaller residual leak of the same shape as v1's. `prep$Z` is gone
#       (replaced by fold_z/fold_zte). Added the monotone squash R2_bounded and the
#       bounded effect sizes; .ccc() now fails together with R2 instead of silently
#       computing on the finite survivors.
#   4 = 2026-08-04: the matched-control NULL changed -- control features are now
#       PSI-independent (09-1 builds per-event feature tables over the FULL cohort
#       row set, not just rows where that event's own PSI was observed), so a
#       control no longer has to cover every focal sample to be eligible. Different
#       control set => different null distribution => different p_emp/q/hits. The
#       FOCAL statistic is unchanged by this: 09s-ridge-screen.R still applies
#       `!is.na(PSI)` at load, so the real fit sees exactly the rows it saw at v3.
#       Bumped anyway, because the row this stamp guards carries the null summary
#       and the p/q derived from it, not only the focal numbers.
#
# A hash of the two function bodies was considered instead (no discipline needed)
# but rejected: it would also invalidate on a comment-only edit, and a false
# invalidation here costs a multi-hour 34k-job rerun.
SCREEN_STAT_VERSION <- 4L


# FEATURE_TABLE_VERSION — provenance stamp for the per-event feature TABLES, the
# exact analogue of SCREEN_STAT_VERSION one layer upstream.
#
# BUMP THIS whenever 09-1-ml-local.R's Phase-1 build changes the CONTENTS of
# feature_table_<id>.csv.gz (its row set, its column set, or any value in it).
# It is baked into the directory name via feature_table_dir_for() below, so a bump
# makes the pipeline write to a fresh directory and the old tables are neither read
# nor overwritten.
#
# Why a version marker and not `rm -rf`: 09-1's Phase-1 loop opens with
# `if (file.exists(feature_table_file)) return(invisible(NULL))`, and there are
# ~34k tables on disk. WITHOUT a bump, a build with changed logic skips every
# single event and exits 0 -- the change is silently a no-op, exactly the failure
# mode SCREEN_STAT_VERSION exists to prevent one stage later. Deleting instead
# would work but throws away ~72 GB with no provenance trail and no way to diff old
# against new while verifying, so the two versions are kept side by side and the
# stale one removed by hand once the new screen is trusted.
#
#   1 = rows = only the samples where THIS event's PSI was observed
#       (`aggregated_dt_filtered[ID == id]`, ~391 of 415 for biotype_filtered).
#   2 = 2026-08-04: rows = the FULL cohort sample set (all 415), PSI/IJC/SJC NA
#       where unobserved. Epigenetic/expression features are filled for every
#       sample, which is what lets a control event be scored on the focal event's
#       samples regardless of where its own PSI happens to be missing.
FEATURE_TABLE_VERSION <- 2L


# Canonical per-event feature-table directory. Single definition, so 09-1 (writer),
# any verification script, and anything reading cfg$feature_table_dir cannot drift
# apart on the version suffix. v1 keeps the historical unsuffixed name so the 72 GB
# of existing tables stay addressable without being moved.
feature_table_dir_for <- function(tf, version = FEATURE_TABLE_VERSION) {
  base <- sprintf("event_feature_tables_%s", tf)
  if (version > 1L) base <- sprintf("%s_v%d", base, version)
  file.path("processed_data", base)
}


# Resolve the Tier-1 rotation count for ONE event, given the cfg field and that
# event's Event Type.
#
# Why this is per-Event-Type: 09s-ridge-screen.R caps usage at
# `min(n_rotations, length(other_ids))`, so the eligible pool is only ever an upper
# bound -- raising eligibility without raising this changes nothing. But the two
# event types need wildly different values. RI has ~1,776 modelable events and a
# matched-control pool of 791-887 (measured), and needs nearly all of it to get its
# empirical-p floor 1/(R+1) below BH's `k*q/m` bound. SE has ~32,370 events and a
# pool of ~15,000, is already comfortably inside its bound at 200, and would cost
# ~10 days if raised to match RI. A single scalar cannot serve both.
#
# Accepts either form:
#   * a length-1 unnamed value -> used for every event (backward compatible with
#     every cfg written before 2026-08-04)
#   * a named vector/list keyed by Event Type, e.g. c(RI = 818L, SE = 200L)
#
# An unrecognised Event Type is a hard ERROR, not a fallback to some default. A
# silent fallback here would reintroduce the precise bug this function exists to
# fix (a new event type quietly screening at 200 rotations and failing FDR for
# arithmetic reasons); the fix is to add the event type to the cfg, the same
# convention splicing_ml's FEATURE_GROUPS uses for a new feature column.
resolve_screen_rotations <- function(spec, event_type, default = 200L) {
  if (is.null(spec) || length(spec) == 0L) {
    return(as.integer(default))
  }
  nms <- names(spec)
  if (is.null(nms) || !any(nzchar(nms))) {
    if (length(spec) != 1L) {
      stop(
        "cfg$screen_rotations is unnamed but has length ", length(spec),
        " -- give it one value per Event Type, e.g. c(RI = 818L, SE = 200L)"
      )
    }
    return(as.integer(spec[[1L]]))
  }
  et <- as.character(event_type)
  if (!length(et) || is.na(et) || !nzchar(et)) {
    stop("resolve_screen_rotations(): event_type is missing/empty")
  }
  if (!et %in% nms) {
    stop(
      "cfg$screen_rotations has no entry for Event Type '", et,
      "' (has: ", paste(nms, collapse = ", "),
      "). Add it rather than defaulting -- a wrong rotation count silently ",
      "decides whether this event type can clear FDR at all."
    )
  }
  as.integer(spec[[et]])
}


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
# .build_fold_z(): construct the confound design for ONE fold, deciding which columns
# survive and what NAs are imputed with using that fold's TRAIN rows only.
#
# Why this cannot be a single shared Z any more: per-fold NA imputation changes the
# VALUES in the design, not merely which columns are kept, so the train and test blocks
# have to be materialised per fold. They are tiny (n <= ~415 rows, a handful of confound
# columns), so the memory cost is irrelevant.
#
# Two leakage-relevant decisions, both train-only:
#   * numeric NAs are filled with the TRAIN mean (was: the mean over ALL n rows, so a
#     held-out group's own values entered the constant used to adjust it -- the same
#     leak shape as v1's global standardisation, smaller in magnitude)
#   * a column constant in THIS fold's train rows is dropped (was: constant over all n
#     rows), which is what caused the rank-deficient QR -> NA out-of-fold residuals
.build_fold_z <- function(confound_df, train_idx, test_idx) {
  ntr <- length(train_idx)
  nte <- length(test_idx)
  Ztr <- matrix(1.0, nrow = ntr, ncol = 1L, dimnames = list(NULL, "(Intercept)"))
  Zte <- matrix(1.0, nrow = nte, ncol = 1L, dimnames = list(NULL, "(Intercept)"))
  notes <- character(0)
  add <- function(vtr, vte, nm) {
    Ztr <<- cbind(Ztr, vtr)
    Zte <<- cbind(Zte, vte)
    colnames(Ztr)[ncol(Ztr)] <<- nm
    colnames(Zte)[ncol(Zte)] <<- nm
  }
  if (!is.null(confound_df) && ncol(confound_df) > 0L) {
    for (nm in names(confound_df)) {
      v <- confound_df[[nm]]
      if (is.numeric(v)) {
        vtr <- v[train_idx]
        vte <- v[test_idx]
        if (all(is.na(vtr))) {
          notes <- c(notes, paste0("confound_all_na_in_train:", nm))
          next
        }
        mu <- mean(vtr, na.rm = TRUE) # TRAIN mean only
        vtr[is.na(vtr)] <- mu
        vte[is.na(vte)] <- mu
        if (stats::sd(vtr) <= 1e-12) {
          # constant in this fold's train rows -> absorbed by the intercept; keeping it
          # would alias the QR and NA out the whole event
          next
        }
        add(vtr, vte, nm)
      } else {
        ch <- ifelse(is.na(v), "NA", as.character(v))
        lv_tr <- sort(unique(ch[train_idx]))
        if (length(lv_tr) < 2L) next # constant in train -> absorbed by the intercept
        # dummies over TRAIN levels, baseline = first. A test row whose level never
        # appears in train gets all-zero dummies, i.e. is folded into the baseline --
        # the only available choice, since no coefficient can exist for a level the fit
        # never saw. Noted rather than silently absorbed.
        for (lvl in lv_tr[-1L]) {
          add(
            as.numeric(ch[train_idx] == lvl),
            as.numeric(ch[test_idx] == lvl),
            paste0(nm, "_", lvl)
          )
        }
        unseen <- base::setdiff(unique(ch[test_idx]), lv_tr)
        if (length(unseen)) {
          notes <- c(notes, paste0(
            "confound_level_unseen_in_train:", nm, "=",
            paste(unseen, collapse = "/")
          ))
        }
      }
    }
  }
  list(Ztr = Ztr, Zte = Zte, notes = notes)
}

prep_event <- function(y, confound_df, groups) {
  n <- length(y)
  group_idx <- split(seq_len(n), groups)
  fold_z <- vector("list", length(group_idx))
  names(fold_z) <- names(group_idx)
  fold_zte <- vector("list", length(group_idx))
  names(fold_zte) <- names(group_idx)
  ytil_train <- vector("list", length(group_idx))
  names(ytil_train) <- names(group_idx)
  ytil_oof <- numeric(n)
  notes <- character(0)

  for (g in names(group_idx)) {
    test_idx <- group_idx[[g]]
    train_idx <- base::setdiff(seq_len(n), test_idx)
    fz <- .build_fold_z(confound_df, train_idx, test_idx)
    notes <- c(notes, fz$notes)
    qrZtr <- qr(fz$Ztr)
    if (qrZtr$rank < ncol(fz$Ztr)) {
      # should be unreachable now that constancy is decided per fold, but a residual
      # collinearity (two confounds identical within this fold's train rows) would still
      # land here, and it must never again pass silently
      notes <- c(notes, "fold_rank_deficient")
    }
    b_y <- qr.coef(qrZtr, y[train_idx])
    b_y[!is.finite(b_y)] <- 0 # aliased columns contribute nothing rather than NA
    fold_z[[g]] <- qrZtr
    fold_zte[[g]] <- fz$Zte
    ytil_train[[g]] <- qr.resid(qrZtr, y[train_idx])
    ytil_oof[test_idx] <- y[test_idx] - fz$Zte %*% b_y
  }

  list(
    n = n, group_idx = group_idx, fold_z = fold_z, fold_zte = fold_zte,
    ytil_train = ytil_train, ytil_oof = ytil_oof, ss_tot = sum(ytil_oof^2),
    notes = unique(notes)
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
    R2 = NA_real_, R2_bounded = NA_real_, CCC = NA_real_, lambda = NA_real_,
    df = NA_real_, n_features = 0L,
    max_abs_oof = NA_real_, max_abs_z = NA_real_, frac_z_gt10 = NA_real_,
    lambda_min = NA_real_
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
  # B3 out-of-support diagnostics. The leave-one-ontology-supergroup-out design means a
  # held-out group can sit far outside the training support (measured: |z| up to 48,404
  # on real data, and >half the held-out rows beyond 10 train SDs in 37.5% of
  # event x feature-set combos). That is the finding the screen exists to expose, so it
  # is recorded rather than suppressed.
  z_max_folds <- numeric(0)
  z_gt10 <- numeric(0)
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

    # how far outside the training support does this held-out group actually sit?
    z_max_folds <- c(z_max_folds, max(abs(Xte)))
    z_gt10 <- c(z_gt10, mean(apply(abs(Xte), 1L, max) > 10))

    # residualise X on Z using this fold's TRAIN-fitted confound coefficients
    # (prep$fold_z[[g]] was fit on train rows only in prep_event) — applied to
    # test rows by explicit matrix projection, never refit on test. The test block
    # comes from prep$fold_zte[[g]], NOT a slice of a shared Z: the design is now
    # built per fold (train-only imputation + train-only constancy), so there is no
    # single Z to slice.
    qrZtr <- prep$fold_z[[g]]
    b_X <- qr.coef(qrZtr, Xtr)
    b_X[!is.finite(b_X)] <- 0 # defensive: aliased column contributes 0, never NA
    MXtr <- qr.resid(qrZtr, Xtr)
    MXte <- Xte - prep$fold_zte[[g]] %*% b_X

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

  r2 <- 1 - sum((prep$ytil_oof - oof)^2) / prep$ss_tot
  list(
    # RAW R2 is retained unchanged and is what p_emp / q are computed from. It is
    # unbounded (measured min -290,429 on real data) but p_emp is RANK-based, so the
    # magnitude is irrelevant to inference.
    R2 = r2,
    # Bounded companion for reporting, tables and effect sizes: R2/(1+|R2|) in (-1,1).
    # Strictly increasing, therefore p_emp computed from it is EXACTLY identical to
    # p_emp from raw R2 (f(null) >= f(real) <=> null >= real; no ties are created).
    # That is the property clamping lacks -- clamping is not injective, it collapses
    # values onto the bound, manufactures ties, and so genuinely moved p (Spearman
    # 0.70 vs raw, 20-23 hit flips across 678 real event x feature-set combos).
    # Bounded values also make `effect` usable again: a mean over (-1,1) cannot be
    # dominated by a single -290,000 null, which was the substance of TODO-2.
    R2_bounded = .squash_r2(r2),
    CCC = .ccc(prep$ytil_oof, oof),
    lambda = mean(lambdas),
    df = mean(dfs),
    n_features = p,
    max_abs_oof = max(abs(oof)),
    max_abs_z = if (length(z_max_folds)) max(z_max_folds) else NA_real_,
    frac_z_gt10 = if (length(z_gt10)) mean(z_gt10) else NA_real_,
    # min lambda over folds: GCV running to the bottom of its grid is the signature of a
    # fold that nearly interpolates its train rows and then extrapolates. Tracks blow-up
    # magnitude monotonically on real data (median 83.1 -> 2.08 -> 1.55 -> 1.04 across
    # |R2| strata), so it is kept as a diagnostic. NB raising lambda was TESTED as a fix
    # and REJECTED: a floor at 1e-2*trace(K)/n_tr binds in 81.8% of the worst stratum yet
    # only moves median R2 -8.18 -> -5.10 and min -290,429 -> -282,831. The blow-up is
    # driven by test-side kernel magnitude, not by lambda being too small.
    lambda_min = if (length(lambdas)) min(lambdas) else NA_real_
  )
}
