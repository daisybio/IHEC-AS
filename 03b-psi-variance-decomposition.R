#!/usr/bin/env Rscript
# PSI variance decomposition — DrEval (Bernett et al. 2026, Nat Commun) Table 1 in splicing form.
#
# QUESTION: before any epigenetic feature is considered, how much of the variation in PSI is
# explained by trivial structure alone — which event it is, which sample it is, which cell type?
#
# WHY IT MATTERS HERE. Tier 1 fits WITHIN an event, so the event mean is absorbed by the intercept:
# the within-event fraction is the entire variance budget Tier 1 could ever explain. splicing_ml
# pools ACROSS events, so it is free to spend the between-event fraction. This decomposition is what
# makes "0 hits" and "AUROC 0.957" provably compatible rather than merely arguable, and it supplies
# the denominator that makes a single event's R2 (e.g. ID 38354's 0.349) interpretable.
#
# Deliberately needs NO epigenetic data — inputs stop at stage 03 plus sample metadata. That is the
# point: the budget is a property of the response, established before any mark is touched.
#
# Fit-free, ~10M rows. Safe on the login node.
#
# STAGE 03b, and the position is load-bearing: every input stops at stage 03 plus sample metadata,
# so this runs BEFORE create_aggregated_dt and touches no epigenetic data at all. That is the
# analysis's own claim -- the variance budget is a property of the RESPONSE, established before any
# mark is considered -- so placing it downstream of the marks would undercut what it is for.
# Promoted from an undeclared prototype 2026-09-01: 12 of paper_numbers' keys come from here (the whole
# of leg A), and Nature requires the code behind a quoted number to be available at submission.

suppressPackageStartupMessages(library(data.table))
setDTthreads(4L)

tf <- Sys.getenv("TRANSCRIPT_FILTER", "biotype_filtered")
ONTO_COL <- "harmonized_sample_ontology_term_high_order_fig1" # the column 01/05/Tier-1 all use

cat(sprintf("=== PSI variance decomposition | transcript_filter = %s ===\n\n", tf))

# ---------------------------------------------------------------- inputs (stage 03 + metadata)
krm <- readRDS(sprintf("processed_data/keep_rows_manual_%s.rds", tf)) # INTEGER event IDs, not logical
stopifnot(is.integer(krm) || is.numeric(krm))
ann <- fread(sprintf("processed_data/event_annotations_dt_%s.csv.gz", tf),
  select = c("ID", "Event Type", "Variability"))
ann <- ann[ID %in% krm]
cat(sprintf("events: %s (keep_rows_manual)\n", format(nrow(ann), big.mark = ",")))
print(ann[, .N, by = .(`Event Type`, Variability)][order(`Event Type`, Variability)])

pl <- fread(sprintf("processed_data/psi_long_dt_%s.csv.gz", tf), select = c("ID", "uuid", "psi"))
pl <- pl[ID %in% krm][!is.na(psi)]
cat(sprintf("\nobserved (event, sample) rows: %s over %d events, %d samples\n",
  format(nrow(pl), big.mark = ","), uniqueN(pl$ID), uniqueN(pl$uuid)))

# uuid -> EpiRR -> ontology
ft <- unique(fread("processed_data/file_table.csv.gz",
  select = c("uuid", "epirr_id_without_version"))[!is.na(uuid)])
meta <- fread("data/IHEC_sample_metadata_harmonization.v1.4_extended.csv",
  select = c("EpiRR", ONTO_COL))
setnames(meta, c("EpiRR", ONTO_COL), c("epirr_v", "ontology"))
meta[, epirr_id_without_version := sub("\\.\\d+$", "", as.character(epirr_v))]
map <- merge(ft, unique(meta[, .(epirr_id_without_version, ontology)]),
  by = "epirr_id_without_version", all.x = TRUE)
pl[map, on = "uuid", ontology := i.ontology]
cat(sprintf("samples with an ontology label: %d of %d (%d distinct labels)\n",
  uniqueN(pl[!is.na(ontology), uuid]), uniqueN(pl$uuid), uniqueN(pl$ontology)))

# NB `i.` cannot prefix a backticked/spaced column name (`i.`Event Type`` is a parse error) --
# alias it in the join table's j first. Documented in CLAUDE.md; this is its third occurrence here.
pl[ann[, .(ID, k_et = `Event Type`, k_vb = Variability)], on = "ID",
  `:=`(event_type = i.k_et, variability = i.k_vb)]

# ------------------------------------------------------------------------- decomposition helpers
# R2 = 1 - SSE/SST for a prediction vector.
.r2 <- function(y, yhat, sst) 1 - sum((y - yhat)^2) / sst

# Leave-one-out prediction for a group mean, closed form: (n*m - y)/(n-1).
# Singleton groups have no LOO prediction; fall back to the stratum grand mean and count them, since
# a silent fallback would quietly inflate the LOO R2 toward the in-sample one.
.loo_group_mean <- function(y, g, grand) {
  d <- data.table(y = y, g = g)
  d[, `:=`(n = .N, m = mean(y)), by = g]
  out <- fifelse(d$n > 1L, (d$n * d$m - d$y) / (d$n - 1L), grand)
  list(pred = out, n_singleton = sum(d$n == 1L))
}

decompose <- function(d, label) {
  y <- d$psi
  grand <- mean(y)
  sst <- sum((y - grand)^2)
  if (!is.finite(sst) || sst <= 0) return(NULL)

  ev <- .loo_group_mean(y, d$ID, grand)
  sa <- .loo_group_mean(y, d$uuid, grand)
  on <- .loo_group_mean(y, d$ontology, grand)

  mu_ev <- d[, mean(psi), by = ID][d, on = "ID", x.V1]
  mu_sa <- d[, mean(psi), by = uuid][d, on = "uuid", x.V1]
  mu_on <- d[, mean(psi), by = ontology][d, on = "ontology", x.V1]

  data.table(
    stratum = label,
    n_rows = length(y), n_events = uniqueN(d$ID), n_samples = uniqueN(d$uuid),
    var_psi = round(sst / (length(y) - 1L), 5),
    # in-sample (descriptive budget)
    r2_global = 0,
    r2_event = round(.r2(y, mu_ev, sst), 4),
    r2_sample = round(.r2(y, mu_sa, sst), 4),
    r2_ontology = round(.r2(y, mu_on, sst), 4),
    r2_additive = round(.r2(y, mu_ev + mu_sa - grand, sst), 4),
    # leave-one-out (honest, over-fit removed)
    r2_event_loo = round(.r2(y, ev$pred, sst), 4),
    r2_sample_loo = round(.r2(y, sa$pred, sst), 4),
    r2_ontology_loo = round(.r2(y, on$pred, sst), 4),
    n_singleton_events = ev$n_singleton
  )
}

# ------------------------------------------------------------------------------------ strata
res <- list()
for (et in c("SE", "RI")) {
  for (vb in c("All", "Low", "High")) {
    d <- if (vb == "All") pl[event_type == et] else pl[event_type == et & variability == vb]
    if (!nrow(d)) next
    res[[paste(et, vb)]] <- decompose(copy(d), sprintf("%s / %s", et, vb))
  }
}
out <- rbindlist(res, fill = TRUE)

cat("\n\n=== IN-SAMPLE (descriptive budget) ===\n")
print(out[, .(stratum, n_rows, n_events, n_samples, var_psi,
  r2_event, r2_sample, r2_ontology, r2_additive)])

cat("\n=== LEAVE-ONE-OUT (over-fit removed; singleton groups fall back to the grand mean) ===\n")
print(out[, .(stratum, r2_event_loo, r2_sample_loo, r2_ontology_loo, n_singleton_events)])

cat("\n=== WHAT IS LEFT FOR EPIGENETICS (within-event budget) ===\n")
cat("1 - r2_event_loo is the fraction of PSI variance that survives removing event identity.\n")
cat("Tier 1 operates entirely inside that fraction.\n\n")
print(out[, .(stratum,
  between_event_pct = round(100 * r2_event_loo, 1),
  within_event_pct = round(100 * (1 - r2_event_loo), 1),
  ontology_pct_of_total = round(100 * r2_ontology_loo, 2))])

# Snakemake owns the path (VARIANCE_OUT is {output}); the sprintf default keeps the script
# runnable standalone. Same contract as 09f-floor-sample.R's FLOOR_OUT.
f <- Sys.getenv(
  "VARIANCE_OUT",
  sprintf("processed_data/psi_variance_decomposition_%s.csv", tf)
)
dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
fwrite(out, f)
cat(sprintf("\nwritten: %s\n", f))

cat("\nCAVEATS\n")
cat("- Unbalanced design (mean 228.9 observed samples of 415 per event), so the components are NOT\n")
cat("  orthogonal: r2_event + r2_sample != r2_additive. Read r2_additive on its own row.\n")
cat("- r2_*_loo is the honest number; the in-sample event mean carries one parameter per event.\n")
cat("- Computed on raw PSI, which has heavy boundary mass (~40% of values exactly at 0 or 1).\n")
cat("  Variance is still variance, but the budget is not a Gaussian one.\n")
