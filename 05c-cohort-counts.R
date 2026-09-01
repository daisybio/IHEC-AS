#!/usr/bin/env Rscript
# Cohort counts -- the MODELLED cohort, as a declared artifact.
#
# WHY THIS IS ITS OWN STAGE. `11-paper-figures.Rmd` used to derive these inline and cache them to a
# path that was not one of its declared outputs, so Snakemake neither tracked nor rebuilt the cache:
# a stale file would have been read silently, and the cohort counts are the numbers most often quoted
# wrongly in this project. It was also the ONE place `11` touched `aggregated_dt` at all, against that
# file's own rule that it must stay fast enough to re-run on every figure tweak.
#
# TWO COHORTS, AND CONFLATING THEM IS THE RECURRING ERROR:
#   modelled   379 epirrs / 415 uuids  (188 mRNA-Seq + 227 total-RNA-Seq, 36 epirrs with both)
#   QC scope   405 epirrs / 441 uuids
# a clean superset -- 26 QC rows are never modelled, 0 modelled uuids lack QC. The 26 are the
# AMED-CREST gap and ALL were mRNA-Seq, which is why total-RNA-Seq is unchanged at 227 while mRNA-Seq
# falls 214 -> 188. Quote the MODELLED counts for anything about what was analysed.
#
# Both scopes are emitted as rows of one table, deliberately: split across two artifacts, the
# 405-vs-379 drift recurs.

suppressPackageStartupMessages(library(data.table))
setDTthreads(2L)

tf <- Sys.getenv("TRANSCRIPT_FILTER", "biotype_filtered")
out <- Sys.getenv("COHORT_OUT", sprintf("processed_data/cohort_counts_%s.csv", tf))

# 3 columns only -- the whole point is not to pay for the 9 M-row x 98-column table.
a <- fread(
  sprintf("processed_data/aggregated_dt_filtered_%s.csv.gz", tf),
  select = c("uuid", "IHEC", "protocol")
)
u <- unique(a, by = "uuid")

modelled <- data.table(
  scope          = "modelled",
  uuids          = uniqueN(a$uuid),
  epirrs         = uniqueN(a$IHEC),
  mrna           = u[protocol == "mRNA-Seq", .N],
  total_rna      = u[protocol == "total-RNA-Seq", .N],
  both_protocols = unique(a[, .(IHEC, protocol)])[, .N, by = IHEC][N > 1, .N],
  source         = sprintf("aggregated_dt_filtered_%s.csv.gz", tf)
)
rm(a, u); gc()

# QC scope, from 01's own output. Emitted beside the modelled row rather than in a separate file so a
# reader cannot pick up one without seeing the other.
qc_file <- "processed_data/qc_summary.csv"
qc <- if (file.exists(qc_file)) {
  q <- fread(qc_file)
  data.table(
    scope          = "qc_metadata",
    uuids          = if ("uuid" %in% names(q)) uniqueN(q$uuid) else NA_integer_,
    epirrs         = uniqueN(q$epirr_id_without_version),
    mrna           = NA_integer_,
    total_rna      = NA_integer_,
    both_protocols = NA_integer_,
    source         = qc_file
  )
} else NULL

res <- rbindlist(list(modelled, qc), fill = TRUE)
print(res)

cat("\n=== the distinction that keeps being lost ===\n")
cat(sprintf("modelled: %d uuids / %d epirrs (%.1f%% total-RNA-Seq)\n",
  modelled$uuids, modelled$epirrs, 100 * modelled$total_rna / modelled$uuids))
if (!is.null(qc)) {
  cat(sprintf("QC scope: %s uuids / %d epirrs -- a SUPERSET; never quote it as the analysed cohort\n",
    format(qc$uuids), qc$epirrs))
}

dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
fwrite(res, out)
cat(sprintf("\nwritten: %s\n", out))
