#!/usr/bin/env Rscript
# Supplementary QC table for the -omics inputs -- the per-sample sheet the manuscript's supplement
# needs, assembled from artifacts that already exist.
#
# WHY THIS IS A SEPARATE STAGE RATHER THAN A CHUNK IN `01-gather-data.Rmd`.
# The obvious place to build this is right after `01`'s line 642, where `merged_qc` is still in scope.
# It is the wrong place. `01-gather-data.Rmd` is a declared input of `rule gather_data`, whose outputs
# (`file_table.csv.gz`, `qc_summary.csv`, `qc_flag_covariates.csv`) get new mtimes on any re-run, and
# that cascades through `create_aggregated_dt`, the 117 GB / 3h37 feature-table build and all 34,146
# screen jobs. Everything needed here is already ON DISK as `01` output, so a downstream reader costs
# nothing and triggers nothing. Same reasoning, same conclusion as `revision/deferred-cosmetic-edits.md`
# reached for the `qc/...txt` header: fix the DELIVERABLE, do not re-run the pipeline for it.
#
# WHY IT IS `05d` AND NOT `01b`. It joins the MODELLED cohort roster, which only exists after `05c`,
# so its true DAG position is after `05c` and the number names that position. It reads no epigenetic
# data and does no modelling.
#
# THE COHORT DISTINCTION IS THE POINT OF THE `modelled` COLUMN.
# `qc_summary.csv` covers 441 RNA-seq uuids / 405 epirrs. Only 415 / 379 were modelled. The QC scope is
# a clean superset -- 26 QC rows are never modelled, 0 modelled uuids lack QC -- and the 26 are the
# AMED-CREST gap, all of them mRNA-Seq. `qc/...txt` prints the header "441 RNA-Seq uuids / 405 epirrs",
# which is correct for the QC pipeline's own scope and WRONG as a description of the analysis; that
# header is deliberately not being edited (it would cascade). So this table has to carry the
# distinction itself, per row, or a reader will take 441 for the analysed cohort.
#
# POLICY A. No sample is dropped for QC. Every flagged sample is retained and `qc_flag_count` (0-8)
# enters the models as a covariate. The `*_potentially_problematic` columns are therefore descriptive,
# not a filter -- that has to be stated in the legend or the table reads as a rejection log.

suppressPackageStartupMessages(library(data.table))
setDTthreads(2L)

tf         <- Sys.getenv("TRANSCRIPT_FILTER", "biotype_filtered")
qc_file    <- Sys.getenv("QC_SUMMARY", "processed_data/qc_summary.csv")
ft_file    <- Sys.getenv("FILE_TABLE", "processed_data/file_table.csv.gz")
roster_file <- Sys.getenv("COHORT_UUIDS", sprintf("processed_data/cohort_uuids_%s.csv", tf))
meta_file  <- Sys.getenv("HARMONIZATION_CSV",
                         "data/IHEC_sample_metadata_harmonization.v1.4_extended.csv")
out        <- Sys.getenv("QC_TABLE_OUT", sprintf("qc/supplementary_qc_table_%s.csv", tf))
thr_out    <- Sys.getenv("QC_THRESHOLDS_OUT",
                         sprintf("qc/supplementary_qc_thresholds_%s.csv", tf))

MARKS <- c("H3K27ac", "H3K4me1", "H3K4me3", "H3K27me3", "H3K36me3", "H3K9me3")
NARROW_MARKS <- c("H3K27ac", "H3K4me3", "H3K4me1")

# `qc_summary.csv` carries 553 columns -- every flagstat/dup/pbc/xcor/jsd field for 6 marks plus their
# input controls. A supplement needs the metrics the thresholds are actually defined on, not all of
# them. These 11 per mark are exactly the ones `01` tests against `chip_thr`, plus library depth.
CHIP_METRICS <- c(
  "potentially_problematic",
  "flagstat_qc:total",            # depth, reported for context; not thresholded
  "flagstat_qc:mapped_pct",
  "dup_qc:dupes_pct",
  "pbc_qc:NRF", "pbc_qc:PBC1", "pbc_qc:PBC2",
  "xcor_score:NSC", "xcor_score:RSC",   # reported for context; not thresholded
  "frip_macs2_qc:rep1.FRiP",
  "jsd_qc:jsd"                          # reported for context; not thresholded
)
RNA_COLS  <- c("rna_potentially_problematic", "total", "fraction_mapped",
               "fraction_intergenic", "fraction_rrna", "fraction_duplicates")
WGBS_COLS <- c("wgbs_potentially_problematic", "Mean_CpG_Coverage",
               "BS_conversion_Rate", "GC_cov_Correlation")

qc <- fread(qc_file)
stopifnot(all(c("epirr_id_without_version", "uuid", "qc_flag_count") %in% names(qc)))

chip_cols <- unlist(lapply(MARKS, function(m) paste(m, CHIP_METRICS, sep = "_")))
missing <- setdiff(c(RNA_COLS, WGBS_COLS, chip_cols), names(qc))
# A hard error, not a warning. A silently absent metric would publish as a column of NA, which reads
# as "not measured" rather than "this script asked for the wrong name".
if (length(missing)) {
  stop("qc_summary.csv is missing ", length(missing), " requested column(s): ",
       paste(utils::head(missing, 10), collapse = ", "))
}

tab <- qc[, c("epirr_id_without_version", "uuid", RNA_COLS, WGBS_COLS, chip_cols,
              "qc_flag_count"), with = FALSE]
setnames(tab, "epirr_id_without_version", "epirr")

# --- protocol, from 01's own file table -------------------------------------------------------
# `experiment_type` is what `05` itself renames to `protocol` (05-create-aggregated-dt.Rmd:397), so
# this is the pipeline's own definition rather than a parallel one.
ft <- fread(ft_file, select = c("uuid", "experiment_type"))
ft <- unique(ft[experiment_type %in% c("mRNA-Seq", "total-RNA-Seq")], by = "uuid")
tab[ft, on = "uuid", protocol := i.experiment_type]

# --- ontology, from the harmonization table ---------------------------------------------------
# `harmonized_sample_ontology_term_high_order_fig1` is the column `05` turns into the `ontology`
# factor (05-create-aggregated-dt.Rmd:402-403). Joined on EpiRR, so it covers the QC-only rows too.
ONTOLOGY_COL <- "harmonized_sample_ontology_term_high_order_fig1"
meta <- fread(meta_file, select = c("EpiRR", ONTOLOGY_COL))
setnames(meta, c("epirr", "ontology"))
# The harmonization table's EpiRR carries a VERSION suffix (`IHECRE00001032.6`); every id on the QC
# side is `epirr_id_without_version` (`IHECRE00000001`). Joining them raw silently matches nothing --
# it produced 441 NA ontologies, i.e. an entirely empty column, not an error. Strip the version, then
# assert the join landed rather than trusting it a second time.
meta[, epirr := sub("\\.[0-9]+$", "", epirr)]
meta <- unique(meta, by = "epirr")
tab[meta, on = "epirr", ontology := i.ontology]

# --- modelled flag ----------------------------------------------------------------------------
# From 05c's roster, never re-derived here: one definition of "the 415".
roster <- fread(roster_file, select = "uuid")
tab[, modelled := uuid %in% roster$uuid]

setcolorder(tab, c("epirr", "uuid", "ontology", "protocol", "modelled", "qc_flag_count"))
setorder(tab, -modelled, epirr, uuid)

# --- checks that would otherwise surface as a wrong supplementary table ------------------------
n_modelled <- tab[modelled == TRUE, .N]
# The superset relation is the claim the `modelled` column makes; assert it rather than trust it.
if (nrow(roster) != n_modelled) {
  stop("roster has ", nrow(roster), " modelled uuids but only ", n_modelled,
       " of them appear in qc_summary.csv -- the QC scope is NOT a superset of the modelled cohort")
}
if (anyNA(tab$protocol)) stop(sum(is.na(tab$protocol)), " uuid(s) have no protocol in ", ft_file)
# A hard stop, not a warning. The first version of this join matched zero rows and only warned, which
# would have shipped a supplementary table with an empty ontology column.
if (anyNA(tab$ontology)) {
  stop(sum(is.na(tab$ontology)), " of ", nrow(tab), " uuid(s) have no ontology in ", meta_file)
}

# --- thresholds, as their own sheet -----------------------------------------------------------
# Mirrors 01-gather-data.Rmd's `rna_thr` (:240), `chip_thr` (:362) and `wgbs_thr` (:533). Duplicated
# deliberately: the alternative is editing `01` to export them, which is the cascade this stage
# exists to avoid. If those lists change, this table is the thing to update -- it is a legend, not an
# input to any computation, so a drift here misdescribes the supplement but cannot corrupt a result.
thr <- rbindlist(list(
  data.table(assay = "RNA-seq", metric = "total reads",            direction = "min", value = 50e6,
             note = "IHEC standard"),
  data.table(assay = "RNA-seq", metric = "fraction_mapped",        direction = "min", value = 0.70,
             note = "ENCODE >=70%"),
  data.table(assay = "RNA-seq", metric = "fraction_intergenic",    direction = "max", value = 0.15,
             note = ""),
  data.table(assay = "RNA-seq", metric = "fraction_rrna",          direction = "max", value = 0.25,
             note = "ENCODE total RNA <25%"),
  data.table(assay = "RNA-seq", metric = "fraction_duplicates",    direction = "max", value = 0.70,
             note = "ENCODE <80%; conservative"),
  data.table(assay = "ChIP-seq", metric = "flagstat mapped_pct",   direction = "min", value = 80,
             note = "ENCODE >80%"),
  data.table(assay = "ChIP-seq", metric = "dup_qc dupes_pct",      direction = "max", value = 0.80,
             note = "ENCODE <80%; stored as a fraction, not a percent"),
  data.table(assay = "ChIP-seq", metric = "pbc_qc NRF",            direction = "min", value = 0.5,
             note = "ENCODE minimum (compliant 0.8, ideal 0.9)"),
  data.table(assay = "ChIP-seq", metric = "pbc_qc PBC1",           direction = "min", value = 0.5,
             note = "ENCODE minimum (compliant 0.8)"),
  data.table(assay = "ChIP-seq", metric = "pbc_qc PBC2",           direction = "min", value = 1,
             note = "ENCODE minimum (compliant 3)"),
  data.table(assay = "ChIP-seq", metric = "FRiP (broad marks)",    direction = "min", value = 0.02,
             note = paste("IHEC minimum;", paste(setdiff(MARKS, NARROW_MARKS), collapse = ", "))),
  data.table(assay = "ChIP-seq", metric = "FRiP (narrow marks)",   direction = "min", value = 0.05,
             note = paste("Landt et al. 2012;", paste(NARROW_MARKS, collapse = ", "))),
  data.table(assay = "WGBS", metric = "Mean_CpG_Coverage",         direction = "min", value = 10,
             note = "IHEC/ENCODE >=10x; matches the meth10 input downstream"),
  data.table(assay = "WGBS", metric = "BS_conversion_Rate",        direction = "min", value = 0.97,
             note = "IHEC >=97%"),
  data.table(assay = "WGBS", metric = "abs(GC_cov_Correlation)",   direction = "max", value = 0.70,
             note = "no community standard; project-set")
))

dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
fwrite(tab, out)
fwrite(thr, thr_out)

cat(sprintf("rows: %d (%d modelled, %d QC-only)   columns: %d\n",
            nrow(tab), n_modelled, nrow(tab) - n_modelled, ncol(tab)))
cat(sprintf("epirrs: %d total, %d modelled\n",
            uniqueN(tab$epirr), uniqueN(tab[modelled == TRUE, epirr])))
cat("protocol (all rows):\n"); print(tab[, .N, by = protocol])
cat("protocol (modelled only -- these are the counts to quote):\n")
print(tab[modelled == TRUE, .N, by = protocol])
cat(sprintf("qc_flag_count range: %d-%d, NA: %d\n",
            min(tab$qc_flag_count), max(tab$qc_flag_count), sum(is.na(tab$qc_flag_count))))
cat("\nPolicy A: no sample dropped for QC; qc_flag_count is a model covariate.\n")
cat(sprintf("written: %s\nwritten: %s\n", out, thr_out))
