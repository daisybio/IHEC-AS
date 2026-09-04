#!/usr/bin/env Rscript
# Atlas landscape summaries -- the small tables the resource figure is drawn from.
#
# WHY THIS STAGE EXISTS AT ALL. `11-paper-figures.Rmd` must stay fast enough to re-run on every
# figure tweak, so it is forbidden from touching anything large; the atlas is 5.1 GB gzipped and
# 1,134,719,490 rows. That rule is not negotiable for a figure, so the reduction happens here, once,
# in its own rule with its own declared outputs, and `11` reads three small CSVs.
#
# WHAT IT DELIBERATELY DOES NOT DO. No modelling, no filtering, no cohort restriction. This describes
# the RESOURCE -- all 1,522 quantified samples and all 745,545 events across five classes -- which is
# the whole point of the figure: the modelled cohort is a small corner of it (415 samples, 39,176
# events, ~0.8% of the PSI values), and that contrast only means something if the denominator is the
# real one.
#
# PSI IS NA WHERE A SAMPLE HAS NO READS ON EITHER SIDE OF THE JUNCTION (incl + skip == 0), so the
# observed FRACTION is a real coverage measurement rather than bookkeeping, and it is reported per
# ontology rather than pooled -- pooling it would hide exactly the cell-type-dependent detectability
# a reuser needs to know about.
#
# MEMORY. Each class is read with `select = c("uuid", "psi")` and released before the next. R's global
# string cache means the 1,522 distinct uuids cost 8 bytes per row, not 37, so SE's 851 M rows land at
# roughly 14 GB rather than 60 -- but SE is still the peak and the rule is sized for it, not for the
# 0.02 GB RI file that makes this look cheap.

suppressPackageStartupMessages({
  library(data.table)
})
setDTthreads(4L)

tf        <- Sys.getenv("TRANSCRIPT_FILTER", "biotype_filtered")
atlas_dir <- Sys.getenv("ATLAS_DIR", file.path("processed_data", "atlas", tf))
ft_file   <- Sys.getenv("FILE_TABLE", "processed_data/file_table.csv.gz")
meta_file <- Sys.getenv("HARMONIZATION_CSV",
                        "data/IHEC_sample_metadata_harmonization.v1.4_extended.csv")
cmd_file  <- Sys.getenv("RMATS_COMMANDS", "data/rmats_split_post_commands.tsv")
roster_f  <- Sys.getenv("COHORT_UUIDS", sprintf("processed_data/cohort_uuids_%s.csv", tf))

out_land  <- Sys.getenv("ATLAS_LANDSCAPE_OUT", sprintf("processed_data/atlas_landscape_%s.csv", tf))
out_hist  <- Sys.getenv("ATLAS_PSI_HIST_OUT",  sprintf("processed_data/atlas_psi_hist_%s.csv", tf))
out_strat <- Sys.getenv("ATLAS_STRATA_OUT",    sprintf("processed_data/atlas_strata_%s.csv", tf))

EVENT_TYPES <- c("SE", "RI", "MXE", "A3SS", "A5SS")
N_BINS <- 50L

# --- sample annotation: uuid -> epirr -> ontology --------------------------------------------
# Both joins are over the FULL 1,522, not the modelled 415. `file_table.csv.gz` is 01's own output
# and covers every RNA-seq uuid; the harmonization table keys on a VERSIONED EpiRR
# (`IHECRE00001032.6`) while every id on this side is unversioned, so the suffix must be stripped --
# joining them raw matches nothing and fails silently as an all-NA column (hit for real in 05d).
ft <- fread(ft_file, select = c("uuid", "epirr_id_without_version", "experiment_type"))
ft <- unique(ft[experiment_type %in% c("mRNA-Seq", "total-RNA-Seq")], by = "uuid")
setnames(ft, "epirr_id_without_version", "epirr")

ONTOLOGY_COL <- "harmonized_sample_ontology_term_high_order_fig1"
meta <- fread(meta_file, select = c("EpiRR", ONTOLOGY_COL))
setnames(meta, c("epirr", "ontology"))
meta[, epirr := sub("\\.[0-9]+$", "", epirr)]
meta <- unique(meta, by = "epirr")
ft[meta, on = "epirr", ontology := i.ontology]

roster <- if (file.exists(roster_f)) fread(roster_f, select = "uuid")$uuid else character(0)
ft[, modelled := uuid %in% roster]

# --- library strata, from the rMATS command table --------------------------------------------
# The 14 runs ARE the batch axis of this resource: rMATS was run once per (paired/single, libtype,
# read length) combination, so a reuser needs them named, not folded away.
runs <- fread(cmd_file)
# One row per (transcript_filter x run): the same 14 physical strata are re-run for each filter, and
# `out_folder` is the only column naming which. Without this filter the table reads 42 runs and 4,566
# samples -- each stratum counted three times -- which would have published a 3x sample count.
# Same predicate as `02-3b-atlas-summary.R:43`, deliberately.
runs <- runs[grepl(tf, out_folder, fixed = TRUE)]
if (!nrow(runs)) stop("no rows in ", cmd_file, " whose out_folder names transcript_filter '", tf, "'")
strata <- runs[, .(
  run_name, paired_end_str, libtype, read_length,
  n_uuids = lengths(strsplit(uuid_string, ",", fixed = TRUE))
)][order(-n_uuids)]
fwrite(strata, out_strat)
cat(sprintf("strata: %d runs, %d uuids total\n", nrow(strata), sum(strata$n_uuids)))

# --- per-class pass ---------------------------------------------------------------------------
land <- list(); hist_l <- list()

for (et in EVENT_TYPES) {
  f <- file.path(atlas_dir, sprintf("%s.csv.gz", et))
  if (!file.exists(f)) {
    warning("missing atlas file, skipping: ", f, call. = FALSE, immediate. = TRUE)
    next
  }
  t0 <- Sys.time()
  dt <- fread(f, select = c("uuid", "psi"))
  dt[ft, on = "uuid", c("ontology", "modelled") := .(i.ontology, i.modelled)]

  # Per (ontology) coverage and location. `frac_observed` is over the FULL grid for that ontology,
  # so it is comparable across classes; the boundary fractions are computed over OBSERVED values
  # only, because an unobserved event is not "at a bound", it is absent.
  land[[et]] <- dt[, {
    obs <- !is.na(psi)
    n_obs <- sum(obs)
    .(
      event_type    = et,
      n_uuids       = uniqueN(uuid),
      n_rows        = .N,
      n_observed    = n_obs,
      frac_observed = n_obs / .N,
      median_psi    = if (n_obs) median(psi[obs]) else NA_real_,
      mean_psi      = if (n_obs) mean(psi[obs]) else NA_real_,
      frac_at_0     = if (n_obs) sum(psi[obs] == 0) / n_obs else NA_real_,
      frac_at_1     = if (n_obs) sum(psi[obs] == 1) / n_obs else NA_real_
    )
  }, by = .(ontology, modelled)]

  # Histogram. The two boundary atoms get their OWN rows rather than being folded into the first and
  # last bin: 40.5% of all PSI values sit exactly at a bound, which is the single most important
  # distributional fact about this resource (it is why a logit transform is unusable here), and a
  # 50-bin histogram that buries them inside [0, 0.02] and [0.98, 1] hides it.
  hist_l[[et]] <- dt[!is.na(psi), {
    at0 <- sum(psi == 0); at1 <- sum(psi == 1)
    mid <- psi[psi > 0 & psi < 1]
    br  <- seq(0, 1, length.out = N_BINS + 1L)
    cnt <- tabulate(cut(mid, breaks = br, include.lowest = TRUE, labels = FALSE), nbins = N_BINS)
    rbind(
      data.table(event_type = et, bin_lo = 0,           bin_hi = 0,          atom = TRUE,  count = at0),
      data.table(event_type = et, bin_lo = br[-(N_BINS + 1L)], bin_hi = br[-1L], atom = FALSE, count = cnt),
      data.table(event_type = et, bin_lo = 1,           bin_hi = 1,          atom = TRUE,  count = at1)
    )
  }]

  cat(sprintf("%-5s %13s rows  %6.1f%% observed  %.1f min\n", et, format(nrow(dt), big.mark = ","),
              100 * mean(!is.na(dt$psi)), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  rm(dt); gc(verbose = FALSE)
}

if (!length(land)) stop("no atlas class files were readable under ", atlas_dir)

landscape <- rbindlist(land)
fwrite(landscape, out_land)
fwrite(rbindlist(hist_l), out_hist)

# --- checks -----------------------------------------------------------------------------------
# The uuid count is the one number a reader will carry away, so assert it rather than print it.
tot_uuids <- landscape[, sum(n_uuids), by = event_type]
if (uniqueN(tot_uuids$V1) != 1L) {
  stop("classes disagree on sample count: ", paste(sprintf("%s=%d", tot_uuids$event_type, tot_uuids$V1),
                                                   collapse = ", "))
}
if (anyNA(landscape$ontology)) {
  stop(landscape[is.na(ontology), sum(n_uuids)], " uuid(s) have no ontology -- check the EpiRR join")
}

cat(sprintf("\nsamples: %d  ontologies: %d  modelled: %d\n",
            tot_uuids$V1[1], uniqueN(landscape$ontology),
            landscape[modelled == TRUE, sum(n_uuids)] / uniqueN(landscape$event_type)))
cat(sprintf("written: %s\nwritten: %s\nwritten: %s\n", out_land, out_hist, out_strat))
