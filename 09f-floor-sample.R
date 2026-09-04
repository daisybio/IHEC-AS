#!/usr/bin/env Rscript
# Detection-floor event sampler -- `checkpoint floor_sample` in the Snakefile.
# See revision/EXECUTION-PLAN.md item 1.
#
# Stratified sample for the rho sweep. Strata are Event Type x Variability x blow-up, and the
# point of running it now is to see the REALISED cell sizes: RI has only 1,846 events, so some
# cells may not reach the target and a detection rate estimated from 12 events has an SE of
# ~14% that must not hide behind an "n=600" headline.
#
# `All` is deliberately NOT a stratum -- it is the pooled row, reported alongside High/Low.
#
# Blow-up gets its own dimension because in v5 it is NO LONGER a small-n phenomenon (median
# n 279 vs 318; v4 was 90 vs 328), so it is roughly independent of n and stratifying on n
# would miss it. Within a cell, events are spread across the n range by ordered systematic
# sampling, which keeps n coverage without needing a fourth dimension.
#
# Blow-up status is taken from the `local` space: that is where the signal is injected
# (fs=local) and, per the rho=1 anchor, the only space with any recovery power at all.
#
# Reads only screen_results.csv.gz + event_annotations_dt. Writes ONE declared output,
# processed_data/event_models/{tf}/floor/floor_events.tsv (FLOOR_OUT), and nothing else --
# in particular it never touches the real screen/ directory.

suppressPackageStartupMessages(library(data.table))
setDTthreads(4L)

tf <- "biotype_filtered"
PER_CELL <- as.integer(Sys.getenv("PER_CELL", "75"))
SEED <- 20260827L

scr <- fread(sprintf("processed_data/event_models/%s/screen_results.csv.gz", tf),
  select = c("ID", "Event Type", "feature_set", "n_samples", "screen_R2", "note"))
scr <- scr[feature_set == "local"]

ann <- fread(sprintf("processed_data/event_annotations_dt_%s.csv.gz", tf),
  select = c("ID", "Variability"))
scr[ann[, .(ID, k_vb = Variability)], on = "ID", Variability := i.k_vb]

# eligible = screened cleanly with a finite statistic; a refused event has nothing to inject into
elig <- scr[(is.na(note) | note == "") & is.finite(screen_R2) & !is.na(Variability)]
elig[, blowup := abs(screen_R2) > 1]

cat(sprintf("eligible events (local, clean, finite): %s of %s\n\n",
  format(nrow(elig), big.mark = ","), format(nrow(scr), big.mark = ",")))

cat("=== available per stratum ===\n")
avail <- elig[, .(available = .N,
  n_min = min(n_samples), n_med = as.integer(median(n_samples)), n_max = max(n_samples)),
  by = .(`Event Type`, Variability, blowup)][order(`Event Type`, Variability, blowup)]
print(avail)

# ordered systematic sampling within a cell: sort by n, take evenly spaced ranks. Spreads the
# sample across the n range deterministically -- no seed dependence for the spread itself.
set.seed(SEED)
pick <- elig[, {
  k <- min(PER_CELL, .N)
  idx <- unique(round(seq(1, .N, length.out = k)))
  o <- order(n_samples, ID)
  .SD[o][idx]
}, by = .(`Event Type`, Variability, blowup)]

cat(sprintf("\n=== realised sample: %d events (target %d per cell x %d cells = %d) ===\n",
  nrow(pick), PER_CELL, nrow(avail), PER_CELL * nrow(avail)))
real <- pick[, .(sampled = .N,
  n_min = min(n_samples), n_med = as.integer(median(n_samples)), n_max = max(n_samples)),
  by = .(`Event Type`, Variability, blowup)][order(`Event Type`, Variability, blowup)]
real[avail, on = c("Event Type", "Variability", "blowup"), available := i.available]
real[, short_by := pmax(0L, PER_CELL - sampled)]
# detection rate SE at p=0.5, the worst case -- the number that says whether a cell is usable
real[, se_at_p50_pct := round(100 * 0.5 / sqrt(sampled), 1)]
print(real[, .(`Event Type`, Variability, blowup, available, sampled, short_by,
  se_at_p50_pct, n_min, n_med, n_max)])

cat("\n=== pooled rows as they will be reported ===\n")
print(pick[, .(sampled = .N), by = .(`Event Type`, Variability)][order(`Event Type`, Variability)])
print(pick[, .(sampled = .N), by = .(`Event Type`)][order(`Event Type`)])

cat(sprintf("\ntotal jobs at 8 rho values: %s\n", format(nrow(pick) * 8L, big.mark = ",")))
cat(sprintf("est. wall time at ~2 min/job, 50 concurrent: %.1f h\n", nrow(pick) * 8 * 2 / 50 / 60))

# FLOOR_OUT is set by the Snakemake rule so it owns the path; the default keeps the
# script runnable by hand.
f <- Sys.getenv("FLOOR_OUT", "processed_data/event_models/biotype_filtered/floor/floor_events.tsv")
dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
fwrite(pick[, .(ID, `Event Type`, Variability, blowup, n_samples)], f, sep = "\t")
cat(sprintf("\nwritten: %s\n", f))
cat("\nCells with sampled << 75 give a wide detection-rate SE -- report realised cell sizes\n")
cat("alongside every floor estimate rather than a pooled n.\n")
