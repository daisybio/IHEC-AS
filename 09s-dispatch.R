## ---------------------------------------------------------------------------
## Tier-1 ridge-screen dispatcher
##
## Fires the fit-free ridge screen as a SLURM array over ALL modelable events.
## Mirrors 09-1's fire-and-forget sbatch convention: returns immediately, the
## array runs async. Aggregation (qvalue → screen_results.csv.gz + tier1_hits)
## is a SEPARATE Snakemake rule (`screen_aggregate` → 09s-aggregate.R), run once
## the array has finished — it guards against partial results.
##
## Chunks the id list into blocks ≤ MaxArraySize (SLURM caps array indices), one
## sbatch per chunk, each with an OFFSET into the shared id file.
##
## CLI:
##   Rscript 09s-dispatch.R <cfg_rds> [throttle]
## Requires 09-1-ml-local.R to have run first (writes cfg + all-events id file).
## ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
cfg_path <- args[1]
throttle <- if (length(args) >= 2L) as.integer(args[2L]) else 50L

cfg <- readRDS(cfg_path)
setwd(cfg$project_dir)

tf <- basename(cfg$event_dir)
all_ids_file <- file.path("processed_data", sprintf("event_glmnet_all_ids_%s.txt", tf))
if (!file.exists(all_ids_file)) {
  stop("Missing all-events id file (run 09-1 build first): ", all_ids_file)
}
n <- length(readLines(all_ids_file))
if (n == 0L) stop("Empty id file: ", all_ids_file)

array_sh <- normalizePath("09s-ridge-screen-array.sh")
cfg_abs <- normalizePath(cfg_path)
proj <- normalizePath(".")
dir.create("event_glmnet_logs", showWarnings = FALSE)

# SLURM MaxArraySize (max array INDEX + 1). Probe scontrol; default 1001.
max_array <- tryCatch({
  cfgtxt <- system("scontrol show config 2>/dev/null", intern = TRUE)
  ln <- grep("MaxArraySize", cfgtxt, value = TRUE)
  if (length(ln)) as.integer(sub(".*=\\s*", "", ln[1])) else 1001L
}, error = function(e) 1001L)
chunk <- max_array # events per array submission

offsets <- seq(0L, n - 1L, by = chunk)
job_ids <- character(0)
for (off in offsets) {
  m <- min(chunk, n - off) # tasks in this chunk
  cmd <- sprintf(
    'sbatch --parsable --array=0-%d%%%d "%s" "%s" "%s" "%s" %d',
    m - 1L, throttle, array_sh, all_ids_file, cfg_abs, proj, off
  )
  jid <- system(cmd, intern = TRUE)
  jid <- sub("[;].*$", "", trimws(jid[length(jid)])) # strip cluster suffix
  job_ids <- c(job_ids, jid)
  message(sprintf("Submitted screen chunk offset=%d size=%d (job %s)", off, m, jid))
}
message(sprintf(
  "Submitted %d screen chunk(s) over %d events (jobs %s). When the array has ",
  length(job_ids), n, paste(job_ids, collapse = ",")
))
message("finished, run `snakemake ... screen_aggregate` (09s-aggregate.R).")
