source("renv/activate.R")
# if (interactive() && Sys.getenv("RSTUDIO") == "") {
#   source(file.path(
#     Sys.getenv(if (.Platform$OS.type == "windows") "USERPROFILE" else "HOME"),
#     ".vscode-R",
#     "init.R"
#   ))
# }
renv::settings$ignored.packages(
  c("cCRE_hits", "hits", "hits_used", "agg_hits", "chromhmm_hits"),
  persist = FALSE
)

# global seed
options(EpiATLAS_AS_SEED = 42L)

# --- pipeline scope + VST gene-expression gate (PLAN §4.12b/§4.12e) ----------
# All overridable via environment variable (Snakemake passes them per rule);
# defaults reproduce standard behaviour for interactive / un-migrated runs.
.epiatlas_env_chr <- function(name, default) {
  v <- Sys.getenv(name, unset = "")
  if (nchar(v) > 0L) v else default
}
.epiatlas_env_num <- function(name, default) {
  v <- Sys.getenv(name, unset = "")
  if (nchar(v) > 0L) as.numeric(v) else default
}
options(
  # default transcript_filter when TRANSCRIPT_FILTER env is unset (interactive)
  EpiATLAS_AS_PRIMARY_FILTER = .epiatlas_env_chr(
    "EpiATLAS_AS_PRIMARY_FILTER", "biotype_filtered"
  ),
  # VST low-expression gate used to build the high-confidence event set in
  # 02-3-rmats-event-filtering.Rmd and asserted in 05-create-aggregated-dt.Rmd
  EpiATLAS_AS_VST_CUTOFF_METHOD = .epiatlas_env_chr(
    "EpiATLAS_AS_VST_CUTOFF_METHOD", "gmm"
  ), # {"gmm","fixed","percentile"}
  EpiATLAS_AS_VST_CUTOFF_VALUE = .epiatlas_env_num(
    "EpiATLAS_AS_VST_CUTOFF_VALUE", NA_real_
  ), # used when method="fixed" or as gmm fallback
  EpiATLAS_AS_VST_CUTOFF_PERCENTILE = .epiatlas_env_num(
    "EpiATLAS_AS_VST_CUTOFF_PERCENTILE", 0.25
  ), # used when method="percentile" or as final fallback
  EpiATLAS_AS_VST_NA_POLICY = .epiatlas_env_chr(
    "EpiATLAS_AS_VST_NA_POLICY", "keep"
  ), # {"keep","drop"}: (event,sample) with no host-gene VST
  EpiATLAS_AS_VST_FLOOR_EPS = .epiatlas_env_num(
    "EpiATLAS_AS_VST_FLOOR_EPS", 0
  ) # extra trim ABOVE the zero-inflation floor atom before the GMM fit. 0 (strict
  # '>' min) removes only the exact structural-zero point mass -> antimode in the
  # true valley (~7.15). >0 eats the real low-expression bump and inflates the
  # cutoff (0.5->7.64, 1.0->8.31); raise only for a deliberately stricter call.
)

# vscode specific libraries
if (interactive()) {
  library(jsonlite)
  library(rlang)
  library(languageserver)
  library(httpgd)
  options(setWidthOnResize = TRUE)
}

# general libraries
# R.utils is loaded but NOT used anywhere in the pipeline (no R.utils:: call or
# insert()/withTimeout()/etc.). On some compute nodes the module-R 4.2.1 renv
# library is platform-incomplete and R.utils is absent — an unconditional
# library(R.utils) then kills .Rprofile at startup, failing every R rule there.
# Guard it so a missing (unused) R.utils only warns.
if (requireNamespace("R.utils", quietly = TRUE)) {
  library(R.utils)
} else {
  warning(".Rprofile: R.utils unavailable — skipping (unused by the pipeline)")
}
library(data.table)
library(pbmcapply)
library(rtracklayer)
library(ggplot2)
library(svglite)
library(grDevices)
# library(ggpubr)
# library(ggrepel)
# library(umap)
# library(pheatmap)
# library(UpSetR)
# library(patchwork)
# library(httr)
# library(glmnet)
# library(ranger)
# library(caret)
# library(RSNNS)
# library(MLmetrics)

main_theme <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(strip.background = element_rect(fill = NA))
}
ggplot2::theme_set(
  main_theme()
)

sample_metadata_file <- "data/IHEC_sample_metadata_harmonization.v1.4_extended.csv"
ontology_column <- "harmonized_sample_ontology_term_high_order_fig1"


slurm_cpus <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = "")
ncores <- if (nchar(slurm_cpus) > 0L) as.integer(slurm_cpus) else 40L
data.table::setDTthreads(ncores)
# Cap the multithreaded OpenBLAS to the allocation. Unset, openblas-pthread
# grabs the whole node's core count (ignoring the SLURM cgroup) → thread
# oversubscription/contention on the matrix-heavy steps (DESeq2 vst transform,
# rowSds diagnostics, big data.table joins). Set at runtime via RhpcBLASctl:
# OPENBLAS_NUM_THREADS is read when the BLAS is dlopen'd, before .Rprofile runs,
# so Sys.setenv() here would be too late. (Snakefile shell.prefix also exports
# it before R starts, for when RhpcBLASctl isn't installed.)
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  try(RhpcBLASctl::blas_set_num_threads(ncores), silent = TRUE)
}
message(sprintf("mc.cores: %d", options()$mc.cores))
options(mc.cores = ncores)
message(sprintf("set mc.cores to: %d", options()$mc.cores))

# pandoc is needed by rmarkdown::render. It is NOT on every SLURM compute node
# (only the login node has /usr/bin/pandoc), and the old single hardcoded
# RStudio path fails there. Probe candidate dirs and point RSTUDIO_PANDOC at the
# first that actually contains a pandoc binary of a sufficient version.
if (!rmarkdown::pandoc_available("1.12.3")) {
  .pandoc_candidates <- unique(c(
    dirname(Sys.which("pandoc")), # pandoc already on PATH (any node)
    file.path(Sys.getenv("CONDA_PREFIX"), "bin"), # conda/mamba env (all nodes)
    "/usr/lib/rstudio-server/bin/quarto/bin/tools", # RStudio Server
    "/usr/lib/rstudio/bin/quarto/bin/tools", # RStudio Desktop
    "/usr/bin" # system
  ))
  for (.p in .pandoc_candidates) {
    if (nzchar(.p) && file.exists(file.path(.p, "pandoc"))) {
      Sys.setenv(RSTUDIO_PANDOC = .p)
      if (rmarkdown::pandoc_available("1.12.3")) break
    }
  }
  if (!rmarkdown::pandoc_available("1.12.3")) {
    warning(
      "pandoc >= 1.12.3 not found on this node. Install it into the env ",
      "(`mamba install -n ihec-as pandoc`) so it is available on all SLURM ",
      "nodes, or module-load it before rendering."
    )
  }
}

data_dir <- "/nfs/data/IHEC/RNAseq"
data_dir2 <- "/nfs/data3/IHEC"
rna_data_dir <- file.path(data_dir, "RNA-Seq")
chip_data_dir <- file.path(data_dir2, "ChIP-Seq")
wgbs_data_dir <- file.path(data_dir, "WGBS")
wgbs_matrices_data_dir <- file.path(data_dir, "WGBS_matrices")
sample_dt_dir <- "sample_dts"

histone_marks <- c(
  "H3K9me3",
  "H3K27me3",
  "H3K27ac",
  "H3K4me1",
  "H3K4me3",
  "H3K36me3"
)
#TODO: make decision on
as_events <- c("SE", "RI", "A3", "A5", "MX") #, "AL", "AF")
to_analyze <- c("SE", "RI")
cor_methods <- c("pearson", "spearman")

control_class <- "excluded"
case_class <- "included"
class_levels <- c(control_class, case_class)

minimum_events <- 25

vicinity <- 5e5

variability_colors <- c(
  "Low" = "#56B4E9",
  "All" = "#999999",
  "High" = "#D55E00"
)

# IHEC IA plot palette (cosmetic). .Rprofile is sourced at R startup — BEFORE
# any Rmd chunk runs (incl. 01's download chunk that fetches this file) — so an
# unconditional read here would crash every stage on a cold checkout. Guard it:
# missing file → warn + NULL palettes (ggplot defaults), never a fatal halt.
# Place data/IHEC_EpiATLAS_IA_colors_Apl01_2024.json (01 downloads it) for the
# real IHEC colors.
.ia_colors_file <- "data/IHEC_EpiATLAS_IA_colors_Apl01_2024.json"
if (file.exists(.ia_colors_file)) {
  ihec_ia_colors <- unlist(
    jsonlite::read_json(.ia_colors_file),
    recursive = FALSE
  )
  sample_hex_colors <- sapply(
    unlist(ihec_ia_colors$fig1_ontology_intermediate_merged, recursive = FALSE),
    function(x) {
      cols <- as.numeric(strsplit(x, ",")[[1]])
      rgb(cols[1], cols[2], cols[3], maxColorValue = 255)
    }
  )
  mark_hex_colors <- sapply(
    unlist(ihec_ia_colors$experiment, recursive = FALSE),
    function(x) {
      cols <- as.numeric(strsplit(x, ",")[[1]])
      rgb(cols[1], cols[2], cols[3], maxColorValue = 255)
    }
  )
  mark_hex_colors <- c(mark_hex_colors, DNAm = mark_hex_colors[["WGBS"]])
} else {
  warning(
    sprintf(
      "%s not found — plot palettes fall back to ggplot defaults. Run 01's download chunk or place the file, then re-run for IHEC colors.",
      .ia_colors_file
    )
  )
  ihec_ia_colors <- NULL
  sample_hex_colors <- NULL
  mark_hex_colors <- NULL
}

plot_dir <- "images/Rplots"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}
my_lambda <- "lambda.1se"
