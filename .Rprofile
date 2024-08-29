source("renv/activate.R")
if (interactive() && Sys.getenv("RSTUDIO") == "") {
  source(file.path(Sys.getenv(if (.Platform$OS.type == "windows") "USERPROFILE" else "HOME"), ".vscode-R", "init.R"))
}
renv::settings$ignored.packages(c("cCRE_hits", "hits", "hits_used", "agg_hits", "chromhmm_hits"), persist = FALSE)
# options(renv.config.pak.enabled = TRUE)

options(
  clustermq.scheduler = "slurm",
  clustermq.template = "clustermq.template.slurm" # if using your own template
)


library(R.utils)
library(data.table)
library(pbmcapply)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(umap)
library(pheatmap)
library(UpSetR)
library(rtracklayer)
library(patchwork)
library(httr)
library(glmnet)
library(ranger)
library(caret)
library(RSNNS)
library(MLmetrics)

ncores <- 40
data.table::setDTthreads(ncores)
message(sprintf("mc.cores: %d", options()$mc.cores))
options(mc.cores = ncores)
message(sprintf("set mc.cores to: %d", options()$mc.cores))
setOption("ignore.interactive", TRUE)

if (!rmarkdown::pandoc_available()) {
  Sys.setenv(RSTUDIO_PANDOC = "/usr/lib/rstudio-server/bin/quarto/bin/tools")
}

data_dir <- "/nfs/data/IHEC/RNAseq"
data_dir2 <- "/nfs/data3/IHEC"
rna_data_dir <- file.path(data_dir, "RNA-Seq")
chip_data_dir <- file.path(data_dir2, "ChIP-Seq")
wgbs_data_dir <- file.path(data_dir, "WGBS")
wgbs_matrices_data_dir <- file.path(data_dir, "WGBS_matrices")
sample_dt_dir <- "sample_dts"

histone_marks <- c("H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K9me3")
as_events <- c("SE", "RI", "AL", "AF", "A3", "A5", "MX")
to_analyze <- c("SE", "RI")
aggregation_functions <- c("median" = median, "mean" = mean, "max" = max, "sum" = sum)
keep_cols <- function(dt, aggregation_method) {
  !grepl(paste(names(aggregation_functions)[names(aggregation_functions) != aggregation_method], collapse = "|"), names(dt))
}
flank_size <- 150
cor_methods <- c("pearson", "spearman")
my_lambda <- "lambda.1se"

control_class <- "excluded"
case_class <- "included"
class_levels <- c(control_class, case_class)

minimum_events <- 25

vicinity <- 5e5

sample_metadata_file <- "IHEC_metadata_harmonization.v1.2.extended.csv"
ontology_column <- "harmonized_sample_ontology_term_high_order_fig1"
variability_colors <- c('Low'="#56B4E9", 'All'="#999999", 'High'="#D55E00")  

plot_dir <- "images/Rplots"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}
publication_plots_dir <- "../65f465a3da21f80541769df9/images"

split_features <- function(dt, feature_col, sep = ";") {
  dt[, (c("mark", "region")) := tstrsplit(gsub("`", "", get(feature_col)), sep, fixed = TRUE, keep = c(1, 2))]
  dt[, mark := as.factor(mark)]

  dt[`Event Type` == "SE" & region == "event_name", region := "SE"]
  dt[`Event Type` == "RI" & region == "event_name", region := "RI"]
  dt[`Event Type` == "SE", region := gsub("other_region", "Intron", region, fixed = TRUE)]
  dt[`Event Type` == "RI", region := gsub("other_region", "Exon", region, fixed = TRUE)]

  dt[is.na(region), region := "Gene-Wide"]
  dt[, region := as.factor(region)]
}
