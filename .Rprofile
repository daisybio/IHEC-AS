source("renv/activate.R")
if (interactive() && Sys.getenv("RSTUDIO") == "") {
  source(file.path(Sys.getenv(if (.Platform$OS.type == "windows") "USERPROFILE" else "HOME"), ".vscode-R", "init.R"))
}
renv::settings$ignored.packages(c("cCRE_hits", "hits", "hits_used", "agg_hits", "chromhmm_hits"), persist = FALSE)

library(R.utils)
library(data.table)
library(pbmcapply)
library(rtracklayer)
library(ggplot2)
library(svglite)
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

sample_metadata_file <- "data/IHEC_metadata_harmonization.v1.3.extended.csv"
ontology_column <- "harmonized_sample_ontology_term_high_order_fig1"


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
cor_methods <- c("pearson", "spearman")

control_class <- "excluded"
case_class <- "included"
class_levels <- c(control_class, case_class)

minimum_events <- 25

vicinity <- 5e5

variability_colors <- c('Low'="#56B4E9", 'All'="#999999", 'High'="#D55E00")  

plot_dir <- "images/Rplots"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}
 my_lambda <- 'lambda.1se'
