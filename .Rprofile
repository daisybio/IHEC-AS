source("renv/activate.R")

library(R.utils)
library(data.table)
library(pbmcapply)
library(rtracklayer)
library(ggplot2)
library(pheatmap)
ncores <- 20
# data.table::setDTthreads(ncores)
message(sprintf('mc.cores: %d', options()$mc.cores))
options(mc.cores = ncores)
message(sprintf('set mc.cores to: %d', options()$mc.cores))
setOption("ignore.interactive", TRUE)

Sys.setenv(RSTUDIO_PANDOC='/usr/lib/rstudio-server/bin/pandoc')

data_dir <- '/nfs/data/IHEC/RNAseq'
data_dir2 <- '/nfs/data3/IHEC'
rna_data_dir <- file.path(data_dir, 'RNA-Seq')
chip_data_dir <- file.path(data_dir2, 'ChIP-Seq')
wgbs_data_dir <- file.path(data_dir, 'WGBS')

as_events <- c('SE', 'RI', 'AL', 'AF', 'A3', 'A5', 'MX')
to_analyze <- c('SE', 'RI')
aggregation_functions <- c('median' = median, 'mean' = mean, 'max' = max, 'sum' = sum)
keep_cols <- function(dt, aggregation_method) {!grepl(paste(names(aggregation_functions)[names(aggregation_functions) != aggregation_method], collapse = '|'), names(dt))}
flank_size <- 150
cor_methods <- c('pearson', 'spearman')
