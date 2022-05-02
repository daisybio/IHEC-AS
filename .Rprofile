source("renv/activate.R")

library(data.table)
library(pbmcapply)
library(rtracklayer)
ncores <- 20
message(sprintf('mc.cores: %d', options()$mc.cores))
options(mc.cores = ncores)
message(sprintf('set mc.cores to: %d', options()$mc.cores))

Sys.setenv(RSTUDIO_PANDOC='/usr/lib/rstudio-server/bin/pandoc')

data_dir <- '/nfs/data/IHEC/RNAseq'
rna_data_dir <- file.path(data_dir, 'RNA-Seq')
chip_data_dir <- file.path(data_dir, 'ChIP-Seq')
wgbs_data_dir <- file.path(data_dir, 'WGBS')

as_events <- c('SE', 'RI')
aggregation_functions <- c('median' = median, 'mean' = mean, 'max' = max)
flank_size <- 150
