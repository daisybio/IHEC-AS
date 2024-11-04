# IHEC-AS

This is the code repository for [Revisiting Evidence for Epigenetic Control of Alternative Splicing](https://doi.org/10.1101/2024.08.30.610315).

Most analyses are done in R using the .Rmd files in this repo. To keep track of R package versions, we use `renv`. To restore this project's versions use `renv::restore()`.

Some analyses are done in other languages. For those, we have a mamba/conda repository that you can restore using `mamba env create -f env.yml`.

Please run the scripts in the order indicated by the number with which the files start. The first notebook is [01-gather-data.Rmd](01-gather-data.Rmd) and the last one is [10-experimental-events.Rmd](10-experimental-events.Rmd). 

You can adjust parameters for your local machine in the [.Rprofile](.Rprofile) file.