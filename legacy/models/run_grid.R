## Run the experiment grid from OR_joint_update2.Rmd outside RStudio.
##
## Sources the model setup (everything before the `run-single` chunk) and then the
## `run-grid` chunk, with the settings written there (grid.outdir, grid.bias,
## grid.cores, grid.sites, ...). Finished runs are skipped, so an interrupted grid
## can be restarted with the same command. Forked workers can deadlock inside
## RStudio, so prefer this for the full grid.
##
## Run from the repo root:  Rscript run_grid.R

pdf(NULL)   # swallow the diagnostic plots drawn by the setup chunks
rfile <- tempfile(fileext = ".R")
knitr::purl("OR_joint_update2.Rmd", output = rfile, quiet = TRUE)
code  <- readLines(rfile)
hdr   <- grep("^## ----", code)                    # chunk headers
rs    <- grep("^## ----run-single", code)
rg    <- grep("^## ----run-grid", code)
stopifnot(length(rs) == 1, length(rg) == 1)
re    <- c(hdr[hdr > rg], length(code) + 1)[1]     # end of the run-grid chunk
source(textConnection(c(code[seq_len(rs - 1)], code[rg:(re - 1)])))
