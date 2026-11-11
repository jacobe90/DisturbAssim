## Load the OR-model data (load-data-or.R), then run the filter on it (run-model-or.R)
## in the same R session. Takes the options of run-model-or.R.
##
## Run from the repository root:
##   Rscript or-model/load-data-and-run-model.R
##   Rscript or-model/load-data-and-run-model.R --sites 300 --R GEDI_het --distprod TRUE

if (!file.exists("or-model/R/config.R")) stop("run from the repository root (the folder holding or-model/)")
t.total <- Sys.time()
message("==== 1/2: load data ====")
source("or-model/load-data-or.R")
message("==== 2/2: run model ====")
source("or-model/run-model-or.R")
say("load and run finished in %s", fmt_secs(as.numeric(difftime(Sys.time(), t.total, units = "secs"))))
