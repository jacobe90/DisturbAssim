## Source the OR-model code and load the inputs cached by or-model/load-data-or.R into the
## global environment, for scripts that run ORJointUpdate themselves.
##
##   source("or-model/R/load-or-model.R")    # from the repository root

source("or-model/R/config.R")
check_repo_root()
suppressPackageStartupMessages(library(BayesianTools))   # VSEM
source("R/disturbance.R")                                 # disturbance.p.rows
source("or-model/R/or-joint-update.R")
source("or-model/R/make-or-args.R")
source("or-model/R/run-grid.R")
source("or-model/R/plots.R")
if (!file.exists(or.files$orInputs))
  stop("no ", or.files$orInputs, ": run  Rscript or-model/load-data-or.R  first", call. = FALSE)
load(or.files$orInputs, envir = globalenv())
