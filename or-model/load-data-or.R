## Load and preprocess everything the OR model needs, and cache it.
##
## Reads the tracked inputs in data/inputs (disturbance prior, site extract, disturbance-
## size mixtures, VSEM calibration), loads daymet for every site (downloading the sites not
## yet in data/cache/daymet.h03v04.RData), and builds
##   rad, radDA                   daily PAR per site (all years / DA years)
##   Yw                           wood-carbon observations per site and year
##   d.prior, dist.priors,        disturbance priors and the forecast parameters of each
##   dist.forecastparams,         disturbance class / size cluster
##   dist.prior.probs, n.dist
##   ICprior                      initial-condition ensemble (VSEM spin-up + recovery)
##   theta, Q                     VSEM parameters and process-error covariance
## in the global environment, and saves them to data/cache/or-model-inputs.RData for
## run-model-or.R. Diagnostic plots go to experiments/or-model/load-data-diagnostics.pdf.
##
## Run from the repository root:
##   Rscript or-model/load-data-or.R                                   # the tracked sites
##   Rscript or-model/load-data-or.R --extract data/inputs/extract.h03v04-n300-seed1.RData
##
## Options:
##   --extract   site extract to use (from select-sites-or.R)  [data/inputs/extract.h03v04.RData]

source("or-model/R/config.R")
check_repo_root()
suppressPackageStartupMessages(library(BayesianTools))   # VSEM
source("or-model/R/data-met.R")
source("or-model/R/data-disturbance.R")
source("or-model/R/data-biomass.R")

or.extract.file <- cli_opt("--extract", or.files$extract)
say("OR model: loading data for the sites in %s", or.extract.file)
t.all <- Sys.time()
nsteps <- 6

run_step(1, nsteps, "tracked inputs (data/inputs)", {
  for (f in c(or.files[c("STprior", "disPrior", "fit")], or.extract.file)) {
    if (!file.exists(f)) stop("missing ", f, " (see data/README.md)")
    message(sprintf("  %s: %s", f, paste(load(f, envir = globalenv()), collapse = ", ")))
  }
  colnames(dextract) <- paste0("t", dtime)
  dextract[is.na(dextract)] <- 0
  d.prior <- clean_disturbance_prior(d.prior)
  Qsoil   <- Qs
  message(sprintf("  %d sites; disturbance product %d-%d, biomass %d-%d, DA years %d-%d",
                  nrow(dextract), min(dtime), max(dtime), min(ctime), max(ctime), min(atime), max(atime)))
})

run_step(2, nsteps, "meteorology (daymet -> PAR)", {
  latlon <- points@coords
  dm  <- load_met(latlon)
  par.out <- met_to_par(dm, atime)
  rad <- par.out$rad; radDA <- par.out$radDA; metTime <- par.out$metTime
  rm(dm)
  message(sprintf("  PAR: %d sites x %d days (DA window: %d days)", nrow(rad), ncol(rad), ncol(radDA)))
})

run_step(3, nsteps, "disturbance priors and forecast parameters", {
  dist.forecastparams <- build_dist_forecastparams(dist.priors)
  dist.prior.probs    <- build_dist_prior_probs(dist.forecastparams, d.prior)
  n.dist <- length(dist.forecastparams)
  NT     <- length(atime)
  message(sprintf("  %d disturbance components: %s", n.dist,
                  paste(sapply(dist.forecastparams, `[[`, "name"), collapse = ", ")))
})

run_step(4, nsteps, "biomass observations", {
  ob <- build_obs(bextract); Yw <- ob$Yw; Rw <- ob$Rw
  AGBprior <- Yw[, paste0("t", atime[1])]   # wood at the start of the DA window, every site
  message(sprintf("  Yw: %d sites x %d years; %.1f%% missing", nrow(Yw), ncol(Yw), 100 * mean(is.na(Yw))))
})

run_step(5, nsteps, "initial-condition prior (VSEM spin-up and recovery trajectories)", {
  ic <- build_ic_prior(theta, rad, AGBprior, nNBG)
  ICprior <- ic$ICprior; age <- ic$age
  ## diagonal process-error covariance (per-pool sds from the calibration); Q[3, 3] is
  ## replaced by 1/tau_g inside ORJointUpdate
  Q <- diag(c(Ql, Qsoil, Qw)^2)
  message(sprintf("  ICprior: %d members; median age %g years", nrow(ICprior), median(age)))
})

run_step(6, nsteps, "diagnostic plots and cache", {
  dir.create("experiments/or-model", recursive = TRUE, showWarnings = FALSE)
  pdf("experiments/or-model/load-data-diagnostics.pdf", width = 8, height = 6)
  plot(metTime, rad[1, ], type = "l", xlab = "date", ylab = "PAR (MJ m-2 day-1)", main = "PAR, site 1")
  hist(age, main = "age of each IC member on its recovery trajectory", xlab = "years since disturbance")
  pairs(ICprior, main = "initial-condition prior", pch = ".")
  invisible(dev.off())
  or.input.objects <- c("dtime", "ctime", "atime", "points", "bextract", "dextract", "d.prior", "dist.priors",
                        "dist.forecastparams", "dist.prior.probs", "n.dist", "NT", "rad", "radDA", "Yw", "Rw",
                        "ICprior", "age", "theta", "Ql", "Qw", "Qs", "Qsoil", "Q", "or.extract.file")
  dir.create(dirname(or.files$orInputs), recursive = TRUE, showWarnings = FALSE)
  save(list = c(or.input.objects, "or.input.objects"), file = or.files$orInputs)
  message(sprintf("  saved %s (%.1f MB)", or.files$orInputs, file.size(or.files$orInputs) / 2^20))
})

or.inputs.loaded <- TRUE
say("OR model data loaded in %s", fmt_secs(as.numeric(difftime(Sys.time(), t.all, units = "secs"))))
