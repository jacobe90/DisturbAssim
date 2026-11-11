## Settings shared by the OR-model scripts: file locations, year axes, unit conversions,
## observation-error parameters and the experiment-grid settings.

## ---- files --------------------------------------------------------------------------
## data/inputs is tracked; data/cache is local (see data/README.md)
or.files <- list(
  STprior   = "data/inputs/STprior.RData",          # d.prior (background disturbance rates)
  extract   = "data/inputs/extract.h03v04.RData",   # points, bextract, cextract, dextract
  disPrior  = "data/inputs/dis.postprocessed.RData",# dist.priors (disturbance-size mixtures)
  fit       = "data/inputs/fit.RData",              # VSEM calibration: theta, Ql, Qw, Qs
  met       = "data/cache/daymet.h03v04.RData",     # daymet per site (missing sites are downloaded)
  ## rasters for select-sites-or.R (local only, see data/README.md)
  disR      = "data/derived/disR.tif",               # disturbance class per year, 1985-2012
  cover     = "data/raw/conus_land_cover_ARD_tile_h03v04.tif",  # land cover, 1990-2017
  biomass   = "data/raw/conus_biomass_ARD_tile_h03v04.tif",     # biomass (Mg/ha), 1990-2017
  orInputs  = "data/cache/or-model-inputs.RData"    # everything load-data-or.R builds
)

## ---- year axes ----------------------------------------------------------------------
## The products cover different years (see legacy/models/StateTransition.Rmd), so every
## lookup goes through calendar year rather than column index:
##   dextract (disturbance product): 1985:2012
##   bextract (biomass / Yw):        1990:2017
##   DA timesteps:                   1990:2017 (timestep t <-> year atime[t])
dtime <- 1985:2012
ctime <- 1990:2017
atime <- 1990:2017

## ---- units and set-up constants ----------------------------------------------------------
B2C    <- 0.5     # biomass to carbon ratio
SW2PAR <- 0.45    # shortwave to PAR
nNBG   <- 200     # years of the post-disturbance recovery trajectories used for the IC prior

## ---- observation error: R = y * slope + intercept ------------------------------------
Rparams.GEDI.het         <- c(0.01331, 0.09761)
Rparams.GEDI.fixed       <- c(0.0000000, 0.8848249)
Rparams.Landtrendr.het   <- c(0.013310, 3.012645)
Rparams.Landtrendr.fixed <- c(0.000, 3.79986)
## observation-error models of the experiment grid: label -> (Rparams, source, model)
Rparams.grid <- list(
  GEDI_het       = list(Rparams = Rparams.GEDI.het,         source = "GEDI",       model = "het"),
  GEDI_con       = list(Rparams = Rparams.GEDI.fixed,       source = "GEDI",       model = "con"),
  LandTrendr_het = list(Rparams = Rparams.Landtrendr.het,   source = "LandTrendr", model = "het"),
  LandTrendr_con = list(Rparams = Rparams.Landtrendr.fixed, source = "LandTrendr", model = "con")
)

## ---- helpers ----------------------------------------------------------------------------
check_repo_root <- function() {
  if (!file.exists("or-model/R/config.R"))
    stop("run the OR-model scripts from the repository root (the folder holding or-model/)", call. = FALSE)
}

## "[12:03:04] message", for the progress printouts
say <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), sprintf(...)))

fmt_secs <- function(s) {
  if (!is.finite(s)) return("--")
  s <- round(s); sprintf("%dh%02dm%02ds", s %/% 3600, (s %% 3600) %/% 60, s %% 60)
}

## value of command-line option `name` (the argument after it), or `default`
cli_opt <- function(name, default = NULL, args = commandArgs(trailingOnly = TRUE)) {
  i <- which(args == name); if (length(i)) args[i[1] + 1] else default
}

## run expr as numbered step i of n, printing its name before and its run time after
run_step <- function(i, n, name, expr) {
  say("step %d/%d: %s ...", i, n, name)
  t0 <- Sys.time()
  val <- force(expr)
  say("step %d/%d: done in %s", i, n, fmt_secs(as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  invisible(val)
}
