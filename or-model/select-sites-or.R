## Select the sites the OR model runs on, and extract their time series.
##
## Ported from the last chunk of legacy/models/StateTransition.Rmd ("select sites of
## interest"). Draws --nsamp cells per disturbance class (0 = none, 10 = other, 20 = cut,
## 40 = fire, 50 = biotic) of the --year disturbance map (raster::sampleStratified), keeps
## those that were evergreen forest (land-cover class 42) the year before, and extracts at
## each one
##   dextract  disturbance class, 1985-2012
##   cextract  land-cover class,  1990-2017
##   bextract  biomass (Mg/ha),   1990-2017
##   points    site coordinates (lon / lat, WGS84)
## About 70% of the drawn cells are forest, so --nsamp 150 gives about 540 sites (the
## tracked data/inputs/extract.h03v04.RData was made this way, with nsamp = 150 and no seed).
##
## Needs the local rasters data/derived/disR.tif and data/raw/conus_{land_cover,biomass}_
## ARD_tile_h03v04.tif (see data/README.md). Then build the model inputs for the new sites:
##   Rscript or-model/load-data-or.R --extract <the file written here>
##
## Run from the repository root:
##   Rscript or-model/select-sites-or.R --nsamp 300 --seed 1
##
## Options (defaults in brackets):
##   --nsamp   cells drawn per disturbance class                                [150]
##   --seed    random seed                                                      [1]
##   --year    year of the disturbance map the classes are drawn from           [2004]
##   --out     output file             [data/inputs/extract.h03v04-n<nsamp>-seed<seed>.RData]

source("or-model/R/config.R")
check_repo_root()
suppressPackageStartupMessages({ library(raster); library(sp) })

nsamp <- as.integer(cli_opt("--nsamp", 150))
seed  <- as.integer(cli_opt("--seed", 1))
year  <- as.integer(cli_opt("--year", 2004))
out   <- cli_opt("--out", sprintf("data/inputs/extract.h03v04-n%d-seed%d.RData", nsamp, seed))
for (f in or.files[c("disR", "cover", "biomass")])
  if (!file.exists(f)) stop("missing raster ", f, " (see data/README.md)", call. = FALSE)
if (!year %in% intersect(dtime, ctime[-1])) stop("--year must be in ", min(ctime) + 1, ":", max(dtime))
say("OR site selection: %d cells per disturbance class in %d, seed %d -> %s", nsamp, year, seed, out)
t.all <- Sys.time()
nsteps <- 4

run_step(1, nsteps, "rasters", {
  disR <- stack(or.files$disR)
  covS <- stack(or.files$cover)
  bioS <- stack(or.files$biomass)
  stopifnot(nlayers(disR) == length(dtime), nlayers(covS) == length(ctime), nlayers(bioS) == length(ctime))
  message(sprintf("  %d x %d cells; disturbance %d layers, cover %d, biomass %d",
                  nrow(disR), ncol(disR), nlayers(disR), nlayers(covS), nlayers(bioS)))
})

run_step(2, nsteps, sprintf("stratified sample of the %d disturbance map", year), {
  set.seed(seed)
  samp  <- raster::sampleStratified(disR[[which(dtime == year)]], size = nsamp, exp = 1000)
  csamp <- covS[[which(ctime == year - 1)]][samp[, 1]]      # cover the year before
  fsamp <- samp[!is.na(csamp) & csamp == 42, , drop = FALSE] # evergreen forest only
  tab   <- table(factor(fsamp[, 2], levels = sort(unique(samp[, 2]))))
  message(sprintf("  drawn %d cells, %d forest: %s", nrow(samp), nrow(fsamp),
                  paste(sprintf("class %s: %d", names(tab), tab), collapse = ", ")))
})

run_step(3, nsteps, sprintf("extract the time series at %d sites", nrow(fsamp)), {
  dextract <- raster::extract(disR, fsamp[, 1])
  cextract <- raster::extract(covS, fsamp[, 1])
  bextract <- raster::extract(bioS, fsamp[, 1])
  coord  <- xyFromCell(bioS, fsamp[, 1])
  points <- spTransform(SpatialPoints(coord, crs(bioS)), CRS("+proj=longlat +datum=WGS84"))
  message(sprintf("  biomass missing at %.1f%% of site-years", 100 * mean(is.na(bextract))))
})

run_step(4, nsteps, "save", {
  site.selection <- list(nsamp = nsamp, seed = seed, year = year, cells = fsamp[, 1], made = Sys.time())
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  save(points, bextract, cextract, dextract, site.selection, file = out)
  message(sprintf("  saved %s", out))
})
say("OR site selection done in %s; next:  Rscript or-model/load-data-or.R --extract %s",
    fmt_secs(as.numeric(difftime(Sys.time(), t.all, units = "secs"))), out)
