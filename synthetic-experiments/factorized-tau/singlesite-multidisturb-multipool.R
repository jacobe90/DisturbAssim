## Factorized-tau version of the single site, multiple disturbance, multi-pool experiment.
##
## Reuses the set-up and filter arguments of ../singlesite-multidisturb-multipool.R (harvest
## in year 10; scenarios undisturbed / harvest / defoliation; three history groups), then runs
## the filter with a factorized tau posterior (a 1-D grid per pool, the other taus held at
## their current estimates; or-joint-update-factorized.R) next to the 3-D grid filter and the
## known-tau filter. Set oracle <- TRUE in the base script to give every member the true
## parameters and driver.
##
## Run from the repository root:
##   Rscript synthetic-experiments/factorized-tau/singlesite-multidisturb-multipool.R

source("synthetic-experiments/factorized-tau/or-joint-update-factorized.R")
source_until("synthetic-experiments/singlesite-multidisturb-multipool.R", "^## tau learned")
outdir <- file.path("experiments/synthetic/factorized-tau",
                    if (oracle) "singlesite-multidisturb-multipool-oracle" else "singlesite-multidisturb-multipool")

## each axis spans process sd 0.01 .. 3.2 with 12 points, as in the base experiment
## (factorized: 36 grid points in all; 3-D grid: 12^3 = 1728)
grid <- list(tauMin = rep(0.1, 3), tauMax = rep(1e4, 3), G = rep(12, 3))
cmp  <- run_comparison(filter_args, grid, tau.known = 1/diag(Q))
compare_report(cmp$fits, cmp$times, Ntrue, No, pool.names, pool.names,
               truth.sd = sqrt(diag(Q)), dist.years = td, outdir = outdir,
               title = "Single site, multiple disturbance, multi-pool (factorized tau vs 3-D grid)",
               events = function(f) cbind(`P(harvest)` = f$pK[, 2], `P(defoliation)` = f$pK[, 3]))
