## Factorized-tau version of the two-site, single disturbance, multi-pool experiment.
##
## Reuses the set-up and filter arguments of ../multisite-singledisturb-multipool.R (site 1
## disturbed in year 10; scenarios none / site 1 / site 2 / both; the three taus shared by
## the two sites), then runs the filter with a factorized tau posterior (a 1-D grid per pool,
## the other taus held at their current estimates; or-joint-update-factorized.R) next to the
## 3-D grid filter and the known-tau filter. Set oracle <- TRUE in the base script to give
## every member the true parameters and driver.
##
## Run from the repository root:
##   Rscript synthetic-experiments/factorized-tau/multisite-singledisturb-multipool.R

source("synthetic-experiments/factorized-tau/or-joint-update-factorized.R")
source_until("synthetic-experiments/multisite-singledisturb-multipool.R", "^## tau learned")
outdir <- file.path("experiments/synthetic/factorized-tau",
                    if (oracle) "multisite-singledisturb-multipool-oracle" else "multisite-singledisturb-multipool")

## each axis spans process sd 0.01 .. 3.2 with 12 points, as in the base experiment
## (factorized: 36 grid points in all; 3-D grid: 12^3 = 1728)
grid <- list(tauMin = rep(0.1, 3), tauMax = rep(1e4, 3), G = rep(12, 3))
cmp  <- run_comparison(filter_args, grid, tau.known = 1/diag(Q))
site.p <- function(f, s) rowSums(f$pK[, sapply(scen.sites, function(ss) s %in% ss), drop = FALSE])
compare_report(cmp$fits, cmp$times, Ntrue, No, state.names, pool.names,
               truth.sd = sqrt(diag(Q)), dist.years = td, outdir = outdir,
               title = "Two sites, single disturbance, multi-pool (factorized tau vs 3-D grid)",
               events = function(f) cbind(`P(site 1 disturbed)` = site.p(f, 1), `P(site 2 disturbed)` = site.p(f, 2)))
