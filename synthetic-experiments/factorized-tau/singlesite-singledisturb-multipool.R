## Factorized-tau version of the single site, single disturbance, multi-pool experiment.
##
## Reuses the set-up of ../singlesite-singledisturb-multipool.R (drivers, parameter
## ensemble, spin-up, proportional harvest in year 10, truth and observations, seed 42),
## then runs the filter with a factorized tau posterior (a 1-D grid per pool, the other
## taus held at their current estimates; or-joint-update-factorized.R) next to the 3-D
## grid filter and the known-tau filter, all through the same implementation.
##
## Filter set-up: scenarios k = 0 undisturbed, 1 harvest; two history groups (undisturbed /
## disturbed at the last step), as in ORJointUpdate; the harvest keeps a fraction mu0 of leaf
## and stem, so their process variance is scaled by mu0^2 under it.
##
## Run from the repository root:
##   Rscript synthetic-experiments/factorized-tau/singlesite-singledisturb-multipool.R

source("synthetic-experiments/factorized-tau/or-joint-update-factorized.R")
source_until("synthetic-experiments/singlesite-singledisturb-multipool.R", "^## ---- EnKF")
outdir <- "experiments/synthetic/factorized-tau/singlesite-singledisturb-multipool"

simulate <- function(x, t, label, mem) {
  tsel <- (t-1)*365 + 1:365
  out <- matrix(0, nrow(x), 3)
  for (n in seq_len(nrow(x))) {
    th <- THETA[mem[n], ]
    th[9:11] <- x[n, ]
    out[n, ] <- VSEM(th, PAR = PAR_driverdata[met[mem[n]], tsel])[365, 2:4]
  }
  out
}

filter_args <- list(
  Y = Observed_states, R = R, Xic = Ensemble_initial.orig, simulate = simulate,
  nMembers = n_ensemble_members, labels = c(1, 1),
  disturb = function(x, k) disturbance.rows(x, mu0, V0, alloc.soil),
  dist.p = rbind(rep(1 - p.dist, NYears), rep(p.dist, NYears)),
  procScale = rbind(c(1, 1, 1), c(mu0[1]^2, 1, mu0[2]^2)),   # leaf, soil, stem
  transition = matrix(1:2, 2, 2, byrow = TRUE),             # group after the step = scenario
  w0 = c(1, 1), poolOf = 1:3, prior = list(type = "loguniform"), N = 100)

## each axis spans process sd 0.01 .. 3.2: 16 points per 1-D grid (factorized: 48 grid
## points in all) and 16^3 = 4096 for the 3-D grid
grid <- list(tauMin = rep(0.1, 3), tauMax = rep(1e4, 3), G = rep(16, 3))
cmp  <- run_comparison(filter_args, grid, tau.known = 1/diag(Q))
compare_report(cmp$fits, cmp$times, True_states, Observed_states, pool.names, pool.names,
               truth.sd = sqrt(diag(Q)), dist.years = which(disturbance_at == 1), outdir = outdir,
               title = "Single site, single disturbance, multi-pool (factorized tau vs 3-D grid)")
