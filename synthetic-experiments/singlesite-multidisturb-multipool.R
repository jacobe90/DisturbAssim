## Single site, multiple disturbance types, multi-pool (VSEM) synthetic experiment.
##
## Migrated from SingleSite_MultiDisturb_MultiPool.Rmd: the same drivers, parameter
## ensemble, spin-up and "true" run (a harvest in year 10), with the JAGS multinomial
## analysis + multi.ens.adj replaced by the generalised ORJointUpdate filter of
## or-joint-update-filter.R. The process error of each pool (leaf, soil, stem) is learned
## on a 3-D tau grid.
##
## Filter set-up:
##   - scenarios k = 0 undisturbed, 1 harvest, 2 defoliation, with the Rmd's prior p.dist
##   - three history groups (undisturbed / harvested / defoliated at the last step), so the
##     analysis is a 3-component mixture per tau grid point
##   - disturbance2 keeps a fraction mu0 of leaf and stem, so their process variance is
##     scaled by mu0^2 under a disturbance; the soil keeps its full process variance
##   - each importance sample takes a random ensemble member's parameters and driver every
##     year. That spread is not persistent, so the filter cannot learn that the truth keeps
##     one member's parameters; a systematic gap between that member and the ensemble shows
##     up as extra process error (leaf and stem sd well above the truth here). Set
##     oracle <- TRUE to give every member the true parameters and driver, which isolates
##     the process-error learning (the learned sd then recovers the truth).
##
## Run from the repository root:  Rscript synthetic-experiments/singlesite-multidisturb-multipool.R

library(BayesianTools)
library(mvnfast)
source("disturbance.R")                              # tobit, disturbance2
source("synthetic-experiments/or-joint-update-filter.R")

oracle <- FALSE   ## TRUE: every ensemble member uses the true parameters and driver
outdir <- file.path("experiments/synthetic", if (oracle) "singlesite-multidisturb-multipool-oracle" else "singlesite-multidisturb-multipool")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
set.seed(42)

## ---- Set-up (as in the Rmd) ------------------------------------------------------------
ne = 500     ## number of ensemble members
NT = 25      ## number of years
nt = 365*NT

## DRIVERS
nd   = 31 ## number of driver ensembles
PAR <- matrix(NA,nd,nt)
for(i in 1:nd){
  PAR[i,] <- VSEMcreatePAR(1:nt)
}

## Parameters
refPars <- VSEMgetDefaults()
refPars$lower[9:11] = 1
np = 50
sigmaP = diag((refPars$upper - refPars$lower)/10)^2 ## assume range is +/- 5 SD
Theta <- as.data.frame(mvtnorm::rmvnorm(np,refPars$best,sigmaP))
colnames(Theta) <- rownames(refPars)
for(i in 1:11){
  Theta[Theta[,i] < refPars$lower[i] , i] <- refPars$lower[i]
  Theta[Theta[,i] > refPars$upper[i] , i] <- refPars$upper[i]
}

## initial conditions
met    = sample.int(nd,ne,replace=TRUE)
parIdx = sample.int(np,ne,replace=TRUE)
Na.orig <- matrix(NA,ne,3)
for(i in 1:ne){
  Na.orig[i,] = VSEM(as.numeric(Theta[parIdx[i],]),PAR = PAR[met[i],])[nt,2:4]
}

Q = diag(c(0.1,0.1,0.1)) ## process error
R = diag(c(.25,10,1))    ## observation error
THETA <- as.matrix(Theta[parIdx,]) ## ensemble parameters

## disturbance
mu0 <- matrix(c(0.1,0.1,
                0.1,0.95),2,2)    ## mean leaf & stem biomass after disturbance as % of original
V1  <- diag(c(0.15,0.15))^2       ## var in post disturbance leaf and stem biomass (one for each type)
V2  <- diag(c(0.15,0.01))^2
V0  <- abind::abind(V1,V2,along=3)
alloc.soil = 0.25                 ## fraction of removed C that goes in the soil
p.dist = c(0.9,0.05,0.05)         ## Multinomial transition probability; first is no change
dist.names = c("undisturbed","harvest","defoliation")
pool.names = c("Cleaf","Csoil","Cstem")

## pseudodata
No = matrix(NA,NT,3)
Ntrue = matrix(NA,NT,3)
true = 43                              ## reference ensemble member
D.true = rep(0,NT)
D.true[10] = 1                         ## prescribe a harvest
refpar = refpar0 = as.numeric(Theta[parIdx[true],])
refpar[9:11] = Na.orig[true,]
for(t in 1:NT){
  tsel <- (t-1)*365 + 1:365
  ref.mu = VSEM(refpar,PAR[met[true],tsel])[365,2:4]
  ref.f  = mvtnorm::rmvnorm(1,ref.mu,Q)
  if(D.true[t]){
    ref.f = disturbance2(ref.f,type=D.true[t])
  }
  ref.f = tobit(ref.f)
  refpar[9:11] = ref.f
  Ntrue[t,] = ref.f
  No[t,] = tobit(mvtnorm::rmvnorm(1,ref.f,R))
}

## ---- Filter --------------------------------------------------------------------------
## disturbance2 applied to every row of X (n x 3), for disturbance type `type`
disturbance2.rows <- function(X, type){
  old.biomass <- X[, c(1,3), drop = FALSE]
  new.biomass <- old.biomass * rep(mu0[,type], each = nrow(X)) +
                 mvtnorm::rmvnorm(nrow(X), c(0,0), V0[,,type])
  new.biomass <- pmin(tobit(new.biomass), old.biomass)
  X[, c(1,3)] <- new.biomass
  X[, 2] <- X[, 2] + rowSums(old.biomass - new.biomass)*alloc.soil
  X
}

simulate <- function(x, t, label, mem) {
  tsel <- (t-1)*365 + 1:365
  out <- matrix(0, nrow(x), 3)
  for (n in seq_len(nrow(x))) {
    th <- THETA[mem[n], ]
    th[9:11] <- x[n, ]
    out[n, ] <- VSEM(th, PAR = PAR[met[mem[n]], tsel])[365, 2:4]
  }
  out
}

if (oracle) { THETA[] <- rep(refpar0, each = ne); met[] <- met[true] }

td <- which(D.true > 0)
filter_args <- list(
  Y = No, R = R, Xic = Na.orig, simulate = simulate, nMembers = ne,
  labels = c(1, 1, 1),                                    # one model for every group
  disturb = disturbance2.rows,
  dist.p = matrix(p.dist, 3, NT),
  procScale = rbind(c(1, 1, 1),                           # leaf, soil, stem
                    c(mu0[1,1]^2, 1, mu0[2,1]^2),
                    c(mu0[1,2]^2, 1, mu0[2,2]^2)),
  transition = matrix(1:3, 3, 3, byrow = TRUE),           # group after the step = scenario
  w0 = c(1, 1, 1), poolOf = 1:3, prior = list(type = "loguniform"),
  N = 1000, keepMix = td)

## tau learned: independent log-uniform priors, each axis spanning process sd 0.01 .. 3.2
or.grid <- list(tauMin = rep(0.1, 3), tauMax = rep(1e4, 3), G = rep(12, 3))
run.time <- system.time(or.fit <- do.call(ORJointUpdateGeneral, c(filter_args, list(grid = or.grid))))
## tau known (= 1/diag(Q))
or.known <- do.call(ORJointUpdateGeneral, c(filter_args, list(
  grid = list(tauMin = 1/diag(Q), tauMax = 1/diag(Q), G = rep(1, 3)))))

## ---- Figures ------------------------------------------------------------------------
yr <- 1:NT
pdf(file.path(outdir, "figures.pdf"), width = 7, height = 6)

## posterior probability of each disturbance type (cf. the Rmd's post.dist plot)
cols <- c("black", "#e41a1c", "#4daf4a")
matplot(yr, or.fit$pK, type = 'l', lty = 1, lwd = 2, col = cols, ylim = c(0, 1), xlab = "year",
        ylab = "posterior probability", main = "Disturbance type at each step")
matpoints(yr, or.known$pK, pch = 19, cex = 0.5, col = cols)
abline(v = td, lty = 2)
legend("right", legend = c(dist.names, "lines: tau learned", "points: tau known"),
       col = c(cols, NA, NA), lwd = c(2, 2, 2, NA, NA), bg = "white")

for (a in 1:3) plot_state(or.fit, a, Ntrue, No, main = pool.names[a], ylab = pool.names[a],
                          dist.years = td, known = or.known)

plot_process_sd(or.fit, sqrt(diag(Q)), pool.names)
plot_dropped_mass(or.fit, "Tau learned: dropped forecast components")

## forecast vs analysis at the disturbance year, leaf vs stem (cf. the Rmd's ensemble-adjustment plot)
xf <- rmix(or.fit$mixF[[as.character(td)]], 1000); xa <- rmix(or.fit$mixA[[as.character(td)]], 1000)
plot(xf[, 1], xf[, 3], cex = 0.5, xlim = range(c(xf[, 1], xa[, 1])), ylim = range(c(xf[, 3], xa[, 3])),
     xlab = pool.names[1], ylab = pool.names[3], main = sprintf("Forecast vs analysis, year %d", td))
points(xa[, 1], xa[, 3], pch = "+", col = adjustcolor("red", 0.5))
points(No[td, 1], No[td, 3], cex = 3, col = adjustcolor("purple", 0.75), pch = 18)
points(Ntrue[td, 1], Ntrue[td, 3], cex = 2, pch = 4, lwd = 2)
legend("topleft", legend = c("forecast", "analysis", "data", "truth"), pch = c(1, 3, 18, 4),
       col = c(1, 2, "purple", 1))

dev.off()

## ---- Save outputs -------------------------------------------------------------------
saveRDS(list(Ntrue = Ntrue, No = No, R = R, Q = Q, D.true = D.true, or.fit = or.fit, or.known = or.known),
        file.path(outdir, "results.rds"))

rmse <- function(x) sqrt(colMeans((x - Ntrue)^2))
fmt  <- function(x, nm = pool.names) paste(sprintf("%s %.3f", nm, x), collapse = ", ")
cat(sprintf("filter run time (tau learned): %.1fs\n", run.time[3]))
cat("analysis-median RMSE: tau learned:", fmt(rmse(or.fit$q[3, , ])), "| tau known:", fmt(rmse(or.known$q[3, , ])), "\n")
cat(sprintf("year %d posterior: %s (tau known %s)\n", td, fmt(or.fit$pK[td, ], dist.names), fmt(or.known$pK[td, ], dist.names)))
cat(sprintf("max P(disturbed) in other years: %.3f\n", max(1 - or.fit$pK[-td, 1])))
cat(sprintf("process sd, sqrt(E[1/tau]) at year %d: %s (truth %.3f each)\n", NT,
            fmt(sqrt(or.fit$procVar[NT, ])), sqrt(Q[1,1])))
cat(dropped_summary(or.fit, "tau learned: "), "\n")
