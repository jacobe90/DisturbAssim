## Multi-site (two sites), single disturbance type, multi-pool (VSEM) synthetic experiment.
##
## Migrated from TwoSite_SingleDisturb_MultiPool.Rmd: the same drivers, parameter
## ensemble, spin-up and "true" run (site 1 disturbed in year 10; both sites share the
## true member's parameters and driver), with the JAGS multinomial analysis +
## multi.ens.adj replaced by the generalised ORJointUpdate filter of
## or-joint-update-filter.R.
##
## Filter set-up:
##   - the state is (Cleaf, Csoil, Cstem) at site 1 then site 2 (d = 6); all six are observed
##   - scenarios k = 0 none, 1 site 1 only, 2 site 2 only, 3 both. Their prior comes from
##     the Rmd's correlated-Bernoulli draw (rcbern, p = 0.025, rho = 0.5), evaluated exactly
##     with the bivariate normal CDF instead of by ensemble counts
##   - four history groups (the scenario at the last step)
##   - the process error is diag(1/tau_leaf, 1/tau_soil, 1/tau_stem) at both sites: one
##     3-D tau grid shared by the two sites, as the truth uses the same Q at both
##   - disturbance() resets leaf and stem, so they get no process variance under a
##     disturbance at that site; the soil keeps its full process variance
##   - each importance sample takes one random ensemble member's parameters and driver for
##     both sites every year, so a systematic gap between the true member and the
##     ensemble shows up as extra process error. oracle <- TRUE gives every member the
##     true parameters and driver.
## The Rmd's EnKF blocks are disabled (if(FALSE)) or broken for two sites, so they are
## not ported.
##
## Run from the repository root:  Rscript synthetic-experiments/multisite-singledisturb-multipool.R

library(BayesianTools)
library(mvnfast)
source("disturbance.R")                              # tobit, disturbance
source("synthetic-experiments/or-joint-update-filter.R")

oracle <- FALSE   ## TRUE: every ensemble member uses the true parameters and driver
outdir <- file.path("experiments/synthetic", if (oracle) "multisite-singledisturb-multipool-oracle" else "multisite-singledisturb-multipool")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
set.seed(42)

## ---- Set-up (as in the Rmd) ------------------------------------------------------------
NT = 25
nt = 365*NT

## drivers
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

## initial conditions (spin up)
ne = 1000
met    = sample.int(nd,ne,replace=TRUE)
parIdx = sample.int(np,ne,replace=TRUE)
Na.orig <- matrix(NA,ne,3)
for(i in 1:ne){
  Na.orig[i,] = VSEM(as.numeric(Theta[parIdx[i],]),PAR = PAR[met[i],])[nt,2:4]
}

nsite = 2
Q = diag(c(0.1,0.1,0.1))            ## process error (each site)
R = diag(c(.5,10,1))                ## observation error
R = as.matrix(Matrix::bdiag(R,R))   ## block diagonal over sites
THETA <- as.matrix(Theta[parIdx,])  ## ensemble parameters (shared by the sites)
p.dist = 0.025
pool.names = c("Cleaf","Csoil","Cstem")
state.names = paste(rep(pool.names, nsite), rep(paste("site", 1:nsite), each = 3))
scen.names = c("none", "site 1", "site 2", "both")

## pseudodata
No = matrix(NA,NT,3*nsite)
Ntrue = matrix(NA,NT,3*nsite)
true = 42                                    ## reference ensemble member
D.true = matrix(0,nsite,NT)
D.true[1,10] = 1                             ## disturbance at site 1 in year 10
refpar = matrix(as.numeric(Theta[parIdx[true],]),nsite,ncol(Theta),byrow = TRUE)
refpar[,9:11] = Na.orig[true+(1:nsite)-1,]   ## different initial biomass at each site
ref.f = list()
for(t in 1:NT){
  tsel <- (t-1)*365 + 1:365
  for(s in 1:nsite){
    ref.mu = VSEM(refpar[s,],PAR[met[true],tsel])[365,2:4]
    ref.f[[s]]  = mvtnorm::rmvnorm(1,ref.mu,Q)
    if(D.true[s,t]){
      ref.f[[s]] = disturbance(ref.f[[s]])
    }
    ref.f[[s]] = tobit(ref.f[[s]])
    refpar[s,9:11] = ref.f[[s]]
  }
  Ntrue[t,] = unlist(ref.f)
  No[t,] = tobit(mvtnorm::rmvnorm(1,Ntrue[t,],R))
}
## initial ensemble: site 2 gets the spin-up ensemble resampled without replacement
Xic <- cbind(Na.orig, Na.orig[sample.int(ne,ne),])

## ---- Filter --------------------------------------------------------------------------
if (oracle) { THETA[] <- rep(as.numeric(Theta[parIdx[true],]), each = ne); met[] <- met[true] }

## scenario prior: rcbern draws x ~ N(0, SIGMA) with SIGMA = rho^|i-j| / (1 - rho^2) and
## flags site s disturbed when pnorm(x_s) > 1 - p, i.e. x_s > qnorm(1 - p)
rho   <- 0.5
SIGMA <- 1/(1-rho^2) * rho^as.matrix(dist(1:nsite, diag = TRUE, upper = TRUE))
cut   <- qnorm(1 - p.dist)
p.both <- mvtnorm::pmvnorm(lower = c(cut, cut), upper = c(Inf, Inf), sigma = SIGMA)[1]
p.one  <- 1 - pnorm(cut, sd = sqrt(SIGMA[1, 1])) - p.both          # one site only (each)
p.scen <- c(1 - 2 * p.one - p.both, p.one, p.one, p.both)

## disturbance (absolute leaf / stem draw, from disturbance.R) applied to every row of X
## (n x 3): defaults mu0 = c(1, 1), V0 = diag(c(0.15, 0.25))^2, alloc.soil = 0.25
mu0 <- c(1,1); V0 <- diag(c(0.15,0.25))^2; alloc.soil <- 0.25
disturbance.rows <- function(X){
  old.biomass <- X[, c(1,3), drop = FALSE]
  new.biomass <- pmin(tobit(mvtnorm::rmvnorm(nrow(X), mu0, V0)), old.biomass)
  X[, c(1,3)] <- new.biomass
  X[, 2] <- X[, 2] + rowSums(old.biomass - new.biomass)*alloc.soil
  X
}
site.cols <- lapply(1:nsite, function(s) 3*(s-1) + 1:3)
scen.sites <- list(integer(0), 1, 2, 1:2)            # disturbed sites of scenario k = 0..3
disturb <- function(x, k) {
  for (s in scen.sites[[k + 1]]) x[, site.cols[[s]]] <- disturbance.rows(x[, site.cols[[s]], drop = FALSE])
  x
}

## both sites forecast with the same member's parameters and driver
simulate <- function(x, t, label, mem) {
  tsel <- (t-1)*365 + 1:365
  out <- matrix(0, nrow(x), 3*nsite)
  for (n in seq_len(nrow(x))) {
    th <- THETA[mem[n], ]
    for (s in 1:nsite) {
      th[9:11] <- x[n, site.cols[[s]]]
      out[n, site.cols[[s]]] <- VSEM(th, PAR = PAR[met[mem[n]], tsel])[365, 2:4]
    }
  }
  out
}

## process variance under each scenario: none for leaf / stem of a disturbed site
procScale <- t(sapply(scen.sites, function(ss) {
  v <- rep(1, 3*nsite)
  for (s in ss) v[site.cols[[s]][c(1, 3)]] <- 0
  v
}))

td <- which(D.true[1, ] > 0)
filter_args <- list(
  Y = No, R = R, Xic = Xic, simulate = simulate, nMembers = ne,
  labels = rep(1, 4), disturb = disturb,
  dist.p = matrix(p.scen, 4, NT),
  procScale = procScale,
  transition = matrix(1:4, 4, 4, byrow = TRUE),           # group after the step = scenario
  w0 = rep(1, 4), poolOf = rep(1:3, nsite), prior = list(type = "loguniform"),
  N = 1000, keepMix = td)

## tau learned: independent log-uniform priors, each axis spanning process sd 0.01 .. 3.2
or.grid <- list(tauMin = rep(0.1, 3), tauMax = rep(1e4, 3), G = rep(12, 3))
run.time <- system.time(or.fit <- do.call(ORJointUpdateGeneral, c(filter_args, list(grid = or.grid))))
## tau known (= 1/diag(Q))
or.known <- do.call(ORJointUpdateGeneral, c(filter_args, list(
  grid = list(tauMin = 1/diag(Q), tauMax = 1/diag(Q), G = rep(1, 3)))))
## P(site s disturbed at each step)
pSite <- function(fit) sapply(1:nsite, function(s) rowSums(fit$pK[, sapply(scen.sites, function(ss) s %in% ss), drop = FALSE]))

## ---- Figures ------------------------------------------------------------------------
yr <- 1:NT
pdf(file.path(outdir, "figures.pdf"), width = 7, height = 6)

cols <- c("#e41a1c", "#377eb8")
matplot(yr, pSite(or.fit), type = 'h', lwd = c(6, 3), lty = 1, col = cols, ylim = c(0, 1), xlab = "year",
        ylab = "P(disturbed)", main = "Posterior probability of disturbance by site")
matpoints(yr, pSite(or.known), pch = 19, cex = 0.6, col = cols)
abline(v = td, lty = 2); abline(h = 1 - p.scen[1], lty = 3)
legend("topright", legend = c("site 1", "site 2", "points: tau known", "prior P(any)"),
       col = c(cols, 1, 1), lwd = c(6, 3, NA, 1), lty = c(1, 1, NA, 3), pch = c(NA, NA, 19, NA), bg = "white")

for (a in 1:(3*nsite)) plot_state(or.fit, a, Ntrue, No, main = state.names[a], ylab = state.names[a],
                                  dist.years = td, known = or.known)

plot_process_sd(or.fit, sqrt(diag(Q)), pool.names)

## forecast vs analysis at the disturbance year, site 1 leaf vs stem (cf. the Rmd's ensemble-adjustment plot)
xf <- rmix(or.fit$mixF[[as.character(td)]], 1000); xa <- rmix(or.fit$mixA[[as.character(td)]], 1000)
plot(xf[, 1], xf[, 3], cex = 0.5, col = adjustcolor(1, 0.4), xlim = range(c(xf[, 1], xa[, 1])),
     ylim = range(c(xf[, 3], xa[, 3])), xlab = state.names[1], ylab = state.names[3],
     main = sprintf("Site 1 forecast vs analysis, year %d", td))
points(xa[, 1], xa[, 3], pch = "+", col = adjustcolor("red", 0.4))
points(No[td, 1], No[td, 3], cex = 3, col = adjustcolor("purple", 0.75), pch = 18)
points(Ntrue[td, 1], Ntrue[td, 3], cex = 2, pch = 4, lwd = 2)
legend("topleft", legend = c("forecast", "analysis", "data", "truth"), pch = c(1, 3, 18, 4),
       col = c(1, 2, "purple", 1))

dev.off()

## ---- Save outputs -------------------------------------------------------------------
saveRDS(list(Ntrue = Ntrue, No = No, R = R, Q = Q, D.true = D.true, p.scen = p.scen,
             or.fit = or.fit, or.known = or.known),
        file.path(outdir, "results.rds"))

rmse <- function(x) sqrt(colMeans((x - Ntrue)^2))
fmt  <- function(x, nm) paste(sprintf("%s %.3f", nm, x), collapse = ", ")
cat(sprintf("scenario prior: %s\n", fmt(p.scen, scen.names)))
cat(sprintf("filter run time (tau learned): %.1fs\n", run.time[3]))
cat("analysis-median RMSE, tau learned:", fmt(rmse(or.fit$q[3, , ]), state.names), "\n")
cat("analysis-median RMSE, tau known:  ", fmt(rmse(or.known$q[3, , ]), state.names), "\n")
cat(sprintf("year %d posterior: %s (tau known %s)\n", td, fmt(or.fit$pK[td, ], scen.names), fmt(or.known$pK[td, ], scen.names)))
cat(sprintf("max P(site 1 disturbed) in other years %.3f; max P(site 2 disturbed) %.3f\n",
            max(pSite(or.fit)[-td, 1]), max(pSite(or.fit)[, 2])))
cat(sprintf("process sd, sqrt(E[1/tau]) at year %d: %s (truth %.3f each)\n", NT,
            fmt(sqrt(or.fit$procVar[NT, ]), pool.names), sqrt(Q[1,1])))
