## Single site, multiple disturbance types with PFT switching, multi-pool (VSEM) synthetic
## experiment.
##
## Migrated from SingleSite_MultiDisturb_MultiPool_MultiPFT.Rmd: the same drivers, two-PFT
## parameter ensemble, spin-up (PFT 1) and "true" run, in which a removal in year 10
## switches the site from PFT 1 (tree-like, 10x stem residence time) to PFT 2. The JAGS
## multinomial analysis + multi.ens.adj is replaced by the generalised ORJointUpdate
## filter of or-joint-update-filter.R, and the process error of each pool (leaf, soil,
## stem) is learned on a 3-D tau grid.
##
## Filter set-up:
##   - scenarios k = 0 undisturbed, 1 removal (-> PFT 2), 2 defoliation (PFT unchanged)
##   - four history groups, each with a PFT that selects the VSEM parameters of its
##     forecast: U1 / U2 = undisturbed at the last step with PFT 1 / 2, D1 = defoliated
##     PFT 1, D2 = disturbed with PFT 2 (removed, or defoliated PFT 2). Every importance
##     sample is run under both PFTs, so the same samples serve every group
##   - disturbance2 keeps a fraction mu0 of leaf and stem, so their process variance is
##     scaled by mu0^2 under a disturbance; the soil keeps its full process variance
##   - each importance sample takes a random ensemble member's parameters and driver every
##     year, so any systematic gap between the true member and the ensemble shows up as
##     extra process error. oracle <- TRUE gives every member the true parameters and driver.
##
## Run from the repository root:  Rscript synthetic-experiments/singlesite-multidisturb-multipool-multipft.R

library(BayesianTools)
library(mvnfast)
source("disturbance.R")                              # tobit, disturbance2
source("synthetic-experiments/or-joint-update-filter.R")

oracle <- FALSE   ## TRUE: every ensemble member uses the true parameters and driver
outdir <- file.path("experiments/synthetic", if (oracle) "singlesite-multidisturb-multipool-multipft-oracle" else "singlesite-multidisturb-multipool-multipft")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
if (!exists("seed", inherits = FALSE)) seed <- 42   ## a runner may set seed first
set.seed(seed)

## ---- Set-up (as in the Rmd) ------------------------------------------------------------
ne = 200     ## number of ensemble members
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
pft1 <- refPars$best; pft1[7] = pft1[7]*10 ## increased stem residence time for tree-like PFT
refPars$upper[7] =  refPars$upper[7]*10
np = 50
sigmaP = diag((refPars$upper - refPars$lower)/10)^2 ## assume range is +/- 5 SD
epsilon <- mvtnorm::rmvnorm(np,rep(0,nrow(sigmaP)),sigmaP)
Theta <- list()
Theta[[1]] <- as.data.frame(t(pft1         + t(epsilon))) ## PFT 1
Theta[[2]] <- as.data.frame(t(refPars$best + t(epsilon))) ## PFT 2
colnames(Theta[[1]]) <- rownames(refPars)
colnames(Theta[[2]]) <- rownames(refPars)
for(i in 1:11){
  for(p in 1:2){
    Theta[[p]][Theta[[p]][,i] < refPars$lower[i] , i] <- refPars$lower[i]
    Theta[[p]][Theta[[p]][,i] > refPars$upper[i] , i] <- refPars$upper[i]
  }
}

## initial conditions (25-year spin-up with PFT 1)
met    = sample.int(nd,ne,replace=TRUE)
parIdx = sample.int(np,ne,replace=TRUE)
Na.orig <- matrix(NA,ne,3)
for(i in 1:ne){
  Na.orig[i,] = VSEM(as.numeric(Theta[[1]][parIdx[i],]),PAR = PAR[met[i],])[nt,2:4]
}

Q = diag(c(0.1,0.1,0.1)) ## process error
R = diag(c(.25,10,.5))   ## observation error
THETA.pft <- list(as.matrix(Theta[[1]][parIdx,]),  ## ensemble parameters of each PFT
                  as.matrix(Theta[[2]][parIdx,]))

## disturbance: mean leaf & stem biomass after disturbance as % of original
mu0 <- matrix(c(0.1,0.05,          ## removal, switch to PFT2
                0.1,0.95),2,2)     ## defoliation, PFT1
V1  <- diag(c(0.15,0.20))^2        ## var in post disturbance leaf and stem biomass (one for each type)
V2  <- diag(c(0.15,0.01))^2
V0  <- abind::abind(V1,V2,along=3)
alloc.soil = 0.25                  ## fraction of removed C that goes in the soil
p.dist = c(0.9,0.05,0.05)          ## Multinomial transition probability; first is no change
dist.pft = c(0,2,0)                ## PFT to assign post-disturbance; 0 = no change
dist.names = c("undisturbed","removal","defoliation")
pool.names = c("Cleaf","Csoil","Cstem")

## pseudodata
No = matrix(NA,NT,3)
Ntrue = matrix(NA,NT,3)
true = 42                              ## reference ensemble member
D.true = rep(0,NT)
D.true[10] = 1                         ## prescribe a removal
refpar = as.numeric(Theta[[1]][parIdx[true],])
refpar[9:11] = Na.orig[true,]
for(t in 1:NT){
  tsel <- (t-1)*365 + 1:365
  ref.mu = VSEM(refpar,PAR[met[true],tsel])[365,2:4]
  ref.f  = mvtnorm::rmvnorm(1,ref.mu,Q)
  if(D.true[t]){
    ref.f = disturbance2(ref.f,type=D.true[t])
    refpar = as.numeric(Theta[[2]][parIdx[true],]) ## switch PFT parameters
  }
  ref.f = tobit(ref.f)
  refpar[9:11] = ref.f
  Ntrue[t,] = ref.f
  No[t,] = tobit(mvtnorm::rmvnorm(1,ref.f,R))
}
pft.true <- ifelse(seq_len(NT) >= which(D.true > 0), 2, 1)

## ---- Filter --------------------------------------------------------------------------
if (oracle) {
  for (p in 1:2) THETA.pft[[p]][] <- rep(as.numeric(Theta[[p]][parIdx[true],]), each = ne)
  met[] <- met[true]
}

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

## forecast with the parameters of PFT `label`
simulate <- function(x, t, label, mem) {
  tsel <- (t-1)*365 + 1:365
  TH  <- THETA.pft[[label]]
  out <- matrix(0, nrow(x), 3)
  for (n in seq_len(nrow(x))) {
    th <- TH[mem[n], ]
    th[9:11] <- x[n, ]
    out[n, ] <- VSEM(th, PAR = PAR[met[mem[n]], tsel])[365, 2:4]
  }
  out
}

## history groups: U1, U2, D1, D2 (see header)
group.names <- c("U1", "U2", "D1", "D2")
group.pft   <- c(1, 2, 1, 2)
transition  <- cbind(c(1, 2, 1, 2),      # undisturbed: stay undisturbed, same PFT
                     c(4, 4, 4, 4),      # removal: PFT 2
                     c(3, 4, 3, 4))      # defoliation: same PFT

td <- which(D.true > 0)
filter_args <- list(
  Y = No, R = R, Xic = Na.orig, simulate = simulate, nMembers = ne,
  labels = group.pft, disturb = disturbance2.rows,
  dist.p = matrix(p.dist, 3, NT),
  procScale = rbind(c(1, 1, 1),                           # leaf, soil, stem
                    c(mu0[1,1]^2, 1, mu0[2,1]^2),
                    c(mu0[1,2]^2, 1, mu0[2,2]^2)),
  transition = transition,
  w0 = c(1, 0, 1, 0),                                     # every member starts as PFT 1
  poolOf = 1:3, prior = list(type = "loguniform"), N = 1000, keepMix = td)

## tau learned: independent log-uniform priors, each axis spanning process sd 0.01 .. 3.2
or.grid <- list(tauMin = rep(0.1, 3), tauMax = rep(1e4, 3), G = rep(12, 3))
run.time <- system.time(or.fit <- do.call(ORJointUpdateGeneral, c(filter_args, list(grid = or.grid))))
## tau known (= 1/diag(Q))
or.known <- do.call(ORJointUpdateGeneral, c(filter_args, list(
  grid = list(tauMin = 1/diag(Q), tauMax = 1/diag(Q), G = rep(1, 3)))))
pPFT2 <- function(fit) rowSums(fit$pGroup[, group.pft == 2, drop = FALSE])

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

## posterior probability of PFT 2 (cf. the Rmd's post.pft plot)
plot(yr, pPFT2(or.fit), type = 'l', lwd = 2, col = "#2a78d6", ylim = c(0, 1), xlab = "year",
     ylab = "P(PFT 2)", main = "PFT after each step")
points(yr, pPFT2(or.known), pch = 19, cex = 0.6)
lines(yr, pft.true - 1, lty = 2)
legend("right", legend = c("tau learned", "tau known", "truth"), lwd = c(2, NA, 1),
       lty = c(1, NA, 2), pch = c(NA, 19, NA), col = c("#2a78d6", 1, 1), bg = "white")

for (a in 1:3) plot_state(or.fit, a, Ntrue, No, main = pool.names[a], ylab = pool.names[a],
                          dist.years = td, known = or.known)

plot_process_sd(or.fit, sqrt(diag(Q)), pool.names)
plot_dropped_mass(or.fit, "Tau learned: dropped forecast components")

## forecast vs analysis at the disturbance year, leaf vs stem
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
saveRDS(list(Ntrue = Ntrue, No = No, R = R, Q = Q, D.true = D.true, pft.true = pft.true,
             or.fit = or.fit, or.known = or.known),
        file.path(outdir, "results.rds"))

rmse <- function(x) sqrt(colMeans((x - Ntrue)^2))
fmt  <- function(x, nm = pool.names) paste(sprintf("%s %.3f", nm, x), collapse = ", ")
cat(sprintf("filter run time (tau learned): %.1fs\n", run.time[3]))
cat("analysis-median RMSE: tau learned:", fmt(rmse(or.fit$q[3, , ])), "| tau known:", fmt(rmse(or.known$q[3, , ])), "\n")
cat(sprintf("year %d posterior: %s (tau known %s)\n", td, fmt(or.fit$pK[td, ], dist.names), fmt(or.known$pK[td, ], dist.names)))
cat(sprintf("max P(disturbed) in other years: %.3f\n", max(1 - or.fit$pK[-td, 1])))
cat(sprintf("P(PFT 2): before year %d max %.3f, from year %d min %.3f (tau known: %.3f, %.3f)\n", td,
            max(pPFT2(or.fit)[1:(td-1)]), td, min(pPFT2(or.fit)[td:NT]),
            max(pPFT2(or.known)[1:(td-1)]), min(pPFT2(or.known)[td:NT])))
cat(sprintf("process sd, sqrt(E[1/tau]) at year %d: %s (truth %.3f each)\n", NT,
            fmt(sqrt(or.fit$procVar[NT, ])), sqrt(Q[1,1])))
cat(dropped_summary(or.fit, "tau learned: "), "\n")
