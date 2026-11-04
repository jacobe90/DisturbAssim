## Single site, single disturbance, multi-pool (VSEM) synthetic experiment.
##
## Migrated from SingleSite_SingleDisturb_MultiPool.Rmd. Reproduces
##   1. the drivers, parameter ensemble, spin-up and "true" run / observations,
##   2. the EnKF without disturbance and the EnKF with disturbance,
##   3. the multinomial (disturbed / undisturbed mixture) filter -- here the joint
##      disturbance + process-precision filter of OR_joint_update2.Rmd (ORJointUpdate)
##      in place of the JAGS analysis + multi.ens.adj,
## and the comparison figures.
##
## Differences from ORJointUpdate, to match the multi-pool experiment:
##   - every pool (Cleaf, Csoil, Cstem) is observed, with a full observation covariance R
##     (ORJointUpdate observes only the wood pool); NA entries are skipped
##   - each importance sample is paired with a random ensemble member's VSEM parameters
##     and PAR driver, so the forecast moments include parameter and driver uncertainty
##   - disturbance removes a percentage of the existing biomass, as disturbance2 in
##     SingleSite_MultiDisturb_MultiPool.Rmd (its harvest type): leaf and stem become
##     N(mu0 * old, V0), capped at the old values, with alloc.soil of the removed C added
##     to the soil. (SingleSite_SingleDisturb_MultiPool.Rmd drew leaf and stem from an
##     absolute N(mu0, V0) instead.)
##   - the process error is diag(1/tau_leaf, 1/tau_soil, 1/tau_stem), all three learned on
##     a 3-D log-spaced product grid (ORJointUpdate learns only the wood pool's tau). Under
##     a disturbance, leaf and stem keep a fraction mu0 of a state that already carries
##     process error, so their process variance is scaled by mu0^2; the soil keeps its full
##     process variance
##   - with thousands of grid points, the per-component 3 x 3 algebra (densities, Kalman
##     updates, proposal draws) is done in closed form across all components at once, and
##     grid points whose posterior weight becomes negligible are dropped
##
## Run from the repository root:  Rscript synthetic-experiments/singlesite-singledisturb-multipool.R
## Figures and results go to `outdir`.

library(BayesianTools)
library(mvnfast)

outdir <- "experiments/synthetic/singlesite-singledisturb-multipool"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
set.seed(42)

## set values of state less than zero to be 0
tobit <- function(x){
  x[x < 0] <- 0
  x
}

### Constants to be used for disturbance function: the harvest of
### SingleSite_MultiDisturb_MultiPool.Rmd (disturbance2, type 1)
mu0 <- c(0.1,0.1)               ## mean leaf & stem biomass after disturbance as % of original
V0  <- diag(c(0.15,0.15))^2     ## var in post disturbance leaf and stem biomass
alloc.soil = 0.25               ## fraction of removed C that goes in the soil

## x = (Cleaf, Csoil, Cstem); keeps a fraction ~ mu0 of leaf and stem (they can only
## decrease) and moves alloc.soil of the removed biomass to the soil (as disturbance2)
disturbance <- function(x){
  old.biomass <- x[c(1,3)]
  new.biomass <- tobit(mvtnorm::rmvnorm(1,mu0*old.biomass,V0)) ## draw disturbed leaf and stem
  check.max <- which(new.biomass > old.biomass)
  if(length(check.max)>0) new.biomass[check.max] = old.biomass[check.max]
  residual = sum(old.biomass-new.biomass)
  x[c(1,3)] <- new.biomass
  x[2] <- x[2] + residual*alloc.soil
  return(x)
}

## disturbance() applied to every row of X (n x 3)
disturbance.rows <- function(X, mu0, V0, alloc.soil){
  old.biomass <- X[, c(1,3), drop = FALSE]
  new.biomass <- old.biomass * rep(mu0, each = nrow(X)) + mvtnorm::rmvnorm(nrow(X), c(0,0), V0)
  new.biomass <- pmin(tobit(new.biomass), old.biomass)
  X[, c(1,3)] <- new.biomass
  X[, 2] <- X[, 2] + rowSums(old.biomass - new.biomass)*alloc.soil
  X
}

## ---- Drivers, parameters, initial conditions (as in the Rmd) -------------------------
NYears = 25
ndays = 365*NYears
ndrivers = 31 ## number of driver ensembles
PAR_driverdata <- matrix(NA,ndrivers,ndays)
for(i in 1:ndrivers){
  PAR_driverdata[i,] <- VSEMcreatePAR(1:ndays)
}

refPars <- VSEMgetDefaults()
refParStateIndices <- 9:11 # Cv, Cs, Cr in the VSEM parameter vector
refPars$lower[refParStateIndices] = 1
nparamsamples = 50
sigmaP = diag((refPars$upper - refPars$lower)/10)^2 ## assume range is +/- 5 SD
Theta_reference <- as.data.frame(mvtnorm::rmvnorm(nparamsamples,refPars$best,sigmaP))
colnames(Theta_reference) <- rownames(refPars)
for(i in 1:11){ # truncate to the parameter bounds
  Theta_reference[Theta_reference[,i] < refPars$lower[i] , i] <- refPars$lower[i]
  Theta_reference[Theta_reference[,i] > refPars$upper[i] , i] <- refPars$upper[i]
}

## spin up
n_ensemble_members = 200
met    = sample.int(ndrivers,n_ensemble_members,replace=TRUE)       # driver of each member
parIdx = sample.int(nparamsamples,n_ensemble_members,replace=TRUE)  # parameters of each member
vsem_state_indices <- 2:4 # Cv, Cs, Cr columns of the VSEM output
Ensemble_initial.orig <- matrix(NA,n_ensemble_members,3)
for(i in 1:n_ensemble_members){
  Ensemble_initial.orig[i,] = VSEM(as.numeric(Theta_reference[parIdx[i],]),PAR = PAR_driverdata[met[i],])[ndays,vsem_state_indices]
}
THETA <- as.matrix(Theta_reference[parIdx,]) ## ensemble parameters [J, 11]

Q = diag(c(0.1,0.1,0.1))  ## process error covariance
R = diag(c(.25,10,1))     ## observation error covariance
p.dist = 0.15
pool.names = c("Cleaf","Csoil","Cstem")

## ---- Pseudodata -------------------------------------------------------------------
Observed_states = matrix(NA,NYears,3)
True_states = matrix(NA,NYears,3)
ref_idx = 42                              ## reference ensemble member
disturbance_at = rep(0,NYears)
disturbance_at[10] = 1                    ## prescribe a disturbance
refpar = as.numeric(Theta_reference[parIdx[ref_idx],])
refpar[refParStateIndices] = Ensemble_initial.orig[ref_idx,]
for(year in 1:NYears){
  sel_days <- (year-1)*365 + 1:365
  simulated_state = VSEM(refpar,PAR_driverdata[met[ref_idx],sel_days])[365,vsem_state_indices]
  state_with_process_error = mvtnorm::rmvnorm(1,simulated_state,Q)
  if(disturbance_at[year]){
    state_with_process_error = disturbance(state_with_process_error)
  }
  state_with_process_error = tobit(state_with_process_error)
  refpar[9:11] = state_with_process_error
  True_states[year,] = state_with_process_error
  Observed_states[year,] = tobit(mvtnorm::rmvnorm(1,state_with_process_error,R))
}

## ---- EnKF, with and without disturbance -------------------------------------------
## Forecast the ensemble one year, optionally disturb members at rate p.dist, then a
## Kalman analysis (H = I) and the PEcAn ensemble adjustment. Returns per-year forecast
## and analysis ensembles.
run_enkf <- function(disturb) {
  X <- Ensemble_initial.orig
  TH <- THETA
  Xf.all <- Xa.all <- vector("list", NYears)
  for(t in 1:NYears){
    tsel <- (t-1)*365 + 1:365
    TH[,9:11] <- X
    Xf <- matrix(NA, n_ensemble_members, 3)
    for(i in 1:n_ensemble_members){
      Xf[i,] = VSEM(TH[i,],PAR = PAR_driverdata[met[i],tsel])[365,vsem_state_indices]
    }
    Xf <- Xf + mvtnorm::rmvnorm(n_ensemble_members,rep(0,3),Q)  ## process error
    if(disturb){
      dd <- which(as.logical(rbinom(n_ensemble_members,1,p.dist)))
      if(length(dd)) Xf[dd,] <- disturbance.rows(Xf[dd,,drop=FALSE],mu0,V0,alloc.soil)
      Xf <- tobit(Xf)
    }
    muf = colMeans(Xf)
    Pf  = var(Xf)

    ## Analysis
    I = H <- diag(rep(1,3))
    K <- Pf %*% t(H) %*% solve(H%*%Pf%*%t(H) + R)
    mua <- muf + K%*%(Observed_states[t,] - muf)
    Pa <- (I-K %*% H)%*%Pf
    X <- as.matrix(PEcAnAssimSequential::adj.ens(Pf = Pf,X = Xf,mu.f = muf,mu.a = mua,Pa = Pa))
    Xf.all[[t]] <- Xf
    Xa.all[[t]] <- X
  }
  list(Xf = Xf.all, Xa = Xa.all)
}
enkf0 <- run_enkf(disturb = FALSE)
enkf  <- run_enkf(disturb = TRUE)

## ---- Multinomial filter: ORJointUpdate ported to multi-pool VSEM -----------------------

# Quantiles of a 1-d Gaussian mixture sum_i w_i N(mu_i, sd_i^2), by root-finding on its
# CDF. Components with negligible weight are dropped first.
mix_quantile <- function(p, w, mu, sd) {
  keep <- w > 1e-12 * sum(w)
  w <- w[keep] / sum(w[keep]); mu <- mu[keep]; sd <- sd[keep]
  lo <- min(mu - 8 * sd); hi <- max(mu + 8 * sd)
  sapply(p, function(pp) uniroot(function(x) sum(w * pnorm(x, mu, sd)) - pp,
                                 c(lo, hi), tol = 1e-8)$root)
}

# Initial tau weights on one grid axis: prior mass of each grid cell (see OR_joint_update2.Rmd)
cell_prior_mass <- function(tauGrid, prior) {
  G   <- length(tauGrid)
  if (G == 1) return(1)
  mid <- (tauGrid[-G] + tauGrid[-1]) / 2
  lo  <- c(tauGrid[1], mid)
  hi  <- c(mid, tauGrid[G])
  cdf <- switch(prior$type,
    gamma      = function(x) pgamma(x, shape = prior$shape, rate = prior$rate),
    loguniform = function(x) log(x))
  m <- cdf(hi) - cdf(lo)
  m / sum(m)
}

## ---- Batched small-matrix algebra ----------------------------------------------------
## The 3-D tau grid gives thousands of components per step, so every per-component
## 3 x 3 operation is done in closed form across all components at once. A batch of
## d x d matrices is stored as an n x d^2 matrix, column-major: entry (a, b) is column
## a + d (b - 1).

# inverse and determinant of a batch of symmetric d x d matrices, d <= 3
inv_sym <- function(S, d) {
  if (d == 1) return(list(inv = 1 / S, det = S[, 1]))
  if (d == 2) {
    det <- S[, 1] * S[, 4] - S[, 2]^2
    return(list(inv = cbind(S[, 4], -S[, 2], -S[, 2], S[, 1]) / det, det = det))
  }
  a <- S[, 1]; b <- S[, 4]; c <- S[, 7]; e <- S[, 5]; f <- S[, 8]; i <- S[, 9]
  A11 <- e * i - f^2; A12 <- c * f - b * i; A13 <- b * f - c * e
  A22 <- a * i - c^2; A23 <- b * c - a * f; A33 <- a * e - b^2
  det <- a * A11 + b * A12 + c * A13
  list(inv = cbind(A11, A12, A13, A12, A22, A23, A13, A23, A33, deparse.level = 0) / det, det = det)
}

# lower Cholesky factor of a batch of symmetric positive-definite 3 x 3 matrices
chol3 <- function(S) {
  L11 <- sqrt(S[, 1]); L21 <- S[, 2] / L11; L31 <- S[, 3] / L11
  L22 <- sqrt(S[, 5] - L21^2); L32 <- (S[, 6] - L21 * L31) / L22
  L33 <- sqrt(S[, 9] - L31^2 - L32^2)
  cbind(L11, L21, L31, 0, L22, L32, 0, 0, L33, deparse.level = 0)
}

# P x N log densities of the samples X (N x 3) under each Gaussian component
# (mu: P x 3, sig: P x 9). The quadratic form (x - m)' A (x - m) = x'Ax - 2 (Am)'x + m'Am
# is expanded so that every (component, sample) pair comes out of two matrix products.
log_dens_batch <- function(X, mu, sig) {
  ctr <- colMeans(X)                          # centre both on a common point first
  X   <- sweep(X, 2, ctr); mu <- sweep(mu, 2, ctr)
  iv  <- inv_sym(sig, 3); A <- iv$inv
  quadX <- cbind(X[, 1]^2, X[, 2]^2, X[, 3]^2, 2 * X[, 1] * X[, 2], 2 * X[, 1] * X[, 3], 2 * X[, 2] * X[, 3])
  Aq    <- A[, c(1, 5, 9, 4, 7, 8)]          # A11 A22 A33 A12 A13 A23
  Am    <- cbind(rowSums(A[, 1:3] * mu), rowSums(A[, 4:6] * mu), rowSums(A[, 7:9] * mu))
  quad  <- tcrossprod(Aq, quadX) - 2 * tcrossprod(Am, X) + rowSums(Am * mu)
  -0.5 * quad - 0.5 * log(iv$det) - 1.5 * log(2 * pi)
}

# Kalman update of a batch of components (mu: n x 3, sig: n x 9) with H = I on the
# observed pools o. Returns the analysis mu / sig and the log likelihood of y[o].
kalman_batch <- function(mu, sig, y, R, o) {
  n  <- nrow(mu); d <- length(o)
  So <- sig[, as.vector(outer(o, 3 * (o - 1), "+")), drop = FALSE] +
        rep(as.vector(R[o, o]), each = n)
  iv <- inv_sym(So, d)
  v  <- matrix(rep(y[o], each = n), n, d) - mu[, o, drop = FALSE]      # innovations
  K  <- array(0, c(n, 3, d))                                           # C[, o] So^-1
  for (a in 1:3) for (b in 1:d)
    K[, a, b] <- rowSums(sig[, a + 3 * (o - 1), drop = FALSE] * iv$inv[, (b - 1) * d + 1:d, drop = FALSE])
  muA <- mu; sigA <- sig
  for (a in 1:3) {
    Ka <- matrix(K[, a, ], n, d)
    muA[, a] <- mu[, a] + rowSums(Ka * v)
    for (b in 1:3) sigA[, a + 3 * (b - 1)] <- sig[, a + 3 * (b - 1)] - rowSums(Ka * sig[, o + 3 * (b - 1), drop = FALSE])
  }
  qf <- 0
  for (b in 1:d) for (c in 1:d) qf <- qf + v[, b] * iv$inv[, b + d * (c - 1)] * v[, c]
  list(mu = muA, sig = sigA, logL = -0.5 * qf - 0.5 * log(iv$det) - 0.5 * d * log(2 * pi))
}

# Moment-matching projection onto 2 Gaussians per tau grid point (from OR_joint_update2.Rmd).
# w: G x 2 x (D+1), mu: G x 2 x (D+1) x 3, sig: G x 2 x (D+1) x 3 x 3. Group 1 is the
# undisturbed components (k = 0), group 2 the disturbed ones (k >= 1).
project_grid <- function(w, mu, sig) {
  G <- dim(w)[1]; C <- prod(dim(w)[-1])        # components c = j + 2k, j fastest
  muC  <- array(mu,  c(G, C, 3))
  sigC <- array(sig, c(G, C, 3, 3))
  w1 <- w[, 1, 1] + w[, 2, 1]
  w2 <- apply(w[, , -1, drop = FALSE], 1, sum)  # not 1 - w1, which can round below 0
  grp <- list(seq_len(C) <= 2, seq_len(C) > 2)  # group 1: k = 0
  wg  <- list(w1, w2)
  out <- list(w = cbind(w1, w2, deparse.level = 0), mu = array(0, c(G, 2, 3)),
              sig = array(0, c(G, 2, 3, 3)))
  for (h in 1:2) {
    A <- matrix(w, G, C) / wg[[h]]              # alpha: within-group weights
    A[, !grp[[h]]] <- 0
    for (a in 1:3) out$mu[, h, a] <- rowSums(A * muC[, , a])
    for (a in 1:3) for (b in 1:3)
      out$sig[, h, a, b] <- rowSums(A * (muC[, , a] * muC[, , b] + sigC[, , a, b])) -
                            out$mu[, h, a] * out$mu[, h, b]
  }
  # a group with zero total weight has undefined moments (0/0 above); give it the other
  # group's moments so it stays finite -- its weight is 0, so it only enters the proposal
  z1 <- !(w1 > 0); z2 <- !(w2 > 0) & !z1
  out$mu[z1, 1, ]    <- out$mu[z1, 2, ];    out$sig[z1, 1, , ] <- out$sig[z1, 2, , ]
  out$mu[z2, 2, ]    <- out$mu[z2, 1, ];    out$sig[z2, 2, , ] <- out$sig[z2, 1, , ]
  out
}

##' Approximate joint disturbance update, multi-pool VSEM with an ensemble of parameters
##' and drivers, learning the process-error variance of every pool
##'
##' The process error is diag(1/tau_leaf, 1/tau_soil, 1/tau_stem). (tau_leaf, tau_soil,
##' tau_stem) is carried on a 3-D log-spaced product grid, flattened to Gt points
##' (leaf axis fastest); everything per grid point is as in ORJointUpdate.
##'
##' @param Y NT x 3 observations of (Cleaf, Csoil, Cstem); NA entries are not assimilated
##' @param R 3 x 3 observation-error covariance
##' @param PAR ndrivers x (NT*365) PAR drivers
##' @param met driver row of each ensemble member
##' @param THETA J x 11 VSEM parameters of each ensemble member (cols 9:11 overwritten)
##' @param Xic J x 3 initial-condition ensemble
##' @param dist disturbance list: p ((D+1) x NT prior), n (# classes), vals (per-class mu0/V0),
##'   alloc.soil
##' @param N importance-sample size
##' @param prior prior for each tau, independent across pools:
##'   list(type = "gamma", shape, rate) or list(type = "loguniform")
##' @param grid list(tauMin, tauMax, G), each a length-3 vector (leaf, soil, stem): range and
##'   number of points of each axis (tauMin == tauMax, G = 1 => that tau is known)
##' @param pruneLog grid points whose log posterior weight falls more than pruneLog below the
##'   maximum are dropped for the rest of the run (their weight is < exp(-pruneLog))
##' @return list with, marginalised over the tau posterior and the mixture components:
##'   mean (NT x 3), q / fq (5 x NT x 3 analysis / forecast quantiles at qs5), pDist,
##'   pit (NT x 3, marginal one-step-ahead PIT of each observed pool), wDrop,
##'   procVar (NT x 3) = E[1/tau | y_{1:t}] per pool, tauAxes, tauMarg (per-pool marginal
##'   posterior on its axis, G_a x NT), nActive (grid points kept at each step), and
##'   mixF / mixA (per-t forecast / analysis mixtures: w, mu, sig) for sampling
qs5 <- c(0.025,0.25,0.5,0.75,0.975)
ORJointUpdateMP <- function(Y, R, PAR, met, THETA, Xic, dist, N = 1000, prior, grid, pruneLog = 30) {
  NT <- nrow(Y)
  D  <- dist$n
  J  <- nrow(THETA)
  tauAxes <- lapply(1:3, function(a) exp(seq(log(grid$tauMin[a]), log(grid$tauMax[a]), length.out = grid$G[a])))
  gridIdx <- as.matrix(expand.grid(lapply(tauAxes, seq_along)))   # Gt x 3 axis indices
  tauMat  <- do.call(cbind, lapply(1:3, function(a) tauAxes[[a]][gridIdx[, a]]))  # Gt x 3 tau values
  logMass <- lapply(tauAxes, function(x) log(cell_prior_mass(x, prior)))
  qs      <- rowSums(do.call(cbind, lapply(1:3, function(a) logMass[[a]][gridIdx[, a]])))
  Gt      <- nrow(tauMat)

  active  <- seq_len(Gt)                       # grid points still carried
  prevW   <- array(0.5, dim = c(Gt, 2))
  prevMu  <- array(rep(colMeans(Xic), each = 2 * Gt), dim = c(Gt, 2, 3))
  prevSig <- array(rep(as.vector(cov(Xic)), each = 2 * Gt), dim = c(Gt, 2, 3, 3))

  hist <- list(
    mean = matrix(NA_real_, NT, 3), q = array(NA_real_, c(5, NT, 3)), fq = array(NA_real_, c(5, NT, 3)),
    pDist = rep(NA_real_, NT), pit = matrix(NA_real_, NT, 3), wDrop = rep(NA_real_, NT),
    procVar = matrix(NA_real_, NT, 3), tauAxes = tauAxes,
    tauMarg = lapply(tauAxes, function(x) matrix(NA_real_, length(x), NT)),
    nActive = rep(NA_integer_, NT), mixF = vector("list", NT), mixA = vector("list", NT)
  )
  ia <- rep(1:3, 3); ib <- rep(1:3, each = 3)   # column (a, b) of a flattened 3 x 3

  for (t in 1:NT) {
    G  <- length(active)
    tg <- tauMat[active, , drop = FALSE]
    P  <- 2L * G
    Cn <- P * (D + 1)                            # flattened components (g fastest, then j, k)
    compMu  <- array(prevMu,  c(P, 3))           # r = (j - 1) * G + g
    compSig <- array(prevSig, c(P, 9))
    hist$nActive[t] <- G

    # sample from the proposal: equal-weight mixture over every (g, j), all pools positive
    compL <- chol3(compSig)
    propSamples <- matrix(0, 0, 3); tries <- 0L
    while (nrow(propSamples) < N) {
      tries <- tries + 1L
      if (tries > 100L) stop(sprintf("t=%d: proposal rejection rate too high", t))
      ri <- sample.int(P, N, replace = TRUE)
      z  <- matrix(rnorm(3 * N), N, 3)
      L  <- compL[ri, , drop = FALSE]
      x  <- compMu[ri, , drop = FALSE] +
            cbind(L[, 1] * z[, 1], L[, 2] * z[, 1] + L[, 5] * z[, 2],
                  L[, 3] * z[, 1] + L[, 6] * z[, 2] + L[, 9] * z[, 3])
      propSamples <- rbind(propSamples, x[rowSums(x > 0) == 3, , drop = FALSE])
    }
    propSamples <- propSamples[1:N, , drop = FALSE]

    # 1. VSEM for each sample, with a random ensemble member's parameters and driver
    tsel <- (t-1)*365 + 1:365
    mem  <- sample.int(J, N, replace = TRUE)
    sim0 <- matrix(0, N, 3)
    for (n in seq_len(N)) {
      th <- THETA[mem[n], ]
      th[9:11] <- propSamples[n, ]
      sim0[n, ] <- VSEM(th, PAR = PAR[met[mem[n]], tsel])[365, 2:4]
    }

    # 2-3. importance weights of each component's samples, in log space. W is reused
    # for each stage to keep a single P x N matrix in memory.
    W       <- log_dens_batch(propSamples, compMu, compSig)      # P x N log densities
    logProp <- matrixStats::colLogSumExps(W) - log(P)           # log proposal density
    W       <- W - rep(logProp, each = P)                        # log importance ratios
    W       <- exp(W - matrixStats::rowMaxs(W))
    W       <- W / rowSums(W)
    cN      <- 1 - rowSums(W^2)

    # components whose weights collapse onto ~1 sample have no estimable forecast moments:
    # unweighted placeholder moments and zero likelihood below (see ORJointUpdate)
    degen <- cN < 1e-8
    W[degen, ] <- 1 / N
    cN[degen]  <- 1 - 1 / N
    dropF <- matrix(degen, G, 2)

    # 4. forecast moments for every (g, j, k), flattened to Cn x 3 / Cn x 9
    muF  <- matrix(0, Cn, 3)
    sigF <- matrix(0, Cn, 9)
    tauP <- rbind(tg, tg)                        # tau of each (g, j) row
    for (k in 0:D) {
      sim <- sim0
      if (k > 0) sim <- disturbance.rows(sim0, dist$vals[[k]]$mu0, dist$vals[[k]]$V0, dist$alloc.soil)
      MU  <- W %*% sim
      ctr <- colMeans(sim)
      Xc  <- sweep(sim, 2, ctr)
      MUc <- sweep(MU, 2, ctr)
      S   <- (W %*% (Xc[, ia] * Xc[, ib]) - MUc[, ia] * MUc[, ib]) / cN   # P x 9
      if (k == 0) {
        S[, c(1, 5, 9)] <- S[, c(1, 5, 9)] + 1 / tauP                   # process variance 1/tau
      } else {
        # the disturbance keeps a fraction mu0 of leaf and stem, which already carry process
        # error, so their variance is scaled by mu0^2; the soil keeps its full variance
        f2 <- c(dist$vals[[k]]$mu0[1]^2, 1, dist$vals[[k]]$mu0[2]^2)
        S[, c(1, 5, 9)] <- S[, c(1, 5, 9)] + sweep(1 / tauP, 2, f2, "*")
      }
      rows <- k * P + seq_len(P)
      muF[rows, ]  <- MU
      sigF[rows, ] <- S
    }

    beta   <- dist$p[, t]
    piPrev <- exp(qs - max(qs)); piPrev <- piPrev / sum(piPrev)          # pi_{t-1}(g)
    hist$wDrop[t] <- sum(piPrev * rowSums(prevW * dropF))
    wPrior <- outer(prevW, beta)                                         # G x 2 x (D+1)
    for (k in 0:D) wPrior[, , k + 1][dropF] <- 0

    # forecast mixture, marginalised over pi_{t-1}
    wF <- as.vector(wPrior * piPrev); wF <- wF / sum(wF)
    for (a in 1:3) hist$fq[, t, a] <- mix_quantile(qs5, wF, muF[, a], sqrt(sigF[, 4 * a - 3]))

    # 5. Kalman update of every component on the observed pools
    o    <- which(!is.na(Y[t, ]))
    muA  <- muF; sigA <- sigF; logL <- numeric(Cn)
    if (length(o)) {
      ka   <- kalman_batch(muF, sigF, Y[t, ], R, o)
      muA  <- ka$mu; sigA <- ka$sig; logL <- ka$logL
      for (a in o) hist$pit[t, a] <- sum(wF * pnorm(Y[t, a], muF[, a], sqrt(sigF[, 4 * a - 3] + R[a, a])))
    }

    # 6. tau posterior on the grid
    lwl  <- log(wPrior) + array(logL, dim(wPrior))
    mz   <- apply(lwl, 1, max)
    logz <- mz + log(apply(exp(lwl - mz), 1, sum))
    if (any(!is.finite(logz))) stop(sprintf("t=%d: zero or non-finite marginal likelihood at %d tau grid points", t, sum(!is.finite(logz))))
    qs <- qs + logz
    ws <- exp(qs - max(qs)); ws <- ws / sum(ws)
    hist$procVar[t, ] <- colSums(ws / tg)
    for (a in 1:3) {
      m <- numeric(length(tauAxes[[a]]))
      s <- rowsum(ws, gridIdx[active, a])
      m[as.integer(rownames(s))] <- s
      hist$tauMarg[[a]][, t] <- m
    }

    # analysis weights
    wA <- exp(lwl - logz)

    # summaries, marginalised over the tau posterior and the components
    wAll <- as.vector(wA * ws)
    hist$pDist[t] <- sum(wAll[-seq_len(P)])
    for (a in 1:3) {
      hist$mean[t, a] <- sum(wAll * muA[, a])
      hist$q[, t, a]  <- mix_quantile(qs5, wAll, muA[, a], sqrt(sigA[, 4 * a - 3]))
    }
    keepF <- wF > 1e-8; keepA <- wAll > 1e-8
    hist$mixF[[t]] <- list(w = wF[keepF], mu = muF[keepF, , drop = FALSE], sig = sigF[keepF, , drop = FALSE])
    hist$mixA[[t]] <- list(w = wAll[keepA], mu = muA[keepA, , drop = FALSE], sig = sigA[keepA, , drop = FALSE])

    # 7. project back onto 2 components per tau grid point
    out <- project_grid(wA, array(muA, c(G, 2, D + 1, 3)), array(sigA, c(G, 2, D + 1, 3, 3)))

    # drop grid points whose posterior weight has become negligible
    keep    <- qs >= max(qs) - pruneLog
    active  <- active[keep]
    qs      <- qs[keep]
    prevW   <- out$w[keep, , drop = FALSE]
    prevMu  <- out$mu[keep, , , drop = FALSE]
    prevSig <- out$sig[keep, , , , drop = FALSE]
  }
  hist
}

## draw n samples from a stored mixture (w, mu: n x 3, sig: n x 9)
rmix <- function(mix, n) {
  ci <- sample.int(length(mix$w), n, replace = TRUE, prob = mix$w)
  t(sapply(ci, function(c) mvnfast::rmvn(1, mix$mu[c, ], matrix(mix$sig[c, ], 3, 3))))
}

or.dist <- list(
  p    = rbind(rep(1 - p.dist, NYears), rep(p.dist, NYears)),
  n    = 1,
  vals = list(list(mu0 = mu0, V0 = V0)),
  alloc.soil = alloc.soil
)
## tau unknown: independent log-uniform priors on a 3-D grid, each axis spanning process
## sd 0.01 .. 3.2 (truth: sd 0.32 for every pool)
or.grid <- list(tauMin = rep(0.1, 3), tauMax = rep(1e4, 3), G = rep(16, 3))
or.fit <- ORJointUpdateMP(Y = Observed_states, R = R, PAR = PAR_driverdata, met = met,
                          THETA = THETA, Xic = Ensemble_initial.orig, dist = or.dist,
                          N = 1000, prior = list(type = "loguniform"), grid = or.grid)
## tau known (= 1/diag(Q)), the setting of the Rmd's multinomial filter
or.known <- ORJointUpdateMP(Y = Observed_states, R = R, PAR = PAR_driverdata, met = met,
                            THETA = THETA, Xic = Ensemble_initial.orig, dist = or.dist,
                            N = 1000, prior = list(type = "loguniform"),
                            grid = list(tauMin = 1/diag(Q), tauMax = 1/diag(Q), G = rep(1, 3)))

## ---- Figures ------------------------------------------------------------------------
ens_q <- function(Xlist) array(sapply(Xlist, function(X) apply(X, 2, quantile, qs5)), c(5, 3, NYears))
an0.stats <- ens_q(enkf0$Xa)   # 5 x pool x year
an.stats  <- ens_q(enkf$Xa)
yr <- 1:NYears
band <- function(x, lo, hi, col) polygon(c(x, rev(x)), c(lo, rev(hi)), border = NA, col = col)
cols <- c(enkf0 = "orange", enkf = "#33a02c", or = "blue")

pdf(file.path(outdir, "figures.pdf"), width = 7, height = 6)

## filter comparison per pool (95% bands, medians)
for (a in 1:3) {
  rng <- range(c(an0.stats[c(1,5), a, ], an.stats[c(1,5), a, ], or.fit$q[c(1,5), , a],
                 Observed_states[, a], True_states[, a]))
  plot(yr, Observed_states[, a], type = 'n', ylim = rng, xlab = "year", ylab = pool.names[a],
       main = sprintf("%s: analysis 95%% intervals", pool.names[a]))
  band(yr, an0.stats[1, a, ], an0.stats[5, a, ], adjustcolor(cols["enkf0"], 0.35))
  band(yr, an.stats[1, a, ],  an.stats[5, a, ],  adjustcolor(cols["enkf"], 0.35))
  band(yr, or.fit$q[1, , a],  or.fit$q[5, , a],  adjustcolor(cols["or"], 0.3))
  lines(yr, an0.stats[3, a, ], col = cols["enkf0"], lwd = 2)
  lines(yr, an.stats[3, a, ],  col = cols["enkf"], lwd = 2)
  lines(yr, or.fit$q[3, , a],  col = cols["or"], lwd = 2)
  lines(yr, True_states[, a], lty = 2, lwd = 2)
  points(yr, Observed_states[, a], pch = 18, col = 2)
  abline(v = which(disturbance_at == 1), col = "grey60")
  legend("topright", legend = c("EnKF, no disturb", "EnKF, disturb", "Multinomial", "truth", "data"),
         col = c(cols, 1, 2), lwd = c(5, 5, 5, 2, NA), lty = c(1, 1, 1, 2, NA), pch = c(NA, NA, NA, NA, 18),
         bg = "white")
}

## P(disturbed)
plot(yr, or.fit$pDist, type = 'h', lwd = 6, col = "#eb6834", ylim = c(0, 1), xlab = "year",
     ylab = "P(disturbed)", main = "Posterior probability of disturbance")
points(yr, or.known$pDist, pch = 19, cex = 0.7)
abline(h = p.dist, lty = 3)
abline(v = which(disturbance_at == 1), lty = 2)
legend("topright", legend = c("tau unknown", "tau known", "prior p.dist", "true disturbance"),
       lty = c(1, NA, 3, 2), lwd = c(6, NA, 1, 1), pch = c(NA, 19, NA, NA), col = c("#eb6834", 1, 1, 1))

## process-error sd posterior of each pool over time (marginals of the 3-D tau grid)
op <- par(mfrow = c(3, 1), mar = c(3.5, 4.5, 2, 1))
for (a in 1:3) {
  tauQ <- apply(or.fit$tauMarg[[a]], 2, function(w) {
    cdf <- cumsum(w); or.fit$tauAxes[[a]][sapply(c(0.025,0.5,0.975), function(p) which(cdf >= p)[1])]
  })
  sdQ <- 1/sqrt(tauQ)
  plot(yr, sdQ[2,], type = 'n', log = "y", ylim = range(1/sqrt(or.fit$tauAxes[[a]])), xlab = "year",
       ylab = "process sd", main = sprintf("%s process-error posterior", pool.names[a]))
  band(yr, sdQ[3,], sdQ[1,], adjustcolor("#2a78d6", 0.3))
  lines(yr, sdQ[2,], lwd = 2, col = "#2a78d6")
  lines(yr, sqrt(or.fit$procVar[, a]), lty = 2, col = "#2a78d6")
  abline(h = sqrt(Q[a,a]), lwd = 2)
  if (a == 1) legend("topright", legend = c("posterior median", "95% interval", "sqrt(E[1/tau])", "truth"),
                     col = c("#2a78d6", adjustcolor("#2a78d6", 0.3), "#2a78d6", 1), lwd = c(2, 10, 1, 2),
                     lty = c(1, 1, 2, 1), bg = "white", cex = 0.8)
}
par(op)

## final-year marginal posterior of each pool's process sd on its grid axis, vs the prior
op <- par(mfrow = c(1, 3), mar = c(4.5, 4, 2.5, 0.5))
for (a in 1:3) {
  sdx <- 1/sqrt(or.fit$tauAxes[[a]])
  pr  <- cell_prior_mass(or.fit$tauAxes[[a]], list(type = "loguniform"))
  plot(sdx, or.fit$tauMarg[[a]][, NYears], type = 'h', lwd = 4, col = "#2a78d6", log = "x",
       ylim = c(0, max(c(or.fit$tauMarg[[a]][, NYears], pr))), xlab = "process sd", ylab = "posterior mass",
       main = sprintf("%s, year %d", pool.names[a], NYears))
  lines(sdx, pr, lty = 3)
  abline(v = sqrt(Q[a,a]), lwd = 2)
}
par(op)

## forecast vs analysis at the disturbance year, leaf vs stem (cf. the Rmd's ensemble-adjustment plot)
td <- which(disturbance_at == 1)
xf <- rmix(or.fit$mixF[[td]], 1000); xa <- rmix(or.fit$mixA[[td]], 1000)
plot(xf[, 1], xf[, 3], cex = 0.5, xlim = range(c(xf[, 1], xa[, 1])), ylim = range(c(xf[, 3], xa[, 3])),
     xlab = pool.names[1], ylab = pool.names[3], main = sprintf("Multinomial filter, year %d", td))
points(xa[, 1], xa[, 3], pch = "+", col = adjustcolor("red", 0.5))
points(Observed_states[td, 1], Observed_states[td, 3], cex = 3, col = adjustcolor("purple", 0.75), pch = 18)
points(True_states[td, 1], True_states[td, 3], cex = 2, pch = 4, lwd = 2)
legend("topleft", legend = c("forecast", "analysis", "data", "truth"), pch = c(1, 3, 18, 4),
       col = c(1, 2, "purple", 1))

dev.off()

## ---- Save outputs -------------------------------------------------------------------
saveRDS(list(True_states = True_states, Observed_states = Observed_states, R = R, Q = Q,
             disturbance_at = disturbance_at, an0.stats = an0.stats, an.stats = an.stats,
             or.fit = or.fit, or.known = or.known),
        file.path(outdir, "results.rds"))

## headline numbers
rmse <- function(x) sqrt(colMeans((x - True_states)^2))
fmt  <- function(x) paste(sprintf("%s %.3f", pool.names, x), collapse = ", ")
cat("analysis-median RMSE by pool\n")
cat("  EnKF no-dist :", fmt(rmse(t(an0.stats[3, , ]))), "\n")
cat("  EnKF dist    :", fmt(rmse(t(an.stats[3, , ]))), "\n")
cat("  multinomial  :", fmt(rmse(or.fit$q[3, , ])), "\n")
cat("  (tau known)  :", fmt(rmse(or.known$q[3, , ])), "\n")
cat(sprintf("P(disturbed) at year %d: %.3f (tau known %.3f); max elsewhere %.3f\n",
            td, or.fit$pDist[td], or.known$pDist[td], max(or.fit$pDist[-td])))
cat(sprintf("process sd, sqrt(E[1/tau]) at year %d: %s (truth %s)\n", NYears,
            fmt(sqrt(or.fit$procVar[NYears, ])), paste(sprintf("%.3f", sqrt(diag(Q))), collapse = ", ")))
