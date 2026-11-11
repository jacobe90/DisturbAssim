## The OR joint-update filter (ORJointUpdate) and its helpers.
##
## Ported from OR_joint_update2.Rmd (now legacy/models/). Needs VSEM (BayesianTools),
## disturbance.p.rows (R/disturbance.R) and mvnfast.

# Helper functions for ORJointUpdate

KalmanUpdate <- function(m, C, H, R, y) {
  S <- H %*% C %*% t(H) + R
  K <- C %*% t(H) %*% solve(S)
  m_hat <- m + K %*% (y - H %*% m)
  C_hat <- (diag(nrow(C)) - K %*% H) %*% C
  L_hat <- dnorm(y, mean = H %*% m, sd = sqrt(S))
  
  return(list(m_hat=m_hat, C_hat=C_hat, L_hat=L_hat))
}

# Quantiles of a 1-d Gaussian mixture sum_i w_i N(mu_i, sd_i^2), by root-finding on its
# CDF. Components with negligible weight are dropped first.
mix_quantile <- function(p, w, mu, sd) {
  keep <- w > 1e-12 * sum(w)
  w <- w[keep] / sum(w[keep]); mu <- mu[keep]; sd <- sd[keep]
  lo <- min(mu - 8 * sd); hi <- max(mu + 8 * sd)
  sapply(p, function(pp) uniroot(function(x) sum(w * pnorm(x, mu, sd)) - pp,
                                 c(lo, hi), tol = 1e-8)$root)
}

# Initial tau weights: grid point g gets the prior mass of its cell,
#   int_{lo_g}^{hi_g} p(tau) dtau,  lo_g = tau_g - (tau_g - tau_{g-1})/2,  hi_g = tau_g + (tau_{g+1} - tau_g)/2,
# i.e. the cell edges are the midpoints to the neighbouring grid points, so the cell
# width follows the grid spacing. The edge cells stop at tau_1 and tau_G and the
# weights are renormalised (the prior is truncated to the grid range).
cell_prior_mass <- function(tauGrid, prior) {
  G   <- length(tauGrid)
  mid <- (tauGrid[-G] + tauGrid[-1]) / 2
  lo  <- c(tauGrid[1], mid)
  hi  <- c(mid, tauGrid[G])
  cdf <- switch(prior$type,
    gamma      = function(x) pgamma(x, shape = prior$shape, rate = prior$rate),
    loguniform = function(x) log(x))
  m <- cdf(hi) - cdf(lo)
  m / sum(m)
}

# moment matching procedure to project down to mixture of two gaussians
project <- function(w, l, mu, sig, D) {
  w1 <- 0
  for (j in 1:2) {
    w1 <- w1 + w[j,1]
  }
  w2 <- 1 - w1
  alpha <- w
  for (j in 1:2) {
    for (k in 0:D) {
      if (k == 0) {
        alpha[j,k+1] <- w[j,k+1] / w1
      } else {
        alpha[j,k+1] <- w[j,k+1] / w2
      }
    }
  }
  mu1 <- alpha[1,1] * mu[1,1,] + alpha[2,1] * mu[2,1,]
  sig1 <- array(0, dim=c(3,3))
  for (j in 1:2) {
    sig1 <- sig1 + alpha[j,1] * (tcrossprod(mu[j,1,]) + sig[j,1,,])
  }
  sig1 <- sig1 - tcrossprod(mu1)
  mu2 <- array(0, dim=c(3))
  for (j in 1:2) {
    for (k in 1:D) {
      mu2 <- mu2 + alpha[j,k+1] * mu[j,k+1,]
    }
  }
  sig2 <- array(0, dim=c(3,3))
  for (j in 1:2) {
    for (k in 1:D) {
      sig2 <- sig2 + alpha[j,k+1] * (tcrossprod(mu[j,k+1,]) + sig[j,k+1,,])
    }
  }
  sig2 <- sig2 - tcrossprod(mu2)
  # a group with zero total weight has undefined moments (0/0 above); give it the other
  # group's moments so it stays finite -- its weight is 0, so it only enters the proposal
  if (!(w1 > 0)) { mu1 <- mu2; sig1 <- sig2 }
  if (!(w2 > 0)) { mu2 <- mu1; sig2 <- sig1 }
  return(list(w1=w1,w2=w2,mu1=mu1,mu2=mu2,sig1=sig1,sig2=sig2))
  
}

# project() applied to every tau grid point g at once. w: G x 2 x (D+1) weights,
# mu: G x 2 x (D+1) x 3, sig: G x 2 x (D+1) x 3 x 3. Group 1 is the undisturbed
# components (k = 0), group 2 the disturbed ones (k >= 1). Returns w (G x 2),
# mu (G x 2 x 3) and sig (G x 2 x 3 x 3), indexed [g, group, ...].
project_grid <- function(w, mu, sig) {
  G <- dim(w)[1]; C <- prod(dim(w)[-1])        # components c = j + 2k, j fastest
  muC  <- array(mu,  c(G, C, 3))
  sigC <- array(sig, c(G, C, 3, 3))
  w1 <- w[, 1, 1] + w[, 2, 1]
  w2 <- 1 - w1
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

##' Approximate joint disturbance update
##'
##' @param gamInit c(a, b) Gamma prior (shape, rate) for tau
##' @param Y NT x 3 observations (wood/AGB pool in column 3; NA => forecast-only step)
##' @param Rparams c(slope, intercept): observation error R = y*slope + intercept
##' @param Q process-error covariance; Q[1,1],Q[2,2] used for foliage/soil, Q[3,3] set to 1/tau_g.
##'   Undisturbed components (k = 0) get Q. Disturbed components (k >= 1) get D_k Q D_k with
##'   D_k = diag(mu0_leaf, 1, mu0_stem), mu0 = dist$vals[[k]]$mu0 (see distProcErr)
##' @param PAR daily PAR driver (length NT*365)
##' @param Xic J x 3 initial-condition ensemble (ICprior)
##' @param THETA VSEM parameter vector/matrix (cols 9:11 overwritten with the IC ensemble)
##' @param dist disturbance list: p (D+1 x NT prior), n (# classes), vals (per-class mu0/V0), alloc.soil
##' @param G number of precision values in the grid
##' @param N ensemble size
##' @param prior prior for tau: list(type = "gamma", shape, rate) or list(type = "loguniform")
##' @param grid list(tauMin, tauMax) for the log-spaced tau grid
##' @param distProcErr TRUE: propagate the process error through the disturbance. Under the
##'   percent-decrease model (disturbance.p: new = mu0 * old + noise) the process noise e of
##'   the undisturbed year passes through the disturbance, so the disturbed leaf and stem keep
##'   mu0 e (variance mu0^2 Q) and the soil keeps all of it. FALSE: disturbed components get
##'   no process error (the behaviour before this option was added).
##' @param bias NULL (no bias correction) or a list configuring an online additive
##'   bias correction b_t of the undisturbed (k = 0) wood-pool forecast mean:
##'     rho    discount factor in (0, 1] on past residuals (default 0.9)
##'     denom  "N": b_{t+1} = S_t / N_t (weighted discounted mean, default);
##'            "t": b_{t+1} = S_t / n_t, n_t = # obs used so far (shrinks b toward 0)
##'     skip   number of initial observed steps not used for b (default 1). The first
##'            forecast starts from the population-wide IC prior (ICprior), so its
##'            residual measures IC mismatch, not VSEM bias.
##'     burnin number of observed steps after `skip` that use w_t = 1 instead of
##'            P(z_t = 0 | y_{1:t}) (guards against bias being absorbed as disturbance)
##'     wMin   b is only updated on steps with P(z_t = 0 | y_{1:t}) >= wMin (after burn-in)
##'   Each observed step, with posterior component weights
##'     tw_{gjk} prop. to pi_{t-1}(g) prevW[g,j] beta_k N(y_t; muF_{gjk,3} + b_t 1{k=0}, sigF_{gjk,33} + R_t),
##'   r_t = sum_{gj} tw_{gj0} (y_t - muF_{gj0,3}) / sum_{gj} tw_{gj0}   (residual vs the RAW forecast),
##'   w_t = sum_{gj} tw_{gj0} = P(z_t = 0 | y_{1:t}),
##'   S_t = rho S_{t-1} + w_t r_t,  N_t = rho N_{t-1} + w_t,  b_{t+1} = S_t / N_t.
##'   b_t is added to muF[, , k = 0, 3] (only the undisturbed components) before the
##'   analysis at step t; b_1 = 0. y_t only enters b from t + 1 onward.
##' @return list of length-NT vectors, each marginalised over the tau grid with the
##'   tau posterior pi_t(g) and over the mixture components (j, k):
##'   mean  = analysis mean of the wood pool
##'   var   = analysis variance of the wood pool (law of total variance)
##'   pDist = posterior probability of disturbance
##'   ci95  = NT x 2 matrix, 2.5% / 97.5% quantiles of the wood-pool analysis mixture
##'   pit   = probability integral transform of y_t under its one-step-ahead predictive
##'           (marginalised over the tau posterior pi_{t-1}(g), the components j and the
##'           disturbance prior beta_k, i.e. before y_t is assimilated); NA without obs
##'   wDrop = forecast weight of components dropped for degenerate importance weights
##'   wDropTau = tau-posterior mass pi_{t-1}(g) of grid points dropped because every
##'           component there was degenerate (zero marginal likelihood)
##'   bias, biasResid, biasW = b_t applied at step t, and r_t, w_t (NA without obs)
##'   predMeanRaw = predictive mean of y_t without the bias shift
##'   predMean0, predVar0, predProcVar0 = one-step-ahead predictive mean / variance (incl.
##'           R) / process-variance part, conditioned on no disturbance (k = 0 only)
ORJointUpdate <- function(Y, Rparams, Q, PAR, Xic, THETA, dist, G = 200, N=100, prior, grid,
                          bias = NULL, distProcErr = TRUE) {
  if (!is.null(bias)) {
    bias <- modifyList(list(rho = 0.9, denom = "N", skip = 1, burnin = 0, wMin = 0), bias)
    bias$denom <- match.arg(bias$denom, c("N", "t"))
  }
  bS <- 0; bN <- 0; bObs <- 0; bHat <- 0   # bias-correction state (S_t, N_t, # obs, b_{t+1})
  NT <- dim(Y)[1]
  D  <- dist$n                      # number of disturbance sub-classes
  icMean <- colMeans(Xic)
  icCov <- cov(Xic)
  
  # process-noise scaling d_k of each scenario: the forecast covariance of scenario k gets
  # D_k Q D_k, D_k = diag(d_k) (row k + 1; Q[3, 3] = 1/tau_g). Undisturbed: d = (1, 1, 1).
  # Disturbed, percent-decrease model: d = (mu0_leaf, 1, mu0_stem), see distProcErr.
  procD <- rbind(c(1, 1, 1), t(vapply(seq_len(D), function(k)
    if (distProcErr) c(dist$vals[[k]]$mu0[1], 1, dist$vals[[k]]$mu0[2]) else c(0, 0, 0),
    numeric(3))))
  Q0 <- Q; Q0[3, 3] <- 0                       # Q without the wood term (1/tau_g added per g)

  # initialize grid
  tauGrid <- exp(seq(log(grid$tauMin), log(grid$tauMax), length.out = G)) # log-spaced grid
  # tauGrid <- seq(grid$tauMin, grid$tauMax, length.out = G)
  
  # initial tau weights: log prior mass of each grid cell (see cell_prior_mass)
  qs <- log(cell_prior_mass(tauGrid, prior))
  prevW <- array(0.5, dim=c(G, 2))
  prevMu <- array(0, dim=c(G, 2, 3))
  prevSig <- array(0, dim=c(G, 2, 3, 3))
  for (i in 1:G) {
    for (j in 1:2) {
      prevMu[i,j,] <- icMean
      prevSig[i,j,,] <- icCov
    }
  }
  
  hist <- list(
    mean  = rep(NA_real_, NT),   # wood-pool analysis mean
    var   = rep(NA_real_, NT),   # wood-pool analysis variance
    ci95  = matrix(NA_real_, NT, 2, dimnames = list(NULL, c("lo", "hi"))),  # 95% interval
    pDist = rep(NA_real_, NT),   # P(disturbed)
    pit   = rep(NA_real_, NT),   # P(y_t' <= y_t | y_{1:t-1}), one-step-ahead PIT
    wDrop = rep(NA_real_, NT),   # forecast weight of degenerate (dropped) components
    wDropTau = rep(0, NT),       # tau-posterior mass of grid points with zero likelihood
    ## one-step-ahead predictive of y_t (NA on steps with no observation):
    predMean    = rep(NA_real_, NT),   # E[y_t | y_{1:t-1}]
    predVar     = rep(NA_real_, NT),   # Var[y_t | y_{1:t-1}], including R
    predProcVar = rep(NA_real_, NT),   # its process-error part: E over components of d_{k,3}^2/tau_g
    predR       = rep(NA_real_, NT),   # its observation-error part R_t
    ## the same predictive conditioned on no disturbance (k = 0 components only):
    predMean0   = rep(NA_real_, NT),   # E[y_t | y_{1:t-1}, z_t = 0]
    predVar0    = rep(NA_real_, NT),   # Var[y_t | y_{1:t-1}, z_t = 0], including R
    predProcVar0 = rep(NA_real_, NT),  # its process-error part E[1/tau_g | y_{1:t-1}, z_t = 0]
    procVar     = rep(NA_real_, NT),   # posterior mean process variance E[1/tau | y_{1:t}]
    predMeanRaw = rep(NA_real_, NT),   # E[y_t | y_{1:t-1}] without the bias shift
    bias        = rep(0, NT),          # b_t added to the k = 0 wood forecast mean at step t
    biasResid   = rep(NA_real_, NT),   # r_t: P(z=0)-weighted residual vs the raw k = 0 forecast
    biasW       = rep(NA_real_, NT)    # w_t = P(z_t = 0 | y_{1:t})
  )
  
  for (t in 1:NT) {
    # sample from proposal distribution
    propSamples <- array(0, dim=c(N, 3))
    n <- 1L; tries <- 0L
    while (n <= N) {
      tries <- tries + 1L
      i <- sample.int(G, 1L); j <- sample.int(2L, 1L)
      x <- mvnfast::rmvn(1, prevMu[i,j,], prevSig[i,j,,])
      if (all(x > 0)) { propSamples[n, ] <- x; n <- n + 1L }
      if (tries > 100L * N) stop(sprintf("t=%d: proposal rejection rate too high", t))
    }
    
    muF <- array(0, dim=c(G, 2, D+1, 3))
    sigF <- array(0, dim=c(G, 2, D+1, 3, 3))
    tsel = (t-1)*365 + 1:365
    ## ---- once per timestep ----

    # 1. VSEM is deterministic in THETA -> identical for every k
    sim0 <- matrix(0, N, 3)
    for (n in seq_len(N)) {
      THETA[9:11] <- propSamples[n, ]
      sim0[n, ]   <- VSEM(THETA, PAR = PAR[tsel])[365, 2:4]
    }
    
    # 2. flatten (g,j) -> r, cache the Cholesky factors
    P        <- 2L * G
    compMu   <- matrix(0, P, 3)
    compChol <- vector("list", P)
    for (j in 1:2) for (g in 1:G) {
      r <- (j - 1L) * G + g
      compMu[r, ]   <- prevMu[g, j, ]
      compChol[[r]] <- chol(prevSig[g, j, , ])
    }
    
    # 3. one P x N density matrix, vectorized over samples
    logDens <- matrix(0, P, N)
    for (r in seq_len(P)) {
      logDens[r, ] <- mvnfast::dmvn(propSamples, compMu[r, ], compChol[[r]],
                                    log = TRUE, isChol = TRUE)
    }
    
    # 4. proposal density and all weights, in log space
    mx      <- apply(logDens, 2, max)
    logProp <- mx + log(colSums(exp(sweep(logDens, 2, mx)))) - log(P)
    logW    <- sweep(logDens, 2, logProp)
    W       <- exp(logW - apply(logW, 1, max))
    W       <- W / rowSums(W)          # P x N, rows sum to 1
    cN      <- 1 - rowSums(W^2)        # P

    # A component whose importance weights collapse onto ~1 sample (cN ~ 0) has no
    # estimable forecast moments (S / cN -> Inf, then NaN likelihoods). This happens to
    # components left far from every proposal sample, e.g. the undisturbed component at
    # high tau right after a clear-cut, whose weight is then negligible. Such components
    # get the unweighted ensemble moments as a finite placeholder and zero likelihood
    # below (so zero weight); hist$wDrop records the forecast weight removed this way.
    degen <- cN < 1e-8
    W[degen, ] <- 1 / N
    cN[degen]  <- 1 - 1 / N
    dropF <- matrix(degen, G, 2)       # r = (j - 1) * G + g  ->  [g, j]
    
    ia <- rep(1:3, 3); ib <- rep(1:3, each = 3)   # column (a, b) of a flattened 3 x 3
    for (k in 0:D) {
      sim <- sim0
      if (k > 0) {
        sim <- disturbance.p.rows(sim0, dist$vals[[k]]$mu0, dist$vals[[k]]$V0, dist$alloc.soil)
      }

      MU <- W %*% sim                              # P x 3, every mean in one matmul

      # every weighted covariance in one matmul (rows of W sum to 1):
      #   S_r = (sum_n W[r,n] x_n x_n' - MU_r MU_r') / cN_r,
      # with x centred on the ensemble mean first to limit cancellation
      ctr <- colMeans(sim)
      Xc  <- sweep(sim, 2, ctr)
      MUc <- sweep(MU, 2, ctr)
      S   <- (W %*% (Xc[, ia] * Xc[, ib]) - MUc[, ia] * MUc[, ib]) / cN   # P x 9
      # the subtraction above can leave a variance slightly negative by round-off when
      # every member maps to (almost) the same value, e.g. a disturbance taking a
      # near-zero wood pool to ~0; clamp the diagonal at 0
      S[, c(1, 5, 9)] <- pmax(S[, c(1, 5, 9)], 0)
      # process error D_k Q D_k (for k = 0 just Q); wood: d_3^2 / tau_g
      d <- procD[k + 1, ]
      S <- sweep(S, 2, as.vector(Q0 * outer(d, d)), "+")
      S[, 9] <- S[, 9] + d[3]^2 / rep(tauGrid, 2)          # S[3, 3]: process var d_3^2 / tau_g
      # r = (j - 1) * G + g runs g fastest, so P x 3 / P x 9 fill [g, j, , ] in order
      muF[, , k + 1, ]    <- MU
      sigF[, , k + 1, , ] <- S
    }

    # bias correction: shift the undisturbed wood forecast mean by b_t (raw copy kept
    # for the residual r_t, which is measured against the uncorrected forecast)
    muF0raw <- muF[, , 1, 3]                       # G x 2
    if (!is.null(bias)) muF[, , 1, 3] <- muF0raw + bHat
    hist$bias[t] <- bHat
    
    # compute analysis
    lA <- array(0, dim=c(G, 2, D+1))
    muA <- array(0, dim=c(G, 2, D+1, 3))
    sigA <- array(0, dim=c(G, 2, D+1, 3, 3))
    if (is.na(Y[t, 3])) {
      # no observation: analysis = forecast, flat likelihood (weights move only by the prior p)
      muA  <- muF
      sigA <- sigF
      lA[] <- 1
    } else {
      R <- Y[t, 3] * Rparams[1] + Rparams[2] # heteroskedastic obs error
      # KalmanUpdate for every (g, j, k) at once. H = (0, 0, 1) and R is scalar, so with
      # s = C33 + R:  K = C[, 3] / s,  m_hat = m + K (y - m3),  C_hat = C - K C[3, ]
      s    <- sigF[, , , 3, 3] + R                                  # G x 2 x (D+1)
      K    <- sigF[, , , , 3] / as.vector(s)                         # G x 2 x (D+1) x 3
      muA  <- muF + K * as.vector(Y[t, 3] - muF[, , , 3])
      for (a in 1:3) for (b in 1:3) sigA[, , , a, b] <- sigF[, , , a, b] - K[, , , a] * sigF[, , , 3, b]
      lA[] <- dnorm(Y[t, 3], mean = muF[, , , 3], sd = sqrt(s))
    }
    for (k in 0:D) lA[, , k + 1][dropF] <- 0   # degenerate forecast components (see above)
    
    # compute tau posterior on grid
    beta <- array(0, dim=c(D+1))
    for (k in 0:D) {
      beta[k+1] = dist$p[k+1,t]
    }

    # PIT of y_t under its one-step-ahead predictive, a mixture over (g, j, k) with
    # weights pi_{t-1}(g) * prevW[g,j] * beta_k (degenerate components dropped) and
    # components N(muF[g,j,k,3], sigF[g,j,k,3,3] + R); sigF[,,k+1,3,3] includes the
    # process variance d_{k,3}^2/tau_g (1/tau_g for k = 0). Left NA on forecast-only steps with no observation.
    piPrev <- exp(qs - max(qs)); piPrev <- piPrev / sum(piPrev)            # pi_{t-1}(g)
    hist$wDrop[t] <- sum(piPrev * rowSums(prevW * dropF))                  # dropped weight
    if (!is.na(Y[t, 3])) {
      wF <- outer(outer(piPrev, rep(1, 2)) * prevW, beta) * as.vector(!dropF)  # G x 2 x (D+1)
      wF <- wF / sum(wF)
      hist$pit[t] <- sum(wF * pnorm(Y[t, 3], muF[, , , 3], sqrt(sigF[, , , 3, 3] + R)))
      hist$predMean[t]    <- sum(wF * muF[, , , 3])
      hist$predMeanRaw[t] <- hist$predMean[t] - bHat * sum(wF[, , 1])
      hist$predVar[t]     <- sum(wF * (sigF[, , , 3, 3] + R + muF[, , , 3]^2)) - hist$predMean[t]^2
      hist$predProcVar[t] <- sum(wF * rep(procD[, 3]^2, each = 2 * G) / tauGrid)  # component k carries d_3^2/tau_g
      hist$predR[t]       <- R
      w0 <- wF[, , 1] / sum(wF[, , 1])                     # k = 0 weights, renormalised
      hist$predMean0[t]    <- sum(w0 * muF[, , 1, 3])
      hist$predVar0[t]     <- sum(w0 * (sigF[, , 1, 3, 3] + R + muF[, , 1, 3]^2)) - hist$predMean0[t]^2
      hist$predProcVar0[t] <- sum(w0 / tauGrid)
    }
    wPrior <- outer(prevW, beta)                     # G x 2 x (D+1): prevW[g,j] * beta_k
    zsum   <- apply(wPrior * lA, 1, sum)             # marginal likelihood of y_t given tau_g
    logz   <- log(zsum)
    # A grid point whose components are ALL degenerate has zero marginal likelihood. This
    # happens at extreme tau that hold negligible posterior mass (e.g. tau ~ 0.1, process
    # sd ~ 3 kg C m-2, when the mode is at tau ~ 500): drop it (pi_t(g) = 0 from here on)
    # rather than stop. Stop only if the dropped mass is not negligible.
    badG <- !is.finite(logz)
    if (any(badG)) {
      hist$wDropTau[t] <- sum(piPrev[badG])
      if (all(badG) || hist$wDropTau[t] > 1e-6)
        stop(sprintf("t=%d: zero or non-finite marginal likelihood at %d tau grid points (tau-posterior mass %.2g)",
                     t, sum(badG), hist$wDropTau[t]))
      logz[badG] <- -Inf
    }
    qs     <- qs + logz
    
    # tau posterior on the grid, pi_t(g)
    m <- max(qs)
    ws <- exp(qs - m) / sum(exp(qs - m))
    hist$procVar[t] <- sum(ws / tauGrid)
    
    # analysis weights
    wA <- wPrior * lA / zsum
    # dropped grid points (pi_t(g) = 0): prior weights as a finite placeholder, so their
    # projected moments stay finite (they still enter the proposal, with zero weight)
    if (any(badG)) wA[badG, , ] <- wPrior[badG, , ]

    # bias update from y_t (applied from t + 1 on). Posterior weights of the undisturbed
    # components, tw_{gj0} = pi_t(g) wA[g, j, 0], which equals
    # pi_{t-1}(g) prevW[g,j] beta_0 N(y_t; muF + b, sigF + R) / normaliser.
    if (!is.null(bias) && !is.na(Y[t, 3])) {
      tw0 <- ws * wA[, , 1]                        # G x 2
      wt  <- sum(tw0)                              # P(z_t = 0 | y_{1:t})
      if (wt > 0) {
        rt <- sum(tw0 * (Y[t, 3] - muF0raw)) / wt
        bObs <- bObs + 1
        burn <- bObs <= bias$skip + bias$burnin
        if (burn) wt <- 1                          # burn-in: full weight
        hist$biasResid[t] <- rt
        hist$biasW[t]     <- wt
        if (bObs > bias$skip && (burn || wt >= bias$wMin)) {
          bS <- bias$rho * bS + wt * rt
          bN <- bias$rho * bN + wt
          bHat <- if (bias$denom == "N") bS / bN else bS / (bObs - bias$skip)
        }
      }
    }
    
    # recompute prev posterior vals (KL projection), for every tau grid point at once
    out     <- project_grid(wA, muA, sigA)
    prevW   <- out$w
    prevMu  <- out$mu
    prevSig <- out$sig
    
    # store analysis values for this tstep, marginalised over the tau grid (weights
    # pi_t(g) = ws) and the components (j, k) (weights wA[g, , ], which sum to 1 for each g)
    pDist.g <- apply(wA[, , -1, drop = FALSE], 1, sum)                  # P_t(g)
    mean.g  <- apply(wA * muA[, , , 3], 1, sum)                           # E[wood | g]
    msq.g   <- apply(wA * (sigA[, , , 3, 3] + muA[, , , 3]^2), 1, sum)   # E[wood^2 | g]
    hist$mean[t]  <- sum(ws * mean.g)
    hist$var[t]   <- sum(ws * msq.g) - hist$mean[t]^2
    hist$pDist[t] <- sum(ws * pDist.g)
    hist$ci95[t, ] <- mix_quantile(c(0.025, 0.975), ws * wA, muA[, , , 3], sqrt(pmax(sigA[, , , 3, 3], 0)))
  }
  
  hist$procScale <- procD   # d_k per scenario (rows k = 0..D; leaf, soil, wood)
  return(hist)
}
