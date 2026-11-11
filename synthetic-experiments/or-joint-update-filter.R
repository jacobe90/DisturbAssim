## Generalised ORJointUpdate filter shared by the multi-pool synthetic experiments.
##
## Extends the filter of singlesite-singledisturb-multipool.R (itself a port of
## ORJointUpdate in or-model/R/or-joint-update.R) along three directions:
##   - any state dimension d (e.g. 3 pools x 2 sites)
##   - any number of disturbance scenarios k = 0..K (k = 0: undisturbed) and of discrete
##     history groups h = 1..H. The analysis at every tau grid point is a mixture of one
##     Gaussian per history group; scenario k sends group h to group transition[h, k+1],
##     and each group carries a model label (e.g. its PFT) that selects the forecast model
##   - the process error is diagonal, element e having variance 1/tau[poolOf[e]], with
##     (tau_1, tau_2, tau_3) learned on a 3-D log-spaced product grid
##
## Usage: source("synthetic-experiments/or-joint-update-filter.R") from the repository root.

source("synthetic-experiments/filter-diagnostics.R")   # apply_drops, plot_dropped_mass

qs5 <- c(0.025, 0.25, 0.5, 0.75, 0.975)

# Quantiles of a 1-d Gaussian mixture sum_i w_i N(mu_i, sd_i^2), by root-finding on its
# CDF. Components with negligible weight are dropped first.
mix_quantile <- function(p, w, mu, sd) {
  keep <- w > 1e-10 * sum(w)
  w <- w[keep] / sum(w[keep]); mu <- mu[keep]; sd <- sd[keep]
  lo <- min(mu - 8 * sd); hi <- max(mu + 8 * sd)
  sapply(p, function(pp) uniroot(function(x) sum(w * pnorm(x, mu, sd)) - pp,
                                 c(lo, hi), tol = 1e-6 * (hi - lo))$root)
}

# Initial tau weights on one grid axis: prior mass of each grid cell (see or-model/R/or-joint-update.R)
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
## A batch of n symmetric d x d matrices is stored as an n x d^2 matrix, column-major:
## entry (i, j) is column i + d (j - 1). Every routine loops over matrix indices and is
## vectorised over the batch, so thousands of components cost a few vector operations.
mi <- function(i, j, d) i + d * (j - 1)

# lower Cholesky factors
chol_batch <- function(S, d) {
  L <- matrix(0, nrow(S), d * d)
  for (j in 1:d) {
    s <- S[, mi(j, j, d)]
    for (k in seq_len(j - 1)) s <- s - L[, mi(j, k, d)]^2
    L[, mi(j, j, d)] <- sqrt(s)
    for (i in seq_len(d - j) + j) {
      s <- S[, mi(i, j, d)]
      for (k in seq_len(j - 1)) s <- s - L[, mi(i, k, d)] * L[, mi(j, k, d)]
      L[, mi(i, j, d)] <- s / L[, mi(j, j, d)]
    }
  }
  L
}

# inverse and log determinant from the Cholesky factors: S^-1 = L^-T L^-1
inv_chol_batch <- function(L, d) {
  Li <- matrix(0, nrow(L), d * d)
  for (i in 1:d) {
    Li[, mi(i, i, d)] <- 1 / L[, mi(i, i, d)]
    for (j in seq_len(i - 1)) {
      s <- 0
      for (k in j:(i - 1)) s <- s + L[, mi(i, k, d)] * Li[, mi(k, j, d)]
      Li[, mi(i, j, d)] <- -s / L[, mi(i, i, d)]
    }
  }
  A <- matrix(0, nrow(L), d * d)
  for (i in 1:d) for (j in i:d) {
    s <- 0
    for (k in j:d) s <- s + Li[, mi(k, i, d)] * Li[, mi(k, j, d)]
    A[, mi(i, j, d)] <- s; A[, mi(j, i, d)] <- s
  }
  list(inv = A, logdet = 2 * rowSums(log(L[, mi(1:d, 1:d, d), drop = FALSE])))
}

# P x N log densities of the samples X (N x d) under each component (mu: P x d,
# sig: P x d^2). The quadratic form (x - m)' A (x - m) = x'Ax - 2 (Am)'x + m'Am is
# expanded so that every (component, sample) pair comes out of two matrix products.
log_dens_batch <- function(X, mu, sig) {
  d   <- ncol(X)
  ctr <- colMeans(X)                          # centre both on a common point first
  X   <- sweep(X, 2, ctr); mu <- sweep(mu, 2, ctr)
  ic  <- inv_chol_batch(chol_batch(sig, d), d); A <- ic$inv
  ut  <- which(upper.tri(diag(d), diag = TRUE), arr.ind = TRUE)
  quadX <- X[, ut[, 1], drop = FALSE] * X[, ut[, 2], drop = FALSE] *
           rep(ifelse(ut[, 1] == ut[, 2], 1, 2), each = nrow(X))
  Aq <- A[, mi(ut[, 1], ut[, 2], d), drop = FALSE]
  Am <- sapply(1:d, function(i) rowSums(A[, mi(i, 1:d, d), drop = FALSE] * mu))
  quad <- tcrossprod(Aq, quadX) - 2 * tcrossprod(Am, X) + rowSums(Am * mu)
  -0.5 * quad - 0.5 * ic$logdet - 0.5 * d * log(2 * pi)
}

# Kalman update of a batch of components (mu: n x d, sig: n x d^2) with H = I on the
# observed elements o. Returns the analysis mu / sig and the log likelihood of y[o].
kalman_batch <- function(mu, sig, y, R, o) {
  n <- nrow(mu); d <- ncol(mu); m <- length(o)
  So <- sig[, as.vector(outer(o, d * (o - 1), "+")), drop = FALSE] + rep(as.vector(R[o, o]), each = n)
  ic <- inv_chol_batch(chol_batch(So, m), m); A <- ic$inv
  v  <- matrix(rep(y[o], each = n), n, m) - mu[, o, drop = FALSE]     # innovations
  muA <- mu; sigA <- sig
  K  <- vector("list", d)                                             # rows of C[, o] So^-1
  for (a in 1:d) K[[a]] <- sapply(1:m, function(b)
    rowSums(sig[, a + d * (o - 1), drop = FALSE] * A[, (b - 1) * m + 1:m, drop = FALSE]))
  for (a in 1:d) {
    Ka <- matrix(K[[a]], n, m)
    muA[, a] <- mu[, a] + rowSums(Ka * v)
    for (b in 1:d) sigA[, mi(a, b, d)] <- sig[, mi(a, b, d)] - rowSums(Ka * sig[, o + d * (b - 1), drop = FALSE])
  }
  qf <- 0
  for (b in 1:m) for (c in 1:m) qf <- qf + v[, b] * A[, b + m * (c - 1)] * v[, c]
  list(mu = muA, sig = sigA, logL = -0.5 * qf - 0.5 * ic$logdet - 0.5 * m * log(2 * pi))
}

# Moment-matching projection of the analysis mixture onto one Gaussian per history group,
# for every tau grid point. w: G x C weights (C = H (K+1) components, sum 1 per row),
# mu: G x C x d, sig: G x C x d^2; grp: new history group of each component. A group with
# zero weight gets the moments of the whole mixture so it stays finite (it only enters
# the proposal).
project_groups <- function(w, mu, sig, grp, H) {
  G <- nrow(w); d <- dim(mu)[3]
  moments <- function(cols) {
    wg <- rowSums(w[, cols, drop = FALSE])
    A  <- w[, cols, drop = FALSE] / wg
    m  <- sapply(1:d, function(a) rowSums(A * mu[, cols, a]))
    m  <- matrix(m, G, d)
    s  <- sapply(1:(d * d), function(ab) {
      a <- (ab - 1) %% d + 1; b <- (ab - 1) %/% d + 1
      rowSums(A * (mu[, cols, a] * mu[, cols, b] + sig[, cols, ab])) - m[, a] * m[, b]
    })
    list(w = wg, mu = m, sig = matrix(s, G, d * d))
  }
  all <- moments(seq_len(ncol(w)))
  out <- list(w = matrix(0, G, H), mu = array(0, c(G, H, d)), sig = array(0, c(G, H, d * d)))
  for (h in 1:H) {
    cols <- which(grp == h)
    mh <- if (length(cols)) moments(cols) else list(w = rep(0, G), mu = all$mu, sig = all$sig)
    z  <- !(mh$w > 0)
    mh$mu[z, ] <- all$mu[z, ]; mh$sig[z, ] <- all$sig[z, ]
    out$w[, h] <- mh$w; out$mu[, h, ] <- mh$mu; out$sig[, h, ] <- mh$sig
  }
  out
}

##' Approximate joint disturbance + process-precision update (generalised ORJointUpdate)
##'
##' @param Y NT x d observations; NA entries are not assimilated
##' @param R d x d observation-error covariance
##' @param Xic J x d initial-condition ensemble
##' @param simulate function(x, t, label, mem) returning the N x d deterministic forecast of
##'   the states x (N x d) over step t, with the model for group label `label` and the
##'   parameters / drivers of ensemble members mem (length N)
##' @param nMembers number of ensemble members mem is drawn from
##' @param labels length-H model label of each history group
##' @param disturb function(x, k) applying disturbance scenario k >= 1 to the forecasts x
##' @param dist.p (K+1) x NT prior probability of each scenario per step (row 1 = undisturbed)
##' @param procScale (K+1) x d multiplier of element e's process variance 1/tau under
##'   scenario k: 1 when undisturbed, 0 where the disturbance resets the state, f^2 where it
##'   keeps a fraction f of a state that already carries process error
##' @param transition H x (K+1) new history group of group h under scenario k
##' @param w0 length-H initial weights of the history groups
##' @param poolOf length-d tau axis (1..3) of each state element
##' @param prior prior for each tau, independent across axes:
##'   list(type = "gamma", shape, rate) or list(type = "loguniform")
##' @param grid list(tauMin, tauMax, G), each length 3: range and number of points per axis
##'   (tauMin == tauMax, G = 1 => that tau is known)
##' @param N importance-sample size
##' @param pruneLog grid points whose log posterior weight falls more than pruneLog below
##'   the maximum are dropped for the rest of the run
##' @param keepMix steps at which to store the full forecast / analysis mixtures
##' @param forecastQuantiles also compute forecast quantiles (fq)
##' @return list with, marginalised over the tau posterior and the components:
##'   mean (NT x d), q / fq (5 x NT x d analysis / forecast quantiles at qs5),
##'   pK (NT x (K+1)) posterior probability of each scenario at that step,
##'   pGroup (NT x H) posterior weight of each history group after the step,
##'   pit (NT x d) marginal one-step-ahead PIT, wDrop (forecast weight of dropped
##'   components, net of restored grid points), wDropIS / wDropCov (forecast weight dropped for
##'   collapsed importance weights / a singular covariance), nNone (grid points left with no
##'   weight, restored), nBad (components dropped for a singular covariance), procVar (NT x 3) = E[1/tau],
##'   tauAxes, tauMarg (per-axis marginal posterior, G_a x NT), tauJoint (per-step joint
##'   posterior on the grid: idx = axis indices of the active grid points, w = their weights),
##'   nActive, and mixF / mixA
##'   (list of w, mu, sig at the keepMix steps)
ORJointUpdateGeneral <- function(Y, R, Xic, simulate, nMembers, labels, disturb, dist.p,
                                 procScale, transition, w0, poolOf, prior, grid, N = 1000,
                                 pruneLog = 30, keepMix = integer(0), forecastQuantiles = TRUE) {
  NT <- nrow(Y); d <- ncol(Y)
  H  <- length(labels); K1 <- nrow(dist.p)
  ia <- rep(1:d, d); ib <- rep(1:d, each = d)  # column (a, b) of a flattened d x d

  tauAxes <- lapply(1:3, function(a) exp(seq(log(grid$tauMin[a]), log(grid$tauMax[a]), length.out = grid$G[a])))
  gridIdx <- as.matrix(expand.grid(lapply(tauAxes, seq_along)))   # Gt x 3 axis indices
  tauMat  <- do.call(cbind, lapply(1:3, function(a) tauAxes[[a]][gridIdx[, a]]))
  logMass <- lapply(tauAxes, function(x) log(cell_prior_mass(x, prior)))
  qs      <- rowSums(do.call(cbind, lapply(1:3, function(a) logMass[[a]][gridIdx[, a]])))
  Gt      <- nrow(tauMat)

  active  <- seq_len(Gt)
  prevW   <- matrix(w0 / sum(w0), Gt, H, byrow = TRUE)
  prevMu  <- array(rep(colMeans(Xic), each = Gt * H), c(Gt, H, d))
  prevSig <- array(rep(as.vector(cov(Xic)), each = Gt * H), c(Gt, H, d * d))
  grpOf   <- as.vector(transition)                # new group of component (h, k), h fastest

  hist <- list(
    mean = matrix(NA_real_, NT, d), q = array(NA_real_, c(5, NT, d)), fq = array(NA_real_, c(5, NT, d)),
    pK = matrix(NA_real_, NT, K1), pGroup = matrix(NA_real_, NT, H), pit = matrix(NA_real_, NT, d),
    wDrop = rep(NA_real_, NT), procVar = matrix(NA_real_, NT, 3), tauAxes = tauAxes,
    tauMarg = lapply(tauAxes, function(x) matrix(NA_real_, length(x), NT)),
    nActive = rep(NA_integer_, NT), nBad = rep(NA_integer_, NT), mixF = list(), mixA = list(),
    wDropIS = rep(NA_real_, NT), wDropCov = rep(NA_real_, NT), nNone = rep(NA_integer_, NT),
    tauJoint = vector("list", NT)
  )

  for (t in 1:NT) {
    G  <- length(active)
    tg <- tauMat[active, , drop = FALSE]
    P  <- G * H                                   # proposal components r = (h - 1) G + g
    Cn <- P * K1                                  # all components, index g + G (h-1) + P k
    compMu  <- array(prevMu,  c(P, d))
    compSig <- array(prevSig, c(P, d * d))
    hist$nActive[t] <- G

    # sample from the proposal: equal-weight mixture over every (g, h), all states positive
    compL <- chol_batch(compSig, d)
    propSamples <- matrix(0, 0, d); tries <- 0L
    while (nrow(propSamples) < N) {
      tries <- tries + 1L
      if (tries > 100L) stop(sprintf("t=%d: proposal rejection rate too high", t))
      ri <- sample.int(P, N, replace = TRUE)
      z  <- matrix(rnorm(d * N), N, d)
      L  <- compL[ri, , drop = FALSE]
      x  <- compMu[ri, , drop = FALSE]
      for (a in 1:d) for (b in 1:a) x[, a] <- x[, a] + L[, mi(a, b, d)] * z[, b]
      propSamples <- rbind(propSamples, x[rowSums(x > 0) == d, , drop = FALSE])
    }
    propSamples <- propSamples[1:N, , drop = FALSE]

    # 1. deterministic forecast of each sample under each model label, then each scenario
    mem  <- sample.int(nMembers, N, replace = TRUE)
    labs <- unique(labels)
    sims <- lapply(setNames(labs, labs), function(lab) {
      s0 <- simulate(propSamples, t, lab, mem)
      c(list(s0), lapply(seq_len(K1 - 1), function(k) disturb(s0, k)))
    })

    # 2-3. importance weights of each component's samples, in log space
    W       <- log_dens_batch(propSamples, compMu, compSig)      # P x N log densities
    logProp <- matrixStats::colLogSumExps(W) - log(P)
    W       <- W - rep(logProp, each = P)
    W       <- exp(W - matrixStats::rowMaxs(W))
    W       <- W / rowSums(W)
    cN      <- 1 - rowSums(W^2)
    degen   <- cN < 1e-8                          # see ORJointUpdate
    W[degen, ] <- 1 / N
    cN[degen]  <- 1 - 1 / N
    dropF <- matrix(degen, G, H)

    # 4. forecast moments of every (g, h, k); muF0 / sigF0 hold the unweighted moments,
    # used as a finite placeholder for components whose weighted covariance is singular
    muF  <- muF0  <- matrix(0, Cn, d)
    sigF <- sigF0 <- matrix(0, Cn, d * d)
    for (lab in labs) {
      hs   <- which(labels == lab)
      rowsP <- as.vector(outer(seq_len(G), (hs - 1) * G, "+"))  # proposal rows of these groups
      Wl  <- W[rowsP, , drop = FALSE]; cNl <- cN[rowsP]
      for (k in 0:(K1 - 1)) {
        sim <- sims[[as.character(lab)]][[k + 1]]
        MU  <- Wl %*% sim
        ctr <- colMeans(sim)
        Xc  <- sweep(sim, 2, ctr)
        MUc <- sweep(MU, 2, ctr)
        S   <- (Wl %*% (Xc[, ia] * Xc[, ib]) - MUc[, ia] * MUc[, ib]) / cNl
        S0  <- matrix(as.vector(cov(sim)), nrow(S), d * d, byrow = TRUE)
        tl  <- tg[rep(seq_len(G), length(hs)), , drop = FALSE]
        for (e in which(procScale[k + 1, ] > 0)) {
          pv <- procScale[k + 1, e] / tl[, poolOf[e]]
          S[, mi(e, e, d)]  <- S[, mi(e, e, d)] + pv
          S0[, mi(e, e, d)] <- S0[, mi(e, e, d)] + pv
        }
        rows <- k * P + rowsP
        muF[rows, ]  <- MU; muF0[rows, ]  <- matrix(ctr, length(rows), d, byrow = TRUE)
        sigF[rows, ] <- S;  sigF0[rows, ] <- S0
      }
    }

    # A component whose importance weights rest on fewer effective samples than there are
    # state elements can have a singular weighted covariance (e.g. leaf / stem of a
    # disturbed site, which get no process variance). Like the collapsed components above,
    # it gets the unweighted moments as a placeholder and zero weight below.
    pivots <- chol_batch(sigF, d)[, mi(1:d, 1:d, d), drop = FALSE]
    badC   <- !(rowSums(is.finite(pivots)) == d) |
              matrixStats::rowMins(pivots^2 / pmax(sigF[, mi(1:d, 1:d, d), drop = FALSE], 1e-300), na.rm = TRUE) < 1e-10
    muF[badC, ] <- muF0[badC, ]; sigF[badC, ] <- sigF0[badC, ]
    rm(muF0, sigF0)

    beta   <- dist.p[, t]
    piPrev <- exp(qs - max(qs)); piPrev <- piPrev / sum(piPrev)
    # drop the collapsed and singular components; a grid point left with no usable
    # component keeps its placeholder moments at their prior weights rather than dropping
    # out (filter-diagnostics.R)
    dr     <- apply_drops(outer(prevW, beta), dropF, badC, piPrev)        # G x H x K1
    wPrior <- dr$w
    hist$wDrop[t]    <- dr$massNet                                       # forecast weight dropped
    hist$wDropIS[t]  <- dr$massIS
    hist$wDropCov[t] <- dr$massCov
    hist$nNone[t]    <- sum(dr$none)
    hist$nBad[t]     <- sum(badC)

    wF <- as.vector(wPrior * piPrev); wF <- wF / sum(wF)
    if (forecastQuantiles)
      for (a in 1:d) hist$fq[, t, a] <- mix_quantile(qs5, wF, muF[, a], sqrt(sigF[, mi(a, a, d)]))

    # 5. Kalman update of every component on the observed elements
    o    <- which(!is.na(Y[t, ]))
    muA  <- muF; sigA <- sigF; logL <- numeric(Cn)
    if (length(o)) {
      ka  <- kalman_batch(muF, sigF, Y[t, ], R, o)
      muA <- ka$mu; sigA <- ka$sig; logL <- ka$logL
      for (a in o) hist$pit[t, a] <- sum(wF * pnorm(Y[t, a], muF[, a], sqrt(sigF[, mi(a, a, d)] + R[a, a])))
    }

    # 6. tau posterior on the grid
    lwl  <- matrix(log(as.vector(wPrior)) + logL, G)                     # G x (H K1)
    logz <- matrixStats::rowLogSumExps(lwl)
    if (any(!is.finite(logz))) stop(sprintf("t=%d: zero or non-finite marginal likelihood at %d tau grid points", t, sum(!is.finite(logz))))
    qs <- qs + logz
    ws <- exp(qs - max(qs)); ws <- ws / sum(ws)
    hist$procVar[t, ] <- colSums(ws / tg)
    hist$tauJoint[[t]] <- list(idx = gridIdx[active, , drop = FALSE], w = ws)
    for (a in 1:3) {
      m <- numeric(length(tauAxes[[a]]))
      s <- rowsum(ws, gridIdx[active, a])
      m[as.integer(rownames(s))] <- s
      hist$tauMarg[[a]][, t] <- m
    }

    # analysis weights and summaries
    wA   <- exp(lwl - logz)                                              # G x (H K1)
    wAll <- as.vector(wA * ws)
    hist$pK[t, ] <- colSums(matrix(colSums(wA * ws), H))
    for (a in 1:d) {
      hist$mean[t, a] <- sum(wAll * muA[, a])
      hist$q[, t, a]  <- mix_quantile(qs5, wAll, muA[, a], sqrt(sigA[, mi(a, a, d)]))
    }
    if (t %in% keepMix) {
      kF <- wF > 1e-8; kA <- wAll > 1e-8
      hist$mixF[[as.character(t)]] <- list(w = wF[kF], mu = muF[kF, , drop = FALSE], sig = sigF[kF, , drop = FALSE])
      hist$mixA[[as.character(t)]] <- list(w = wAll[kA], mu = muA[kA, , drop = FALSE], sig = sigA[kA, , drop = FALSE])
    }

    # 7. project back onto one Gaussian per history group, per tau grid point
    out <- project_groups(wA, array(muA, c(G, H * K1, d)), array(sigA, c(G, H * K1, d * d)), grpOf, H)
    hist$pGroup[t, ] <- colSums(out$w * ws)

    keep    <- qs >= max(qs) - pruneLog
    active  <- active[keep]
    qs      <- qs[keep]
    prevW   <- out$w[keep, , drop = FALSE]
    prevMu  <- out$mu[keep, , , drop = FALSE]
    prevSig <- out$sig[keep, , , drop = FALSE]
  }
  hist
}

## draw n samples from a stored mixture (w, mu: m x d, sig: m x d^2)
rmix <- function(mix, n) {
  d  <- ncol(mix$mu)
  ci <- sample.int(length(mix$w), n, replace = TRUE, prob = mix$w)
  t(sapply(ci, function(c) mvnfast::rmvn(1, mix$mu[c, ], matrix(mix$sig[c, ], d, d))))
}

## per-axis process-sd posterior over time (one panel per tau axis)
plot_process_sd <- function(fit, truth, axis.names, yr = seq_len(ncol(fit$tauMarg[[1]]))) {
  band <- function(x, lo, hi, col) polygon(c(x, rev(x)), c(lo, rev(hi)), border = NA, col = col)
  op <- par(mfrow = c(3, 1), mar = c(3.5, 4.5, 2, 1)); on.exit(par(op))
  for (a in 1:3) {
    tauQ <- apply(fit$tauMarg[[a]], 2, function(w) {
      cdf <- cumsum(w); fit$tauAxes[[a]][sapply(c(0.025, 0.5, 0.975), function(p) which(cdf >= p)[1])]
    })
    sdQ <- 1 / sqrt(tauQ)
    plot(yr, sdQ[2, ], type = 'n', log = "y", ylim = range(1 / sqrt(fit$tauAxes[[a]])), xlab = "year",
         ylab = "process sd", main = sprintf("%s process-error posterior", axis.names[a]))
    band(yr, sdQ[3, ], sdQ[1, ], adjustcolor("#2a78d6", 0.3))
    lines(yr, sdQ[2, ], lwd = 2, col = "#2a78d6")
    lines(yr, sqrt(fit$procVar[, a]), lty = 2, col = "#2a78d6")
    abline(h = truth[a], lwd = 2)
    if (a == 1) legend("topright", legend = c("posterior median", "95% interval", "sqrt(E[1/tau])", "truth"),
                       col = c("#2a78d6", adjustcolor("#2a78d6", 0.3), "#2a78d6", 1), lwd = c(2, 10, 1, 2),
                       lty = c(1, 1, 2, 1), bg = "white", cex = 0.8)
  }
}

## analysis 95% band, median, truth and observations of one state element
plot_state <- function(fit, a, truth, obs, main, ylab, dist.years = integer(0), known = NULL) {
  yr <- seq_len(nrow(truth))
  rng <- range(c(fit$q[c(1, 5), , a], obs[, a], truth[, a]), na.rm = TRUE)
  plot(yr, obs[, a], type = 'n', ylim = rng, xlab = "year", ylab = ylab, main = main)
  polygon(c(yr, rev(yr)), c(fit$q[1, , a], rev(fit$q[5, , a])), border = NA, col = adjustcolor("blue", 0.3))
  lines(yr, fit$q[3, , a], col = "blue", lwd = 2)
  if (!is.null(known)) lines(yr, known$q[3, , a], col = "orange", lwd = 2, lty = 3)
  lines(yr, truth[, a], lty = 2, lwd = 2)
  points(yr, obs[, a], pch = 18, col = 2)
  abline(v = dist.years, col = "grey60")
  legend("topright", legend = c("analysis 95%", "median (tau learned)", if (!is.null(known)) "median (tau known)", "truth", "data"),
         col = c(adjustcolor("blue", 0.3), "blue", if (!is.null(known)) "orange", 1, 2),
         lwd = c(10, 2, if (!is.null(known)) 2, 2, NA), lty = c(1, 1, if (!is.null(known)) 3, 2, NA),
         pch = c(NA, NA, if (!is.null(known)) NA, NA, 18), bg = "white", cex = 0.8)
}
