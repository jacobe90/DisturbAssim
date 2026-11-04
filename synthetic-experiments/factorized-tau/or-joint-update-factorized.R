## ORJointUpdate with a factorized tau posterior, shared by the factorized-tau experiments.
##
## The filter of ../or-joint-update-filter.R carries (tau_1, tau_2, tau_3) on a 3-D product
## grid: G_1 G_2 G_3 grid points, each with its own state mixture. Here each tau gets its own
## 1-D grid instead. Block a holds G_a grid points: tau_a takes its grid values, and the other
## two taus are held at their current posterior estimates. Each block is a filter of its own
## (state mixtures per grid point, weights updated with that grid point's marginal likelihood)
## whose weights give the posterior of tau_a, so the run carries G_1 + G_2 + G_3 grid points.
##
## The plug-in estimate of tau_b used by the other blocks is 1 / E[1/tau_b] (the posterior
## mean of the process variance), from block b's posterior after the previous step. State
## summaries (quantiles, scenario probabilities, ...) average the three blocks' posteriors
## with equal weight.
##
## Usage, from the repository root:
##   source("synthetic-experiments/factorized-tau/or-joint-update-factorized.R")
## which also loads ../or-joint-update-filter.R (the batched algebra and ORJointUpdateGeneral,
## the 3-D grid version used for comparison).

source("synthetic-experiments/or-joint-update-filter.R")

##' Approximate joint disturbance + process-precision update with a factorized tau posterior
##'
##' Arguments as for ORJointUpdateGeneral (see ../or-joint-update-filter.R); `grid$G[a]` is
##' the size of tau_a's 1-D grid.
##' @param plugin estimate of each tau held fixed by the other blocks: "meanvar" (1 / E[1/tau],
##'   the default), "median" or "mean"
##' @return as ORJointUpdateGeneral (tauMarg is each block's own 1-D posterior), plus
##'   tauHat (NT x 3, the plug-in values used at each step) and pKaxis (NT x (K+1) x 3, the
##'   scenario probabilities of each block)
ORJointUpdateFactorized <- function(Y, R, Xic, simulate, nMembers, labels, disturb, dist.p,
                                    procScale, transition, w0, poolOf, prior, grid, N = 1000,
                                    pruneLog = 30, keepMix = integer(0), forecastQuantiles = TRUE,
                                    plugin = c("meanvar", "median", "mean")) {
  plugin <- match.arg(plugin)
  NT <- nrow(Y); d <- ncol(Y)
  H  <- length(labels); K1 <- nrow(dist.p)
  ia <- rep(1:d, d); ib <- rep(1:d, each = d)  # column (a, b) of a flattened d x d

  tauAxes <- lapply(1:3, function(a) exp(seq(log(grid$tauMin[a]), log(grid$tauMax[a]), length.out = grid$G[a])))
  logMass <- lapply(tauAxes, function(x) log(cell_prior_mass(x, prior)))
  axisOf  <- rep(1:3, lengths(tauAxes))         # block (tau axis) of each grid point
  idxOf   <- unlist(lapply(tauAxes, seq_along)) # its index on that axis
  Gt      <- length(axisOf)
  qs      <- unlist(logMass)                    # log posterior weight within its block

  ## plug-in estimate of one tau from weights w on its axis
  estimate <- function(a, w) {
    x <- tauAxes[[a]]; w <- w / sum(w)
    switch(plugin, meanvar = 1 / sum(w / x), mean = sum(w * x),
           median = x[which(cumsum(w) >= 0.5)[1]])
  }
  tauHat <- sapply(1:3, function(a) estimate(a, exp(logMass[[a]])))
  tau_rows <- function(rows) {                  # tau vector of each grid point
    tm <- matrix(tauHat, length(rows), 3, byrow = TRUE)
    tm[cbind(seq_along(rows), axisOf[rows])] <- tauAxes_val(rows)
    tm
  }
  tauAxes_val <- function(rows) mapply(function(a, i) tauAxes[[a]][i], axisOf[rows], idxOf[rows])

  active  <- seq_len(Gt)
  prevW   <- matrix(w0 / sum(w0), Gt, H, byrow = TRUE)
  prevMu  <- array(rep(colMeans(Xic), each = Gt * H), c(Gt, H, d))
  prevSig <- array(rep(as.vector(cov(Xic)), each = Gt * H), c(Gt, H, d * d))
  grpOf   <- as.vector(transition)              # new group of component (h, k), h fastest

  hist <- list(
    mean = matrix(NA_real_, NT, d), q = array(NA_real_, c(5, NT, d)), fq = array(NA_real_, c(5, NT, d)),
    pK = matrix(NA_real_, NT, K1), pKaxis = array(NA_real_, c(NT, K1, 3)),
    pGroup = matrix(NA_real_, NT, H), pit = matrix(NA_real_, NT, d),
    wDrop = rep(NA_real_, NT), procVar = matrix(NA_real_, NT, 3), tauHat = matrix(NA_real_, NT, 3),
    tauAxes = tauAxes, tauMarg = lapply(tauAxes, function(x) matrix(NA_real_, length(x), NT)),
    nActive = rep(NA_integer_, NT), nBad = rep(NA_integer_, NT), mixF = list(), mixA = list()
  )

  ## weights of the active grid points: normalised within each block, then 1/3 per block,
  ## so that every quantity marginalised over the grid averages the three blocks
  block_weights <- function(qs, ax) {
    w <- numeric(length(qs))
    for (a in unique(ax)) { s <- ax == a; e <- exp(qs[s] - max(qs[s])); w[s] <- e / sum(e) }
    w / length(unique(ax))
  }

  for (t in 1:NT) {
    G  <- length(active)
    ax <- axisOf[active]
    tg <- tau_rows(active)
    hist$tauHat[t, ] <- tauHat
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

    # 4. forecast moments of every (g, h, k), with each grid point's own tau vector
    muF  <- muF0  <- matrix(0, Cn, d)
    sigF <- sigF0 <- matrix(0, Cn, d * d)
    for (lab in labs) {
      hs   <- which(labels == lab)
      rowsP <- as.vector(outer(seq_len(G), (hs - 1) * G, "+"))
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

    # components with a singular weighted covariance: placeholder moments, zero weight
    pivots <- chol_batch(sigF, d)[, mi(1:d, 1:d, d), drop = FALSE]
    badC   <- !(rowSums(is.finite(pivots)) == d) |
              matrixStats::rowMins(pivots^2 / pmax(sigF[, mi(1:d, 1:d, d), drop = FALSE], 1e-300), na.rm = TRUE) < 1e-10
    muF[badC, ] <- muF0[badC, ]; sigF[badC, ] <- sigF0[badC, ]
    rm(muF0, sigF0)

    beta   <- dist.p[, t]
    piPrev <- block_weights(qs, ax)
    wPrior <- wPrior0 <- outer(prevW, beta)                              # G x H x K1
    for (k in 1:K1) wPrior[, , k][dropF] <- 0
    wPrior[badC] <- 0
    none <- rowSums(matrix(wPrior, G)) == 0
    if (any(none)) wPrior[none, , ] <- wPrior0[none, , ]
    hist$wDrop[t] <- 1 - sum(as.vector(wPrior) * piPrev)
    hist$nBad[t]  <- sum(badC)

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

    # 6. each block's tau posterior on its own 1-D grid
    lwl  <- matrix(log(as.vector(wPrior)) + logL, G)                     # G x (H K1)
    logz <- matrixStats::rowLogSumExps(lwl)
    if (any(!is.finite(logz))) stop(sprintf("t=%d: zero or non-finite marginal likelihood at %d tau grid points", t, sum(!is.finite(logz))))
    qs <- qs + logz
    ws <- block_weights(qs, ax)                                          # sums to 1 over blocks
    for (a in 1:3) {
      s <- ax == a
      m <- numeric(length(tauAxes[[a]])); m[idxOf[active[s]]] <- ws[s] / sum(ws[s])
      hist$tauMarg[[a]][, t] <- m
      hist$procVar[t, a]     <- sum(m / tauAxes[[a]])
    }

    # analysis weights and summaries (the three blocks averaged)
    wA   <- exp(lwl - logz)                                              # G x (H K1)
    wAll <- as.vector(wA * ws)
    hist$pK[t, ] <- colSums(matrix(colSums(wA * ws), H))
    for (a in 1:3) {
      s <- ax == a
      hist$pKaxis[t, , a] <- colSums(matrix(colSums(wA[s, , drop = FALSE] * ws[s] / sum(ws[s])), H))
    }
    for (a in 1:d) {
      hist$mean[t, a] <- sum(wAll * muA[, a])
      hist$q[, t, a]  <- mix_quantile(qs5, wAll, muA[, a], sqrt(sigA[, mi(a, a, d)]))
    }
    if (t %in% keepMix) {
      kF <- wF > 1e-8; kA <- wAll > 1e-8
      hist$mixF[[as.character(t)]] <- list(w = wF[kF], mu = muF[kF, , drop = FALSE], sig = sigF[kF, , drop = FALSE])
      hist$mixA[[as.character(t)]] <- list(w = wAll[kA], mu = muA[kA, , drop = FALSE], sig = sigA[kA, , drop = FALSE])
    }

    # 7. project back onto one Gaussian per history group, per grid point
    out <- project_groups(wA, array(muA, c(G, H * K1, d)), array(sigA, c(G, H * K1, d * d)), grpOf, H)
    hist$pGroup[t, ] <- colSums(out$w * ws)

    # plug-in estimates for the next step, from each block's posterior
    tauHat <- sapply(1:3, function(a) estimate(a, hist$tauMarg[[a]][, t]))

    # prune within each block
    keep <- qs >= ave(qs, ax, FUN = max) - pruneLog
    active  <- active[keep]
    qs      <- qs[keep]
    prevW   <- out$w[keep, , drop = FALSE]
    prevMu  <- out$mu[keep, , , drop = FALSE]
    prevSig <- out$sig[keep, , , drop = FALSE]
  }
  hist
}

## ---- Comparison report shared by the experiment scripts ---------------------------------
## fits: named list of filter outputs (e.g. factorized, joint, known); times: run time of each
## (seconds); events: function(fit) -> NT x E matrix of named event probabilities to compare
## (default P(disturbed)); state.names / pool.names label the state elements and tau axes.
compare_report <- function(fits, times, truth, obs, state.names, pool.names, truth.sd, dist.years,
                           outdir, title, events = function(fit) cbind(`P(disturbed)` = 1 - fit$pK[, 1])) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  NT <- nrow(truth); d <- ncol(truth); yr <- 1:NT
  cols <- c(factorized = "#d95f02", joint = "#1b78b4", known = "grey35")[names(fits)]
  lty  <- c(factorized = 1, joint = 1, known = 3)[names(fits)]
  ev   <- lapply(fits, events)
  rmse <- sapply(fits, function(f) sqrt(colMeans((f$q[3, , ] - truth)^2)))
  rownames(rmse) <- state.names

  pdf(file.path(outdir, "figures.pdf"), width = 7, height = 6)
  ## events
  E <- ncol(ev[[1]])
  op <- par(mfrow = c(E, 1), mar = c(3.6, 4.4, 2.4, 1), oma = c(0, 0, 2, 0))
  for (e in 1:E) {
    plot(yr, ev[[1]][, e], type = 'n', ylim = c(0, 1), xlab = "year", ylab = "probability",
         main = colnames(ev[[1]])[e])
    abline(v = dist.years, col = "grey70")
    for (m in names(fits)) lines(yr, ev[[m]][, e], type = 'o', pch = 19, cex = 0.5, col = cols[m], lty = lty[m], lwd = 1.8)
    if (e == 1) legend("topright", legend = names(fits), col = cols, lty = lty, lwd = 2, bty = "n", cex = 0.85)
  }
  mtext(sprintf("%s: posterior event probabilities", title), outer = TRUE, font = 2)
  par(op)

  ## states, three per page: factorized and joint 95% bands, medians, truth, data
  for (p0 in seq(1, d, by = 3)) {
    op <- par(mfrow = c(3, 1), mar = c(3.6, 4.4, 2.4, 1))
    for (a in p0:min(d, p0 + 2)) {
      rng <- range(c(sapply(fits, function(f) f$q[c(1, 5), , a]), obs[, a], truth[, a]), na.rm = TRUE)
      plot(yr, obs[, a], type = 'n', ylim = rng, xlab = "year", ylab = state.names[a], main = state.names[a])
      for (m in intersect(c("joint", "factorized"), names(fits)))
        polygon(c(yr, rev(yr)), c(fits[[m]]$q[1, , a], rev(fits[[m]]$q[5, , a])), border = NA, col = adjustcolor(cols[m], 0.22))
      for (m in names(fits)) lines(yr, fits[[m]]$q[3, , a], col = cols[m], lty = lty[m], lwd = 1.8)
      lines(yr, truth[, a], lty = 2, lwd = 1.5); points(yr, obs[, a], pch = 18, col = 2)
      abline(v = dist.years, col = "grey70")
      if (a == p0) legend("topright", bty = "n", cex = 0.75, ncol = 2,
                          legend = c(paste(names(fits), "median"), "95%: factorized / joint", "truth", "data"),
                          col = c(cols, "grey50", 1, 2), lty = c(lty, 1, 2, NA), lwd = c(rep(1.8, length(fits)), 8, 1.5, NA),
                          pch = c(rep(NA, length(fits)), NA, NA, 18))
    }
    par(op)
  }

  ## process sd posterior of each tau over time: factorized vs joint marginals
  learned <- intersect(c("factorized", "joint"), names(fits))
  op <- par(mfrow = c(3, 1), mar = c(3.6, 4.4, 2.4, 1))
  for (a in 1:3) {
    q <- lapply(fits[learned], function(f) 1 / sqrt(apply(f$tauMarg[[a]], 2, function(w) {
      cdf <- cumsum(w); f$tauAxes[[a]][sapply(c(0.025, 0.5, 0.975), function(p) which(cdf >= p)[1])]
    })))
    plot(yr, q[[1]][2, ], type = 'n', log = "y", ylim = range(1 / sqrt(fits[[learned[1]]]$tauAxes[[a]])),
         xlab = "year", ylab = "process sd", main = sprintf("%s process sd (95%% interval, median)", pool.names[a]))
    for (m in learned) {
      polygon(c(yr, rev(yr)), c(q[[m]][3, ], rev(q[[m]][1, ])), border = NA, col = adjustcolor(cols[m], 0.22))
      lines(yr, q[[m]][2, ], col = cols[m], lwd = 2)
    }
    if (!is.null(fits$factorized$tauHat)) lines(yr, 1 / sqrt(fits$factorized$tauHat[, a]), col = cols["factorized"], lty = 3)
    abline(h = truth.sd[a], lwd = 1.5)
    if (a == 1) legend("topright", bty = "n", cex = 0.8, legend = c(learned, "factorized plug-in", "truth"),
                       col = c(cols[learned], cols["factorized"], 1), lwd = c(2, 2, 1, 1.5), lty = c(1, 1, 3, 1))
  }
  par(op)

  ## final-year marginal posterior of each tau, factorized vs joint
  op <- par(mfrow = c(1, 3), mar = c(4.4, 4, 2.6, 0.6), oma = c(0, 0, 2, 0))
  for (a in 1:3) {
    ymax <- max(sapply(fits[learned], function(f) max(f$tauMarg[[a]][, NT])))
    plot(NA, xlim = range(1 / sqrt(fits[[learned[1]]]$tauAxes[[a]])), ylim = c(0, ymax), log = "x",
         xlab = "process sd", ylab = "posterior mass", main = pool.names[a])
    for (m in learned) lines(1 / sqrt(fits[[m]]$tauAxes[[a]]), fits[[m]]$tauMarg[[a]][, NT], type = 'o', pch = 19,
                             cex = 0.5, col = cols[m], lwd = 1.8)
    abline(v = truth.sd[a], lwd = 1.5)
    if (a == 1) legend("topleft", bty = "n", cex = 0.8, legend = learned, col = cols[learned], lwd = 2)
  }
  mtext(sprintf("Marginal tau posterior, year %d", NT), outer = TRUE, font = 2)
  par(op)
  dev.off()

  ## summary
  fmt <- function(x, nm) paste(sprintf("%s %.3f", nm, x), collapse = ", ")
  other <- setdiff(yr, dist.years)
  lines.out <- c(title, "", "run time (s): " %+% fmt(times, names(times)), "",
                 "analysis-median RMSE:", sapply(names(fits), function(m) sprintf("  %-10s %s", m, fmt(rmse[, m], state.names))), "")
  for (e in 1:E) {
    lines.out <- c(lines.out, sprintf("%s:", colnames(ev[[1]])[e]),
      sapply(names(fits), function(m) sprintf("  %-10s at year(s) %s: %s | max elsewhere %.3f",
             m, paste(dist.years, collapse = ","), paste(sprintf("%.3f", ev[[m]][dist.years, e]), collapse = ","),
             max(ev[[m]][other, e]))))
    if (all(c("factorized", "joint") %in% names(fits)))
      lines.out <- c(lines.out, sprintf("  max |factorized - joint| over years: %.3f", max(abs(ev$factorized[, e] - ev$joint[, e]))))
  }
  lines.out <- c(lines.out, "", sprintf("process sd, sqrt(E[1/tau]) at year %d (truth %s):", NT, paste(sprintf("%.3f", truth.sd), collapse = ", ")),
                 sapply(learned, function(m) sprintf("  %-10s %s", m, fmt(sqrt(fits[[m]]$procVar[NT, ]), pool.names))))
  writeLines(lines.out, file.path(outdir, "summary.txt"))
  cat(lines.out, sep = "\n")
  saveRDS(list(fits = lapply(fits, function(f) f[setdiff(names(f), c("mixF", "mixA"))]), times = times,
               truth = truth, obs = obs, rmse = rmse, events = ev), file.path(outdir, "results.rds"))
  invisible(list(rmse = rmse, events = ev))
}
`%+%` <- function(a, b) paste0(a, b)

## ---- Experiment helpers ------------------------------------------------------------------
## Evaluate `file` up to (not including) the first line matching `pattern`, in `envir`. The
## factorized experiments use it to reuse an experiment's set-up (drivers, ensemble, truth,
## filter arguments) exactly, so both filters see the same data.
source_until <- function(file, pattern, envir = parent.frame()) {
  lines <- readLines(file)
  n <- grep(pattern, lines)[1]
  if (is.na(n)) stop(sprintf("%s: no line matching '%s'", file, pattern))
  eval(parse(text = lines[seq_len(n - 1)], keep.source = FALSE), envir)
  invisible(n - 1)
}

## Run the factorized filter, and for comparison the 3-D grid filter (ORJointUpdateGeneral)
## and the known-tau filter, on the same arguments. Each run starts from the same RNG
## state, so differences between them come from the method, not the importance samples.
run_comparison <- function(args, grid, tau.known, joint = TRUE, seed = 1) {
  run <- function(fun, g) { set.seed(seed); tm <- system.time(f <- do.call(fun, c(args, list(grid = g)))); list(fit = f, time = unname(tm[3])) }
  out <- list(factorized = run(ORJointUpdateFactorized, grid))
  if (joint) out$joint <- run(ORJointUpdateGeneral, grid)
  out$known <- run(ORJointUpdateGeneral, list(tauMin = tau.known, tauMax = tau.known, G = rep(1, 3)))
  list(fits = lapply(out, `[[`, "fit"), times = sapply(out, `[[`, "time"))
}
