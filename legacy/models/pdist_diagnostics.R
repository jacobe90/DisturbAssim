## ---------------------------------------------------------------------------
## Diagnosing high P(disturbed | gSamp) at timesteps with no visible event.
##
## Notation (matching the write-up):
##   beta[k,t]   prior disturbance probabilities, dist$p[k+1, t]
##               beta[0,t] = 1 - p_t          (undisturbed)
##               beta[k,t] = p_t * pi_k       (class k, 1 <= k <= D)
##               sum_{k=0}^{D} beta[k,t] = 1
##   betaPlus[t] = sum_{k=1}^{D} beta[k,t] = p_t   (prior prob of ANY disturbance)
##   gamma[j]    = wPrev[t,g,j] / sum_j wPrev[t,g,j], the projected prior branch
##               weights (j = 1 <-> z_{t-1} = 0, j = 2 <-> z_{t-1} = 1)
##   l[j,k]      = lA[t,g,j,k+1] = N(y_t; H mu_f^{j,k}, H Sig_f^{j,k} H' + R)
##
## Analysis weights:  alpha[j,k] = gamma[j] beta[k] l[j,k] / Z_t.
##
## Posterior odds of disturbance in the gSamp slice:
##   O_t = P(z_t = 1 | tau*, D_t) / P(z_t = 0 | tau*, D_t)
##       = ( sum_j sum_{k>=1} alpha[j,k] ) / ( sum_j alpha[j,0] )
## Z_t cancels, giving the exact factorization
##   O_t  = (betaPlus / beta[0]) * BF_t
##   BF_t = ( sum_j gamma[j] sum_{k>=1} (beta[k]/betaPlus) l[j,k] )
##          / ( sum_j gamma[j] l[j,0] )
##
## Componentwise, for the modal disturbed (j_d,k_d) and modal undisturbed
## (j_u,0) components, the exact additive decomposition is
##   log alpha[j_d,k_d] - log alpha[j_u,0]
##     = [log gamma[j_d] - log gamma[j_u]]        (prior branch term)
##     + [log beta[k_d]  - log beta[0]]           (prior class term)
##     + [log l[j_d,k_d] - log l[j_u,0]]          (likelihood term)
##
## Array convention: k is 1-based in the stored arrays, class label 0 sits at
## index 1. innovMean / innovCov hold recycled scalars; read [1] and [1,1].
## ---------------------------------------------------------------------------

HL <- c(10, 17:20, 22:28)          # timesteps under investigation


## ---- extraction -------------------------------------------------------------
branch_frame <- function(hist, Rparams) {
  NT <- length(hist$gSamp)
  D1 <- dim(hist$wA)[4]

  rows <- expand.grid(k = seq_len(D1), j = 1:2, t = seq_len(NT))
  rows <- rows[, c("t", "j", "k")]

  out <- data.frame(rows,
                    klab  = rows$k - 1L,
                    g     = hist$gSamp[rows$t],
                    obs   = hist$obs[rows$t],
                    yhat  = NA_real_, vF = NA_real_, Rt = NA_real_,
                    sd    = NA_real_, z = NA_real_,
                    lA    = NA_real_, wA = NA_real_,
                    beta  = NA_real_, wPrev = NA_real_)

  for (i in seq_len(nrow(out))) {
    t <- out$t[i]; j <- out$j[i]; k <- out$k[i]; g <- out$g[i]
    if (is.na(g)) next
    out$yhat[i]  <- hist$innovMean[t, g, j, k, 1]
    out$vF[i]    <- hist$innovCov[t, g, j, k, 1, 1]
    out$Rt[i]    <- hist$obs[t] * Rparams[1] + Rparams[2]
    out$sd[i]    <- sqrt(out$vF[i] + out$Rt[i])
    out$z[i]     <- (out$obs[i] - out$yhat[i]) / out$sd[i]
    out$lA[i]    <- hist$lA[t, g, j, k]
    out$wA[i]    <- hist$wA[t, g, j, k]
    out$beta[i]  <- hist$beta[t, k]
    out$wPrev[i] <- hist$wPrev[t, g, j]
  }

  ## posterior weights renormalized within the gSamp slice
  out$wAn <- ave(out$wA, out$t, FUN = function(v) v / sum(v))

  ## gamma_j: prior branch weights normalized over j only
  out$gamma <- NA_real_
  for (t in unique(out$t)) {
    idx <- out$t == t
    w   <- tapply(out$wPrev[idx], out$j[idx], function(v) v[1])
    out$gamma[idx] <- (w / sum(w))[as.character(out$j[idx])]
  }
  out
}


odds_frame <- function(hist, bf) {
  NT <- length(hist$gSamp)

  res <- data.frame(
    t = seq_len(NT), obs = hist$obs,
    pDist = NA_real_, pDistPrior = NA_real_,
    logOddsPost = NA_real_, logOddsPrior = NA_real_, logBF = NA_real_,
    likShare = NA_real_, gamma2 = NA_real_,
    z0 = NA_real_, zWin = NA_real_,
    kWin = NA_integer_, jWin = NA_integer_, jUnd = NA_integer_,
    dGamma = NA_real_, dBeta = NA_real_, dLik = NA_real_, dTotal = NA_real_)

  for (t in seq_len(NT)) {
    d <- bf[bf$t == t, ]
    if (all(is.na(d$wA))) next

    ## -- posterior / prior odds ----------------------------------------------
    pd <- sum(d$wAn[d$k > 1])
    b  <- hist$beta[t, ]
    bp <- sum(b[-1]) / sum(b)

    res$pDist[t]        <- pd
    res$pDistPrior[t]   <- bp
    res$logOddsPost[t]  <- log(pd) - log1p(-pd)
    res$logOddsPrior[t] <- log(sum(b[-1])) - log(b[1])
    res$logBF[t]        <- res$logOddsPost[t] - res$logOddsPrior[t]

    ## -- marginal likelihood share of the disturbed branches -----------------
    gam <- tapply(d$gamma, d$j, function(v) v[1])
    num <- sum(sapply(1:2, function(j)
      gam[j] * sum((b[-1] / sum(b[-1])) * d$lA[d$j == j & d$k > 1])))
    den <- sum(sapply(1:2, function(j) gam[j] * d$lA[d$j == j & d$k == 1]))
    res$likShare[t] <- num / (num + den)
    res$gamma2[t]   <- gam[2]

    ## -- modal components and the additive log-weight gap --------------------
    du <- d[d$k == 1, ];  u <- du[which.max(du$wA), ]
    dd <- d[d$k >  1, ];  w <- dd[which.max(dd$wA), ]

    res$jUnd[t] <- u$j; res$jWin[t] <- w$j; res$kWin[t] <- w$klab
    res$z0[t]   <- u$z; res$zWin[t] <- w$z

    res$dGamma[t] <- log(w$gamma) - log(u$gamma)
    res$dBeta[t]  <- log(w$beta)  - log(u$beta)
    res$dLik[t]   <- log(w$lA)    - log(u$lA)
    res$dTotal[t] <- res$dGamma[t] + res$dBeta[t] + res$dLik[t]
  }
  res
}


## ---- plotting helpers -------------------------------------------------------
.page <- function(right = 13) {
  par(mfrow = c(1, 1), mar = c(4.5, 5, 4, right),
      cex.main = 1.25, cex.lab = 1.1, cex.axis = 1.0)
}

.legend_out <- function(...) {
  op <- par(xpd = NA); on.exit(par(op))
  u <- par("usr")
  legend(x = u[2] + 0.03 * (u[2] - u[1]), y = u[4], bty = "n", cex = 0.95, ...)
}

.shade <- function(hl, col = rgb(1, 0, 0, 0.10)) {
  if (!length(hl)) return(invisible(NULL))
  runs <- split(hl, cumsum(c(1, diff(hl) != 1)))
  for (r in runs) rect(min(r) - 0.5, par("usr")[3], max(r) + 0.5, par("usr")[4],
                       col = col, border = NA)
}

## stacked bars that handle mixed signs (positive up, negative down)
.signed_stack <- function(x, M, cols, hl, width = 0.75, ylab = "", main = "") {
  pos <- pmax(M, 0); neg <- pmin(M, 0)
  pos[!is.finite(pos)] <- 0; neg[!is.finite(neg)] <- 0
  ylim <- range(c(0, rowSums(pos), rowSums(neg)))
  plot(range(x) + c(-1, 1), ylim, type = "n", xlab = "t", ylab = ylab, main = main)
  .shade(hl)
  abline(h = 0, col = "grey40")
  for (i in seq_along(x)) {
    b <- 0
    for (m in seq_len(ncol(M))) {
      h <- pos[i, m]
      if (h > 0) {
        rect(x[i] - width / 2, b, x[i] + width / 2, b + h, col = cols[m], border = NA)
        b <- b + h
      }
    }
    b <- 0
    for (m in seq_len(ncol(M))) {
      h <- neg[i, m]
      if (h < 0) {
        rect(x[i] - width / 2, b + h, x[i] + width / 2, b, col = cols[m], border = NA)
        b <- b + h
      }
    }
  }
  box()
}


## ---- page 1 -----------------------------------------------------------------
plot_obs_analysis <- function(hist, od, hl = HL) {
  .page()
  plot(od$t, hist$obs, type = "n",
       ylim = range(c(hist$lo, hist$hi, hist$obs), na.rm = TRUE),
       xlab = "t", ylab = "wood / AGB pool",
       main = "Observation and analysis of the modal component at gSamp")
  .shade(hl)
  polygon(c(od$t, rev(od$t)), c(hist$lo, rev(hist$hi)),
          col = rgb(0, 0, 1, 0.15), border = NA)
  lines(od$t, hist$mean, col = "blue", lwd = 2)
  lines(od$t, hist$obs, type = "o", pch = 16, cex = 0.8)
  .legend_out(legend = c("observation", "analysis mean", "95% band", "flagged t"),
              col = c("black", "blue", rgb(0, 0, 1, 0.3), rgb(1, 0, 0, 0.3)),
              lty = c(1, 1, NA, NA), pch = c(16, NA, 15, 15),
              lwd = c(1, 2, NA, NA))
}


## ---- page 2 -----------------------------------------------------------------
## Posterior disturbance probability marginalized over tau, rather than read off
## the single sampled slot. With ws[t,g] = hist$tauProbs[t,g] approximating
## p(tau_g | D_t),
##
##   P(z_t = 1 | D_t) = sum_g p(tau_g | D_t) P(z_t = 1 | tau_g, D_t)
##                    = sum_g ws[t,g] * sum_{j} sum_{k>=1} alpha[t,g,j,k]
##
## The band is +/- 1 sd of P(z_t = 1 | tau_g, D_t) under the same weights, so it
## shows how much the answer depends on tau. The gSamp line is what every other
## page in this file conditions on; where it departs from the average, the
## single-slot diagnostics are reading an unrepresentative tau.
pdist_by_tau <- function(hist) {
  NT <- dim(hist$wA)[1]; G <- dim(hist$wA)[2]
  pg <- matrix(NA_real_, NT, G)
  for (t in seq_len(NT)) for (g in seq_len(G)) {
    W <- hist$wA[t, g, , ]
    s <- sum(W)
    if (is.finite(s) && s > 0) pg[t, g] <- sum(W[, -1]) / s
  }
  pg
}

plot_pdist_tau_averaged <- function(hist, od, hl = HL) {
  .page()
  NT <- dim(hist$wA)[1]
  t  <- seq_len(NT)
  pg <- pdist_by_tau(hist)

  if (is.null(hist$tauProbs)) {
    warning("hist$tauProbs missing; falling back to a uniform tau weighting")
    ws <- matrix(1 / ncol(pg), NT, ncol(pg))
  } else {
    ws <- hist$tauProbs
    ws <- ws / rowSums(ws, na.rm = TRUE)
  }
  ws[is.na(pg)] <- NA

  pbar <- rowSums(ws * pg, na.rm = TRUE) / rowSums(ws, na.rm = TRUE)
  psd  <- sqrt(rowSums(ws * (pg - pbar)^2, na.rm = TRUE) /
                 rowSums(ws, na.rm = TRUE))
  atG  <- pg[cbind(t, hist$gSamp)]

  plot(t, pbar, type = "n", ylim = c(0, 1), xlab = "t",
       ylab = "P(z_t = 1 | D_t)",
       main = "Posterior disturbance probability, marginalized over tau")
  .shade(hl)
  abline(h = 0.5, col = "grey60", lty = 3)

  ok <- is.finite(pbar) & is.finite(psd)
  polygon(c(t[ok], rev(t[ok])),
          c(pmin(pbar[ok] + psd[ok], 1), rev(pmax(pbar[ok] - psd[ok], 0))),
          col = rgb(0, 0, 0, 0.12), border = NA)
  lines(t, pbar, lwd = 2.5, type = "o", pch = 16, cex = 0.8)
  lines(t, atG, col = "steelblue3", lwd = 2, lty = 2)
  lines(t, od$pDistPrior, col = "darkorange", lwd = 2)

  .legend_out(legend = c("tau-averaged posterior",
                         "+/- 1 sd across tau",
                         "at gSamp(t) only",
                         "prior  betaPlus = p_t"),
              col = c("black", rgb(0, 0, 0, 0.3), "steelblue3", "darkorange"),
              lty = c(1, NA, 2, 1), pch = c(16, 15, NA, NA),
              lwd = c(2.5, NA, 2, 2))

  invisible(data.frame(t = t, pTauAvg = pbar, pTauSd = psd, pAtGSamp = atG))
}


## ---- page 3 -----------------------------------------------------------------
plot_odds_decomposition <- function(od, hl = HL) {
  .page()
  yl <- range(c(od$logOddsPost, od$logOddsPrior, od$logBF), finite = TRUE)
  plot(od$t, od$logOddsPost, type = "n", ylim = yl, xlab = "t",
       ylab = "log odds",
       main = "log O_t  =  log(betaPlus / beta_0)  +  log BF_t")
  .shade(hl)
  abline(h = 0, col = "grey40", lty = 3)
  lines(od$t, od$logOddsPost,  type = "o", pch = 16, cex = 0.8, lwd = 2)
  lines(od$t, od$logOddsPrior, col = "darkorange", lty = 2, lwd = 2)
  lines(od$t, od$logBF,        col = "purple", lwd = 2)
  .legend_out(legend = c("posterior log odds", "prior log odds", "log Bayes factor"),
              col = c("black", "darkorange", "purple"),
              lty = c(1, 2, 1), lwd = 2)
}


## ---- page 4 -----------------------------------------------------------------
plot_innovation_z <- function(bf, hl = HL) {
  .page()
  d <- bf[!is.na(bf$z), ]
  plot(d$t, d$z, type = "n", xlab = "t", ylab = "standardized innovation z",
       main = "z = (y_t - H mu_f) / sqrt(H Sig_f H' + R), best branch per class")
  .shade(hl)
  abline(h = c(-2, 0, 2), lty = c(3, 1, 3), col = "grey50")
  ks <- sort(unique(d$klab))
  for (i in seq_along(ks)) {
    dk <- aggregate(z ~ t, d[d$klab == ks[i], ],
                    function(v) v[which.min(abs(v))])
    lines(dk$t, dk$z, col = i, type = "o", pch = 16, cex = 0.7,
          lwd = if (ks[i] == 0) 3 else 1.2)
  }
  .legend_out(legend = c(paste0("k = ", ks), "|z| = 2"),
              col = c(seq_along(ks), "grey50"),
              lty = c(rep(1, length(ks)), 3),
              lwd = c(ifelse(ks == 0, 3, 1.2), 1))
}


## ---- page 5 -----------------------------------------------------------------
plot_predictive_sd <- function(bf, hl = HL) {
  .page()
  d <- aggregate(sd ~ t + klab, bf, mean)
  plot(d$t, d$sd, type = "n", xlab = "t", ylab = "predictive sd",
       main = "Predictive sd  sqrt(H Sig_f H' + R)  by class")
  .shade(hl)
  ks <- sort(unique(d$klab))
  for (i in seq_along(ks)) {
    dk <- d[d$klab == ks[i], ]
    lines(dk$t, dk$sd, col = i, type = "o", pch = 16, cex = 0.7,
          lwd = if (ks[i] == 0) 3 else 1.2)
  }
  .legend_out(legend = paste0("k = ", ks), col = seq_along(ks),
              lty = 1, lwd = ifelse(ks == 0, 3, 1.2))
}


## ---- page 6 -----------------------------------------------------------------
plot_class_composition <- function(bf, hl = HL) {
  .page()
  D1 <- max(bf$k)
  M  <- t(sapply(seq_len(D1), function(k)
    tapply(bf$wAn[bf$k == k], bf$t[bf$k == k], sum)))
  barplot(M, col = seq_len(D1), border = NA, space = 0,
          xlab = "t", ylab = "posterior weight",
          main = "Posterior mass by class at gSamp")
  .legend_out(legend = paste0("k = ", seq_len(D1) - 1L), fill = seq_len(D1),
              border = NA)
}


## ---- page 7 -----------------------------------------------------------------
plot_factor_shares <- function(od, hl = HL) {
  .page()
  plot(od$t, od$pDist, type = "n", ylim = c(0, 1), xlab = "t",
       ylab = "share attributable to the disturbed branches",
       main = "Relative pull of each factor toward disturbance")
  .shade(hl)
  abline(h = 0.5, col = "grey60", lty = 3)
  lines(od$t, od$pDistPrior, col = "darkorange",  lwd = 2)
  lines(od$t, od$likShare,   col = "purple",      lwd = 2)
  lines(od$t, od$gamma2,     col = "forestgreen", lwd = 2)
  lines(od$t, od$pDist, lwd = 2.5, lty = 2)
  .legend_out(legend = c("beta:       betaPlus",
                         "likelihood: weighted l share",
                         "wPrev:      gamma_2",
                         "resulting posterior"),
              col = c("darkorange", "purple", "forestgreen", "black"),
              lty = c(1, 1, 1, 2), lwd = 2)
}


## ---- page 8 -----------------------------------------------------------------
plot_weight_decomposition <- function(od, hl = HL) {
  .page()
  M <- cbind(od$dGamma, od$dBeta, od$dLik)
  cols <- c("forestgreen", "darkorange", "purple")
  .signed_stack(od$t, M, cols, hl,
                ylab = "log alpha(disturbed) - log alpha(undisturbed)",
                main = "Which factor gives the modal disturbed component its weight")
  lines(od$t, od$dTotal, lwd = 2.5)
  points(od$t, od$dTotal, pch = 16, cex = 0.7)
  .legend_out(legend = c("wPrev term  d log gamma",
                         "beta term   d log beta",
                         "lik. term   d log l",
                         "total gap"),
              fill = c(cols, NA), border = NA,
              col = c(NA, NA, NA, "black"), lty = c(NA, NA, NA, 1),
              lwd = c(NA, NA, NA, 2))
}


## ---- page 9 -----------------------------------------------------------------
plot_relative_influence <- function(od, hl = HL) {
  .page()
  M   <- abs(cbind(od$dGamma, od$dBeta, od$dLik))
  M[!is.finite(M)] <- NA
  tot <- rowSums(M)
  P   <- M / ifelse(is.finite(tot) & tot > 0, tot, NA)
  cols <- c("forestgreen", "darkorange", "purple")
  barplot(t(P), col = cols, border = NA, space = 0, names.arg = od$t,
          xlab = "t", ylab = "share of |log-weight gap|",
          main = "Relative influence of wPrev, beta and likelihood")
  .legend_out(legend = c("wPrev", "beta", "likelihood"), fill = cols,
              border = NA)
}


## ---- page 10 ----------------------------------------------------------------
## Every (j,k) innovation series, no selection. Blue = undisturbed class k = 0,
## red = disturbed classes k >= 1; solid/circles = j = 1 (z_{t-1} = 0),
## dashed/triangles = j = 2 (z_{t-1} = 1). Opacity of each segment is set by the
## posterior weight alpha[j,k] of that component, so faint lines are branches
## the filter did not believe and saturated lines are the ones carrying the mass.
plot_innovation_z_weighted <- function(bf, hl = HL, amin = 0.07, pow = 0.5) {
  .page()
  d <- bf[!is.na(bf$z) & is.finite(bf$z), ]
  wmax <- max(d$wAn, na.rm = TRUE)

  plot(range(d$t), range(d$z), type = "n", xlab = "t",
       ylab = "standardized innovation z",
       main = "z per component, shaded by posterior weight alpha[j,k]")
  .shade(hl)
  abline(h = c(-2, 0, 2), lty = c(3, 1, 3), col = "grey50")

  keys <- unique(d[, c("j", "k")])
  keys <- keys[order(keys$k, keys$j), ]

  for (r in seq_len(nrow(keys))) {
    j <- keys$j[r]; k <- keys$k[r]
    dk <- d[d$j == j & d$k == k, ]
    dk <- dk[order(dk$t), ]
    if (!nrow(dk)) next

    base <- if (k == 1) c(0, 0, 1) else c(1, 0, 0)     # blue vs red
    lt   <- if (j == 1) 1 else 2
    ph   <- if (j == 1) 16 else 17

    a <- amin + (1 - amin) * (pmax(dk$wAn, 0) / wmax)^pow
    n <- nrow(dk)

    if (n >= 2) {
      ok <- which(diff(dk$t) == 1)                      # only consecutive t
      if (length(ok)) {
        aseg <- pmax(a[ok], a[ok + 1])
        segments(dk$t[ok], dk$z[ok], dk$t[ok + 1], dk$z[ok + 1],
                 col = rgb(base[1], base[2], base[3], aseg),
                 lty = lt, lwd = 2.2)
      }
    }
    points(dk$t, dk$z, pch = ph, cex = 0.65,
           col = rgb(base[1], base[2], base[3], a))
  }

  ramp <- c(1, 0.5, 0.2, 0.05)
  .legend_out(legend = c("undisturbed  k = 0",
                         "disturbed    k >= 1",
                         "",
                         "j = 1  (z_{t-1} = 0)",
                         "j = 2  (z_{t-1} = 1)",
                         "",
                         "opacity ~ alpha[j,k]:",
                         paste0("  alpha/alpha_max = ", ramp)),
              col = c(rgb(0, 0, 1), rgb(1, 0, 0), NA,
                      "grey30", "grey30", NA, NA,
                      rgb(0, 0, 0, amin + (1 - amin) * ramp^pow)),
              lty = c(1, 1, NA, 1, 2, NA, NA, rep(1, length(ramp))),
              pch = c(NA, NA, NA, 16, 17, NA, NA, rep(NA, length(ramp))),
              lwd = c(2.2, 2.2, NA, 1.5, 1.5, NA, NA, rep(4, length(ramp))))
}


## ---- page 11 ----------------------------------------------------------------
## Split the likelihood term of the log-weight gap into its only two sources.
## With s = sqrt(H Sig_f H' + R), log l = -0.5 log(2 pi) - log s - z^2 / 2, so
## for the modal disturbed (d) and modal undisturbed (u) components
##
##   d log l  =  log(s_u / s_d)        <- spread term
##             + 0.5 * (z_u^2 - z_d^2) <- centering term
##
## The -0.5 log(2 pi) cancels. A positive total can therefore come from the
## disturbed forecast being better centred OR from it being narrower; these are
## different failure modes and this page tells them apart.
loglik_split <- function(hist, bf) {
  NT <- length(hist$gSamp)
  res <- data.frame(t = seq_len(NT), dLik = NA_real_,
                    spread = NA_real_, center = NA_real_, resid = NA_real_,
                    zu = NA_real_, zd = NA_real_, su = NA_real_, sdd = NA_real_,
                    ju = NA_integer_, jd = NA_integer_, kd = NA_integer_)

  for (t in seq_len(NT)) {
    d <- bf[bf$t == t, ]
    if (all(is.na(d$wA))) next
    du <- d[d$k == 1, ]; u <- du[which.max(du$wA), ]
    dd <- d[d$k >  1, ]; w <- dd[which.max(dd$wA), ]

    res$spread[t] <- log(u$sd) - log(w$sd)
    res$center[t] <- 0.5 * (u$z^2 - w$z^2)
    res$dLik[t]   <- log(w$lA) - log(u$lA)
    res$resid[t]  <- res$dLik[t] - res$spread[t] - res$center[t]

    res$zu[t] <- u$z;  res$zd[t]  <- w$z
    res$su[t] <- u$sd; res$sdd[t] <- w$sd
    res$ju[t] <- u$j;  res$jd[t]  <- w$j; res$kd[t] <- w$klab
  }
  res
}


plot_loglik_split <- function(ls, hl = HL, hist = NULL) {
  .page()
  M <- cbind(ls$spread, ls$center)
  cols <- c("steelblue3", "purple")
  .signed_stack(ls$t, M, cols, hl,
                ylab = "contribution to  d log l",
                main = "Likelihood term: spread advantage vs centring advantage")
  lines(ls$t, ls$dLik, lwd = 2.5)
  points(ls$t, ls$dLik, pch = 16, cex = 0.7)

  ## mark timesteps where the modal component identity changes
  chg <- which(c(FALSE, diff(ls$jd) != 0 | diff(ls$kd) != 0 |
                   diff(ls$ju) != 0))
  if (length(chg)) abline(v = ls$t[chg] - 0.5, col = "grey45", lty = 3)

  ## only k = 0 receives Q(tau), so the spread term should track 1 / tau*
  leg <- c("spread   log(s_u / s_d)",
           "centring 0.5 (z_u^2 - z_d^2)",
           "total    d log l",
           "modal identity change")
  lcol <- c(NA, NA, "black", "grey45")
  llty <- c(NA, NA, 1, 3)
  llwd <- c(NA, NA, 2.5, 1)
  lfil <- c(cols, NA, NA)

  if (!is.null(hist)) {
    par(new = TRUE)
    plot(ls$t, 1 / hist$tau, type = "l", col = "darkgreen", lwd = 2, lty = 4,
         axes = FALSE, xlab = "", ylab = "",
         xlim = range(ls$t) + c(-1, 1))
    axis(4, col = "darkgreen", col.axis = "darkgreen")
    mtext("1 / tau*", side = 4, line = 2.5, col = "darkgreen", cex = 1.1)
    leg  <- c(leg, "1 / tau* (right axis)")
    lcol <- c(lcol, "darkgreen"); llty <- c(llty, 4); llwd <- c(llwd, 2)
    lfil <- c(lfil, NA)
  }

  .legend_out(legend = leg, fill = lfil, border = NA,
              col = lcol, lty = llty, lwd = llwd)
}


## ---- page 12 ----------------------------------------------------------------
## wPrev[t, g, 2] summarized across the tau grid: mean over g at each t, with
## +/- 1 sd error bars.
plot_wprev2_across_grid <- function(hist, hl = HL) {
  .page()
  W  <- hist$wPrev[, , 2]
  t  <- seq_len(nrow(W))
  mu <- apply(W, 1, mean, na.rm = TRUE)
  sg <- apply(W, 1, sd,   na.rm = TRUE)

  plot(t, mu, type = "n", ylim = range(c(mu - sg, mu + sg), na.rm = TRUE),
       xlab = "t", ylab = "wPrev[t, g, 2]",
       main = "wPrev[t, g, 2]: mean and +/- 1 sd across the tau grid")
  .shade(hl)
  arrows(t, mu - sg, t, mu + sg, angle = 90, code = 3, length = 0.05, lwd = 2)
  lines(t, mu, lwd = 2)
  points(t, mu, pch = 16, cex = 0.9)

  .legend_out(legend = c("mean over g", "+/- 1 sd over g"),
              col = "black", lty = 1, pch = c(16, NA), lwd = 2)

  invisible(data.frame(t = t, mean = mu, sd = sg))
}


## ---- page 13 ----------------------------------------------------------------
## Which branch owns the modal UNDISTURBED component?
##
##   alpha[j,0] = gamma_j * beta_0 * l[j,0] / Z,
##
## so beta_0 and Z cancel and (2,0) is modal exactly when
##
##   l[2,0] / l[1,0]  >  gamma_1 / gamma_2.
##
## The left side is how much better the disturbed branch's undisturbed-class
## forecast fits the observation; the right side is the weight deficit it has to
## overcome. Both are summarized across the grid as geometric means (these are
## ratios, so averaging in log space is the meaningful operation), weighted by
## hist$tauProbs when available. Where the solid line rises above the dashed
## one, the modal pair straddles branches and the dGamma term on page 8 turns
## positive -- see the note there about that decomposition.
undist_branch_ratio <- function(hist) {
  NT <- dim(hist$lA)[1]; G <- dim(hist$lA)[2]

  lr <- log(hist$lA[, , 2, 1])   - log(hist$lA[, , 1, 1])   # log l[2,0]/l[1,0]
  gr <- log(hist$wPrev[, , 1])   - log(hist$wPrev[, , 2])   # log gamma_1/gamma_2
  lr[!is.finite(lr)] <- NA; gr[!is.finite(gr)] <- NA

  ws <- if (is.null(hist$tauProbs)) matrix(1 / G, NT, G) else hist$tauProbs
  ws <- ws / rowSums(ws, na.rm = TRUE)

  wmean <- function(X) {
    w <- ws; w[is.na(X)] <- NA
    rowSums(w * X, na.rm = TRUE) / rowSums(w, na.rm = TRUE)
  }
  wsd <- function(X, mu) {
    w <- ws; w[is.na(X)] <- NA
    sqrt(rowSums(w * (X - mu)^2, na.rm = TRUE) / rowSums(w, na.rm = TRUE))
  }

  lrm <- wmean(lr); grm <- wmean(gr)
  data.frame(t = seq_len(NT),
             logLR = lrm, logLRsd = wsd(lr, lrm),
             logThr = grm, logThrSd = wsd(gr, grm),
             logLRat = lr[cbind(seq_len(NT), hist$gSamp)],
             logThrAt = gr[cbind(seq_len(NT), hist$gSamp)])
}


plot_undist_branch_ratio <- function(hist, hl = HL, logy = TRUE) {
  .page()
  d <- undist_branch_ratio(hist)
  t <- d$t

  yl <- range(exp(c(d$logLR - d$logLRsd, d$logLR + d$logLRsd,
                    d$logThr, d$logLRat)), finite = TRUE)
  plot(t, exp(d$logLR), type = "n", ylim = yl, log = if (logy) "y" else "",
       xlab = "t", ylab = "ratio",
       main = "l[2,0] / l[1,0] against the threshold gamma_1 / gamma_2")
  .shade(hl)
  abline(h = 1, col = "grey60", lty = 3)

  ok <- is.finite(d$logLR) & is.finite(d$logLRsd)
  polygon(c(t[ok], rev(t[ok])),
          c(exp(d$logLR + d$logLRsd)[ok], rev(exp(d$logLR - d$logLRsd)[ok])),
          col = rgb(0, 0, 0, 0.12), border = NA)
  lines(t, exp(d$logLR),  lwd = 2.5, type = "o", pch = 16, cex = 0.8)
  lines(t, exp(d$logLRat), col = "steelblue3", lwd = 1.5, lty = 2)
  lines(t, exp(d$logThr), col = "firebrick", lwd = 2.5, lty = 2)

  win <- which(d$logLR > d$logThr)
  if (length(win)) points(win, exp(d$logLR[win]), pch = 1, cex = 2.2,
                          col = "firebrick", lwd = 2)

  .legend_out(legend = c("l[2,0] / l[1,0]  (geom. mean over g)",
                         "+/- 1 sd across g",
                         "same, at gSamp(t)",
                         "threshold  gamma_1 / gamma_2",
                         "(2,0) is modal undisturbed"),
              col = c("black", rgb(0, 0, 0, 0.3), "steelblue3", "firebrick",
                      "firebrick"),
              lty = c(1, NA, 2, 2, NA), pch = c(16, 15, NA, NA, 1),
              lwd = c(2.5, NA, 1.5, 2.5, 2))

  invisible(d)
}


## ---- standalone: tau posterior at one timestep ------------------------------
## hist$tauProbs[t, ] is the normalized grid posterior ws from step t, and
## hist$tauGrid holds the grid it lives on. Falls back to an explicit tauGrid or
## to rebuilding one from the `grid` list for runs made before tauGrid was
## stored. Draws one page per requested timestep.
plot_tau_posterior <- function(hist, tstep = 10, grid = NULL, tauGrid = NULL,
                               logx = TRUE) {
  G <- ncol(hist$tauProbs)
  if (is.null(tauGrid)) tauGrid <- hist[["tauGrid"]]
  if (is.null(grid))    grid    <- hist[["grid"]]
  if (is.null(tauGrid) && is.list(grid) &&
      !is.null(grid$tauMin) && !is.null(grid$tauMax)) {
    tauGrid <- seq(grid$tauMin, grid$tauMax, length.out = G)
  }
  if (is.null(tauGrid)) {
    stop("no tau grid found: store hist$tauGrid, or pass grid = <list> / tauGrid = <numeric>")
  }
  stopifnot(length(tauGrid) == G)

  for (tt in tstep) {
    .page()
    ws <- hist$tauProbs[tt, ]
    ws <- ws / sum(ws, na.rm = TRUE)

    edge <- max(1L, round(0.05 * G))
    mLo  <- sum(ws[seq_len(edge)], na.rm = TRUE)
    mHi  <- sum(ws[seq.int(G - edge + 1L, G)], na.rm = TRUE)
    pMean <- sum(ws * tauGrid, na.rm = TRUE)
    pMode <- tauGrid[which.max(ws)]

    plot(tauGrid, ws, type = "h", lwd = 2, col = "grey55",
         log = if (logx) "x" else "",
         xlab = "tau", ylab = sprintf("p(tau_g | D_%d)", tt),
         main = sprintf("Posterior over the tau grid at t = %d", tt))
    lines(tauGrid, ws, col = "black", lwd = 2)
    points(tauGrid, ws, pch = 16, cex = 0.5)

    abline(v = hist$tau[tt], col = "red",       lty = 2, lwd = 2)
    abline(v = pMean,        col = "blue",      lty = 3, lwd = 2)
    abline(v = pMode,        col = "darkgreen", lty = 4, lwd = 2)

    ## second axis in process-variance units
    at <- axTicks(1)
    axis(3, at = at, labels = signif(1 / at, 2))
    mtext("1 / tau  (process variance)", side = 3, line = 2.3, cex = 1.05)

    .legend_out(legend = c(
      "grid posterior",
      sprintf("sampled tau* = %.3g", hist$tau[tt]),
      sprintf("posterior mean = %.3g", pMean),
      sprintf("posterior mode = %.3g", pMode),
      "",
      sprintf("mass in lowest %d pts: %.3f", edge, mLo),
      sprintf("mass in highest %d pts: %.3f", edge, mHi)),
      col = c("black", "red", "blue", "darkgreen", NA, NA, NA),
      lty = c(1, 2, 3, 4, NA, NA, NA),
      lwd = c(2, 2, 2, 2, NA, NA, NA))
  }

  invisible(data.frame(tau = tauGrid, p = hist$tauProbs[tstep[1], ]))
}


## ---- driver -----------------------------------------------------------------
diagnose_pdist <- function(hist, Rparams, hl = HL, file = NULL,
                           tauSteps = c(hl[1], length(hist$obs)),
                           grid = NULL, tauGrid = NULL,
                           width = 14, height = 8) {
  bf <- branch_frame(hist, Rparams)
  od <- odds_frame(hist, bf)

  if (!is.null(file)) {
    pdf(file, width = width, height = height, onefile = TRUE)
    on.exit(dev.off(), add = TRUE)
  } else {
    dev.new(width = width, height = height)
  }
  op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)

  plot_obs_analysis(hist, od, hl)
  plot_pdist_tau_averaged(hist, od, hl)
  plot_odds_decomposition(od, hl)
  plot_innovation_z(bf, hl)
  plot_innovation_z_weighted(bf, hl)
  plot_predictive_sd(bf, hl)
  plot_class_composition(bf, hl)
  plot_factor_shares(od, hl)
  plot_weight_decomposition(od, hl)
  plot_relative_influence(od, hl)

  ls <- loglik_split(hist, bf)
  plot_loglik_split(ls, hl, hist)
  plot_wprev2_across_grid(hist, hl)
  plot_undist_branch_ratio(hist, hl)

  tauSteps <- sort(unique(tauSteps[tauSteps >= 1 &
                                     tauSteps <= length(hist$obs)]))
  if (length(tauSteps)) {
    tryCatch(plot_tau_posterior(hist, tauSteps, grid = grid, tauGrid = tauGrid),
             error = function(e)
               warning("tau posterior page skipped: ", conditionMessage(e),
                       call. = FALSE))
  }

  cols <- c("t", "obs", "pDist", "pDistPrior", "logOddsPrior", "logBF",
            "dGamma", "dBeta", "dLik", "dTotal", "z0", "zWin", "kWin", "jWin")
  sub <- merge(od[, cols], ls[, c("t", "spread", "center", "su", "sdd",
                                  "ju", "jd", "kd")], by = "t")
  sub <- sub[sub$t %in% hl, ]
  num <- sapply(sub, is.numeric); sub[num] <- lapply(sub[num], round, 3)
  print(sub, row.names = FALSE)

  invisible(list(branch = bf, odds = od, lik = ls))
}


## Usage:
##   res <- diagnose_pdist(hist, Rparams = c(slope, intercept),
##                         file = "pdist_diagnostics.pdf")
