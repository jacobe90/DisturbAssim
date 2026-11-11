## Dropped-component bookkeeping and its figure, shared by every ORJointUpdate-style filter
## (singlesite-singledisturb-multipool.R, or-joint-update-filter.R, and the factorized-tau
## filter).
##
## Before the Kalman update a forecast component (g, h, k) is dropped (prior weight set to 0)
## if
##   - its importance weights collapsed onto ~1 sample (cN < 1e-8; dropF, per (g, h)), or
##   - its weighted forecast covariance is singular (badC, per (g, h, k)).
## A grid point left with no component at all would get a zero marginal likelihood, so its
## components are restored at their prior weights (on placeholder moments) instead.

##' Apply the drops to the prior component weights and measure them
##'
##' @param wPrior0 G x H x K1 prior weights of every component (before dropping)
##' @param dropF G x H logical, collapsed importance weights
##' @param badC logical over the G H K1 components (g fastest), singular covariance; or NULL
##' @param piPrev weight of each grid point (sums to 1)
##' @return list(w = weights after dropping and restoring, massIS = forecast probability mass
##'   dropped for collapsed importance weights, massCov = further mass dropped for a singular
##'   covariance, none = logical per grid point, TRUE if it was left with no mass and restored,
##'   massNet = mass dropped after restoring)
apply_drops <- function(wPrior0, dropF, badC, piPrev) {
  G <- dim(wPrior0)[1]; K1 <- dim(wPrior0)[3]
  wIS <- wPrior0
  for (k in 1:K1) wIS[, , k][dropF] <- 0
  w <- wIS
  if (!is.null(badC)) w[badC] <- 0
  none <- rowSums(matrix(w, G)) == 0
  tot  <- sum(wPrior0 * piPrev)
  massIS  <- sum((wPrior0 - wIS) * piPrev) / tot
  massCov <- sum((wIS - w) * piPrev) / tot
  if (any(none)) w[none, , ] <- wPrior0[none, , ]
  list(w = w, none = none, massIS = massIS, massCov = massCov,
       massNet = 1 - sum(w * piPrev) / tot)
}

##' Dropped forecast probability mass and empty grid points over time
##'
##' Top: the share of the forecast probability (over the tau posterior, history groups and
##' scenarios) put on dropped components at each step, split by cause, on a log scale; zero
##' is drawn at the bottom of the axis. Bottom: the number of grid points that had no mass
##' left after dropping (and were restored), against the number of active grid points.
##' @param fit filter output with wDropIS, wDropCov, wDrop, nNone, nActive (and nBad)
plot_dropped_mass <- function(fit, main = "Dropped forecast components", yr = seq_along(fit$wDrop)) {
  lo   <- 1e-12
  cl   <- function(x) pmax(x, lo)
  cols <- c(IS = "#d95f02", cov = "#1b78b4", net = "grey20")
  op <- par(mfrow = c(2, 1), mar = c(3.8, 5.6, 3.4, 1), oma = c(0, 0, 2, 0), mgp = c(2.4, 0.7, 0))
  on.exit(par(op))

  tot <- fit$wDropIS + fit$wDropCov
  plot(yr, cl(tot), type = 'n', log = "y", ylim = c(lo, 1), yaxt = "n", xlab = "year",
       ylab = "", main = "Forecast probability mass dropped")
  mtext(sprintf("max %.2g in one step; %d components dropped for a singular covariance in all",
                max(tot), sum(fit$nBad)), side = 3, line = 0.3, cex = 0.75)
  mtext("probability mass", side = 2, line = 4.2)
  ticks <- 10^seq(-12, 0, by = 3)
  axis(2, at = ticks, labels = c("<= 1e-12", sprintf("1e%d", seq(-9, -3, by = 3)), "1"), las = 1, cex.axis = 0.75)
  abline(h = ticks, col = "grey90")
  lines(yr, cl(fit$wDropIS),  type = 'o', pch = 19, cex = 0.6, col = cols["IS"],  lwd = 1.8)
  lines(yr, cl(fit$wDropCov), type = 'o', pch = 17, cex = 0.6, col = cols["cov"], lwd = 1.8)
  lines(yr, cl(fit$wDrop),    type = 'o', pch = 1,  cex = 0.8, col = cols["net"], lwd = 1, lty = 2)
  legend("topleft", bty = "n", cex = 0.7,
         legend = c("collapsed importance weights", "singular covariance", "net, after restoring empty grid points"),
         col = cols, pch = c(19, 17, 1), lty = c(1, 1, 2), lwd = c(1.8, 1.8, 1))

  nn <- fit$nNone
  top <- max(1, nn)
  b <- barplot(nn, names.arg = yr, col = "#b2182b", border = NA, ylim = c(0, top * 1.6), yaxt = "n",
               xlab = "year", ylab = "grid points", main = "Grid points with zero mass after dropping")
  axis(2, at = seq(0, top, by = max(1, ceiling(top / 4))), las = 1, cex.axis = 0.75)
  mtext(sprintf("%d in all, in %d of %d steps (labels: zero-mass / active grid points)",
                sum(nn), sum(nn > 0), length(nn)), side = 3, line = 0.3, cex = 0.75)
  if (any(nn > 0))
    text(b[nn > 0], nn[nn > 0] + top * 0.05, labels = sprintf("%d/%d", nn, fit$nActive)[nn > 0],
         srt = 90, adj = 0, cex = 0.65, xpd = TRUE)
  mtext(main, outer = TRUE, font = 2)
}

## one-line summary of the drops of a fit, for the text summaries
dropped_summary <- function(fit, name = "") {
  sprintf("%sdropped mass: max %.2g (importance %.2g, covariance %.2g); components with singular covariance %d; zero-mass grid points %d in %d steps",
          name, max(fit$wDropIS + fit$wDropCov), max(fit$wDropIS), max(fit$wDropCov),
          sum(fit$nBad), sum(fit$nNone), sum(fit$nNone > 0))
}
