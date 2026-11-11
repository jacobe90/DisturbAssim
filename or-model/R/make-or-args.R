## Per-site inputs to ORJointUpdate.
##
## Ported from OR_joint_update2.Rmd (now legacy/models/). Uses the globals created by
## or-model/load-data-or.R: n.dist, NT, d.prior, dist.priors, dist.forecastparams,
## dist.prior.probs, dextract, Yw, radDA, dtime, ctime, atime.

## ---- Per-site argument builder ---------------------------------------------
## make_or_args(site, useDistProduct) assembles the site-specific inputs to
## ORJointUpdate: the (n.dist+1) x NT prior probability matrix p.dist (row 1 =
## undisturbed, rows 2.. = disturbance sub-classes), the bundled disturb.params,
## the NT x 3 observation matrix Y (wood-pool obs in column 3), and the PAR
## driver for that site. Relies on the site-independent globals above plus
## d.prior, dist.priors, dist.forecastparams, dextract, Yw and rad.
## tauPrior selects the prior on the process precision tau:
##   "gamma"      : Gamma(shape = 0.1, rate = 2e-4), which spreads mass over every
##                  decade of the grid (mean tau = 500)
##   "loguniform" : p(tau) proportional to 1/tau (flat in log tau)
make_or_args <- function(site, useDistProduct = TRUE, tauPrior = c("gamma", "loguniform")) {
  tauPrior <- match.arg(tauPrior)
  ## Initialize every year with the background prior (same column for all t).
  p.dist <- array(NA, dim = c(n.dist + 1, NT))
  for (t in 1:NT) {
    p.dist[1:(n.dist+1), t] <- as.numeric(c(d.prior[1], unlist(dist.prior.probs)))
  }

  ## If using the disturbance product, override the prior in years where a
  ## disturbance was detected at this site: replace that class's prior with the
  ## user's-accuracy estimate (Kennedy et al. 2015 table 6). Then give each
  ## component a share of the disturbance probability proportional to its weight,
  ## and the disturbed values to sum to 1 if they sum to greater than 1, and else
  ## set the undisturbed value to 1 - p(disturbed).
  ## within-class mixture weights (mclust proportions; sum to 1 within each class)
  weights <- unlist(sapply(dist.priors, function(x) x$weights))

  if (useDistProduct) {
    distidx.toname = list("10"="OTHER", "20"="CUT", "40"="FIRE", "50"="PEST")
    users_accuracy = c(1-d.prior[1], .77, .92, .77, .78, .78) ## 0 = none, 10 = other, 20 = cut, 30 = road, 40 = fire, 50 = biotic

    ## product years inside the DA window and their DA timesteps
    dyrs <- intersect(dtime, atime)
    dprod <- dextract[site, paste0("t", dyrs)]
    for (d in seq(names(dist.priors))) {
      d.name <- names(dist.priors)[d]
      has.dist <- which(dprod == d.name)                 # index into dyrs
      if (length(has.dist) > 0) {
        is.dist <- which(names(dist.prior.probs) == distidx.toname[[d.name]]) + 1
        p.dist[is.dist, match(dyrs[has.dist], atime)] <- users_accuracy[dprod[has.dist]/10 + 1]
      }
    }
  }

  ## Build the per-year categorical prior for BOTH the product and no-product paths:
  ## split each class's probability across its clusters via `weights` (so multi-cluster
  ## classes are not over-counted), then set the undisturbed prior to the remaining mass
  ## (= 1 - product accuracy in flagged years, = background otherwise). Applying this on
  ## the no-product path too means an UNflagged year has the same prior with or without
  ## the product -- the product only changes flagged years. Previously this lived inside
  ## the useDistProduct branch, so the no-product path kept an un-split, over-counted
  ## disturbance background that inflated probD relative to the product on unflagged sites.
  for (j in 1:ncol(p.dist)) {
    p.dist[2:(n.dist+1), j] <- p.dist[2:(n.dist+1), j] * weights   # split class prob across clusters
    dist.mass <- sum(p.dist[2:(n.dist+1), j])
    if (dist.mass > 1) {                                           # guard: cap total disturbed mass at 1
      p.dist[2:(n.dist+1), j] <- p.dist[2:(n.dist+1), j] / dist.mass
      dist.mass <- 1
    }
    p.dist[1, j] <- 1 - dist.mass                                  # undisturbed = remaining mass
  }

  ## Bundle the disturbance prior + per-class forecast params.
  disturb.params <- list(
    p = p.dist,
    n = n.dist,
    vals = dist.forecastparams,
    alloc.soil = 0.05
  )

  ## NT x 3 observation matrix with the wood-pool obs in column 3.
  ## Rows stay aligned with DA timesteps; years without an observation stay NA
  ## (forecast-only step in ORJointUpdate) rather than being dropped, which would
  ## shift Y out of step with p.dist and PAR.
  Y <- matrix(NA, NT, 3)
  oyrs <- intersect(ctime, atime)
  Y[match(oyrs, atime), 3] <- Yw[site, paste0("t", oyrs)]

  ## tau prior and grid
  prior <- switch(tauPrior,
    gamma      = list(type = "gamma", shape = 0.1, rate = 2e-4),
    loguniform = list(type = "loguniform"))
  grid <- list(tauMin = 0.1, tauMax = 10^4)
  list(Y = Y, disturb.params = disturb.params, PAR = radDA[site, ], site = site, prior = prior, grid = grid)
}
