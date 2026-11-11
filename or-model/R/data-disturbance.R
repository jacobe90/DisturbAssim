## Disturbance inputs: background disturbance-type prior and the per-class forecast
## parameters of the percent-decrease disturbance model (R/disturbance.R, disturbance.p).
##
## Ported from OR_joint_update2.Rmd (now legacy/models/).

## Background prior from STprior.RData: fold the <NA> class into "UN" (undisturbed)
clean_disturbance_prior <- function(d.prior) {
  d.prior["UN"] <- d.prior["UN"] + d.prior[6]   # removes <NA>
  d.prior[-6]
}

## One forecast component per disturbance class and size cluster (mclust fits in
## dist.priors). mu0 = (leaf, stem) mean fraction retained after the disturbance (leaf
## fixed at 0.2, stem from the cluster mean); V0 = their variances.
build_dist_forecastparams <- function(dist.priors) {
  dist_name <- c("0" = "none", "10" = "other", "20" = "cut", "30" = "road", "40" = "fire", "50" = "biotic")
  out <- list()
  for (dclass in names(dist.priors)) {
    dname <- as.character(dist_name[dclass])
    for (i in seq_along(dist.priors[[dclass]]$means)) {
      out[[length(out) + 1]] <- list(mu0  = c(0.2, dist.priors[[dclass]]$means[i]),
                                     V0   = diag(c(0.04, dist.priors[[dclass]]$vars[i])),
                                     name = paste0(dname, i))
    }
  }
  out
}

## Background disturbance-type probability per forecast component, by matching the
## component's name to one of "other"/"cut"/"fire"/"biotic".
build_dist_prior_probs <- function(dist.forecastparams, d.prior) {
  distname.toidx <- list("other" = 2, "cut" = 3, "fire" = 4, "biotic" = 5)
  sapply(dist.forecastparams, function(x) {
    for (d in names(distname.toidx)) if (grepl(d, x$name)) return(d.prior[distname.toidx[[d]]])
  })
}
