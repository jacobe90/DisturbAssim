## Biomass observations and the initial-condition prior.
##
## Ported from OR_joint_update2.Rmd (now legacy/models/).

## Wood carbon observations (kg C m-2), sites x years 1990:2017, from the LandTrendr
## biomass extract (Mg/ha). Also the observation error from Kennedy et al. 2018, Environ.
## Res. Lett. 13 025004, Fig. 2 RMSE (not used by ORJointUpdate, which takes Rparams).
build_obs <- function(bextract) {
  Yw <- bextract * B2C * 1000 * 1e-4   # Mg/ha -> kg C m-2
  colnames(Yw) <- paste0("t", ctime)
  list(Yw = Yw, Rw = (84 * B2C * 1000 * 1e-4)^2)
}

## Initial-condition prior: one ensemble member per site.
##   1. spin-up: VSEM run over the concatenated PAR of every site; ne states sampled
##      from its second half
##   2. recovery: from each spin-up state with leaf and wood cut to 1%, VSEM run for nNBG
##      years (pools at the end of every year)
##   3. each member's age is the year of its recovery trajectory whose wood is closest to
##      that site's 1990 biomass observation; its leaf and soil are taken from there and its
##      wood is the observation
build_ic_prior <- function(theta, rad, AGBprior, nNBG) {
  ne <- length(AGBprior)
  PARlong <- as.vector(t(rad))           # PAR series of every site, one after the other
  message(sprintf("  spin-up: VSEM over %d days", length(PARlong)))
  spin   <- VSEM(theta, PARlong)
  thin   <- seq(nrow(spin) / 2, nrow(spin), length = ne)
  spinIC <- spin[thin, 2:4]
  message(sprintf("  recovery trajectories: %d members x %d years", ne, nNBG))
  NBG <- array(NA, c(ne, nNBG, 3))
  t0 <- Sys.time()
  for (i in 1:ne) {
    param <- theta
    param[9:11] <- spinIC[i, ] * c(0.01, 1, 0.01)   # leaf and wood to 1% (a disturbance)
    raw <- VSEM(param, PARlong[1:(nNBG * 365.25)])
    NBG[i, , ] <- raw[round(seq(365, nNBG * 365.25, length = nNBG)), 2:4]   # end of each year
    if (i %% 100 == 0 || i == ne) {
      el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      message(sprintf("  recovery trajectories: %d / %d (elapsed %s, ETA %s)", i, ne,
                      fmt_secs(el), fmt_secs(el / i * (ne - i))))
    }
  }
  ICprior <- matrix(NA, ne, 3)
  age <- rep(NA, ne)
  for (i in 1:ne) {
    age[i] <- which.min((AGBprior[i] - NBG[i, , 3])^2)
    ICprior[i, ] <- c(NBG[i, age[i], 1:2], AGBprior[i])
  }
  colnames(ICprior) <- c("Cv", "Cs", "CR")
  list(ICprior = ICprior, age = age)
}
