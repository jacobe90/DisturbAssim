## Single site, single disturbance, single pool (logistic LAI) synthetic experiment.
##
## Migrated from SingleSite_SingleDisturb.Rmd. Reproduces
##   1. the "true" run and observations,
##   2. the blind (free-run) forecast without and with disturbance,
##   3. the standard EnKF without and with disturbance,
##   4. the multinomial (disturbed / undisturbed mixture) filter -- here the joint
##      disturbance + process-precision filter of or-model/R/or-joint-update.R (ORJointUpdate),
##      ported from the 3-pool VSEM model to the 1-pool logistic LAI model,
## and the comparison figure of the three filters.
##
## The logistic model has a single scalar process error, n_t ~ N(f(n_{t-1}), 1/tau),
## with f(n) = n + r n (1 - n/K). As in ORJointUpdate, tau is unknown and carried on
## a log-spaced grid; disturbance replaces the state with a draw from N(mu0, V0).
##
## Run from the repository root:  Rscript synthetic-experiments/singlesite-singledisturbance-singlepool.R
## Figures and results go to `outdir` (git-ignored).

outdir <- "experiments/synthetic/singlesite-singledisturbance-singlepool"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

## ---- "True" run ---------------------------------------------------------------
## Same parameters and RNG sequence as SingleSite_SingleDisturb.Rmd, so n.true and
## n.obs are identical to the Rmd's (and to data$Y in SSSD.RData).
set.seed(1066)
p.dist = 0.005   # background disturbance probability per step
n0     = 1       # mean post-disturbance LAI
n0.sd  = 0.5     # sd of post-disturbance LAI
nIC    = 6
r      = .1
K      = 6
n.sd   = 0.25    # process error sd (tau = 1/n.sd^2 = 16)
NT     = 100
time   = 1:NT
xlim   = c(5,25)
n = rep(nIC,NT)
d = rep(0,NT)
d[10] = 1  ## prescribe a disturbance
for(t in 2:NT){
  n[t] = n[t-1] + r*n[t-1]*(1-n[t-1]/K) + rnorm(1,0,n.sd)
  if(d[t]) n[t] = abs(rnorm(1,n0,n0.sd))
}
obs.sd = 0.75
n.true = n
n.obs  = rnorm(NT,n.true,obs.sd)

logistic <- function(x) x + r*x*(1-x/K)

## ---- Shared helpers -----------------------------------------------------------
NE   = 10000                       # ensemble size (blind runs and EnKF)
qs5  = c(0.025,0.25,0.5,0.75,0.975)

## blind forecast: NE free runs from the lognormal IC, optionally with disturbance
blind_forecast <- function(disturb) {
  ne = matrix(rlnorm(NE,log(nIC),n.sd),NE,NT)
  for(t in 2:NT){
    ne[,t] = abs(logistic(ne[,t-1]) + rnorm(NE,0,n.sd))
    if(disturb){
      dd = as.logical(rbinom(NE,1,p.dist))
      ne[dd,t] = abs(rnorm(sum(dd),n0,n0.sd))
    }
  }
  apply(ne,2,quantile,c(0.025,0.5,0.975))
}

## standard EnKF (Gaussian-resampled analysis ensemble), optionally with disturbance
## in the forecast ensemble. Returns forecast and analysis quantiles (5 x NT).
## NOTE: the Rmd used obs.sd (not obs.sd^2) as the observation variance in the Kalman
## gain; that is fixed here, so the EnKF bands are slightly narrower than in the Rmd.
run_enkf <- function(disturb) {
  ens = matrix(NA,NT*2+1,NE)
  ens[1,] = rlnorm(NE,log(nIC),n.sd)
  Xae = Pae = Xfe = Pfe = rep(NA,NT)
  for(i in 1:NT){
    ## approximate forecast with Normal
    k = 2*(i-1)+1
    Xfe[i] = mean(ens[k,])
    Pfe[i] = var(ens[k,])

    ## analysis step: combine previous forecast with observed data
    Kg = Pfe[i]/(Pfe[i] + obs.sd^2)
    Xae[i] = if(!is.na(n.obs[i])) Xfe[i] + Kg*(n.obs[i] - Xfe[i]) else Xfe[i]
    Pae[i] = (1-Kg)*Pfe[i]

    ## generate new ensemble
    ens[k+1,] = rnorm(NE,Xae[i],sqrt(Pae[i]))

    ## forecast step: predict to next step from current
    ens[k+2,] = logistic(ens[k+1,]) + rnorm(NE,0,n.sd)
    if(disturb){
      disturbed = as.logical(rbinom(NE,1,p.dist))
      ens[k+2,disturbed] = abs(rnorm(sum(disturbed),n0,n0.sd)) ## disturbance
    }
  }
  list(fx = apply(ens[2*(1:NT)-1,],1,quantile,qs5),
       an = apply(ens[2*(1:NT),],1,quantile,qs5))
}

## ---- Blind forecasts and EnKF ---------------------------------------------------
n0.stats = blind_forecast(disturb = FALSE)
enkf0    = run_enkf(disturb = FALSE)
fx0.stats = enkf0$fx; an0.stats = enkf0$an

n.stats  = blind_forecast(disturb = TRUE)
enkf     = run_enkf(disturb = TRUE)
fx.stats = enkf$fx; an.stats = enkf$an

## ---- Multinomial filter: ORJointUpdate ported to the logistic LAI model ----------

# Quantiles of a 1-d Gaussian mixture sum_i w_i N(mu_i, sd_i^2), by root-finding on its
# CDF. Components with negligible weight are dropped first.
mix_quantile <- function(p, w, mu, sd) {
  keep <- w > 1e-12 * sum(w)
  w <- w[keep] / sum(w[keep]); mu <- mu[keep]; sd <- sd[keep]
  lo <- min(mu - 8 * sd); hi <- max(mu + 8 * sd)
  sapply(p, function(pp) uniroot(function(x) sum(w * pnorm(x, mu, sd)) - pp,
                                 c(lo, hi), tol = 1e-8)$root)
}

# Initial tau weights: grid point g gets the prior mass of its cell, with the cell
# edges at the midpoints to the neighbouring grid points (see or-model/R/or-joint-update.R).
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

# Moment-matching projection of the analysis mixture onto two Gaussians per tau grid
# point: group 1 = undisturbed components (k = 0), group 2 = disturbed (k >= 1).
# w, mu, sig: G x 2 x (D+1). Returns w, mu, sig as G x 2 matrices [g, group].
project_grid <- function(w, mu, sig) {
  G   <- dim(w)[1]
  grp <- list(1, seq_len(dim(w)[3])[-1])
  out <- list(w = matrix(0, G, 2), mu = matrix(0, G, 2), sig = matrix(0, G, 2))
  for (h in 1:2) {
    wh  <- w[, , grp[[h]], drop = FALSE]
    wg  <- apply(wh, 1, sum)
    out$w[, h]   <- wg
    out$mu[, h]  <- apply(wh * mu[, , grp[[h]], drop = FALSE], 1, sum) / wg
    out$sig[, h] <- apply(wh * (sig[, , grp[[h]], drop = FALSE] + mu[, , grp[[h]], drop = FALSE]^2), 1, sum) / wg -
                    out$mu[, h]^2
  }
  # a group with zero total weight has undefined moments (0/0 above); give it the other
  # group's moments so it stays finite -- its weight is 0, so it only enters the proposal
  z1 <- !(out$w[, 1] > 0); z2 <- !(out$w[, 2] > 0) & !z1
  out$mu[z1, 1] <- out$mu[z1, 2]; out$sig[z1, 1] <- out$sig[z1, 2]
  out$mu[z2, 2] <- out$mu[z2, 1]; out$sig[z2, 2] <- out$sig[z2, 1]
  out
}

##' Approximate joint disturbance + process-precision update, logistic LAI model
##'
##' Port of ORJointUpdate (or-model/R/or-joint-update.R) to a 1-pool state. Per tau grid point g
##' the analysis is a 2-component Gaussian mixture (j = undisturbed / disturbed history).
##' Forecast: importance-sample the previous analysis mixture, push the samples through
##' the logistic model and take weighted moments for every (g, j); the undisturbed branch
##' (k = 0) adds process variance 1/tau_g, a disturbed branch (k >= 1) replaces the state
##' with N(mu0_k, V0_k), which does not depend on the prior state so its moments are exact.
##' Analysis: Kalman update of every (g, j, k), likelihood-weighted, then moment-matched
##' back onto 2 components per g; the tau posterior is updated on the grid.
##'
##' @param Y length-NT observations (NA => forecast-only step)
##' @param R observation-error variance
##' @param Xic initial-condition ensemble (vector)
##' @param model deterministic one-step model f(x), vectorised
##' @param dist disturbance list: p ((D+1) x NT prior; row 1 = undisturbed), n (# classes),
##'   vals (per-class list(mu0, V0))
##' @param G number of tau grid points
##' @param N importance-sample size
##' @param prior prior for tau: list(type = "gamma", shape, rate) or list(type = "loguniform")
##' @param grid list(tauMin, tauMax) for the log-spaced tau grid (tauMin == tauMax => known tau)
##' @param forecastFirst if FALSE, Xic is the prior at t = 1 and step 1 is analysis-only
##'   (the alignment of the EnKF in SingleSite_SingleDisturb.Rmd); if TRUE, Xic is the
##'   state at t = 0 and step 1 forecasts first, as in ORJointUpdate
##' @return list of length-NT summaries, marginalised over the tau posterior and components:
##'   mean, var, q (5 x NT analysis quantiles at qs5), fq (forecast quantiles), pDist,
##'   pit, wDrop, procVar = E[1/tau | y_{1:t}], tauPost (G x NT), tauGrid, and mix
##'   (per-t analysis mixture: w, mu, sd) for density plots
ORJointUpdateLogistic <- function(Y, R, Xic, model, dist, G = 200, N = 1000, prior, grid,
                                  forecastFirst = FALSE) {
  NT <- length(Y)
  D  <- dist$n
  tauGrid <- exp(seq(log(grid$tauMin), log(grid$tauMax), length.out = G))

  qs      <- log(cell_prior_mass(tauGrid, prior))
  prevW   <- matrix(0.5, G, 2)
  prevMu  <- matrix(mean(Xic), G, 2)
  prevSig <- matrix(var(Xic), G, 2)

  hist <- list(
    mean = rep(NA_real_, NT), var = rep(NA_real_, NT),
    q = matrix(NA_real_, 5, NT), fq = matrix(NA_real_, 5, NT),
    pDist = rep(NA_real_, NT), pit = rep(NA_real_, NT), wDrop = rep(NA_real_, NT),
    procVar = rep(NA_real_, NT), tauPost = matrix(NA_real_, G, NT), tauGrid = tauGrid,
    mix = vector("list", NT)
  )

  for (t in 1:NT) {
    muF  <- array(0, dim = c(G, 2, D+1))
    sigF <- array(0, dim = c(G, 2, D+1))
    dropF <- matrix(FALSE, G, 2)

    if (t == 1 && !forecastFirst) {
      # the IC is the time-1 prior: no propagation, no process error
      muF[, , 1]  <- prevMu
      sigF[, , 1] <- prevSig
    } else {
      # sample from the proposal: equal-weight mixture over every (g, j), positive LAI only
      propSamples <- numeric(0); tries <- 0L
      while (length(propSamples) < N) {
        tries <- tries + 1L
        if (tries > 100L) stop(sprintf("t=%d: proposal rejection rate too high", t))
        i <- sample.int(G, N, replace = TRUE); j <- sample.int(2L, N, replace = TRUE)
        x <- rnorm(N, prevMu[cbind(i, j)], sqrt(prevSig[cbind(i, j)]))
        propSamples <- c(propSamples, x[x > 0])
      }
      propSamples <- propSamples[1:N]
      sim0 <- model(propSamples)

      # P x N log density of every sample under every flattened component r = (j-1)*G + g
      compMu  <- as.vector(prevMu)
      compSd  <- sqrt(as.vector(prevSig))
      logDens <- dnorm(outer(-compMu, propSamples, "+") / compSd, log = TRUE) - log(compSd)

      # importance weights of each component's samples, in log space
      mx      <- apply(logDens, 2, max)
      logProp <- mx + log(colSums(exp(sweep(logDens, 2, mx)))) - log(2 * G)
      logW    <- sweep(logDens, 2, logProp)
      W       <- exp(logW - apply(logW, 1, max))
      W       <- W / rowSums(W)
      cN      <- 1 - rowSums(W^2)

      # components whose weights collapse onto ~1 sample have no estimable forecast
      # moments: unweighted placeholder moments and zero likelihood below
      degen <- cN < 1e-8
      W[degen, ] <- 1 / N
      cN[degen]  <- 1 - 1 / N
      dropF <- matrix(degen, G, 2)

      ctr <- mean(sim0)
      MU  <- as.vector(W %*% sim0)
      S   <- (as.vector(W %*% (sim0 - ctr)^2) - (MU - ctr)^2) / cN
      muF[, , 1]  <- MU
      sigF[, , 1] <- S + 1 / rep(tauGrid, 2)        # process variance 1/tau_g
    }
    for (k in seq_len(D)) {                         # disturbance replaces the state
      muF[, , k + 1]  <- dist$vals[[k]]$mu0
      sigF[, , k + 1] <- dist$vals[[k]]$V0
    }

    beta    <- dist$p[, t]
    piPrev  <- exp(qs - max(qs)); piPrev <- piPrev / sum(piPrev)          # pi_{t-1}(g)
    wPrior  <- outer(prevW, beta)                                         # G x 2 x (D+1)
    wPrior[, , 1][dropF] <- 0
    hist$wDrop[t] <- sum(piPrev * prevW * dropF) * beta[1]

    # forecast (one-step-ahead) mixture, marginalised over pi_{t-1}
    wF <- wPrior * piPrev; wF <- wF / sum(wF)
    hist$fq[, t] <- mix_quantile(qs5, wF, muF, sqrt(sigF))

    if (is.na(Y[t])) {
      muA <- muF; sigA <- sigF; lA <- array(1, dim(muF))
    } else {
      s   <- sigF + R
      Kg  <- sigF / s
      muA <- muF + Kg * (Y[t] - muF)
      sigA <- sigF * (1 - Kg)
      lA  <- dnorm(Y[t], muF, sqrt(s))
      hist$pit[t] <- sum(wF * pnorm(Y[t], muF, sqrt(s)))
    }

    # tau posterior on the grid
    zsum <- apply(wPrior * lA, 1, sum)
    logz <- log(zsum)
    if (any(!is.finite(logz))) stop(sprintf("t=%d: zero or non-finite marginal likelihood at %d tau grid points", t, sum(!is.finite(logz))))
    qs <- qs + logz
    ws <- exp(qs - max(qs)); ws <- ws / sum(ws)
    hist$tauPost[, t] <- ws
    hist$procVar[t]   <- sum(ws / tauGrid)

    # analysis weights
    wA <- wPrior * lA / zsum

    # summaries, marginalised over the tau posterior and components
    wAll <- wA * ws
    hist$mean[t]  <- sum(wAll * muA)
    hist$var[t]   <- sum(wAll * (sigA + muA^2)) - hist$mean[t]^2
    hist$pDist[t] <- sum(wAll[, , -1])
    hist$q[, t]   <- mix_quantile(qs5, wAll, muA, sqrt(sigA))
    keep <- wAll > 1e-10
    hist$mix[[t]] <- list(w = wAll[keep], mu = muA[keep], sd = sqrt(sigA[keep]))

    # project back onto 2 components per tau grid point
    out     <- project_grid(wA, muA, sigA)
    prevW   <- out$w
    prevMu  <- out$mu
    prevSig <- out$sig
  }
  hist
}

## disturbance prior: background p.dist every step except t = 1 (the IC step, where
## the truth cannot be disturbed)
or.dist <- list(
  p    = rbind(c(1, rep(1 - p.dist, NT - 1)), c(0, rep(p.dist, NT - 1))),
  n    = 1,
  vals = list(list(mu0 = n0, V0 = n0.sd^2))
)
or.Xic <- rlnorm(NE, log(nIC), n.sd)

## tau unknown: log-uniform prior on a grid spanning process sd 0.01 .. 3.2
or.fit <- ORJointUpdateLogistic(Y = n.obs, R = obs.sd^2, Xic = or.Xic, model = logistic,
                                dist = or.dist, G = 200, N = 1000,
                                prior = list(type = "loguniform"),
                                grid = list(tauMin = 0.1, tauMax = 1e4))
## tau known (= 1/n.sd^2), the setting of the Rmd's multinomial filter
or.known <- ORJointUpdateLogistic(Y = n.obs, R = obs.sd^2, Xic = or.Xic, model = logistic,
                                  dist = or.dist, G = 1, N = 1000,
                                  prior = list(type = "loguniform"),
                                  grid = list(tauMin = 1/n.sd^2, tauMax = 1/n.sd^2))

CI.fN = or.fit$q
Dbar  = or.fit$pDist

## ---- Figures ----------------------------------------------------------------------
band <- function(x, lo, hi, col) polygon(c(x, rev(x)), c(lo, rev(hi)), border = NA, col = col)
col.alpha <- function(col, a) adjustcolor(col, alpha.f = a)
leg <- function() legend("topleft",legend=c("True","Obs","95%","IQR","Model"),lty=c(1,0,1,1,1),
                         col=c(1,2,"lightgreen","green4","lightblue"),lwd=c(1,NA,15,15,15),
                         cex=1.2,pch=c(NA,18,NA,NA,NA))
filter_panel <- function(stats, free, yrng, main, extra = NULL) {
  plot(time,stats[5,],ylim=yrng,xlim=xlim,type='n',ylab="LAI",cex.lab=1.5,main=main)
  band(time,free[1,],free[3,],'lightblue')       ## free run 95%
  band(time,stats[1,],stats[5,],'lightgreen')     ## filter 95%
  band(time,stats[2,],stats[4,],'green4')         ## filter IQR
  if(!is.null(extra)) band(time,extra[2,],extra[4,],'orange')
  lines(time,free[2,],lwd=2,col="blue4")
  lines(time,n.true,lwd=2)
  points(time,n.obs,pch=18,col=2)
  leg()
}

pdf(file.path(outdir, "figures.pdf"), width = 7, height = 6)

## truth and obs
plot(time,n.true,type='l',ylim=range(n.true),xlim=xlim,ylab="LAI",lwd=3,main="Truth")
points(time,n.obs,pch=18,col=2,cex=1.5)
legend("topright",legend=c("truth","obs"),lwd=c(2,NA),pch=c(NA,18),col=1:2)

## blind forecast w/o disturbance
yrng0 = c(0,8)
plot(time,n.obs,type='n',ylim=yrng0,xlim=xlim,ylab="LAI",main="Blind forecast, no disturbance")
band(time,n0.stats[1,],n0.stats[3,],'lightblue')
lines(time,n0.stats[2,],lwd=2,col="blue4")
points(time,n.obs,pch=18,col=2,cex=1.5)

## standard EnKF w/o disturbance
filter_panel(fx0.stats, n0.stats, yrng0, "EnKF, no disturbance: forecast")
filter_panel(an0.stats, n0.stats, yrng0, "EnKF, no disturbance: analysis")

## blind forecast w/ disturbance
yrng = c(0, max(n.stats))
plot(time,n.obs,type='n',ylim=yrng,xlim=xlim,ylab="LAI",main="Blind forecast, with disturbance")
band(time,n.stats[1,],n.stats[3,],'lightblue')
lines(time,n.stats[2,],lwd=2,col="blue4")
points(time,n.obs,pch=18,col=2,cex=1.5)

## standard EnKF + disturbance (orange = analysis IQR without disturbance)
filter_panel(fx.stats, n.stats, yrng, "EnKF + disturbance: forecast")
filter_panel(an.stats, n.stats, yrng, "EnKF + disturbance: analysis", extra = an0.stats)

## multinomial (OR joint) filter
filter_panel(or.fit$fq, n.stats, yrng, "Multinomial filter (tau unknown): forecast")
filter_panel(or.fit$q,  n.stats, yrng, "Multinomial filter (tau unknown): analysis")

## P(disturbed)
plot(time,or.fit$pDist,type='h',lwd=4,col="#eb6834",xlim=xlim,ylim=c(0,1),
     ylab="P(disturbed)",main="Posterior probability of disturbance")
points(time,or.known$pDist,pch=19,cex=0.6)
abline(v=which(d==1),lty=2)
legend("topright",legend=c("tau unknown","tau known","true disturbance"),
       lty=c(1,NA,2),lwd=c(4,NA,1),pch=c(NA,19,NA),col=c("#eb6834",1,1))

## process-error sd: posterior over time vs truth
tauQ <- apply(or.fit$tauPost, 2, function(w) {
  cdf <- cumsum(w); or.fit$tauGrid[sapply(c(0.025,0.5,0.975), function(p) which(cdf >= p)[1])]
})
sdQ <- 1/sqrt(tauQ)
plot(time,sdQ[2,],type='n',log="y",ylim=range(sdQ),ylab="process error sd (1/sqrt(tau))",
     main="Process-error posterior")
band(time,sdQ[3,],sdQ[1,],col.alpha("#2a78d6",0.3))
lines(time,sdQ[2,],lwd=2,col="#2a78d6")
lines(time,sqrt(or.fit$procVar),lty=2,col="#2a78d6")
abline(h=n.sd,lwd=2)
legend("topright",legend=c("posterior median","95% interval","sqrt(E[1/tau])","truth"),
       col=c("#2a78d6",col.alpha("#2a78d6",0.3),"#2a78d6",1),lwd=c(2,10,1,2),lty=c(1,1,2,1))

## analysis densities around the disturbance (mixtures, not Gaussian)
xseq <- seq(-0.5, 8, length.out = 400)
tt   <- xlim[1]:xlim[2]
dens <- sapply(tt, function(t) { m <- or.fit$mix[[t]]
  colSums(m$w * dnorm(outer(m$mu, xseq, function(mu, x) x - mu) / m$sd) / m$sd) })
sc <- 0.9 / max(dens)
plot(NA,xlim=range(xseq),ylim=range(tt)+c(0,1),xlab="LAI",ylab="Time",
     main="Multinomial filter analysis densities")
for(i in rev(seq_along(tt))) polygon(c(xseq,rev(xseq)), tt[i] + c(dens[,i]*sc, rep(0,length(xseq))),
                                     col=col.alpha("#0c4c8a",0.6), border="white")
points(n.true[tt], tt, pch=4, col=2)

## comparison of the three filters (95% / IQR as in the Rmd)
for(lev in list(c(1,5), c(2,4))){
  plot(time,n.obs,ylim=range(c(CI.fN[lev,], n.obs, an.stats[lev,])),xlim=xlim,type='n',
       ylab="LAI",xlab="time",main=sprintf("Filter comparison (%s)", if(lev[1]==1) "95%" else "IQR"))
  band(time,an0.stats[lev[1],],an0.stats[lev[2],],'#FFBF6C')
  band(time,an.stats[lev[1],],an.stats[lev[2],],col.alpha('#b2df8a',0.7))
  band(time,CI.fN[lev[1],],CI.fN[lev[2],],col.alpha("blue",0.3))
  lines(an0.stats[3,],col="orange",lwd=2)
  lines(an.stats[3,],col="#b2df8a",lwd=2)
  lines(CI.fN[3,],col="blue",lwd=2)
  lines(Dbar,lwd=2)
  lines(time,n.true,lty=2)
  points(n.obs,pch=18,col=2)
  legend("topright",legend=c("EnKF, no disturb","EnKF, disturb","Multinomial","data","P(dist)","truth"),
         lwd=c(5,5,5,5,5,1),lty=c(1,1,1,1,1,2),col=c("orange","#b2df8a","blue","red","black","black"))
}

dev.off()

## ---- Save outputs -------------------------------------------------------------
data = list(Y = n.obs)
saveRDS(list(data=data, n.true=n.true, d=d, time=time, NT=NT, xlim=xlim,
             n0.stats=n0.stats, n.stats=n.stats,
             fx0.stats=fx0.stats, an0.stats=an0.stats, fx.stats=fx.stats, an.stats=an.stats,
             CI.fN=CI.fN, Dbar=Dbar, or.fit=or.fit, or.known=or.known),
        file.path(outdir, "results.rds"))

## headline numbers
rmse <- function(x) sqrt(mean((x - n.true)^2))
cat(sprintf("analysis-median RMSE: EnKF no-dist %.3f | EnKF dist %.3f | multinomial %.3f (tau known %.3f)\n",
            rmse(an0.stats[3,]), rmse(an.stats[3,]), rmse(or.fit$q[3,]), rmse(or.known$q[3,])))
cat(sprintf("P(disturbed) at t=10: %.3f (tau known %.3f); max elsewhere %.3f\n",
            or.fit$pDist[10], or.known$pDist[10], max(or.fit$pDist[-10])))
cat(sprintf("process sd: posterior sqrt(E[1/tau]) at t=NT %.3f (truth %.3f)\n",
            sqrt(or.fit$procVar[NT]), n.sd))
