tobit <- function(x){
  x[x < 0] <- 0
  x
}

tobit2 <- function(x,upper){
  x[x < 0] <- 0
  x[x>upper] <- upper
  x
}


disturbance <- function(x,                         ## x = Bleaf, Bsoil, Bstem
                        mu0 = c(1,1),              ## mean leaf & stem biomass after disturbance
                        V0 = diag(c(0.15,0.25))^2, ## var in post disturbance leaf and stem biomass, 
                        alloc.soil = 0.25){        ## fraction of removed C that goes in the soil
  old.biomass <- x[c(1,3)]
  new.biomass <- tobit(mvtnorm::rmvnorm(1,mu0,V0)) ## draw disturbed leaf and stem
  check.max <- which(new.biomass > old.biomass)
  if(length(check.max)>0) new.biomass[check.max] = old.biomass[check.max]
  residual = sum(old.biomass-new.biomass)
  x[c(1,3)] <- new.biomass
  x[2] <- x[2] + residual*alloc.soil
  removal <- residual*(1-alloc.soil)
  return(x)
}

disturbance.p <- function(x,                         ## x = Bleaf, Bsoil, Bstem
                        mu0 = c(1,1),              ## mean leaf & stem % retention after disturbance
                        V0 = diag(c(0.15,0.25))^2, ## var in post disturbance leaf and stem biomass, 
                        alloc.soil = 0.25){        ## fraction of removed C that goes in the soil
  old.biomass <- x[c(1,3)]
  new.biomass <- tobit(mvtnorm::rmvnorm(1,mu0*old.biomass,tcrossprod(old.biomass,old.biomass)*V0)) ## draw disturbed leaf and stem
  check.max <- which(new.biomass > old.biomass)
  if(length(check.max)>0) new.biomass[check.max] = old.biomass[check.max] # constrain new.biomass[i] <= old.biomass[i]
  residual = sum(old.biomass-new.biomass)
  x[c(1,3)] <- new.biomass # update state with disturbed biomass
  x[2] <- x[2] + residual*alloc.soil # some of the biomass gets transported to the soil
  removal <- residual*(1-alloc.soil)
  return(x)
}

## disturbance.p applied to every row of X (n x 3: Bleaf, Bsoil, Bstem) in one call.
## For diagonal V0 the draw covariance tcrossprod(old)*V0 is diagonal, so its square
## root is old*sqrt(diag(V0)); the rnorm draws are consumed in the same order as n
## sequential disturbance.p calls, giving the same result for the same RNG state.
## Non-diagonal V0 falls back to the row-by-row loop.
disturbance.p.rows <- function(X, mu0 = c(1,1), V0 = diag(c(0.15,0.25))^2, alloc.soil = 0.25){
  if (any(V0[upper.tri(V0)] != 0)) {
    for (n in seq_len(nrow(X))) X[n, ] <- disturbance.p(X[n, ], mu0, V0, alloc.soil)
    return(X)
  }
  old.biomass <- X[, c(1,3), drop = FALSE]
  z <- matrix(rnorm(2 * nrow(X)), ncol = 2, byrow = TRUE)
  new.biomass <- old.biomass * rep(mu0, each = nrow(X)) + abs(old.biomass) * z * rep(sqrt(diag(V0)), each = nrow(X))
  new.biomass <- pmin(tobit(new.biomass), old.biomass)   # 0 <= new.biomass <= old.biomass
  residual <- rowSums(old.biomass - new.biomass)
  X[, c(1,3)] <- new.biomass
  X[, 2] <- X[, 2] + residual*alloc.soil
  X
}

disturbance2 <- function(x,type){ ## x = Bleaf, Bsoil, Bstem
  old.biomass <- x[c(1,3)]  # select leaf/stem biomass
  Mu0 <- mu0[,type] * old.biomass                         ## mean for disturbed biomass (proportional to old biomass)
  new.biomass <- tobit(mvtnorm::rmvnorm(1,Mu0,V0[,,type])) ## draw disturbed leaf and stem from distribution for type
  check.max <- which(new.biomass > old.biomass)
  if(length(check.max)>0) new.biomass[check.max] = old.biomass[check.max]   # constrain new.biomass[i] <= old.biomass[i]
  residual = sum(old.biomass-new.biomass)
  x[c(1,3)] <- new.biomass
  x[2] <- x[2] + residual*alloc.soil  # add some of the biomass lost back to the soil
  removal <- residual*(1-alloc.soil)  # not used?
  return(x)
}

# random correlated bernoulli
rcbern <- function(p,rho,mu){
  
  ## calculate cov
  D = as.matrix(dist(seq_along(mu),diag = TRUE,upper = TRUE))
  SIGMA <- 1/(1-rho^2)*rho^D # TODO: why this covariance matrix?
  
  ## draw rMVN
  x = mvtnorm::rmvnorm(1,rep(0,length(p)),SIGMA)
  
  ## inverse distribution transform
  y = pnorm(x) # normal distribution function
  qbinom(y,1,p) # map into binomial quantiles
  
}