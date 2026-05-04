##' @title multi.analytical ens.adj
##' @name  multi.analytica.ens.adj
##' @author Jacob Epstein 
##' 
##' 
##' @param Xf Dataframe or matrix of forecast state variables for different ensembles. [ne x 3]
##' @param cf vector assigning forecast disturbance class to each ensemble member [ne x 1]
##' @param mu.f A vector with forecast mean estimates of state variables. [nc (num disturbance classes) x nS (num state variables)]
##' @param Pf  A vector of cov matrices of forecast state variables.  [nc (num disturbance classes) x nS x nS]
##' @param Pp.mu.a Posterior predictive forecast means of state variables. [nc x nS]
##' @param Pp.Pa Posterior predictive covariance matrices of state variables. [nc x nS x nS]
##' @param Pp.w Posterior frequencies of each class, or posterior component weights [nc x 1]
##' @param uc (optional) A vector with integer entries corresponding to each disturbance class [1 x nc] 
##' 
##' @return Returns a matrix of adjusted analysis mean estimates of state variables and class assignments
##' @export

multi.analytical.ens.adj <- function(Xf, cf, mu.f, Pf, Pp.mu.a, Pp.Pa, Pp.w, uc=NULL){
  
  ## reassign classes
  if (is.null(uc)) {
    uc = sort(unique(cf)) # unique class values sorted
  }
  ff = table_by(cf, uc)/length(cf) ## class frequency in the forecast
  df = pmax(Pp.w-ff,0) ## difference in frequency
  df = df/sum(df)    ## posterior reassignment frequency
  cA = cf            ## class assigned in the Analysis
  for(i in seq_along(uc)){
    if(Pp.w[i] < ff[i]){   ## if a class decreases in frequency in the analysis
      sel.c = which(cf == uc[i])
      cA[sel.c] = sample(uc,length(sel.c),prob=df,replace = TRUE) # sample new classes according to reassignment frequencies
    }
  }
  
  Z <- Xf*0 # [ne x 3]
  
  
  ## rescale foreacast ensemble Xf to be multivariate normal with 0 mean, identity covariance
  
  for(k in seq_along(uc)){   ## loop over disturbance classes
    sel.c = which(cf == uc[k])
    
    ## SVD of forecast covariances
    S_f  <- svd(Pf[k,,])
    L_f  <- S_f$d
    V_f  <- S_f$v
    
    ## normalize
    for(i in sel.c){
      Z[i,] <- 1/sqrt(L_f) * t(V_f)%*%(Xf[i,]-mu.f[k,])
    }
    Z[is.na(Z)]<-0
    Z[is.infinite(Z)] <- 0
  }
  
  ### ANALYSIS
  # rescale the ensemble members from Z-space (zero mean, identity cov)
  #.  into posterior space (posterior mean, posterior sample covariance)
  #   and do this on a disturbance class-by-class basis
  X_a <- Xf*0 # [ne x 3] (the datapoints in posterior space)
  for(k in seq_along(uc)){
    sel.c = which(cA == uc[k])
    if(length(sel.c) == 0) next
    
    S_a  <- svd(Pp.Pa[k,,])
    L_a  <- S_a$d
    V_a  <- S_a$v
    
    ## analysis ensemble 
    for(i in sel.c){
      # we decomposed Pa - then it's putting it back together but with a different Z which comes from the likelihood of that ens    
      X_a[i,] <- V_a %*%diag(sqrt(L_a))%*%Z[i,] + Pp.mu.a[k,]
    }
  }
  
  return(cbind(cA,X_a))
  
}