##' Compute Normal Quantile
##' @param p quantile value
##' @param mean mean of normal distribution
##' @param variance variance of normal distribution
##' @returns pth quantile of normal distribution with specified mean / variance
normal_quantile <- function(p, mean, variance) {
  sd <- sqrt(variance)
  mean + sd * qnorm(p)
}