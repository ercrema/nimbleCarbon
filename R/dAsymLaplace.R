
#' @title Asymmetric Laplace Distribution
#' @aliases dAsymLaplace, rAsymLaplace 
#' @name  dAsymLaplace
#' @description  Density and random generation for a Asymmetric Laplace Distribution.
#' @param x  value to be computed.
#' @param mu location parameter.
#' @param sigma scale parameter.
#' @param tau asymmetry parameter.
#' @param log TRUE or 1 to return log probability. FALSE or 0 to return probability.
#' @param n number of random draws. Currently only n = 1 is supported, but the argument exists for standardization of "r" functions.
#' @author Enrico Crema
NULL
#' @rdname dAsymLaplace
#' @import nimble
#' @export

dAsymLaplace <- nimbleFunction(
  run = function(x = double(0),mu = double(0), sigma = double(0), tau=double(0), log=integer(0)) {
    returnType(double(0))
    if (x<mu)
    {
      Prob=(tau*(1-tau)/sigma)*exp((1-tau)*(x-mu)/sigma)
    } 
    if (x>=mu) 
    {
      Prob=(tau*(1-tau)/sigma)*exp(-(tau)*(x-mu)/sigma)
    }
    if(log) {
      return(log(Prob))
    } else {
      return(Prob)
    }
  })

#' @rdname dAsymLaplace
#' @export 
 
rAsymLaplace <- nimbleFunction(
  run = function(n = integer(0),mu = double(0), sigma = double(0), tau=double(0)) {
    returnType(double(0))
    u = runif(1)
    if (u < tau)
    {
      res =  mu + ((sigma * log(q/tau))/(1 - tau))
    }
    if (u >= tau)
    {
      res =  mu - ((sigma * log((1 - q)/(1 - tau)))/tau)
    }
    return(res)
  })
