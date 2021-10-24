
#' @title Asymmetric Laplace Distribution
#' @aliases dAsymLaplace, rAsymLaplace, qAsymLaplace, pAsymLaplace 
#' @name  dAsymLaplace
#' @description  Density, distribution function, quantile function, and random generation for a Asymmetric Laplace Distribution.
#' @param x  value to be computed.
#' @param q  quantile to be computed.
#' @param p  probability to be computed.
#' @param mu location parameter.
#' @param sigma scale parameter.
#' @param tau asymmetry parameter.
#' @param log,log.p TRUE or 1 to return log probability. FALSE or 0 to return probability.
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x] otherwise, P[X > x].
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
      res =  mu + ((sigma * log(u/tau))/(1 - tau))
    }
    if (u >= tau)
    {
      res =  mu - ((sigma * log((1 - u)/(1 - tau)))/tau)
    }
    return(res)
  })

#' @rdname dAsymLaplace
#' @export

pAsymLaplace <- nimbleFunction(
  run = function(q = double(0), mu = double(0), sigma = double(0),
		 tau = double(0), lower.tail = integer(0, default = 1),
		 log.p = integer(0, default = 0)){
    returnType(double(0))
    if (q < mu)
    {
	    res = tau * exp((1 - tau) * (q - mu)/sigma)
    }
	   
    if (q >= mu)
    {

	    res = 1 - (1 - tau) * exp(-(tau) * (q - mu)/sigma)
    }

    if(!lower.tail) { 
	    res = 1 - res
	    if(log.p) return(log(res))
	    else return(res)
    } else {
	    if(log.p) return(log(res))
	    else return(res)
    }
  })


#' @rdname dAsymLaplace
#' @export

qAsymLaplace <- nimbleFunction(
  run = function(p = double(0), mu = double(0), sigma = double(0),
		 tau = double(0), lower.tail = integer(0, default = 1),
		 log.p = integer(0, default = 0)){
    returnType(double(0))

    if (!lower.tail) {p = 1 - p}
    if (p < tau)
    {
	    res =  mu + ((sigma * log(p/tau))/(1 - tau))
    }

    if (p >= tau)
    {
	    res =  mu - ((sigma * log((1 - p)/(1 - tau)))/tau)
    }


    if(log.p) return(log(res))
	    else return(res)
  })

