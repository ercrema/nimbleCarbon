#' @title Exponential Growth Model
#' @aliases dExponentialGrowth rExponentialGrowth 
#' @name  dExponentialGrowth
#' @description  Density and random generation of an exponential growth model distribution.
#' @param x vector of calendar years (in BP).
#' @param a lower (earliest) limit of the distribution (in BP).
#' @param b upper (latest) limit of the distribution (in BP).
#' @param r intrinsic growth rate.
#' @param log TRUE or 1 to return log probability. FALSE or 0 to return probability.
#' @param n number of random draws. Currently only n = 1 is supported, but the argument exists for standardization of "r" functions.
#' @return For \code{dExponentialGrowth}: the probability (or likelihood) or log probability of an observed date x (in Cal BP). For \code{rExponentialGrowth} a simulated date in Cal BP. 
#' @author Enrico Crema

NULL
#' @examples
#' p = list(r=0.002)
#' modelPlot(model = dExponentialGrowth,a=6000,b=4000,params=p,alpha = 1)
#' @rdname dExponentialGrowth
#' @import nimble
#' @export

dExponentialGrowth=nimbleFunction(
  run = function(x = integer(0),a=double(0),b=double(0),r=double(0), log = integer(0)) {
    returnType(double(0))
    t = 1:(abs(b-a)+1)
    n = numeric(abs(b-a)+1)
    tfinal = abs(b-a)+1
    for (i in 1:tfinal)
    {
      n[i] = (1+r)^t[i]
    }
    p = n/sum(n)
    # This last bit would be the same for any model
    logProb = dcat(a-x+1,prob=p,log=TRUE)
    if(log) {
      return(logProb)
    } else {
      return(exp(logProb))
    }
  })   

#' @rdname dExponentialGrowth
#' @export
rExponentialGrowth = nimbleFunction(
  run = function(n=integer(0),a=double(0),b=double(0),r=double(0)) {
    returnType(double(0))
    t = 1:(abs(b-a)+1)
    pop = numeric(abs(b-a))
    tfinal = abs(b-a)+1
    for (i in 1:tfinal)
    {
      pop[i] = (1+r)^t[i]
    }
    p = pop/sum(pop)
    res=a-rcat(n=1,prob=p)+1
    return(res)
  })
