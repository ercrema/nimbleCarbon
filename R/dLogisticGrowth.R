#' @title Logistic Growth Model
#' @aliases dLogisticGrowth rLogisticGrowth 
#' @name  dLogisticGrowth
#' @description  Density and random generation of a logistic growth model distribution.
#' @param x vector of calendar years (in BP).
#' @param a lower (earliest) limit of the distribution (in BP).
#' @param b upper (latest) limit of the distribution (in BP).
#' @param k initial proportion of the carrying capacity (must be between 0 and 1).
#' @param r intrinsic growth rate.
#' @param log TRUE or 1 to return log probability. FALSE or 0 to return probability.
#' @param n number of random draws. Currently only n = 1 is supported, but the argument exists for standardization of "r" functions.
#' @return For \code{dLogisticGrowth}: the probability (or likelihood) or log probability of an observed date x (in Cal BP). For \code{rLogisticGrowth} a simulated date in Cal BP. 
#' @author Enrico Crema
NULL
#' @examples
#' p = list(k=0.01,r=0.007)
#' modelPlot(model = dLogisticGrowth,a=6000,b=4000,params=p,alpha = 1)
#' @rdname dLogisticGrowth
#' @import nimble
#' @export


dLogisticGrowth=nimbleFunction(
  run = function(x = integer(0),a=double(0),b=double(0),k=double(0),r=double(0), log = integer(0)) {
    returnType(double(0))
    t = 1:(abs(b-a)+1)
    n = 1/(1+((1-k)/k)*exp(-r*t))
    p = n/sum(n)
    # This last bit would be the same for any model
    logProb = dcat(a-x+1,prob=p,log=TRUE)
    if(log) {
      return(logProb)
    } else {
      return(exp(logProb))
    }
  })   

#' @rdname dLogisticGrowth
#' @export
rLogisticGrowth = nimbleFunction(
  run = function(n=integer(0),a=double(0),b=double(0),k=double(0),r=double(0)) {
    returnType(double(0))
    t = 1:(abs(b-a)+1)
    pop = 1/(1+((1-k)/k)*exp(-r*t))
    p = pop/sum(pop)
    res=a-rcat(n=1,prob=p)+1
    return(res)
  })