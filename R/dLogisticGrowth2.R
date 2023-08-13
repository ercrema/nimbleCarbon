#' @title Logistic Growth Model (parametrisation with inflection point)
#' @aliases dLogisticGrowth2 rLogisticGrowth2 
#' @name  dLogisticGrowth2
#' @description  Density and random generation of a logistic growth model distribution with alternative parametrisation.
#' @param x vector of calendar years (in BP).
#' @param a lower (earliest) limit of the distribution (in BP).
#' @param b upper (latest) limit of the distribution (in BP).
#' @param m inflection point in BP (i.e. point with the highest growth rate).
#' @param r intrinsic growth rate.
#' @param log TRUE or 1 to return log probability. FALSE or 0 to return probability.
#' @param n number of random draws. Currently only n = 1 is supported, but the argument exists for standardization of "r" functions.
#' @details This is an alternative parametrisation of \code{dLogisiticGrowth}, where \code{m} is equal to \code{log(-(k-1)/k)/r}.
#' @return For \code{dLogisticGrowth2}: the probability (or likelihood) or log probability of an observed date x (in Cal BP). For \code{rLogisticGrowth2} a simulated date in Cal BP. 
#' @author Enrico Crema
NULL
#' @examples
#' p = list(m=4500,r=0.007)
#' modelPlot(model = dLogisticGrowth2,a=6000,b=4000,params=p,alpha = 1)
#' @rdname dLogisticGrowth2
#' @import nimble
#' @export


dLogisticGrowth2=nimbleFunction(
  run = function(x = integer(0),a=double(0),b=double(0),m=double(0),r=double(0), log = integer(0)) {
    returnType(double(0))
    t = 1:(abs(b-a)+1)
    m2 = a-m
    n = 1/(1+exp(-r*(t-m2)))
    p = n/sum(n)
    # This last bit would be the same for any model
    logProb = dcat(a-x+1,prob=p,log=TRUE)
    if(log) {
      return(logProb)
    } else {
      return(exp(logProb))
    }
  })   

#' @rdname dLogisticGrowth2
#' @export
rLogisticGrowth2 = nimbleFunction(
  run = function(n=integer(0),a=double(0),b=double(0),m=double(0),r=double(0)) {
    returnType(double(0))
    t = 1:(abs(b-a)+1)
    m2 = a-m
    pop = 1/(1+exp(-r*(t-m2)))
    p = pop/sum(pop)
    res=a-rcat(n=1,prob=p)+1
    return(res)
  })
