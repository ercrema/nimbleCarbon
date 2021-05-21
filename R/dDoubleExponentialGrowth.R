#' @title Double Exponential Growth Model
#' @aliases dDoubleExponentialGrowth rDoubleExponentialGrowth 
#' @name  dDoubleExponentialGrowth
#' @description  Density and random generation of an exponential growth model distribution.
#' @param x vector of calendar years (in BP).
#' @param a lower (earliest) limit of the distribution (in BP).
#' @param b upper (latest) limit of the distribution (in BP).
#' @param r1 growth rate before change point mu.
#' @param r2 growth rate after change point mu.
#' @param mu change point (in BP).
#' @param log TRUE or 1 to return log probability. FALSE or 0 to return probability.
#' @param n number of random draws. Currently only n = 1 is supported, but the argument exists for standardization of "r" functions.
#' @return For \code{dDoubleExponentialGrowth}: the probability (or likelihood) or log probability of an observed date x (in Cal BP). For \code{rDoubleExponentialGrowth} a simulated date in Cal BP. 
#' @author Enrico Crema

NULL
#' @examples
#' p = list(r1=0.003,r2=-0.001,mu=5200)
#' modelPlot(model = dDoubleExponentialGrowth,a=6000,b=4000,params=p,alpha = 1)
#' @rdname dDoubleExponentialGrowth
#' @import nimble
#' @export


dDoubleExponentialGrowth=nimbleFunction(
  run = function(x = integer(0),a=double(0),b=double(0),r1=double(0),r2=double(0),mu=double(0), log = integer(0)) {
    returnType(double(0))
    mu = round(mu)
    t1 = 1:(abs(mu-a))
    t2 = 1:(abs(b-mu)+1)
    t1final = abs(mu-a)
    t2final = abs(b-mu)+1
    n1 = numeric(abs(mu-a))
    n2 = numeric(abs(b-mu)+1)
    for (i in 1:t1final)
    {
      n1[i] = (1+r1)^t1[i]
    }
    for (i in 1:t2final)
    {
      n2[i] = ((1+r1)^abs(mu-a))  * (1+r2)^t2[i]
    }
    n = c(n1,n2)
    p = n/sum(n)
    # This last bit would be the same for any model
    logProb = dcat(a-x+1,prob=p,log=TRUE)
    if(log) {
      return(logProb)
    } else {
      return(exp(logProb))
    }
  })   

#' @rdname dDoubleExponentialGrowth
#' @export
rDoubleExponentialGrowth = nimbleFunction(
  run = function(n=integer(0),a=double(0),b=double(0),r1=double(0),r2=double(0),mu=double(0)) {
    returnType(double(0))
    mu = round(mu)
    t1 = 1:(abs(mu-a))
    t2 = 1:(abs(b-mu)+1)
    t1final = abs(mu-a)
    t2final = abs(b-mu)+1
    pop1 = numeric(abs(mu-a))
    pop2 = numeric(abs(b-mu)+1)
    for (i in 1:t1final)
    {
      pop1[i] = (1+r1)^t1[i]
    }
    for (i in 1:t2final)
    {
      pop2[i] = ((1+r1)^abs(mu-a))  * (1+r2)^t2[i]
    }
    pop = c(pop1,pop2)
    p = pop/sum(pop)
    res=a-rcat(n=1,prob=p)+1
    return(res)
  })