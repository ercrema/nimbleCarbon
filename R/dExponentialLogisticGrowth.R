#' @title Exponential-Logistic Growth Model
#' @aliases dExponentialLogisticGrowth rExponentialLogisticGrowth 
#' @name  dExponentialLogisticGrowth
#' @description  Density and random generation of a exponential-logistic growth model distribution.
#' @param x vector of calendar years (in BP).
#' @param a lower (earliest) limit of the distribution (in BP).
#' @param b upper (latest) limit of the distribution (in BP).
#' @param k initial proportion of the carrying capacity (must be between 0 and 1).
#' @param r1 growth rate of the exponential phase.
#' @param mu change point (in BP).
#' @param r2 growth rate of logistic phase.
#' @param log TRUE or 1 to return log probability. FALSE or 0 to return probability.
#' @param n number of random draws. Currently only n = 1 is supported, but the argument exists for standardization of "r" functions.
#' @return For \code{dExponentialLogisticGrowth}: the probability (or likelihood) or log probability of an observed date x (in Cal BP). For \code{rExponentialLogisticGrowth} a simulated date in Cal BP. 
#' @author Enrico Crema

NULL
#' @examples
#' p = list(r1=-0.001,r2=0.01,mu=5200,k=0.2)
#' modelPlot(model = dExponentialLogisticGrowth,a=6000,b=4000,params=p,alpha = 1)
#' @rdname dExponentialLogisticGrowth
#' @import nimble
#' @export

dExponentialLogisticGrowth=nimbleFunction(
  run = function(x = integer(0),a=double(0),b=double(0),k=double(0), r1=double(0),r2=double(0),mu=double(0), log = integer(0)) {
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
      n1[i] = k*(1+r1)^t1[i]
    }
    for (i in 1:t2final)
    {
      n2[i] = 1/(1+((1-(k*(1+r1)^abs(mu-a)))/(k*(1+r1)^abs(mu-a)))*exp(-r2*t2[i]))
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

#' @rdname dExponentialLogisticGrowth
#' @export
rExponentialLogisticGrowth=nimbleFunction(
  run = function(n=integer(0),a=double(0),b=double(0),k=double(0), r1=double(0),r2=double(0),mu=double(0)) {
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
      pop1[i] = k*(1+r1)^t1[i]
    }
    for (i in 1:t2final)
    {
      pop2[i] = 1/(1+((1-(k*(1+r1)^abs(mu-a)))/(k*(1+r1)^abs(mu-a)))*exp(-r2*t2[i]))
    }
    pop = c(pop1,pop2)
    p = pop/sum(pop)
    res=a-rcat(n=1,prob=p)+1
    return(res)
  })   