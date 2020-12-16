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
#' @param n number of random draws, each returning a vector of length len. Currently only n = 1 is supported, but the argument exists for standardization of "r" functions.
#' @author Enrico Crema

NULL

#' @rdname dDoubleExponentialGrowth
#' @import nimble
#' @export


# a=7000
# b=3000
# mu=6000
# nt0 = 1
# nt1 =666666
# t =  1:(abs(b-a)+1)
# t1 = 1:(abs(mu-a))
# t2 = 1:(abs(b-mu)+1)
# r1=0.004
# r2=-0.004
# n = c(nt0*(1+r1)^c(t1),(nt0*(1+r1)^abs(mu-a))*(1+r2)^t2)
# n1 = c(nt1*(1+r1)^c(t1),(nt1*(1+r1)^abs(mu-a))*(1+r2)^t2)
# plot(a:b,n,type='l',xlim=c(a,b))
# abline(v=mu,lty=2)
# all.equal(n/sum(n),n1/sum(n1))

#dDoubleExponentialGrowth(2000,a=3400,b=1850,r1=0.01,r2=0.003,mu=2000,log=TRUE)

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