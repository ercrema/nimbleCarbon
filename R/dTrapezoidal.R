#' @title Trapezoidal Distribution
#' @aliases dTrapezoidal rTrapezoidal 
#' @name  dTrapezoidal
#' @description  Density and random generation of an Trapezoidal distribution.
#' @param x A calendar year (in BP).
#' @param a lower (earliest) limit of the distribution (in BP).
#' @param m1 lower mode (in BP)
#' @param m2 upper mode (in BP).
#' @param b upper (latest) limit of the distribution (in BP).
#' @param log TRUE or 1 to return log probability. FALSE or 0 to return probability.
#' @param n number of random draws. Currently only n = 1 is supported, but the argument exists for standardization of "r" functions.
#' @author Enrico Crema

NULL
#' @examples
#' a=7000
#' b=6700
#' c=4000
#' d=3000
#' x=5400
#' modelPlot(dTrapezoidal,a=7000,b=5000,params=c(m1=6000,m2=5300),alpha=1,col=1)
#'

#' @rdname dTrapezoidal
#' @import nimble
#' @export

dTrapezoidal=nimbleFunction(
  run = function(x = integer(0),a=double(0),m1=double(0),m2=double(0),b=double(0), log = integer(0)) {
    returnType(double(0))
    #normalizing_constant = 2/(a + b - d - c)
    normalizing_constant = 2/(a + m1 - b - m2)
    dens = 0
      #if ((d <= x[i]) & (x[i] < c))
      if ((b <= x) & (x < m2))
      {
        #dens[i] = normalizing_constant * (x[i]-d)/(c-d);
        dens = normalizing_constant * (x-b)/(m2-b);
      }
      #if ((c <= x[i]) & (x[i] < b))
      if ((m2 <= x) & (x < m1))
      {
        dens = normalizing_constant
      }
      #if ((b <= x[i]) & (x[i] <= a))
      if ((m1 <= x) & (x <= a))
      {
        #dens[i] = normalizing_constant * (a-x[i])/(a-b);
        dens = normalizing_constant * (a-x)/(a-m1);
      }
    if(log) {
      return(log(dens))
    } else {
      return(dens)
    }
  })   

#' @rdname dTrapezoidal
#' @export
#' 
rTrapezoidal=nimbleFunction(
  run = function(n = integer(0),a=double(0),m1=double(0),m2=double(0),b=double(0)) {
    returnType(double(0))
    
    p = runif(1);
    
    #pi1 = (b-a)/(c+d-a-b)
    #pi1 = (c-d)/(b+a-d-c)
    pi1 = (m2-b)/(m1+a-b-m2)
    #pi2 = (2*c -2*b)/(c+d-a-b)
    #pi2 = (2*b -2*c)/(b+a-d-c)
    pi2 = (2*m1 -2*m2)/(m1+a-b-m2)
    #pi3 = (d-c)/(c+d-a-b)  
    #pi3 = (a-b)/(b+a-d-c)
    pi3 = (a-m1)/(m1+a-b-m2)
    q = NaN
    
    if ((0 <= p) & (p <= pi1))
    {
      #q = sqrt(p)*sqrt(c+d-a-b)*sqrt(b-a)+a
      #q = sqrt(p)*sqrt(b+a-d-c)*sqrt(c-d)+d
      q = sqrt(p)*sqrt(m1+a-b-m2)*sqrt(m2-b)+b
    }
    if ((pi1 < p) & (p <= (1 - pi3)))
    {
      #q = b + (((p - pi1) / (pi2)) * (c - b))
      #q = c + (((p - pi1) / (pi2)) * (b - c))
      q = m2 + (((p - pi1) / (pi2)) * (m1 - m2))
    }
    if (((1 - pi3) < p) & (p <= 1))
    {
      #q = d + (-d * sqrt(1-p) * sqrt(c+d-a-b)+c*sqrt(1-p)*sqrt(c+d-a-b))/(sqrt(d-c))
      #q = a + (-a * sqrt(1-p) * sqrt(b+a-d-c)+b*sqrt(1-p)*sqrt(b+a-d-c))/(sqrt(a-b))
      q = a + (-a * sqrt(1-p) * sqrt(m1+a-b-m2)+m1*sqrt(1-p)*sqrt(m1+a-b-m2))/(sqrt(a-m1))
    }
    return(q)
  })

#' Variant of dTrapezoidal for modelPlot
#' @param x A calendar year (in BP).
#' @param a lower (earliest) limit of the distribution (in BP).
#' @param m1 lower mode (in BP)
#' @param m2 upper mode (in BP).
#' @param b upper (latest) limit of the distribution (in BP).
#' @noRd

dTrapezoidal_modelPlot<-function(x,a,m1,m2,b,log)
{
  normalizing_constant = 2/(a + m1 - b - m2)
  n = length(x)
  dens=numeric(length=n)
  
  for (i in 1:n)
  {
    #if ((d <= x[i]) & (x[i] < c))
    if ((b <= x[i]) & (x[i] < m2))
    {
      #dens[i] = normalizing_constant * (x[i]-d)/(c-d);
      dens[i] = normalizing_constant * (x[i]-b)/(m2-b);
    }
    #if ((c <= x[i]) & (x[i] < b))
    if ((m2 <= x[i]) & (x[i] < m1))
    {
      dens[i] = normalizing_constant
    }
    #if ((b <= x[i]) & (x[i] <= a))
    if ((m1 <= x[i]) & (x[i] <= a))
    {
      #dens[i] = normalizing_constant * (a-x[i])/(a-b);
      dens[i] = normalizing_constant * (a-x[i])/(a-m1);
    }
    if ((x[i]>a) | (x[i]<b))
    {
      dens[i] = 0 
    }
  }
  return(dens)
}
