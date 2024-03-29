# This file registers all distributions when the package is loaded.
.onAttach <- function(libname, pkgname) {
  
  packageStartupMessage("Loading nimbleCarbon. \nRegistering distrubutions and compiling functions")
  suppressMessages({
    registerDistributions(list(
      dLogisticGrowth = list(
        BUGSdist = "dLogisticGrowth(a,b,k,r)",
        Rdist = "dLogisticGrowth(a,b,k,r)",
        pqAvail = FALSE
      )))
    registerDistributions(list(
      dLogisticGrowth2 = list(
        BUGSdist = "dLogisticGrowth2(a,b,m,r)",
        Rdist = "dLogisticGrowth2(a,b,m,r)",
        pqAvail = FALSE
      )))
    registerDistributions(list(
      dExponentialGrowth = list(
        BUGSdist = "dExponentialGrowth(a,b,r)",
        Rdist = "dExponentialGrowth(a,b,r)",
        pqAvail = FALSE
      )))
    registerDistributions(list(
      dDoubleExponentialGrowth = list(
        BUGSdist = "dDoubleExponentialGrowth(a,b,r1,r2,mu)",
        Rdist = "dDoubleExponentialGrowth(a,b,r1,r2,mu)",
        pqAvail = FALSE
      )))
    registerDistributions(list(
      dExponentialLogisticGrowth = list(
        BUGSdist = "dExponentialLogisticGrowth(a,b,k,r1,r2,mu)",
        Rdist = "dExponentialLogisticGrowth(a,b,k,r1,r2,mu)",
        pqAvail = FALSE
      )))
	registerDistributions(list(
	  dLogisticExponentialGrowth= list(
        BUGSdist = "dLogisticExponentialGrowth(a,b,r1,r2,k,mu)",
        Rdist = "dLogisticExponentialGrowth(a,b,r1,r2,k,mu)",
        pqAvail = FALSE
      )))
	registerDistributions(list(
	  dTrapezoidal= list(
	    BUGSdist = "dTrapezoidal(a,b,m1,m2)",
	    Rdist = "dTrapezoidal(a,b,m1,m2)",
	    pqAvail = FALSE
	  )))
	registerDistributions(list(
	  dAsymLaplace= list(
	    BUGSdist = "dAsymLaplace(mu,sigma,tau)",
	    Rdist = "dAsymLaplace(mu,sigma,tau)",
	    pqAvail = TRUE
	  )))
  })
}
