if(getRversion() >= "2.15.1")  utils::globalVariables(c("i"))

#' @title SPD-based Posterior Predictive Check
#'
#' @description Generates SPDs from posterior samples.
#' @param x a vector of observed uncalibrated radiocarbon ages.
#' @param errors a vector of standard deviations corresponding to each estimated radiocarbon age.
#' @param calCurve character string naming a calibration curve already provided with the rcarbon package (currently 'intcal20','intcal13','intcal13nhpine16','shcal20','shcal13','shcal13shkauri16',''marine13','marine20').
#' @param model growth model
#' @param a lower (earliest) limit of the distribution (in BP).
#' @param b upper (latest) limit of the distribution (in BP).
#' @param params list of vectors containing model parameters. The names attribute of each vector should match growth model parameters.
#' @param nsim number of SPDs to be generated.  Default is the length of the parameter vectors supplied in the argument \code{params}.
#' @param method method for the creation of random dates from the fitted model. Either 'uncalsample' or 'calsample'.
#' @param datenormalised a logical variable indicating whether dates should be normalised to sum to unity or not. Default is TRUE.
#' @param spdnormalised a logical variable indicating whether the total probability mass of the SPD is normalised to sum to unity for both observed and simulated data. Default is TRUE.
#' @param ncores number of cores used for for parallel execution. Default is 1.
#' @param verbose a logical variable indicating whether extra information on progress should be reported. Default is TRUE.
#' @return An object of class \code{spdppc} with the following elements
#' \itemize{
#' \item{\code{obs}} {A data.frame containing the years (in Cal BP) and the corresponding summed probability in the observed data.}
#' \item{\code{spdmat}} {A matrix containing the summed probability distribution of the simulated data.}
#' }
#' @import grDevices
#' @import graphics
#' @import doSNOW 
#' @import snow
#' @import foreach 
#' @import rcarbon
#' @export 

postPredSPD = function(x,errors,calCurve,model,a,b,params,nsim,method=NULL,spdnormalised=TRUE,datenormalised=TRUE,ncores=1,verbose=TRUE)
{
  #Sanity checks ####
  if (length(unique(unlist(lapply(params,length))))>1)
  {
    print('Number of supplied sample parameters should be the same.')
    stop()
  }
  
  if (unique(unlist(lapply(params,length)))<nsim)
  {
    print(paste0('The argument nsim is larger than the supplied length of parameter values. Running with nsim=',unique(unlist(lapply(params,length)))))
    stop()
  }
  
  if (ncores>1&!requireNamespace("doSNOW", quietly=TRUE)){	
    warning("the doSNOW package is required for multi-core processing; ncores has been set to 1")
    ncores=1
  }	
  
  if (!any(method%in%c("uncalsample","calsample"))|is.null(method))
  {
    stop("The 'method' argument must be either 'uncalsample' or 'calsample'")
  }
  
  
  
  # General Settings ####
  #randomly sample nsim parameter values
  index = sample(nsim,size=nsim)
  params = lapply(params,function(x,i){x[i]},i=index)
  ndates = length(x)
  calBP = a:b
  
  #Calibrate Observed Data and Generate SPD ####
  if(verbose){'Generate observed SPD'}
  calibrated.dates=calibrate(x,errors,calCurves = calCurve,verbose=verbose, normalised=datenormalised)
  obs.spd = spd(calibrated.dates,timeRange=c(a,b),spdnormalised = spdnormalised,verbose=verbose)

  
  
  
  
  #Iterate through posterior ####
  
  # Setup foreach loop:
  opts = NULL
  if (verbose)
  {
    if (ncores>1){ print(paste("Running in parallel on ",getDoParWorkers()," workers...",sep=""))}
    pb <- txtProgressBar(min=0, max=nsim, style=3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
  }
  
  `%dofun%` <- `%do%`
  if (ncores>1){
    cl <- snow::makeCluster(ncores)
    registerDoSNOW(cl)
    `%dofun%` <- `%dopar%`
  }
  
  # Loop starts here
  predMat <- foreach (i=1:nsim,.packages = 'rcarbon',.options.snow = opts,.combine='cbind') %dofun% {
    if (verbose & ncores==1) {setTxtProgressBar(pb, i)}
    
    #Create Probability Model
    PrDens=do.call(model,args=c(list(x=a:b,a=a,b=b),lapply(params,function(x,i){x[[i]]},i=i),list(log=FALSE)))
  
    #Uncalibrate
    uncalGrid=uncalibrate(as.CalGrid(data.frame(calBP=calBP,PrDens=PrDens)),verbose=FALSE)
    
    #Sample Dates
    if (method=='uncalsample')
    {
      x.tmp=sample(uncalGrid$CRA,replace=TRUE,size=ndates,prob=uncalGrid$PrDens)
    }
    if (method=='calsample')
    {
      x.tmp=sample(uncalGrid$CRA,replace=TRUE,size=ndates,prob=uncalGrid$Raw)
    }
    
    #Calibrate
    x.tmp.calibrated = calibrate(x.tmp,errors = sample(errors,size=ndates,replace=TRUE),calCurves = calCurve, normalised=datenormalised, verbose=FALSE)
    
    #Make SPD
    x.tmp.spd = spd(x.tmp.calibrated,timeRange = c(a,b),spdnormalised = spdnormalised,verbose=FALSE)
    return(x.tmp.spd$grid$PrDens)
  }
  if (ncores>1){
    snow::stopCluster(cl)
  }
  if (verbose){ close(pb) }
  
  out = list(obs=obs.spd$grid,simmatrix=predMat)
  class(out) <- c("spdppc",class(out))
  return(out)
}
