#' @title Calculate Agreement Indices.
#' @description  Computes OxCal-style (Bronk-Ramsey 1995) individual and overall agreement index for evaluating model consistency.
#' @param CRA vector of C14 ages.
#' @param CRAError vector of C14 errors associated with \code{CRA}. 
#' @param calCurve character string naming a calibration curve, one between 'intcal20','intcal13','shcal20','shcal13','marine13' and 'marine20'. 
#' @param theta a Matrix containing the posterior samples of each date. 
#' @param verbose a logical variable indicating whether extra information on progress should be reported. Default is TRUE.
#' @return a list containing the individual and overall agreement indices.
#' @references Bronk-Ramsey, C. (1995). Radiocarbon Calibration and Analysis of Stratigraphy: The OxCal Program. Radiocarbon, 37(2), 425–430.
#' @import rcarbon
#' @export 

agreementIndex = function(CRA,CRAError,calCurve='intcal20',theta,verbose=TRUE)
{
  p = calibrate(CRA,CRAError,calCurves = calCurve,verbose=verbose)
  agreement = numeric(length=length(p))
  if (verbose)
  {
    print('Computing agreement indices...')
    pb <- txtProgressBar(min=0, max=length(p), style=3)
  }
  
  for (i in 1:length(p))
  {
    if (verbose)
    {
      setTxtProgressBar(pb, i)
    }
    p.tmp = p$grids[[i]]
    p.prime.samples = round(theta[,i])
    p.prime.samples=factor(p.prime.samples,levels=max(p.prime.samples):min(p.prime.samples))
    p.prime.probs = table(p.prime.samples)/sum(table(p.prime.samples)) 
    p.prime.t = as.numeric(names(p.prime.probs))
    p.t = p.tmp$calBP 
    p.prime = data.frame(calBP=p.prime.t,PrDens=p.prime.probs)
    xx=merge(x=p.tmp,y=p.prime,by.x='calBP',by.y='calBP',all=TRUE)
    xx$PrDens.Freq[which(is.na(xx$PrDens.Freq))]=0
    xx$PrDens[which(is.na(xx$PrDens))]=0
    agreement[i]=sum(xx$PrDens*xx$PrDens.Freq)/  sum(xx$PrDens*xx$PrDens)
  }
  if (verbose){close(pb)}
  overall.agreement = (prod(agreement))^(1/sqrt(length(p)))
  return(list(agreement=agreement*100,overall.agreement=overall.agreement*100))
}
