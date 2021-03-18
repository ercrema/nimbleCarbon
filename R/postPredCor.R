#' @title Calculates correlation between observed and posterior generated SPD.
#' @description  Computes the correlation between observed SPDs and posterior generated SPD from the output of \code{postPredSPD()} function as an heuristic of the goodness-of-fit of the model. 
#' @param x An object of class \code{spdppc}.
#' @param method a character string indicating which correlation coefficient is to be computed. One of "pearson" (default), "kendall", or "spearman": can be abbreviated.
#' @returns A vector of correlation values. 
#' @export 

postPredCor= function(x,method='pearson')
{
  if (!c('spdppc')%in%class(x))
  {
    stop('x must be of class "spdppc"')
  }
  xx = cbind(x$obs$PrDens,x$simmatrix)
  res = as.numeric(cor(xx,method=method)[1,-1])
  return(res)
}