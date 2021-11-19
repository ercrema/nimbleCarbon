#' @title WAIC-based model comparison
#'
#' @description  Compute delta WAIC and WAIC weights for model comparison.
#' @param ... MCMC output from either \code{\link{nimbleMCMC}} or \code{\link{runMCMC}} functions in the \code{nimble} R package. Note that in argument \code{WAIC} should be set to TRUE.
#' @return A table containing WAIC, delta WAIC, and WAIC weights.
#' @export 


compare.models = function(...)
{
  models <- list(...)
  L = length(models)
  model.names <- match.call()
  model.names <- as.character(model.names)[2:(L + 1)]
  waic=unlist(lapply(models,function(x){x$WAIC$WAIC}))
  if (is.null(waic)){stop("MCMC samples do not contain WAIC values. Rerun 'nimbleMCMC()' setting the 'WAIC' argument to TRUE.")}
  deltaWAIC = waic - min(waic)
  w = exp(-0.5*deltaWAIC)/sum(exp(-0.5*deltaWAIC))
  res=data.frame(WAIC=waic,deltaWAIC=deltaWAIC,Weights=w,row.names=model.names)
  res = res[order(res$deltaWAIC),]
  return(res)
}
