#' @title Plot Marginal Posterior Distribution
#'
#' @description Plot marginal posterior distribution highlighting user-defined higher posterior density interval.
#' @param x Posterior samples
#' @param prob Highest posterior density interval. Default is 0.9.
#' @param bw The smoothing bandwidth to be used. See \code{\link{density}} for details. Default is "SJ".
#' @param hpd.col Fill colour for the highest density interval. Default is 'lightblue'. Ignored when \code{HPD} is set to FALSE.
#' @param line.col Line color for the density plot. Default is 'darkgrey'.
#' @param HPD Whether the highest posterior density interval is highlighted or not. Default is TRUE.
#' @param show.hpd.val Whether the highest posterior density interval is displayed as subtitle. Default is TRUE.
#' @param rnd Integer indicating the number of decimal places to be used in the reporting of the highest posterior density interval.
#' @param ... other graphical parameters.
#' @return None.
#' @import grDevices
#' @import graphics
#' @import utils
#' @import rcarbon
#' @import coda
#' @export 

postHPDplot = function(x,prob=0.9,bw='SJ',hpd.col='lightblue',line.col='darkgrey',rnd=3,HPD=TRUE,show.hpd.val=TRUE,...)
{
  x.lo = HPDinterval(mcmc(x),prob = prob)[1]
  x.hi = HPDinterval(mcmc(x),prob = prob)[2]
  dens=density(x,bw = bw)
  plot(dens$x,dens$y,type='n',...)
  hpdi.x = dens$x[which(dens$x>=x.lo&dens$x<=x.hi)]
  hpdi.y = dens$y[which(dens$x>=x.lo&dens$x<=x.hi)]
  if (HPD){polygon(x=c(hpdi.x,rev(hpdi.x)),y=c(hpdi.y,rep(0,length(hpdi.y))),border=NA,col=hpd.col)}
  lines(dens,col=line.col)
  if (show.hpd.val){title(sub=paste0(prob*100,'%HPDI:',round(x.lo,rnd),'~',round(x.hi,rnd)))}
}

