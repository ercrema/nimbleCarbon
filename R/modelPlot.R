#' @title Plot Growth Models
#'
#' @description Plots growth models based on user provided parameters for prior and posterior predictive checks. 
#' @param model growth model.
#' @param a lower (earliest) limit of the distribution (in BP).
#' @param b upper (latest) limit of the distribution (in BP).
#' @param params a \link{list} of vectors containing model parameters. The names attribute of each vector should match growth model parameters.
#' @param type either a 'spaghetti' plot or a quantile based 'envelope' plot. Default is 'spaghetti'.
#' @param nsample number of samples to be used. Default is the length of the parameter vectors supplied in the argument \code{params}.
#' @param interval quantile interval used for the envelope plot. Ignored when type is set to 'spaghetti'.
#' @param calendar  either \code{'BP'} or \code{'BCAD'}. Indicate whether the calibrated date should be displayed in BP or BC/AD. Default is  \code{'BP'}.
#' @param col fill color for the quantile envelope (when \code{type=='envelope'}) or line colour (when \code{type=='spaghetti'}).
#' @param alpha transparency value for each line in the spaghetti plot. Ignored when type is set to 'envelope'. Default is 0.1.
#' @param ylim the y limits of the plot.
#' @param xlim the x limits of the plot (in Cal BP).
#' @param add whether or not the new graphic should be added to an existing plot. 
#' @param ... additional arguments affecting the plot
#' @return None.
#' @examples
#' \donttest{
#' params = list(k=runif(100,0.01,0.02),r=runif(100,0.003,0.004))
#' modelPlot(model=dLogisticGrowth,a=5000,b=2000,params=params,type=c('spaghetti'),alpha=0.5)
#' }
#' @export 
modelPlot = function(model,a,b,params,type=c('spaghetti'),nsample=NULL,interval=0.9,calendar='BP',col='lightgrey',alpha=0.1,ylim=NULL,xlim=NULL,add=FALSE,...)
{
  #Check provided model is supported
  modelName <- as.character(substitute(model))
  if (!modelName%in%c('dLogisticGrowth','dLogisticExponentialGrowth','dExponentialGrowth','dDoubleExponentialGrowth','dExponentialLogisticGrowth','dTrapezoidal'))
  {
    stop(paste0(modelName,' is currently not supported'))
  }
  
  #If dTrapezoidal consider using alternative function
  if (modelName=='dTrapezoidal')
  {
    model = dTrapezoidal_modelPlot
  }
  #Check length parameters
  if (length(unique(unlist(lapply(params,length))))>1)
  {
    stop('Number of supplied sample parameters should be the same.')
  }
  supplied.nsample = unique(unlist(lapply(params,length)))
  if (!is.null(nsample)){nsample.index=sample(supplied.nsample,size=nsample)}
  if (is.null(nsample)){nsample=supplied.nsample;nsample.index=1:nsample}
  mat = matrix(NA,nrow=length(a:b),ncol=nsample)
  
  for (i in 1:nsample)
  {
    mat[,i] = do.call(model,args=c(list(x=a:b,a=a,b=b),lapply(params,function(x,i){x[[i]]},i=nsample.index[i]),list(log=FALSE)))
  }
  
  #Setting calendar and xlim
  if (calendar=="BP"){
    plotyears <- a:b
    xlabel <- "Years cal BP"
    if (!is.null(xlim)){xlim=sort(xlim,T)}
    if (is.null(xlim)){xlim <- c(max(plotyears),min(plotyears))}
  } else if (calendar=="BCAD"){
    plotyears <- BPtoBCAD(a:b)
    xlabel <- "Years BC/AD"
    if (all(range(plotyears)<0)){xlabel <- "Years BC"}
    if (all(range(plotyears)>0)){xlabel <- "Years AD"}
    if (!is.null(xlim)){xlim=BPtoBCAD(sort(xlim,T))}
    if (is.null(xlim)){xlim <- c(min(plotyears),max(plotyears))}
  } else {
    stop("Unknown calendar type")
  }
  
  
  if (type=='envelope')
  {
    lo=apply(mat,1,quantile,prob=0+(1-interval)/2)
    hi=apply(mat,1,quantile,prob=1-(1-interval)/2)
    median = apply(mat,1,median)
    if (is.null(ylim)){ylim=c(0,max(hi))}
    if(!add){plot(plotyears, median, xlim=xlim, ylim=ylim, type="n", col="white", ylab='Probability', xlab=xlabel, xaxt="n",...)}
    polygon(c(plotyears,rev(plotyears)),c(lo,rev(hi)),col=col,border=NA)
    lines(plotyears, median,lwd=2)
  }
  
  if (type=='spaghetti')
  {
    if (is.null(ylim)){ylim=c(0,max(mat))}
    if(!add){plot(0, 0, xlim=xlim, ylim=ylim, type="n", col="white", ylab='Probability', xlab=xlabel, xaxt="n",...)}
    apply(mat,2,lines,x=plotyears,col=rgb(col2rgb(col)[1]/255,col2rgb(col)[2]/255,col2rgb(col)[3]/255,alpha))
  }
  
  if (calendar=="BP"){
    rr <- range(pretty(xlim))    
    if(!add){axis(side=1,at=seq(rr[2],rr[1],-100),labels=NA,tck = -.01)
      axis(side=1,at=pretty(xlim),labels=abs(pretty(xlim)))}
  } else if (calendar=="BCAD"){
    yy <-  xlim
    rr <- range(pretty(xlim))    
    prettyTicks <- seq(rr[1],rr[2],+100)
    prettyTicks[which(prettyTicks>=0)] <-  prettyTicks[which(prettyTicks>=0)]-1
    axis(side=1,at=prettyTicks, labels=NA,tck = -.01)
    py <- pretty(yy)
    pyShown <- py
    if (any(pyShown==0)){pyShown[which(pyShown==0)]=1}
    py[which(py>1)] <-  py[which(py>1)]-1
    if(!add){axis(side=1,at=py,labels=abs(pyShown))}
  }
}