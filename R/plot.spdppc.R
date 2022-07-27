#' @title Plot SPD-based Posterior Predictive Check
#'
#' @description Plots \code{spdppc} class object for SPD-based Posterior Predictive Check.
#' @param x An \code{spdppc} class object.
#' @param type Either a 'spaghetti' plot or a quantile based envelope plot. Default is 'envelope'.
#' @param interval Quantile interval used for the envelope plot. Ignored when \code{type} is set to 'spaghetti'. Default is 0.90.
#' @param nsample Number of samples to be displayed in the 'spaghetti' plot. Default is the total number of simulations supplied in the 'spdppc' class object, ignored when \code{type} is set to 'envelope'.
#' @param obs.lwd Line width of the observed SPD. Default is 1.5.
#' @param obs.col Line colour of the observed SPD. Default is 'black'.
#' @param sim.col Line colour of simulated SPDs. Default is 'lightgrey', ignored when \code{type} is set to 'envelope'.
#' @param alpha Transparency value for each line in the spaghetti plot. Default is 1, ignored when \code{type} is set to 'envelope'.
#' @param envelope.col Fill colour of the simulation envelope. Default is 'lightgrey', ignored when \code{type} is set to 'envelope.'spaghetti'.
#' @param positive.col Fill colour for the area with positive deviation from the simulation envelope.  Default is 'red', ignored when \code{type} is set to 'spaghetti'.
#' @param negative.col Fill colour for the area with positive deviation from the simulation envelope.  Default is 'blue', ignored when \code{type} is set to 'spaghetti'.
#' @param xlab a label for the x axis. Default is 'Years cal BP','Years BC/AD','Years BC', or 'Years AD' depending on data range and settings of \code{calendar}.  
#' @param ylab a label for the y axis. Default is 'Probability'.
#' @param calendar  Either \code{'BP'} or \code{'BCAD'}. Indicate whether the calibrated date should be displayed in BP or BC/AD. Default is  \code{'BP'}.
#' @param ... Additional arguments affecting the plot
#' @return None.
#' @method plot spdppc
#' @export  

plot.spdppc = function(x, type='envelope', nsample=NULL, interval=0.90, obs.lwd=1.5,obs.col='black',sim.col='lightgrey',alpha=1,envelope.col='lightgrey',positive.col='red',negative.col='blue',calendar='BP', xlab=NULL, ylab=NULL, ...)
{
  if (!type%in%c('spaghetti','envelope')) {stop("The argument 'type' should be either 'spaghetti' or 'envelope'.")}
  if (is.null(nsample)) {nsample = ncol(x$simmatrix)}
  if (type=='spaghetti' & nsample > ncol(x$simmatrix))
  {
    warning(paste0('nsample large than the number of posterior simulations. Running with nsample=',ncol(x$simmatrix)))
    nsample = ncol(x$simmatrix)
  }
  #Setting y Label
  ylabel  <- ifelse(is.null(ylab),"Probability",ylab)
  #Setting calendar
  if (calendar=="BP"){
    plotyears <- x$obs$calBP
    xlabel <- ifelse(is.null(xlab),"Years cal BP",xlab)
#     xlabel <- "Years cal BP"
    xlim <- c(max(plotyears),min(plotyears))
  } else if (calendar=="BCAD"){
    plotyears <- BPtoBCAD(x$obs$calBP)
    xlabel <- ifelse(is.null(xlab),"Years BC/AD",xlab)
#     xlabel <- "Years BC/AD"
    if (all(range(plotyears)<0))
	    {
		    xlabel <- ifelse(is.null(xlab),"Years BC",xlab)
# 		    xlabel <- "Years BC"
	    }
    if (all(range(plotyears)>0))
	    {
		    xlabel <- ifelse(is.null(xlab),"Years AD",xlab)
# 		    xlabel <- "Years AD"
	    }
    xlim <- c(min(plotyears),max(plotyears))
  } else {
    stop("Unknown calendar type")
  }
  
  
  
  if (type=='envelope')
  {
    lo = apply(x$simmatrix,1,quantile,prob=(1-interval)/2)
    hi = apply(x$simmatrix,1,quantile,prob=1-(1-interval)/2)
    obs = x$obs$PrDens
    ylim = c(0,max(hi,x$obs$PrDens))
    
    # Boom and Bust Handling ####
    booms <- which(obs>hi)
    busts <- which(obs<lo)
    baseline <- rep(NA,length(obs))
    colpts = rep('grey',length(obs))
    colpts[booms] = 'red'
    colpts[busts] = 'blue'
    
    boomPlot <- baseline
    if (length(booms)>0){ boomPlot[booms]=obs[booms] }
    bustPlot <- baseline
    if (length(busts)>0){ bustPlot[busts]=obs[busts] }           
    
    boomBlocks <- vector("list")
    counter <- 0
    state <- "off"
    for (i in 1:length(boomPlot)){
      if (!is.na(boomPlot[i])&state=="off"){
        counter <- counter+1
        boomBlocks <- c(boomBlocks,vector("list",1))
        boomBlocks[[counter]] <- vector("list",2)
        boomBlocks[[counter]][[1]] <- boomPlot[i]
        boomBlocks[[counter]][[2]] <- plotyears[i]
        state <- "on"
      }
      if (state=="on"){
        if (!is.na(boomPlot[i])){
          boomBlocks[[counter]][[1]] <- c(boomBlocks[[counter]][[1]],boomPlot[i])
          boomBlocks[[counter]][[2]] <- c(boomBlocks[[counter]][[2]],plotyears[i])
        }
        if (is.na(boomPlot[i])){
          state <- "off"
        }
      }   
    }
    
    bustBlocks <- vector("list")
    counter <- 0
    state <- "off"
    for (i in 1:length(bustPlot)){
      if (!is.na(bustPlot[i])&state=="off"){
        counter <- counter+1
        bustBlocks <- c(bustBlocks,vector("list",1))
        bustBlocks[[counter]] <- vector("list",2)
        bustBlocks[[counter]][[1]] <- bustPlot[i]
        bustBlocks[[counter]][[2]] <- plotyears[i]
        state <- "on"
      }
      if (state=="on"){
        if (!is.na(bustPlot[i])){
          bustBlocks[[counter]][[1]] <- c(bustBlocks[[counter]][[1]],bustPlot[i])
          bustBlocks[[counter]][[2]] <- c(bustBlocks[[counter]][[2]],plotyears[i])
        }
        if (is.na(bustPlot[i])){
          state <- "off"
        }
      }   
    }
    
    plot(0, 0, xlim=xlim, ylim=ylim, type="n", col="white", ylab=ylabel, xlab=xlabel, xaxt="n", ...)
    
    polygon(c(plotyears,rev(plotyears)),c(lo,rev(hi)),col=envelope.col,border=NA)
    
    if (length(booms)>0){
      for (i in 1:length(boomBlocks)){
        bbb = unique(boomBlocks[[i]][[2]])
        index = which(plotyears%in%bbb)
        polygon(c(bbb,rev(bbb)),c(x$obs$PrDens[index],rev(hi[index])),border=NA,col=positive.col)
      }  
    }
    
    if (length(busts)>0){
      for (i in 1:length(bustBlocks)){
        bbb = unique(bustBlocks[[i]][[2]])
        index = which(plotyears%in%bbb)
        polygon(c(bbb,rev(bbb)),c(x$obs$PrDens[index],rev(lo[index])),border=NA,col=negative.col)
      }  
    }
    lines(plotyears,x$obs$PrDens,lwd=obs.lwd,col=obs.col)
  } 
  
  if (type=='spaghetti')
  {
    simmat = x$simmatrix[,sample(1:ncol(x$simmatrix),size=nsample)]
    ylim=c(0,max(x$obs$PrDens,simmat))
    plot(0, 0, xlim=xlim, ylim=ylim, type="n", col="white", ylab='Probability', xlab=xlabel, xaxt="n", ...)
    sim.col=col2rgb(sim.col)
    sim.col = rgb(sim.col[1,]/255,sim.col[2,]/255,sim.col[3,]/255,alpha=alpha)
    apply(simmat,2,lines,x=plotyears,col=sim.col)
    lines(plotyears,x$obs$PrDens,lwd=obs.lwd,col=obs.col)
  }

  if (calendar=="BP"){
    rr <- range(pretty(plotyears))    
    axis(side=1,at=seq(rr[2],rr[1],-100),labels=NA,tck = -.01)
    axis(side=1,at=pretty(plotyears),labels=abs(pretty(plotyears)))
  } else if (calendar=="BCAD"){
    yy <-  plotyears
    rr <- range(pretty(yy))    
    prettyTicks <- seq(rr[1],rr[2],+100)
    prettyTicks[which(prettyTicks>=0)] <-  prettyTicks[which(prettyTicks>=0)]-1
    axis(side=1,at=prettyTicks, labels=NA,tck = -.01)
    py <- pretty(yy)
    pyShown <- py
    if (any(pyShown==0)){pyShown[which(pyShown==0)]=1}
    py[which(py>1)] <-  py[which(py>1)]-1
    axis(side=1,at=py,labels=abs(pyShown))
  }
}
