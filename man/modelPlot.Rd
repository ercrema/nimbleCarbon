% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modelPlot.R
\name{modelPlot}
\alias{modelPlot}
\title{Plot Growth Models}
\usage{
modelPlot(
  model,
  a,
  b,
  params,
  type = c("spaghetti"),
  nsample = NULL,
  interval = 0.9,
  calendar = "BP",
  col = "lightgrey",
  alpha = 0.1,
  ylim = NULL,
  xlim = NULL,
  xlab = NULL,
  ylab = NULL,
  add = FALSE,
  lwd = 1,
  ...
)
}
\arguments{
\item{model}{growth model.}

\item{a}{lower (earliest) limit of the distribution (in BP).}

\item{b}{upper (latest) limit of the distribution (in BP).}

\item{params}{a \link{list} of vectors containing model parameters. The names attribute of each vector should match growth model parameters.}

\item{type}{either a 'spaghetti' plot or a quantile based 'envelope' plot. Default is 'spaghetti'.}

\item{nsample}{number of samples to be used. Default is the length of the parameter vectors supplied in the argument \code{params}.}

\item{interval}{quantile interval used for the envelope plot. Ignored when type is set to 'spaghetti'.}

\item{calendar}{either \code{'BP'} or \code{'BCAD'}. Indicate whether the calibrated date should be displayed in BP or BC/AD. Default is  \code{'BP'}.}

\item{col}{fill color for the quantile envelope (when \code{type=='envelope'}) or line colour (when \code{type=='spaghetti'}).}

\item{alpha}{transparency value for each line in the spaghetti plot or the fill color in the 'envelope' plot.  Default is 1.}

\item{ylim}{the y limits of the plot.}

\item{xlim}{the x limits of the plot (in Cal BP).}

\item{xlab}{a label for the x axis. Default is 'Years cal BP','Years BC/AD','Years BC', or 'Years AD' depending on data range and settings of \code{calendar}.}

\item{ylab}{a label for the y axis. Default is 'Probability'.}

\item{add}{whether or not the new graphic should be added to an existing plot.}

\item{lwd}{line width. Default is 1.}

\item{...}{additional arguments affecting the plot}
}
\value{
None.
}
\description{
Plots growth models based on user provided parameters for prior and posterior predictive checks.
}
\examples{
\donttest{
params = list(k=runif(100,0.01,0.02),r=runif(100,0.003,0.004))
modelPlot(model=dLogisticGrowth,a=5000,b=2000,params=params,type=c('spaghetti'),alpha=0.5)
}
}
