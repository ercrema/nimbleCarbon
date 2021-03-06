% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/postPredCor.R
\name{postPredCor}
\alias{postPredCor}
\title{Calculates correlation between observed and posterior generated SPD.}
\usage{
postPredCor(x, method = "pearson")
}
\arguments{
\item{x}{An object of class \code{spdppc}.}

\item{method}{a character string indicating which correlation coefficient is to be computed. One of "pearson" (default), "kendall", or "spearman": can be abbreviated.}
}
\value{
A vector of correlation values.
}
\description{
Computes the correlation between observed SPDs and posterior generated SPD from the output of \code{postPredSPD()} function as an heuristic of the goodness-of-fit of the model.
}
