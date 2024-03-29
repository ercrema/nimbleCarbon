% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dAsymLaplace.R
\name{dAsymLaplace}
\alias{dAsymLaplace}
\alias{dAsymLaplace,}
\alias{rAsymLaplace,}
\alias{qAsymLaplace,}
\alias{pAsymLaplace}
\alias{rAsymLaplace}
\alias{qAsymLaplace}
\title{Asymmetric Laplace Distribution}
\usage{
dAsymLaplace(x, mu, sigma, tau, log)

rAsymLaplace(n, mu, sigma, tau)

pAsymLaplace(q, mu, sigma, tau, lower.tail = 1, log.p = 0)

qAsymLaplace(p, mu, sigma, tau, lower.tail = 1, log.p = 0)
}
\arguments{
\item{x}{value to be computed.}

\item{mu}{location parameter.}

\item{sigma}{scale parameter.}

\item{tau}{asymmetry parameter.}

\item{log, log.p}{TRUE or 1 to return log probability. FALSE or 0 to return probability.}

\item{n}{number of random draws. Currently only n = 1 is supported, but the argument exists for standardization of "r" functions.}

\item{q}{quantile to be computed.}

\item{lower.tail}{logical; if TRUE (default), probabilities are P[X <= x] otherwise, P[X > x].}

\item{p}{probability to be computed.}
}
\description{
Density, distribution function, quantile function, and random generation for a Asymmetric Laplace Distribution.
}
\author{
Enrico Crema
}
