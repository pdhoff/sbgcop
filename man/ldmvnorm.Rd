% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ldmvnorm.R
\name{ldmvnorm}
\alias{ldmvnorm}
\title{Log Multivariate Normal Density}
\usage{
ldmvnorm(Y, S)
}
\arguments{
\item{Y}{an n x p matrix}

\item{S}{a p x p positive definite matrix}
}
\value{
A real number.
}
\description{
Computes the log of the multivariate normal density
}
\details{
This function computes the log density of the data matrix \code{Y} under the
model that the rows are independent samples from a mean-zero multivariate
normal distribution with covariance matrix \code{S}.
}
\examples{

Y<-matrix(rnorm(9*7),9,7) 
ldmvnorm(Y,diag(7))



}
\author{
Peter Hoff
}
\keyword{distribution}
\keyword{multivariate}
