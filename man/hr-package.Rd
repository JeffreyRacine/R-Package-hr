\name{hr-package}
\alias{hr-package}
\alias{hr}
\docType{package}
\title{
\packageTitle{hr}
}
\description{
\packageDescription{hr}
}
\details{
The DESCRIPTION file:
\packageDESCRIPTION{hr}
\packageIndices{hr}
\code{vignette("bmmaur", package = "hr")}
}
\author{
\packageAuthor{hr}

Maintainer: \packageMaintainer{hr}
}
\references{
Racine (2017), \dQuote{Bootstrap Mallows Model Averaging Unit Root Inference,} Technical Report.
}

\keyword{ package }
\seealso{
 \code{\link[tseries]{pp.test}}
}
\examples{
library(tseries)

## Generate a series with a unit root, MA errors

set.seed(42)
x <- arima.sim(list(order=c(0,1,1),ma=-0.8),n=100)

## Compute the test

hr.test(x,B=99)

## The Phillips-Perron test rejects 100% of the time for this DGP

pp.test(x)
}
