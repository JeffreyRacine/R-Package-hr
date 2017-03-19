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
## Generate a series with a unit root, MA errors

set.seed(42)
T <- 100
x <- arima.sim(list(order=c(0,1,1),ma=-0.8),n=T)

## Compute the test

hr.test(x,B=99)

## The Phillips-Perron test rejects 100% of the time for this DGP

suppressPackageStartupMessages(library("tseries"))

pp.test(x)

## The ADF test with BIC model selection incorrectly rejects the null a large 
## fraction of the time for this DGP

suppressPackageStartupMessages(library("CADFtest"))

CADFtest(x,max.lag.y=round(12*(T/100)^0.25),criterion="BIC",type="drift")
}
