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
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
\packageAuthor{hr}

Maintainer: \packageMaintainer{hr}
}
\references{
~~ Literature or other references for background information ~~
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
~~ the R documentation directory ~~
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
require(tseries)
## Generate a series with a unit root
set.seed(42)
x <- arima.sim(list(order=c(0,1,1),ma=-0.8),n=100)
hr.test(x,B=99)
## The Phillips-Perron test rejects 100% of the time for this DGP
pp.test(x)
}
