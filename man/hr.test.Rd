\name{hr.test}

\alias{hr.test}

\title{Bootstrap Model Averaged Unit Root Test}

\description{This function takes a univariate time series object and performs a nonparametric bootstrap model average test for the null hypothesis that x has a unit root.}

\usage{
hr.test(x = NULL, 
        K.vec = NULL, 
        random.seed = 42, 
        B = 399, 
        alpha = 0.05, 
        trend = TRUE, 
        verbose = TRUE, 
        q = c(0.005, 0.01, 0.025, 0.05, 0.95, 0.975, 0.99, 0.995))
}

\arguments{
  \item{x}{a univariate time series (if not a time series it will be converted into one)}
  \item{K.vec}{optional vector of lags for the augmented Dickey-Fuller candidate models}
  \item{random.seed}{seed for random number generator}
  \item{B}{number of bootstrap replications}
  \item{alpha}{size of test}
  \item{trend}{whether to include a deterministic trend in the candidate models or not}
  \item{verbose}{whether to display progress or not}
  \item{q}{vector of probabilities for the quantiles returned}
}

\details{

Given a set of candidate estimating equation models that include i) a model with no constant and first lag of the series, ii) a model with a constant and first lag of the series, iii) a model with a constant, trend and first lag of the series (when trend=TRUE) and iv) models with constant, trend, first lag of the series and up to k lags of the differenced series (when trend=TRUE), the unit root test statistics obtained when regressing the difference of a series on the above covariates are averaged over the candidate models using a Mallows criterion. 

When \code{trend=TRUE} there are 3+k candidate models that we average over, otherwise there are 2+k candidate models. By default Schwert's ad-hoc rule is used to determine the maximum lag k to be used is adopted (\code{round(12*(n/100)^0.25)} where n is the length of the series). By way of illustration, for n=100 there would be 15 candidate models when a trend is included, 14 otherwise.

The null distribution of the model average unit root statistic is obtained via a model-free bootstrap constructed from residuals of the differenced series. The bootstrap applied to these residuals is a geometric bootstrap with automatic block length selection, and the null is imposed by taking the cumulative sum (\code{cumsum()}) of these residuals starting from the first value of the series. See \code{vignette("bmmaur", package = "hr")} for details.

Simulations reveal that the test's performance is quite balanced regardless of how large the number of candidate models may be, unlike its peers, and it also has more balanced size than its peers in cases known to confound existing tests. To the best of our knowledge, this is the first such unit root test based on model averaging.

We use the function \code{\link[np]{b.star}} from the \pkg{np} package for the automatic dependent bootstrap procedure, \code{\link[quadprog]{solve.QP}} from the \pkg{quadprog} package to solve the quadratic program for the Mallow model average weights, \code{\link[boot]{tsboot}} from the \pkg{boot} package to generate the geometric bootstrap series, and the \code{\link[fUnitRoots]{adfTest}} function from the \pkg{fUnitRoots} package to generate the statistics that are averaged.

}

\value{
\item{tau}{model averaged test statistic}
\item{tau.alpha.low}{lower critical value}
\item{tau.alpha.up}{upper critical value}
\item{decision}{outcome of test}
\item{reject}{1 if reject, 0 otherwise}
\item{quantiles}{vector of nonparametric null quantiles}
\item{alpha}{size of test}
\item{trend}{whether to include a deterministic trend or not}
\item{mma.weights}{vector of Mallows model average weights}
\item{tau.boot}{vector of bootstrapped test statistics generated under the null of a unit root}
\item{e.block.length}{expected block length for geometric bootstrap returned by \code{\link[np]{b.star}}}
\item{adf.lags}{vector of lags used for the augmented models}
}

\references{
Racine (2017), \dQuote{Bootstrap Mallows Model Averaging Unit Root Inference,} Technical Report.
}

\author{
Jeffrey S. Racine
}

\note{
For a preliminary technical report on the methodology run \code{vignette("bmmaur", package = "hr")}. This function is in beta mode until further notice.
}

\seealso{
\code{\link[tseries]{pp.test}},  \code{\link[stats]{PP.test}}, \code{\link[tseries]{adf.test}}, \code{\link[fUnitRoots]{adfTest}}, \code{\link[fUnitRoots]{unitrootTest}}, \code{\link[urca]{ur.df}}
}

\examples{
library(tseries)

## Generate a series that contains a unit root and has moving average errors

set.seed(42)
x <- arima.sim(list(order=c(0,1,1),ma=-0.8),n=100)

## Compute the test

hr.test(x,B=99)

## The Phillips-Perron test incorrectly rejects the null 100% of the time for 
## this DGP

pp.test(x)
}

\keyword{ nonparametric }
