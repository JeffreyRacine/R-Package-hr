\name{hr.test}

\alias{hr.test}

\title{Bootstrap Model Averaged Unit Root Test}

\description{Computes a nonparametric bootstrap model average test for the null hypothesis that x has a unit root against stationary and explosive alternatives.}

\usage{
hr.test(x = NULL, 
        K.vec = NULL, 
        random.seed = 42, 
        B = 399, 
        alpha = 0.05, 
        trend = TRUE, 
        verbose = TRUE, 
        q = c(0.005,0.01,0.025,0.05,0.95,0.975,0.99,0.995))
}

\arguments{
  \item{x}{a univariate time series (\code{numeric} will be converted automatically)}
  \item{K.vec}{optional vector of lags for the (augmented) Dickey-Fuller estimating equations}
  \item{random.seed}{random seed for the procedure}
  \item{B}{number of bootstrap replications}
  \item{alpha}{size of the test}
  \item{trend}{whether or not to include a deterministic trend in the estimating equations}
  \item{verbose}{whether or not to display progress}
  \item{q}{vector of probabilities for the null quantile values returned}
}

\details{

Typical usages are
\preformatted{
    hr.test(x)
    hr.test(x,B=99)
    hr.test(x,K.vec=1:4)
}

Take a set of (augmented) Dickey-Fuller estimating equations with covariates i) the first lag of x, ii) a constant and the first lag of x, iii) a constant, a deterministic trend and the first lag of x (when trend=TRUE) and iv) a constant, a deterministic trend, the first lag of x and 1,2,...,k lags of the differenced x yielding k such models (when trend=TRUE, otherwise no deterministic trend is included). Next, take the unit root test statistics from each of the estimating equations (in each estimating equation the first difference of x is regressed on the covariates). Finally, average these statistics using a Mallows criterion  (Hansen (2007)). 

The null distribution of the test statistic is constructed using a model-free nonparametric bootstrap applied to the first difference of x similar to that outlined in Swensen (2003). A stationary geometric bootstrap (Politis and Romano (1994)) with automatic expected block length selection (Patton, Politis and White (2009)) is employed, and the null is imposed by taking the cumulative sum of the bootstrap values initialized to the first value of x; see \code{vignette("bmmaur", package = "hr")} for further details.

When \code{trend=TRUE} there are 3+k estimating equations, otherwise 2+k. By default, Schwert's (1989) ad-hoc rule is used to determine the maximum lag k to be used (\code{round(12*(T/100)^0.25)} where T is the length of x). By way of illustration, for T=100 there would be 15 estimating equations when a trend is included, 14 otherwise.

Simulations reveal that the test's size and power is balanced regardless of how large the number of estimating equations may be, unlike its peers, particularly in cases known to confound existing tests. To the best of our knowledge, this is the first such unit root test to be proposed that is based on model averaging.

We use the function \code{\link[np]{b.star}} from the \pkg{np} package for the automatic dependent bootstrap procedure, \code{\link[quadprog]{solve.QP}} from the \pkg{quadprog} package to solve the quadratic program for the Mallows model average weights (Hansen (2007)), \code{\link[boot]{tsboot}} from the \pkg{boot} package to generate the stationary geometric bootstrap resample, and the \code{\link[fUnitRoots]{adfTest}} function from the \pkg{fUnitRoots} package to generate the statistics that are averaged.

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
\item{e.block.length}{expected block length for the stationary geometric bootstrap returned by \code{\link[np]{b.star}}}
\item{adf.lags}{vector of lags used for the augmented models}
}

\references{

Dickey, D. A. & Fuller, W. A. (1979), \dQuote{Distribution of the estimators for autoregressive time series with a unit root,} Journal of the American Statistical Association 74(366), 427-431.

Hansen, B. E. (2007), \dQuote{Least squares model averaging,} Econometrica 75, 1175-1189.

Patton, A. and D.N. Politis and H. White (2009), \dQuote{CORRECTION TO "Automatic block-length selection for the dependent bootstrap" by D. Politis and H. White}, Econometric Reviews 28(4), 372-375.

Phillips, P. C. B. & Perron, P. (1988), \dQuote{Testing for a unit root in time series regression,} Biometrika 75(2), 335-346.

Politis, D.N. and J.P. Romano (1994), \dQuote{Limit theorems for weakly dependent Hilbert space valued random variables with applications to the stationary bootstrap}, Statistica Sinica 4, 461-476.

Politis, D.N. and H. White (2004), \dQuote{Automatic block-length selection for the dependent bootstrap}, Econometric Reviews 23(1), 53-70.

Racine, J. S. (2017), \dQuote{Bootstrap Mallows Model Averaging Unit Root Inference,} Technical Report.

Schwert, G. W. (1989), \dQuote{Tests for unit roots: A Monte Carlo investigation,} Journal of Business & Economic Statistics 7(2), 147-159.

Swensen, A. R. (2003), \dQuote{Bootstrapping unit root tests for integrated processes,} Journal of Time Series Analysis 24(1), 99-126.

}

\author{
Jeffrey S. Racine
}

\note{
For a draft of the methodology run \code{vignette("bmmaur", package = "hr")}. This function and related material is in beta mode until further notice.
}

\seealso{
\code{\link[tseries]{pp.test}},  \code{\link[stats]{PP.test}}, \code{\link[tseries]{adf.test}}, \code{\link[fUnitRoots]{adfTest}}, \code{\link[fUnitRoots]{unitrootTest}}, \code{\link[urca]{ur.df}}
}

\examples{
library(tseries)

## Generate a series that contains a unit root and has moving average errors
## with a large negative coefficient (this can confound existing tests)

set.seed(42)
x <- arima.sim(list(order=c(0,1,1),ma=-0.8),n=100)

## Compute the test

hr.test(x,B=99)

## The Phillips-Perron test incorrectly rejects the null 100% of the time for 
## this DGP

pp.test(x)
}

\keyword{ nonparametric }
