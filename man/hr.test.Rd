\name{hr.test}

\alias{hr.test}

\title{Hansen-Racine Bootstrap Model Averaged Unit Root Test}

\description{Computes a nonparametric bootstrap model average test for the null hypothesis that x has a unit root against stationary, explosive, or both alternatives.}

\usage{
hr.test(x = NULL, 
        adf.type = c("c","ct","nc","all"),
        alpha = 0.05, 
        alternative = c("both","stationary","explosive"),
        B = 399, 
        boot.method = c("geom","fixed","iid"),
        df.type = c("nc","nccct","ncc","nct","cct","c","ct","none"),
        group.start = 4,
        group.by = 4,        
        lag.vec = NULL, 
        method = c("jma","mma"),
        quantile.vec = c(0.005,0.01,0.025,0.05,0.95,0.975,0.99,0.995),
        random.seed = 42, 
        S = 12,
        verbose = TRUE)
}

\arguments{
  \item{x}{a univariate time series (\code{numeric} will be converted automatically)}
  \item{adf.type}{whether to include in the estimating equations Augmented Dickey-Fuller (ADF) models with no constant (nc), constant (c) or constant and trend (ct), or to include all three}
  \item{alpha}{size of the test}
  \item{alternative}{alternative hypothesis (null hypothesis is that x has a unit root)}
  \item{B}{number of bootstrap replications}
  \item{boot.method}{whether to use \code{sim="geom"} or \code{sim="fixed"} in the call to \code{\link[boot]{tsboot}} (\code{sim="iid"} generates a naive bootstrap)}
  \item{df.type}{whether to include in the estimating equations a Dickey-Fuller (DF) model with no constant (nc), constant(c), or constant and a trend (ct), two DF models with no constant and constant (ncc) or no constant and constant and trend (nct) or constant and constant and trend (cct), or three DF models with no constant and constant and constant and trend (nccct), or to include no DF models (none)}
  \item{group.start}{lag for first starting group for Hansen's (2014) grouping of covariates}
  \item{group.by}{group size for Hansen's (2014) grouping of covariates}
  \item{lag.vec}{optional vector of lags for the ADF estimating equations (default uses a modification of Schwert's (1989) ad hoc rule)}
  \item{method}{whether to use Hansen's (2007) Mallows model averaging criterion or Hansen and Racine's (2012) jackknife model averaging criterion}
  \item{quantile.vec}{vector of probabilities for the null quantile values returned}
  \item{random.seed}{random seed for the procedure}
  \item{S}{the constant in Schwert's (1989) ad-hoc formula for maximum lag selection}
  \item{verbose}{whether or not to display progress}
}

\details{

Typical usages are
\preformatted{
    hr.test(x)
    hr.test(x,B=999)
    hr.test(x,lag.vec=c(4,8,12))
    hr.test(x,S=2)
    hr.test(x,method="mma")
    hr.test(x,df.type="none",adf.type="nc")}

The default involves taking a set of DF and ADF estimating equations (ADF contain k > 0 lags of the differenced series, DF contain 0 lags) with covariates i) the first lag of x, ii) a constant and the first lag of x, and iii) a constant, the first lag of x and 1,4,8,...,k lags of the differenced x yielding k/4 such models (Hansen (2014)). Next, take the unit root test statistics from each of the estimating equations (in each estimating equation the first difference of x is regressed on the covariates). Finally, average these statistics using either a Mallows criterion (Hansen (2007)) or a jackknife criterion (Hansen and Racine (2012)). 

The null distribution of the model averaged test statistic is constructed using a model-free nonparametric bootstrap applied to the first difference of x similar to that outlined in Swensen (2003). A stationary geometric bootstrap (Politis and Romano (1994)) with automatic expected block length selection (Patton, Politis and White (2009)) is employed, and the null is imposed by taking the cumulative sum of the bootstrap values initialized to the first value of x; see \code{vignette("bmmaur", package = "hr")} for further details.

By default, Schwert's (1989) ad-hoc rule is used to determine the maximum lag k to be used (\code{round(S*(T/100)^0.25)}) where T is the length of x). By way of illustration, for T=100 there would be 10 estimating equations using the default settings (lags 1-8 for the augmented equations (S=8, T=100), no lag and no constant/no lag and constant (ncc)). Simulations indicate that the choice of constant S=12, though reasonable, involves more candidate models that necessary, hence the default S=8.

Simulations reveal that the test's size and power is balanced regardless of how large the number of estimating equations may be, unlike its peers, particularly in cases known to confound existing tests. The test displays the lowest size distortions in such cases among the range of tests considered (Phillips-Perron, Ng-Perron with modified AIC model selection, ADF with BIC model selection). Of these three peers, the Ng-Perron test exhibits the smallest size distortions in cases known to confound existing tests. Not only does our model average test have the lowest size distortions in such cases among all tests considered, but it also exhibits higher power than Ng-Perron's approach. The Phillips-Perron and ADF tests can fail to distinguish stationary from non-stationary cases even for large samples of data (see the example below) and therefore cannot be recommended for general application. To the best of our knowledge, this is the first such unit root test to be proposed that is based on model averaging.

We use the function \code{\link[np]{b.star}} from the \pkg{np} package for the automatic dependent bootstrap procedure, \code{\link[quadprog]{solve.QP}} from the \pkg{quadprog} package to solve the quadratic program for the model average weights), \code{\link[boot]{tsboot}} from the \pkg{boot} package to generate the stationary geometric bootstrap resample, and the \code{\link[fUnitRoots]{adfTest}} function from the \pkg{fUnitRoots} package to generate the statistics that are averaged.

}

\value{

A list with class "hrtest" containing the following components:

\item{tau}{model averaged test statistic}
\item{tau.alpha.low}{lower critical value}
\item{tau.alpha.up}{upper critical value}
\item{decision}{outcome of test}
\item{reject}{1 if reject, 0 otherwise}
\item{quantiles}{vector of nonparametric null quantiles}
\item{alpha}{size of test}
\item{ma.weights}{vector of model average weights}
\item{tau.boot}{sorted vector of bootstrapped test statistics generated under the null of a unit root}
\item{e.block.length}{expected block length for the stationary geometric bootstrap returned by \code{\link[np]{b.star}}}
\item{adf.lags}{vector of lags used for the augmented models}
}

\references{
Dickey, D. A. & Fuller, W. A. (1979), \dQuote{Distribution of the estimators for autoregressive time series with a unit root,} Journal of the American Statistical Association 74(366), 427-431.

Hansen, B. E. (2007), \dQuote{Least squares model averaging,} Econometrica 75, 1175-1189.

Hansen, B. E. (2014), \dQuote{Model averaging, asymptotic risk, and regressor groups,} Quantitative Economics 5(3), 495-530.

Hansen, B. E. & Racine, J. S. (2012), \dQuote{Jackknife model averaging,} Journal of Econometrics 167(1), 38-46.

Ng, S. & Perron, P. (2001), \dQuote{Lag length selection and the construction of unit root tests with good size and power,} Econometrica 69(6), 1519-1554.

Palm, F. C., Smeekes, S. & Urbain, J. P. (2008), \dQuote{Bootstrap unit-root tests: Comparison and extensions,} Journal of Time Series Analysis 29(2), 371-401.

Patton, A. & D. N. Politis & H. White (2009), \dQuote{CORRECTION TO "Automatic block-length selection for the dependent bootstrap" by D. Politis and H. White}, Econometric Reviews 28(4), 372-375.

Phillips, P. C. B. & Perron, P. (1988), \dQuote{Testing for a unit root in time series regression,} Biometrika 75(2), 335-346.

Politis, D. N. & J. P. Romano (1994), \dQuote{Limit theorems for weakly dependent Hilbert space valued random variables with applications to the stationary bootstrap}, Statistica Sinica 4, 461-476.

Politis, D. N. & H. White (2004), \dQuote{Automatic block-length selection for the dependent bootstrap}, Econometric Reviews 23(1), 53-70.

Racine, J. S. (2017), \dQuote{Bootstrap Model Averaging Unit Root Inference,} Technical Report.

Said S. E., Dickey D. A. (1984), \dQuote{Testing for Unit Roots in Autoregressive-Moving Average Models of Unknown Order,} Biometrika 71, 599-607.

Schwert, G. W. (1989), \dQuote{Tests for unit roots: A Monte Carlo investigation,} Journal of Business & Economic Statistics 7(2), 147-159.

Swensen, A. R. (2003), \dQuote{Bootstrapping unit root tests for integrated processes,} Journal of Time Series Analysis 24(1), 99-126.
}

\author{
Jeffrey S. Racine
}

\note{
This package is in beta status until further notice - proceed accordingly.
}

\seealso{
\code{\link[tseries]{pp.test}},  \code{\link[stats]{PP.test}}, \code{\link[tseries]{adf.test}}, \code{\link[fUnitRoots]{adfTest}}, \code{\link[fUnitRoots]{unitrootTest}}, \code{\link[urca]{ur.df}}, \code{\link[CADFtest]{CADFtest}}
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

\keyword{ nonparametric }
