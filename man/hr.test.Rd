\name{hr.test}
\alias{hr.test}
\title{Hansen and Racine's Nonparametric Model Averaged Unit Root Test}
\description{This function takes a univariate time series object and performs Hansen and Racine's (2017) nonparametric bootstrap frequentist model average unit root test.}
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
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a univariate time series}
  \item{K.vec}{optional vector of lags for the augmented Dickey-Fuller candidate models}
  \item{random.seed}{seed for random number generator}
  \item{B}{number of bootstrap replications}
  \item{alpha}{size of test}
  \item{trend}{whether to include a deterministic trend in the candidate models or not}
  \item{verbose}{whether to display progress or not}
  \item{q}{vector of probabilities for the quantiles returned}
  }
\details{
We use the function \code{\link[np]{b.star}} from the \pkg{np} package for the automatic dependent bootstrap procedure, \code{\link[quadprog]{solve.QP}} from the \pkg{quadprog} package, \code{\link[boot]{tsboot}} from the \pkg{boot} package, and the \code{\link[fUnitRoots]{adfTest}} function from the \pkg{fUnitRoots} package. These packages must be installed.
}
\value{
\item{tau }{model averaged test statistic}
\item{decision }{outcome of test}
\item{reject }{1 if reject, 0 otherwise}
\item{quantiles }{vector of nonparametric null quantiles}
\item{alpha }{size of test}
\item{trend }{whether to include a deterministic trend or not}
\item{mma.weights }{vector of Mallows model average weights}
\item{tau.boot }{vector of bootstrapped test statistics generated under the null of a unit root}
\item{e.block.length }{expected block length for geometric bootstrap returned by \code{\link{b.star}}}
\item{adf.lags }{vector of lags used for the augmented models}
}
\references{
Hansen and Racine (2017)
}
\author{
Jeffrey S. Racine
}
\note{
This function is in beta mode until further notice.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Generate a series with a unit root
set.seed(42)
x <- ts(cumsum(rnorm(100)))
hr.test(x,B=99)
}

\keyword{ ~nonparametric }
