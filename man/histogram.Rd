% This is a template R document file
% Distributed under GPL 3 or later

\name{histogram}
\alias{histogram}
\alias{print.hist}
\alias{plot.hist}
\alias{lines.hist}

\title{Draw histogram based on data with rounding errors}
\description{
Draw histogram based on data with rounding errors.
}
\usage{
  histogram(x, w, nclass, binwidth, lb,
      range.x, freq=FALSE, truncate = TRUE)
}
\arguments{
  \item{x}{a numeric value or a vector.}
  \item{w}{A vector of weights for the values of \code{x}}
  \item{nclass}{Number of classes}
  \item{binwidth}{bin width.}
  \item{lb}{Lower bound of the first bin.}
  \item{range.x}{The range of \code{x} where 
    the histogram will be constructed.}
  \item{freq}{Whether \code{w} (if given) are frequencies.}
  \item{truncate}{Whether or not to truncate the data fall
    beyond the \code{range.x}.}
}
\details{
  If \code{nclass} is missing and the data are not prebinned, the 
  cross validation score based binwidth selector will be used
  by default.  If the data is pre-binned, the histogram will be 
  constructed by using the relative frequencies in \code{weights},
  or constructed using \code{table} if \code{weights} is missing.
}
\value{
  A class \code{hist} is returned, where a component \code{plot}
  is an R object "histogram" for visualization.
}
\references{
  Wang, B. and Wertelecki, W. (2012) 
  Density Estimation for Data With Rounding Errors.
  \emph{Computational Statistics and Data Analysis}, 
  (in press), doi: 10.1016/j.csda.2012.02.016. 
}

\author{
  B. Wang  \email{bwang@southalabama.edu}
}
\seealso{
 \code{\link{weighting}}.
}
\examples{

 ofc = round(rnorm(1000,34.5,1.5))
 (out = histogram(ofc,lb=min(ofc)-0.5, binwidth=1))
 plot(out)
 x0 = seq(28,42, length=100)
 lines(dnorm(x0,34.5,1.5)~x0, col=2)

}
\keyword{smooth}

