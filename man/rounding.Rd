\name{rounding}
\alias{rounding}
\alias{print.bdata}
\alias{hist.bdata}
\alias{bdata}
\title{
  Data rounding.
}
\description{
  To bin raw data and prepare data for finite mixture model fitting.
}
\usage{
  rounding(x, scale=1, method='nearest')
}
\arguments{
  \item{x}{Sample data before rounding.}
  \item{scale}{The scale of data.  Default: scale=1.}
  \item{method}{Method to round the data: up/down, or to the nearest integers.}
}

\details{ Raw data need to be prepared to have specific formats to be analyzed
 using the FMMBD functions.  In summary, data should have three columns: 
x = center of the bins, 
widths = width of the bins, and 
counts = counts (frequencies) of the bins.  

If data are prepared manually, default scale=1.  One may need to transform data first.

}

\references{
  Wang, B. and Wertelecki, W. (2012) 
  Density Estimation for Data With Rounding Errors.
  \emph{Computational Statistics and Data Analysis}, 
  (in press), doi: 10.1016/j.csda.2012.02.016. 
}

\seealso{
\code{\link{bfmm}}
}

\examples{

data(birth)
(ofc = rounding(birth$Head))
(bwt = rounding(birth$Weight, scale=100))

}
\keyword{data}
\keyword{models}
