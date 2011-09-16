\name{binning}
\alias{binning}
\alias{print.bdata}
\alias{hist.bdata}
\alias{bdata}
\title{
  Data binning.
}
\description{
  To bin raw data and prepare data for finite mixture model fitting.
}
\usage{
  binning(x, scale=1, method='nearest')
}
\arguments{
  \item{x}{Sample data before binning.}
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
  Wang, B. and Wertelecki, W. (2011) Density estimation from data with rounding errors, \emph{Computational Statistics and Data Analysis}, (submitted)
}

\seealso{
\code{\link{bfmm}}
}

\examples{

data(birth)
(ofc = binning(birth$Head))
(bwt = binning(birth$Weight, scale=100))

}
\keyword{data}
\keyword{models}
