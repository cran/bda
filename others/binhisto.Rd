% This is a template R document file
% Distributed under GPL 3 or later

\name{binhist}
\alias{binhist}
\title{Draw histogram based on binned data}
\description{
Draw histogram based on binned data.
}
\usage{
binhist(x, just="center",scale=1.0)
}
\arguments{
  \item{x}{a numeric value or a vector.}
  \item{just}{How the data were binned: 'center' (default), 'left' or 'right'.}
  \item{scale}{Data can be rescaled.  Default: 1 (no rescaling).}
}
\details{
	Draw histogram using the binned data.  The data will not be further binned.
}
\value{
If 'scale'=100, the data will be rescaled by dividing all values by 100.  For exam the birth weight are usually rounded to 100mg.
}

\author{
  B. Wang  \email{bwang@jaguar1.usouthal.edu}
}
\seealso{
 \code{\link{edf}}
}
\examples{

data(birth)
binhist(birth$Weight,scale=100)
lines(smkde(birth$Weight,scale=100),col=4)
lines(histospline(birth$Weight,scale=100),col=3)
out = binem(birth$Weight,scale=100)
lines(out,col=2)

}
\keyword{stats}

