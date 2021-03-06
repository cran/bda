% This is a template R document file
% Distributed under GPL 3 or later

\name{ZipfPlot}
\alias{ZipfPlot}
\alias{ZipfPlot.default}
\alias{ZipfPlot.FSD}
\alias{ZipfPlot.bdata}
\alias{ZipfPlot.histogram}

\title{Draw Zipf Plot}

\description{ 
  Draw Zipf Plot.
}

\usage{
  ZipfPlot(x, x0, plot=FALSE,plot.new=TRUE, weights,...)
}
\arguments{
  \item{x}{data: two vectors.}
  \item{x0}{low bound to filter data.}
  \item{plot}{Draw Zipf plot if \code{TRUE}}
  \item{plot.new}{whether draw a new plot.}
  \item{weights}{Compute weighted least squares line
  if \code{weights} is given.}
  \item{...}{plotting parameters.}
}

\value{None}

\references{
Wang, B. (2020)
A Zipf-plot based normalization method for high-throughput RNA-Seq data.
PLoS ONE, (in press).
}

\examples{
data(LCL)
names(LCL)
x <- LCL$p47
y <- LCL$p107
outx <- ZipfPlot(x)
plot(outx,type='l')
outy <- ZipfPlot(y)
lines(outy,col=2)

out2 <- Zipf.Normalize(x,y)
outy2 <- ZipfPlot(out2$y)
lines(outy2,col=4)

  }
\keyword{stats}

