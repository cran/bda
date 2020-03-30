% This is a template R document file
% Distributed under GPL 3 or later

\name{fit.GBP}
\alias{fit.GBP}
\alias{dGBP}
\alias{pGBP}

\title{Fitting Mixture Model of Generalized Beta and Pareto}
\description{
 To fit a mixture model of generalize beta and Pareto to grouped data.
}
\usage{
  fit.GBP(x,breaks)
}
\arguments{
  \item{x}{'x' can be a vector or a matrix, or a 'histogram'.}
  \item{breaks}{a matrix with two columns if 'x' is a matrix. Otherwise, it is a vector. Can be missing if 'x' is a vector and 'x' will be grouped uisng the default parameters with \code{hist}. }
}

\value{
 \item{pars}{estimated parameters.}
}

\examples{
data(FSD)
x <- as.numeric(FirmSize[nrow(FirmSize),])
brks.size <- c(0,4.5,9.5,19.5,49.5,99.5,249.5,499.5,
               999.5,2499.5,4999.5, 9999.5,Inf)
xhist1 <- binning(counts=x,breaks=brks.size)

(out <- fit.GBP(x,brks.size))
(out <- fit.GBP(xhist1))

plot(xhist1,xlim=c(0,110))
x0 <- seq(0,110,length=1000)
f0 <- dGBP(x0,out)
lines(f0~x0, col=2,lwd=2)

ZipfPlot(xhist1,plot=TRUE)
F0 <- pGBP(brks.size,out)
lines(log(1-F0)~log(brks.size), col=2)
}

\keyword{distribution}

