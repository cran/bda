% This is a template R document file
% Distributed under GPL 3 or later

\name{fit.GB}
\alias{fit.GB}
\alias{print.GB}
\alias{plot.GB}
\alias{lines.GB}

\encoding{UTF-8}

\title{Fitting GB}

\description{
 To fit a generalized beta distribution to raw/binned data.
}

\usage{
 fit.GB(x, lbound, ubound)
}

\arguments{
  \item{x}{A vector of raw data, or a histogram or binned data.}

  \item{lbound,ubound}{lower and upper bound for the support of the
  density. The bounds could be finite values, or positive or negative
  infinity.}

}


\examples{

  x <- rnorm(100, 34.5,1.5)
  hofc <- binning(x, bw=1)
  lmd7 <- fit.GB(hofc)
  plot(lmd7)
  gof(lmd7)


}
\keyword{distribution}
\keyword{stats}

