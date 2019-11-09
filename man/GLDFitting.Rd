% This is a template R document file
% Distributed under GPL 3 or later

\name{fit.GLD.FMKL}
\alias{fit.GLD.FMKL}
\alias{fit.GLD}
\alias{print.FMKL}
\alias{plot.FMKL}
\alias{lines.FMKL}

\encoding{UTF-8}

\title{Fitting FMKL GLD}

\description{
 To fit a FMKL GLD to raw/binned data.
}

\usage{
 fit.GLD.FMKL(x, qtl, qtl.levels, lbound, ubound,
              percentile='exact', mle=FALSE)
 fit.GLD(x, lbound, ubound, method='chisquare')
}

\arguments{
  \item{x}{A vector of raw data, or a histogram or binned data.}
  \item{qtl,qtl.levels}{vectors. Five quantiles and quantile levels
  to fit GLD using method of percentiles. If missing, the exact or
  approximate quantiles will be determined depending on the data.}
  \item{percentile}{Use the exact percentiles (\code{exact}) 
  or approxiated values (\code{approximate}).}

  \item{mle}{Logical. To find the MLE or not.}

  \item{lbound,ubound}{lower and upper bound for the support of the
  density. The bounds could be finite values, or positive or negative
  infinity.}

  \item{method}{Method for goodness-of-fit test.}

}


\examples{
  data(hhi)
  hmob <- binning(counts=hhi$mob, breaks=hhi$breaks)
  lmd5 <- fit.GLD.FMKL(hmob)
  lmd6 <-  fit.GLD.FMKL(hmob, mle=TRUE)
  plot(lmd5)
  lines(lmd6, col=4)
  ## GOP example (handbook) -- Hahn & Sapiro (1967)
  ## KS-GLD based on original data: (0.0345, 0.00009604, 0.87, 4.92)
  ## Table 3.6-1
  breaks <- c(-Inf, seq(0.015, length=10, by=0.005), Inf)
  counts <- c(1,9,30,44,58,45,29,17,9,4,4)
  rho.mid <- c(0.0325, 0.0250, 0.667, 0.600)
  rho.unif <- c(0.03352, 0.02531, 0.7786, 0.5009)
  ## histogram for chi-square test
  ## KS = 0.0225, p-value = 0.999.  Chi=0.5176, p-value=0.7720
  breaks <- c(-Inf, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, Inf)
  counts <- c(40,44,58,45, 29,17,17)


}
\keyword{distribution}
\keyword{stats}

