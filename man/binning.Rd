% This is a template R document file
% Distributed under GPL 3 or later

\name{binning}
\alias{binning}
\alias{binning.default}
\alias{binning.histogram}
\alias{plot.bdata}
\alias{print.bdata}
\alias{quantile.bdata}
%\alias{print.bdata2}
%\alias{percentile}

\encoding{UTF-8}

\title{Data Binning}

\description{
  To bin a univariate data set in to a consecutive bins.
}

\usage{
  binning(x, counts, nclass, breaks, bw)
}

\arguments{
  \item{x}{A vector of raw data. 'NA' values will be automatically 
  removed.}
  
  \item{counts}{Frequencies or counts of observations in different 
  classes (bins)}
  
  \item{nclass}{Number of classes}

  \item{breaks}{The break points for data binning.}
  
  \item{bw}{Bin (class) width.}

}

\details{
  
  Data binning can be done in the following cases:  
  (1) both \code{x} and \code{breaks} are specified;  
  (2) both \code{counts} and \code{breaks} are specified;
  (3) both \code{x} and \code{bw} are specified.  
      If (3a) \code{breaks} is of length 1, it gives the lower limit 
      of the first bin/class, and a sequence of equal-width bins will 
      be defined. 
      If (3b) \code{breaks} is missing, a random starting point will 
      be chosen.
 
 \code{x} and \code{counts} cannot be both given.

 \code{nclass} will not be used if \code{x} is missing. 

}

\examples{

  y <- c(10, 21, 56,79,114,122,110,85,85,61,47,49,47,44,31,20,11,4,4)
  x <- 14.5 + c(0:length(y))
  out1 <- binning(counts=y, breaks=x)
  plot(out1)

}
\keyword{distribution}
\keyword{stats}

