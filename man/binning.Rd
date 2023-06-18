% This is a template R document file
% Distributed under GPL 3 or later

\name{binning}
\alias{binning}
\alias{binning.default}
\alias{binning.histogram}
\alias{plot.bdata}
\alias{print.bdata}
\alias{.bdaConnect}

\encoding{UTF-8}

\title{Data Binning}

\description{
  To bin a univariate data set in to a consecutive bins.
}

\usage{
  binning(x, counts, breaks,lower.limit, upper.limit)
}

\arguments{
  \item{x}{A vector of raw data. 'NA' values will be automatically 
  removed.}
  
  \item{counts}{Frequencies or counts of observations in different 
  classes (bins)}
  
  \item{breaks}{The break points for data binning.}
  
  \item{lower.limit, upper.limit}{The lower and upper limits of the bins.}

}

\details{
  
To create a 'bdata' object. If 'x' is given, a histogram will be
created. Otherwise, create a histogram-type data using 'counts' and 
'breaks' (or class limits with 'lower.limit' and/or 'upper.limit').
}

\value{
 \item{ll}{lower limits}
 \item{ul}{upper limits}
 \item{freq}{frequencies}
 \item{xhist}{histogram}
 \item{xZipf}{Zipf plot}
}

\examples{
y <- c(10, 21, 56,79,114,122,110,85,85,61,47,49,47,44,31,20,11,4,4)
x <- 14.5 + c(0:length(y))
out1 <- binning(counts=y, breaks=x)
plot(out1)

z = rnorm(100, 34.5, 1.6)
out1 <- binning(z)
plot(out1)

data(FSD)

x <- as.numeric(FirmAge[38,]);
age <- c(0,1:6,11,16,21,26,38);
y <- binning(counts=x, lower.limit=age)
plot(y)
plot(y, type="Zipf")

x <- as.numeric(FirmSize[38,]);
names(FirmSize)
ll <- c(1,5,10,20,50,100,250,500,1000,2500,5000,10000);
ul <- c(4,9,19,49,99,249,499,999,2499,4999,9999,Inf)
y <- binning(counts=x, lower.limit=ll,upper.limit=ul)
plot(y)
plot(y, type="Zipf")
}

\keyword{distribution}
\keyword{stats}

