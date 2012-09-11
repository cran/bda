% This is a template R document file
% Distributed under GPL 3 or later

\name{smkde}
\alias{smkde}
\alias{print.bde}
\alias{plot.bde}
\title{Fit smoothed KDE to binned data.}
\description{
     Fit smoothed KDE to binned data via an iterative algorithm.}
\usage{
smkde(x,f,bw,gridsize=512L,na.rm=TRUE,just="center",binned=FALSE, scale=1.0)
}
\arguments{
  \item{x}{vector of raw data or distinct center information of binned data.}
  \item{f}{A vector of frequencies.  Can not be missing if 'binned=TRUE'.}
  \item{bw}{Smooth parameter.}
  \item{gridsize}{The size of grid points where the density 
  willbe evaluated.}
  \item{na.rm}{Default "TRUE".  Automatically remove missing values.  
  Otherwise, stop and warn.}
  \item{just}{Character showing how data are binned.  "center", "left"
   or "right".}
  \item{binned}{Logic value. Default value: 'FALSE'.}
  \item{scale}{Numeric.  Default=1.0.}
}
\author{
 Blower, Gordon and Kelsall, Julia, (2002). 
 Nonlinear kernel density estimation for binned data: convergence in entropy. 
 \emph{Bernoulli}, 
 8(4), 423-449.
}
\seealso{
 \code{\link{bfmm}}.
}
\examples{

p = c(.4,.6)
mu = c(34.5,40)
s=c(1.5,1.3)
x = rmixnorm(300,p,mu,s)
y = round(x,0)
out = smkde(y)

}
\keyword{smooth}

