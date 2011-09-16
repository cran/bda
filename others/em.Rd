% This is a template R document file
% Distributed under GPL 3 or later

\name{binem}
\alias{binem}
\alias{plot.em}
\alias{print.em}
\title{An E-M algorithm to fit normal mixture based on binned data.}
\description{
	To fit a normal mixture based on data that are already binned.
}
\usage{
  binem(x,f,m=2,mu=NULL,na.rm=TRUE,just="center",binned=FALSE,
  scale=1.0, iter=100, modsel="AIC")
}
\arguments{
  \item{x}{vector of binned data; or can be a data frame of the binned data.}
  \item{f}{vector of frequencies if the data is not binned.}
  \item{m}{number of mixture components, default 2.}
  \item{mu}{A vector to sepcify the center of different component.  If 'mu' is specified, 'm' will be ignored.}	  
  \item{na.rm}{Option to deal with missing value.  Default: TRUE -- to remove missing values.} 
  \item{just}{Justify the center of the bins for data that are rounded up ('right'), rounded down ('left'), or rounded to the nearest integers ('center').}
  \item{binned}{whether x is binned data or not.  Default: FALSE.}
  \item{scale}{A positive value showing the scale of the measures in 'x'.  Default 1 (data are rounded to the nearest integers). If the data are rounded to 100mg, scale=100. }
  \item{iter}{Iteration number in selecting initial values.  Default 100.}
  \item{modsel}{Specify model selection criteria.  Options: 'AIC', 'AICc' or 'BIC'.}
}
\author{
  B. Wang  \email{bwang@jaguar1.usouthal.edu}
}
\seealso{
\code{\link{dmixnorm}}, \code{\link{qmixnorm}},\code{\link{rmixnorm}}.
}
\examples{
x = c(rnorm(100,34.5,1.5), rnorm(200,40,1.39))
y = round(x,0)
(out = binem(y))
p = out$Para$Prop
mu = out$Para$Mean
s = out$Para$Std.Dev
x0 = seq(min(y),max(y),length=100)
plot(dmixnorm(x0,p,mu,s)~x0,type='l')
(out2 = smkde(y))
lines(out2,col=2)
lines(density(y),col=3)
lines(density(x),col=4)
p0 = c(1/3,2/3)
mu0 = c(34.5,40)
s0 = c(1.5,1.39)
lines(dmixnorm(x0,p0,mu0,s0)~x0,lwd=2,col=2)



}
\keyword{stats}

