\name{weighting}
\alias{weighting}
\alias{weighting.default}
\encoding{UTF-8}
\title{
 To weight the data
}
\description{
  To create an R object by weighting the data with a vector of weights.  
}
\usage{
 weighting(x,w,freq=FALSE,na.rm=TRUE, type, method,...)
}
\arguments{
\item{x}{
  vector of observations from the distribution whose density is to
  be estimated.  Missing values are not allowed.
}
\item{w}{
 The weights of \code{x}.  The weight \code{w_i} of any 
 observation \code{x_i} should
 be non-negative.  If \code{x_i=0}, \code{x_i} will be removed
 from the analysis.  
 For survival data, \code{w} is a vector of censoring status: 
 0 for censored data and 1 for uncensored data.
}
\item{freq}{
  An indicator showing whether \code{w} is a vector of frequecies
  (counts) or weights.
}
\item{na.rm}{
logical flag: if \code{TRUE}, \code{NA} values will be ignored; 
otherwise, the program will be halted with error information.
}
\item{type}{Data type: \code{rc} for random right-censoring.  
If missing, weighted by \code{w}.
}
\item{method}{method to compute the weights.  Not needed if 
\code{type} is missing.}
\item{...}{controls}
}

\value{
a list containing the following components:

\item{x}{
vector of sorted \code{x} values at which the estimate was computed.
}
\item{y}{
vector of weights estimates
at the corresponding \code{x}.
}
\item{n}{
length of \code{x}
}
\item{size}{
sample size (of the original data).
}
\item{totMass}{
total mass.  needed for data with censoring.
}
\item{pars}{
Reserved.
}

}

\author{
  B. Wang  \email{bwang@southalabama.edu}
}

\references{
Wang, B. and Wang, X-F. (2007). 
"Bandwidth Selection for Weighted Kernel Density Estimation".
}
\seealso{
  \code{\link{lprde}}.
}
\examples{

 ## time to relapse of patients in a clinical trial of 6-MP 
 ## against a placebo (Survival Analysis, Klein and Moeschberger 
 ## (1950) 2nd edition). 

 x <- c(10,7,32,23,22,6,16,34,32,25,11,20,19,6,17,35,6,13,9,6,10)
 w <- c(1,1,0,1,1,1,1,0,0,0,0,0,0,1,0,0,1,1,0,0,0)
 (out1 <- weighting(x,w,type="rc",method="Nelson"))
 (out2 <- weighting(x,w,type="rc",method="Kaplan"))

}
\keyword{smooth}
