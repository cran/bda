% This is a template R document file
% Distributed under GPL 3 or later

\name{histo}
\alias{histo}
\alias{print.hist}
\alias{plot.hist}
\alias{lines.hist}

\title{Draw histogram based on data with rounding errors}
\description{
Draw histogram based on data with rounding errors.
}
\usage{
histo(x, m=NULL, alpha=0.05, binned=FALSE,just="center",scale=1.0)
}
\arguments{
  \item{x}{a numeric value or a vector.}
  \item{m}{number of bins.}
  \item{alpha}{Significance level.  Default: 0.05.}
  \item{binned}{Have the data been binned.  Default: FALSE.}
  \item{just}{How the data were binned: 'center' (default), 'left' or 'right'.}
  \item{scale}{Data can be rescaled.  Default: 1 (no rescaling).}
}
\details{
The cross validation score will be return as 'Jh'. 
}
\value{
If 'scale'=100, the data will be rescaled by dividing all values by 100.  For exam the birth weight are usually rounded to 100mg.
}

\author{
  B. Wang  \email{bwang@jaguar1.usouthal.edu}
}
\seealso{
 \code{\link{edf}}, \code{\link{pcb}}.
}
\examples{

x = rnorm(10000)
m=2:100
Jh = NULL
for(i in m){
  out=histo(x,i)
  Jh=c(Jh,out$Jh)
}
plot(Jh~m,type='l')

j = m[which(Jh==min(Jh))]
(out1=histo(x,j))
plot(out1)
lines(dnorm(sort(x))~sort(x))

}
\keyword{stats}

