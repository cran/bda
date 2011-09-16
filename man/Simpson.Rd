% This is a template R document file
% Distributed under GPL 3 or later

\name{rsimpson}
\alias{rsimpson}
\alias{dsimpson}
\title{The simpson (claw) distribution.}
\description{
     Compute the density, and generate a random sample of a simpson (the claw) distribution.}
\usage{
 rsimpson(n)
 dsimpson(x)
}
\arguments{
  \item{n}{number of observations. If 'length(n) > 1', the length is
          taken to be the number required.} 
  \item{x}{vector of quantiles.}
}
\author{
  B. Wang  \email{bwang@jaguar1.usouthal.edu}
}
\seealso{
  \code{\link{smkde}},\code{\link{histospline}}.
}
\examples{

x0 = seq(-3,3,length=100)
x = rsimpson(3000)
y = round(x)
out = histo(y)
plot(out)
lines(smkde(y),col=3)
lines(histospline(y),col=4)
lines(density(x))
lines(dsimpson(x0),col=2)

}
\keyword{stats}

