% This is a template R document file
% Distributed under GPL 3 or later

\name{mle.gamma}
\alias{mle.gamma}
\title{Compute the MLEs of a gamma distribution.}
\description{
     Compute the MLEs of a gamma distribution.
}
\usage{
  mle.gamma(x)
}
\arguments{
  \item{x}{vector of quantiles.}
}
\author{
  B. Wang  \email{bwang@jaguar1.usouthal.edu}
}
\seealso{
  \code{\link{smkde}},\code{\link{histospline}}.
}
\examples{

x = rgamma(100,10,2)
mle.gamma(x)

}
\keyword{stats}

