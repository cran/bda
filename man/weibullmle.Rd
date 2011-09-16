% This is a template R document file
% Distributed under GPL 3 or later

\name{mle.weibull}
\alias{mle.weibull}
\title{Compute the MLEs of a weibull distribution.}
\description{
     Compute the MLEs of a weibull distribution.
}
\usage{
  mle.weibull(x)
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

x = rweibull(100,10,2)
mle.weibull(x)

}
\keyword{stats}

