% This is a template R document file
% Distributed under GPL 3 or later

\name{gof.test}
\alias{gof}
\alias{gof.em}
\alias{gof.bde}
\alias{gof.test}
\alias{gof.default}
\title{To perform goodness-of-fit test.}
\description{
  To perform goodness-of-fit test.
}
\usage{
  gof.test(object, x,...)
}
\arguments{
  \item{x}{A sample data.}
  \item{object}{An R object obtained from fitting a 
        distribution to a data set.
  	Supported classes include 'gamma','norm','weibull','bde', 'em'.}
  \item{...}{Controls}
}
\author{
  B. Wang  \email{bwang@southalabama.edu}
}
\seealso{
 \code{\link{mediation.test},\link{perm.test}}.
}
\examples{

x = rnorm(1000,34.5,1.5)
y = round(x,0)
##  fit by histospline: class=bde


}
\keyword{htest}

