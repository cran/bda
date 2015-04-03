\name{ppower}
\alias{ppower}

\title{
 Compute the power of an ANOVA-type test for proportions
}

\encoding{UTF-8}

\description{
 To compute the power of an ANOVA-type test for proportions.  
 It can be used for the sample size estimation.
}

\usage{
  ppower(effect.size, ngroup, n0, n1, range.p, alpha=0.05, gridsize=512)
}
\arguments{
  \item{effect.size}{Difference between the two proportions}
  \item{ngroup}{Levels of the treatment}
  \item{n0, n1}{Sample sizes.}
  \item{range.p}{Approximate range of the proportion of control group.}
  \item{alpha}{Significance level: experimentwise error rate.}
  \item{gridsize}{A (fine) grid size to compute the power when no (less) 
  prior knowledge about \code{p0} can be used.}
}


\value{ The power of the test.}

\references{
  Wang, Bin. (TBA).
}

\keyword{htest}
