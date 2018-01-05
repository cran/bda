% This is a template R document file
% Distributed under GPL 3 or later

\name{metamean.test}
\alias{metamean.test}
\alias{plot.MetaMeans}
\alias{print.MetaMeans}
\alias{lines.MetaMeans}


\title{Test Equal Means}
\description{
To test the equality of two means based on two sequences of sample means.
}
\usage{
  metamean.test(x,n,group,s,mu, sigma,years,
    alternative = c("two.sided", "less", "greater"),
    conf.level=0.95, sizable=FALSE,legend.pos=2)
}
\arguments{
  \item{x,n,s}{Vectors of sample means, sample sizes, and sample standard
  deviations.\code{s} can be missing.}
  \item{mu, sigma}{The vectors of the mean and standard deviation values for
  the populations or much larger sample.}
  \item{years}{The years of the collected data.  If missing, replace with a
  sequence of integers starting with \code{1}.}
  \item{group}{A vector showing the group IDs.  Must have two levels.}
  \item{alternative}{Type of test.}
  \item{conf.level}{Confidence level if the means.}
  \item{sizable}{A logical variable to control the size of the point symbols.}
  \item{legend.pos}{Position of the legend: 1=top right, 2=top left,
  3=bottom left, 4=bottom right, 0=no legend.}
}
\author{
  B. Wang  \email{bwang@southalabama.edu}
}
\details{
 To be updated.
}
\keyword{stats}

