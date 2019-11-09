% This is a template R document file
% Distributed under GPL 3 or later

\name{compare}
\alias{compare}
\alias{compare.default}
\alias{compare.matrix}
\alias{compare.data.frame}
\alias{compare.NGS.lognormal}
\alias{print.comp}

\title{Compare two or more samples}
\description{
An automated hypothesis testing procedure.
}
\usage{
  compare(y, x, var.name, alternative = c("two.sided", "less", "greater"),
          paired = FALSE, conf.level = 0.95, ...)
}
\arguments{
  \item{x,y}{Two vector. 'NA' values will be automatically removed.}
  \item{var.name}{Names of the two variables being compared for \code{x}
  and \code{y}, respectively.}
  \item{alternative,paired,conf.level,...}{parameters for stat inferences.}
}
\author{
  B. Wang  \email{bwang@southalabama.edu}
}
\details{
 Will be updated later.
}
\keyword{stats}

