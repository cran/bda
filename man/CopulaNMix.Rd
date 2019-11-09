% This is a template R document file
% Distributed under GPL 3 or later

\name{fit.NGS.nmix.copula}
\alias{fit.nmix.copula}
\alias{fit.mlnorm.copula}
\alias{fit.NGS.nmix.copula}
\alias{fit.GLD.copula}
\alias{fit.GLD2bin}
\alias{print.copula}
\alias{print.copula2}
\alias{persp.copula}

\title{Fit a Bivariate Copula Based on Normal Mixture Model for NGS data}
\description{
   Fit bivariate copula for normal mixture models}
\usage{
  fit.NGS.nmix.copula(x,y,mle.large=FALSE,gridsize=101)
  fit.nmix.copula(x,y,mle.large=FALSE,gridsize=101)
  fit.GLD.copula(x,y,mle.large=FALSE,gridsize=101)
  fit.mlnorm.copula(x,y,mle.large=FALSE,gridsize=101,k=2)
}
\arguments{
  \item{x,y}{Two Next-Generation Sequencing data with same structures}
  \item{mle.large}{Default \code{,mle.large=FALSE}, don't use the big
  two-way histogram to find MLE.}
  \item{gridsize}{The size of the grid in the xy-plane to evaluate the
  joint density.}
  \item{k}{number of components}
}

\author{
  B. Wang  \email{bwang@southalabama.edu}
}

\details{ the two profiles need to be sent as parameters.  We don't
  need to fit the normal mixture outside the function.  But we need to
  fit normal mixture model for NGS data withinin the function to get
  the histograms, and the counts within two 2x2 tables. The output
  include 1) the two fitted model, and 2) estimated parameter Psi.  In
  addition, we show 3) the Psi's used to search for the optimal Psi,
  and 4) Chi-square test statistics, and 5) the correponding p-values.
  In addition, we also find 6) the MLE of Psi, and 7) its Chi-square
  test statistic and 8) p-value.

}

\references{
  to be updated
}
\examples{
 # To be updated.
}
\keyword{smooth}

