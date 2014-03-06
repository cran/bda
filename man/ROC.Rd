% This is a template R document file
% Distributed under GPL 3 or later

\name{roc}
\alias{roc}
\alias{print.broc}
\alias{plot.broc}
\title{Receiver operating characteristic (ROC)}
\description{
 To compute and draw a ROC curve.
}

\usage{
  roc(y, x, z=NULL, event="less")
}
\arguments{
  \item{y}{Logical. Actual status:  event or not.}
  \item{x}{Numberical vector based on which the diagnosis is performed.}
  \item{z}{Logical. A covariate for diagnosis in addition to 'x'.} 
  \item{event}{Default: \code{less}.  To classify a case to be an 
       \code{event} if \code{x} is less than a \code{cutoff}.  
       Otherwise, a case will be classified as an event of 
       if \code{x} is greater than a \code{cutoff}.
       }
}

\value{
  The TPR or sensitivity, and FPR or (1-specificity) will be returned 
  together with the cutoffs.

  The area under the ROC curve can be computed using \code{summary} function.

  The ROC curve can be ploted using \code{plot} function.
}

\author{
  B. Wang  \email{bwang@southalabama.edu}
}

\examples{

x <- rnorm(100)
y <- rnorm(100)<0
out = roc(y,x)

}
\keyword{stats}

