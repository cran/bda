% This is a template R document file
% Distributed under GPL 3 or later

\name{blinding.est}
\alias{blinding.est}
\encoding{UTF-8}

\title{Estimating Effect Size using Blinding Survey Data}

\description{
  To estimate the effect size from an RCT with blinding survey data.
}

\usage{
  blinding.est(x, group, guess, type='cpe')
}

\arguments{
  \item{x}{outcome variable, continuous, missing values not allowed.}
  \item{group}{arm/group assignment; must have two levels. Missing values not allowed.}
  \item{guess}{Responses from blinding survey, can be two
  (0=sham/1=active treatment) or three (0=sham/1=active treatment/2=I
  don't know) levels.}
  \item{type}{three options: "simple" - ignoring blinding survey
  result. "adjusted" - unblinding adjusted estimate using 'guess' as a
  covariate in a regression model. "cpe" - change-point approach.}
}

\details{
 Point estimate, bootstrapping estimate, with/without multiple
 imputation if missing values exist in 'guess'.
 TBA  
}

\value{
  TBA
}

\examples{
 x <- sort(rnorm(20,6,2))
 y1 <- 1 + 1.5 * x[1:10] + rnorm(10)
 y2 <- 10 + 0.5 * x[11:20] + rnorm(10)
 y <- c(y1,y2)
 plot(y~x)
 abline(a=10,b=0.5, col='red',lty=2)
 abline(a=1,b=1.5, col='blue',lty=2)
}

\keyword{stats}

