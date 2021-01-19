% This is a template R document file
% Distributed under GPL 3 or later

\name{bootPRO}
\alias{bootPRO}

\title{Effectiveness Evaluation based on PROs with bootstrapping method}

\description{
 Effectiveness Evaluation based on PROs with bootstrapping method.
}

\usage{
 bootPRO(x,type="relative",MCID,iter=999,conf.level=0.95)
}
\arguments{
  \item{x}{Data frame of repeated PRO measures. First column gives the
 treatment groups: 'group' or 'treat'; Repeated measures listed after
 column 1.}
  \item{type}{use absolute changes ("absolute") or relative (percent)
 changes ("relative").}
  \item{MCID}{A positive value to define responders.}
  \item{iter}{number of iterations used by the bootstrap algorithm.}
  \item{conf.level}{Confidence level of the bootstrapping confidence interval.}
}

\value{NONE.}

\references{
To be updated.
}

\examples{
  data(Pain)
  x <- pain[,c(2,3:9,11:17)]
  grp <- rep("treat",nrow(x))
  grp[x[,1]==0] <- "control"
  x[,1] <- grp
  #bootPRO(x,type='mean')
  #bootPRO(x,type='mean',MCID=1)
  #bootPRO(x,type='mean',MCID=1.5)
  #bootPRO(x,type='mean',MCID=2)

  #bootPRO(x,type='rel')
  #bootPRO(x,type='relative',MCID=.2)
  #bootPRO(x,type='relative',MCID=.3)
  ##bootPRO(x,type='relative',MCID=.5)
 }
\keyword{smooth}

