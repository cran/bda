\name{pro.test}
\alias{pro.test}

\encoding{UTF-8}
\title{
 Test effectiveness based on PROs
}
\description{
  Tests for effectiveness evaluations based on PROs.
}
\usage{
 pro.test(x,y,group,cutoff,class.names,conf.level,x.range,delta)
}
\arguments{
\item{x,y}{
  vector of PROs at T0 and T1.
}
\item{group}{
 Group assignment: control or treatment.
}
\item{cutoff}{
  Class boundaries to define states. Works only when 'x' and 'y' are
  numeric.
}
\item{class.names}{
  Labels/names for the states/classes.
}
\item{conf.level}{
  confidence level.
}
\item{x.range}{
  Range of the scores for 'x' and 'y'.
}
\item{delta}{
  an offset added to the class boundaries for sensitivity analyses.
}
}
\value{
 To be added.
}

\details{
 To be added.
}

\references{
 To be added.
}
\examples{
states <- c("low", "moderate", "high")
x0 <- sample(states, size=100, replace=TRUE)
x1 <- sample(states, size=100, replace=TRUE)
grp <- c(rep("control",50),rep("treatment",50))
pro.test(x=x0,y=x1,group=grp)
}
\keyword{test}
