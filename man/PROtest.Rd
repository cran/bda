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
 pro.test(x,y,group,cutoff,x.range,type)
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
\item{x.range}{
  Range of the scores for 'x' and 'y'.
}
\item{type}{
  Data (grouping/binning) type: 'vas', 'nrs', 'wbf'.
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
