% This is a template R document file
% Distributed under GPL 3 or later

\name{getbdata}
\alias{getbdata}
\title{Birth data anaylysis.}
\description{
     To compute the density function of various variables for birth data.
}
\usage{
 getbdata(x,type='head',gender="male",location="rc",filter=FALSE)
}
\arguments{
  \item{x}{Data frame of the birth data.}
  \item{type}{The variable whose density is to be computed: length, weight, head (or ofc), chest, mother.s.age, gestation (or ga).}
  \item{gender}{Gender of the births based on which the density is to be computed: "male" (default),"female".}
  \item{location}{Location of the births based on which the density to be computed: "r","u","zz","bz","sz","rc".}
  \item{filter}{Filter the data? Default: FALSE.  If "TRUE", the data will be filtered such that the analysis will be performed based on those with GA>=38, singleton livebirth, not zero head, mother's age not missing.}
}
\author{
  B. Wang  \email{bwang@jaguar1.usouthal.edu}
}
\seealso{
\code{\link{birth}},\code{\link{smkde}},\code{\link{histospline}}.
}
\examples{

data(birth)
getbdata(birth,type='head')

}
\keyword{stats}

