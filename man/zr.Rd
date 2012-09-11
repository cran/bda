\name{zr}
\docType{data}
\alias{zr}
\usage{data(zr)}
\title{occipitofrontal head circumference data}
\description{
OFC data for singleton live births with gestational age at least 38 weeks in city "R" and county "Z".
}
\format{A data frame with three variables.
  \tabular{lll}{
    \code{Sex}   \tab character  \tab 'male' or 'female' \cr
    \code{Region}  \tab character  \tab city "R" or county "Z". \cr
    \code{Head}   \tab numeric  \tab head size. \cr
  }
}
\references{
  Wang, B. and Wertelecki, W. (2012) 
  Density Estimation for Data With Rounding Errors.
  \emph{Computational Statistics and Data Analysis}, 
  (in press), doi: 10.1016/j.csda.2012.02.016. 
}

\examples{

 data(zr)
 male = zr$Sex=="male"
 zar = zr$Region=="Zarichne"

 out1 = bootkde(zr$Head[male&zar])

}
\keyword{datasets}
