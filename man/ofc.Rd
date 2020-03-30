\name{ofc}
\docType{data}
\alias{ofc}
\title{occipitofrontal head circumference data}
\description{
OFC data for singleton live births with gestational age at least 38 weeks.
}
\format{A data frame with 2019 observations on 4 variables.
  \tabular{lll}{
    \code{Year}  \tab numeric  \tab 2006 -- 2009 \cr
    \code{Sex}   \tab character  \tab 'male' or 'female' \cr
    \code{Gestation}  \tab numeric  \tab Gestational age (in weeks). \cr
    \code{Head}   \tab numeric  \tab head size. \cr
  }
}
\references{
Wang, B and Wertelecki, W, (2013)
Computational Statistics and Data Analysis,
65: 4-12.
}

\examples{
 data(ofc)
 head(ofc)
}

\keyword{datasets}
