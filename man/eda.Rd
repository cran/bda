% This is a template R document file
% Distributed under GPL 3 or later

\name{eda}
\alias{eda}
\alias{eda.default}
\alias{eda.data.frame}
\alias{eda.matrix}
\alias{print.desc}
\alias{plot.desc}
\title{Exploratory Data Analysis (EDA)}
\description{
A graphical method for exploratory data analysis.
}
\usage{
  eda(x,plot=FALSE)
}
\arguments{
  \item{x}{A vector of numeric values. 'NA' values will be automatically removed.}
  \item{plot}{Show distribution and others using graph(s).}
}
\author{
  B. Wang  \email{bwang@southalabama.edu}
}
\details{
  This function is modified based on an S-Plus function written by Dr. Jiayang Sun 
  (Course ST 426. Statistics Department, Case Western Reserve University).
.
}
\references{
 J. W. Tukey. (1977). Exploratory Data Analysis. Addison-Wesley, New York.
}
\examples{
 x = rnorm(100)
 eda(x)
}
\keyword{stats}

