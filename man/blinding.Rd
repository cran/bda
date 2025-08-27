% This is a template R document file
% Distributed under GPL 3 or later

\name{blinding}
\alias{blinding}
\encoding{UTF-8}

\title{Changing Point Estimation for Blinding Data}

\description{
  To find the changing point.
}

\usage{
  blinding(x, y)
}

\arguments{
  \item{x,y}{paired vectors - sorted by 'x'.}
}

\details{
 using the likelihood ratio test with two linear fits  
}

\value{
  the changing point at 'x'.
}

\examples{
 x <- sort(rnorm(20,6,2))
 y1 <- 1 + 1.5 * x[1:10] + rnorm(10)
 y2 <- 10 + 0.5 * x[11:20] + rnorm(10)
 y <- c(y1,y2)
 plot(y~x)
 abline(a=10,b=0.5, col='red',lty=2)
 abline(a=1,b=1.5, col='blue',lty=2)
 blinding(x,y)
}

\keyword{distribution}
\keyword{stats}

