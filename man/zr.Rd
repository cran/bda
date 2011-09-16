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
\seealso{
 \code{\link{birth}},
 \code{\link{ofc}}.
}
\examples{
 data(zr)
 male = zr$Sex=="male"
 zar = zr$Region=="Zarichne"
 sum(male&zar)
 sum(!male&zar)
 sum(male&!zar)
 sum(!male&!zar)

 out1 = bootkde(zr$Head[male&zar])
 out2 = bootkde(zr$Head[!male&zar])
 out3 = bootkde(zr$Head[male&!zar])
 out4 = bootkde(zr$Head[!male&!zar])

 postscript(file='zr.eps', paper='letter')
 plot(out3, lwd=2, lty=1, ylab="Density", xlab="OFC",
   main="", xlim=c(29,40), ylim=c(0,.35))
 lines(out4, lwd=2, lty=2)
 lines(out1, lwd=4, lty=1)
 lines(out2, lwd=4, lty=2)
 legend(40,.35,xjust=1, yjust=1,
 	legend=c(
	 "Rivne City in non-Polissia (males)",
	 "Rivne City in non-Polissia (females)",
	 "Zarichne County in Polissia (males)",
	 "Zarichne County in Polissia (females)"),
        lwd=c(2,2,4,4), lty=c(1,2,1,2))
 dev.off()

}
\keyword{datasets}
