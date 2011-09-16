% This is a template R document file
% Distributed under GPL 3 or later

\name{gof}
\alias{gof}
\alias{gof.em}
\alias{gof.bde}
\alias{gof.test}
\alias{gof.default}
\title{To perform goodness-of-fit test.}
\description{
  To perform goodness-of-fit test.
}
\usage{
  gof(object, x,...)
}
\arguments{
  \item{x}{A sample data.}
  \item{object}{An R object obtained from fitting a 
        distribution to a data set.
  	Supported classes include 'gamma','norm','weibull','bde', 'em'.}
  \item{...}{Controls}
}
\author{
  B. Wang  \email{bwang@jaguar1.usouthal.edu}
}
\seealso{
 \code{\link{mediation.test},\link{perm.test}}.
}
\examples{

x = rnorm(1000,34.5,1.5)
y = round(x,0)
##  fit by histospline: class=bde
out1 = histospline(y)
gof(out1,x)
gof(out1,y)
out2 = smkde(y)
gof(out2,x)
gof(out2,y)

## fit by em-alg: class=em
out3 = bfmm(y)
##gof(out3,x)
##gof(out3,y)
##  Naive estimates
gof("normal",x)
gof("weibull",x)
gof("gamma",x)

x0 = seq(min(x),max(x),length=100)
f0 = dnorm(x0,34.5,1.5)
plot(f0~x0,col=2,lwd=2,type='l')
lines(out1, col=3)
lines(out2, col=4)
lines(density(out3), col=5)
legend(min(x),max(f0),xjust=0,yjust=1,
  legend=c("True","Histospline","Smoothed KDE","Bin-EM"),
  col=c(2,3,4,5), lwd=c(2,1,1,1))

}
\keyword{stats}

