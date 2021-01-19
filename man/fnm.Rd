% This is a template R document file
% Distributed under GPL 3 or later

\name{fnm}
\alias{fnm}
\alias{tkde}

\title{Distribution of Two Finite Gaussian Mixtures}

\description{ 

  To compute the values of the density and distribution functions of
  two finite Gaussian mixture models.
}

\usage{
  fnm(p1,p2,mu1,mu2,sig1,sig2,from,to)
}
\arguments{
  \item{p1,p2}{mixing coefficients.}
  \item{mu1,mu2}{vectors of the mean values of the Gaussian components.}
  \item{sig1,sig2}{vectors of the SD values of the Gaussian components.}
  \item{from,to}{to specify the range of data.} 
}

\value{
  Return the densities ('y') and probabilities ('Fx')
  over a grid of 'x'.
}

\examples{
data(Pain)

group <- pain$treat
x <- pain$recall0
y <- pain$recall1
#out <- tkde(x,y,group)
out <- tkde(x,y,group,type='percent')
plot(out$risk,type='l')
abline(h=1,col='gray')

plot(out$responder,type='l',ylim=c(-.28,.1))
lines(out$resp$ll~out$resp$x,lty=2,col=1+(out$resp$p<0.05))
lines(out$resp$ul~out$resp$x,lty=2,col=1+(out$resp$p<0.05))
abline(h=0,col='gray')

plot(out$g2,type='l')
lines(out$g1,col=2)

}
\keyword{distribution}

