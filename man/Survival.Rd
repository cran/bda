\name{Survival}
\alias{Survival}
\alias{Hazard}
\alias{hazard}
\alias{print.surv}
\alias{weibull.mle}
\encoding{UTF-8}
\title{
 Functions for survival data analysis
}
\description{
 To compute the survival function S(t), harzard function h(t) 
 and cumulative harzard function H(t). 
}
\usage{
 Survival(x, censor, method="Nelson",from, to,gridsize)
 Hazard(x, censor, method="Nelson",from, to,gridsize)
 hazard(x, censor, method="Nelson",from, to,gridsize)
}
\arguments{
  \item{x}{life time data}
  \item{censor}{Status of \code{x}: \code{0} censored, \code{1} not censored.}
  \item{method}{Estimation method: "Aalen-Nelson" or "Kaplan-Meier" method.}
  \item{from,to,gridsize}{start point, end point and size of a fine 
  grid where the EDF will be evaluated.}

}

\author{
  B. Wang  \email{bwang@southalabama.edu}
}

\examples{

par(mfrow=c(1,1))
 n <- 30
 shape <- 1.5
 scale <- .8; 
 x <- rweibull(n,shape,scale)
 s <- rep(1, n)
 sele <- runif(n) > 0.9
 if(sum(sele)>0) s[sele] <- 0

# weibull.mle(x,s)

 x0 <- seq(0, max(x), length=100)
 f0 <- dweibull(x0, shape, scale)
 F0 <- pweibull(x0, shape, scale)
 S0 <- 1-F0
 S.km <- Survival(x, s, method="KM")
 S.na <- Survival(x, s, method="Nelson")
 S.wkde <- Survival(x,s, method='wkde')
 S.weibull <- Survival(x, s, method="weibull")
 S.emae <- Survival(x, s, method="mae.exp")
 S.wmae <- Survival(x, s, method="mae.weibull")

 plot(S0~x0, type='l', col=2, lty=2,ylim=c(0,1))
 lines(S.km, col=4)
 lines(S.na, col=1)
 lines(S.wkde, col=3, lwd=2)
 lines(S.weibull, col=5, lwd=3)
 lines(S.emae, col=2, lwd=3)
 lines(S.wmae, col=4, lwd=3)

 x <- c(10,7,32,23,22,6,16,34,32,25,11,20,19,6,17,35,6,13,9,6,10)
 w <- c(1,1,0,1,1,1,1,0,0,0,0,0,0,1,0,0,1,1,0,0,0)
 S.km <- Survival(x, w, method="KM")
 S.na <- Survival(x, w, method="Nelson")
 S.wkde <- Survival(x,w, method='wkde')
 S.weibull <- Survival(x, w, method="weibull")
 S.emae <- Survival(x, w, method="mae.exp")
 S.wmae <- Survival(x, w, method="mae.weibull")

 plot(S.km, type='l')
 lines(S.wkde, col=2)
 lines(S.na, col=4)
 lines(S.weibull, col=5, lwd=3)
 lines(S.wmae, col=3, lwd=3)
 lines(S.emae, col=6, lwd=3)


# x <- c(10,7,32,23,22,6,16,34,32,25,11,20,19,6,17,35,6,13,9,6,10)
# w <- c(1,1,0,1,1,1,1,0,0,0,0,0,0,1,0,0,1,1,0,0,0)
# H.na <- Hazard(x, w, method="Nelson")
# y <- H.na$y; x <- H.na$x; v <- H.na$var
# y <- diff(c(0,y));
# v <- diff(c(0,v))
# out1 <- lpsmooth(y,x,v,kernel='epan')    
# out2 <- lpsmooth(y,x,v,kernel='unif')    
# out3 <- lpsmooth(y,x,v,kernel='biweight')    
# plot(out1, type='l')
# lines(out2, col=2)
# lines(out3, col=4)
# out4 <- hazard(x,w,method='wkde')
# lines(out4, col=2, lwd=2)

}
\keyword{survival}
