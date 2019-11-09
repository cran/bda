% This is a template R document file
% Distributed under GPL 3 or later

\name{ZipfPlot}
\alias{ZipfPlot}
\alias{Zipf.Normalize}
\alias{ZPlot}

\title{The Zipf Plot}

\description{ 
  Draw Zipf Plot.
}

\usage{
  ZipfPlot(x, x0, plot.new=TRUE, fitted=TRUE,weights,...)
  Zipf.Normalize(x, y, cutoff=6,optim=FALSE, method)
  ZPlot(x, plot.new=TRUE, ...)
}
\arguments{
  \item{x,y}{data: two vectors.}
  \item{x0}{low bound to filter data.}
  \item{plot.new}{whether draw a new plot.}
  \item{weights}{Compute weighted least squares line
  if \code{weights} is given.}
  \item{fitted}{Draw the fitted straight line.}
  \item{cutoff}{a large enought value such that the values larger
  than the \code{cutoff} (approximately) follows a power law distribution.}
  \item{optim}{Find the optimal normalization parameters if TRUE}
  \item{method}{use both power transformation and scalingby default.
  If 'scaling' is specified, skip power transformation. }
  \item{...}{plotting parameters.}
}
\author{
  B. Wang  \email{bwang@southalabama.edu}
}

\examples{

data(LCL)
head(LCL)

x <- LCL[,1]; y <- LCL[,2]
## Zipf plots
out1 <- ZipfPlot(x)
out2 <- ZipfPlot(y)
plot(out1, type='l', ylim=c(0,9), lwd=2,
     xlab="log(X)", ylab="log(Rank)")
lines(out2, lty=2,lwd=2,col=4)

## Normalize 'y' and use 'x' as baseline 
out <- Zipf.Normalize(x,y) # OLS estimates
out3 <- ZipfPlot(out$y)
## Zipf plot of the normalized profile
lines(out3, lty=3,lwd=2, col=2)

## Show fitted least squares lines
out$coef
abline(a=out$coef[1,1],b=out$coef[1,2], lty=1)
abline(a=out$coef[2,1],b=out$coef[2,2], lty=2,col=4)

## search for optimal estimates
### not run --slow
#out4 <- Zipf.Normalize(x,y,optim=TRUE) 
#persp(out4$mat.optim, theta=90, phi=30,
#      xlab="Gamma", ylab="Sigma", zlab="Dmax")
#plot(out1, type='l', ylim=c(0,9), lwd=2,
#     xlab="log(X)", ylab="log(Rank)")
#lines(out2, lty=2,lwd=2,col=4)
#out5 <- ZipfPlot(out4$y)
#lines(out5, lty=3,lwd=2, col=2)


## Linear scaling normalization
out6 <- Zipf.Normalize(x,y,method='scaling') 

plot(out1, type='l', ylim=c(0,9), lwd=2,
     xlab="log(X)", ylab="log(Rank)")
lines(out2, lty=2,lwd=2,col=4)
out7 <- ZipfPlot(out6$y)
lines(out7, lty=3,lwd=2, col=2)

### not run --slow
#out8 <- Zipf.Normalize(x,y,method='scaling',optim=TRUE) 

#plot(out1, type='l', ylim=c(0,9), lwd=2,
#     xlab="log(X)", ylab="log(Rank)")
#lines(out2, lty=2,lwd=2,col=4)
#out9 <- ZipfPlot(out8$y)
#lines(out9, lty=3,lwd=2, col=3)

}
\keyword{stats}

