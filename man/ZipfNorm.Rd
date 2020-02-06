% This is a template R document file
% Distributed under GPL 3 or later

\name{Zipf.Normalize}
\alias{Zipf.Normalize}

\title{Zipf Normalization}

\description{ 
 Zipf plot based normalization.
}

\usage{
  Zipf.Normalize(x, y, cutoff=6,optim=FALSE, method)
}
\arguments{
  \item{x,y}{data: two vectors.}
  \item{cutoff}{a large enought value such that the values larger
  than the \code{cutoff} (approximately) follows a power law distribution.}
  \item{optim}{Find the optimal normalization parameters if TRUE}
  \item{method}{use both power transformation and scalingby default.
  If 'scaling' is specified, skip power transformation. }
}

\value{
\item{x}{reference profile (not normalized)}
\item{y}{normalized profile}
\item{scaler}{Linear rescaling normalization parameter estimate}
\item{power}{power transformation parameter estimate}
\item{scaler.optim}{Optimized estimate of the linear rescaling parameter}
\item{power.optim}{Optimzed estimate of the power transformation parameter.}
\item{mat.optim}{A matrix of the objective function values generated to find the optimal estimates.}
\item{coef}{Coefficient table to display the estimates.}
}

\references{
Wang, B. (2020)
A Zipf-plot based normalization method for high-throughput RNA-Seq data.
PLoS ONE, (Submitted).
}

\examples{
data(LCL)
names(LCL)
x <- LCL$p47
y <- LCL$p107
outx <- ZipfPlot(x)
plot(outx,type='l')
outy <- ZipfPlot(y)
lines(outy,col=2)

out2 <- Zipf.Normalize(x,y)
outy2 <- ZipfPlot(out2$y)
lines(outy2,col=4)

  }
\keyword{stats}

