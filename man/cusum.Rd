\name{cusum}
\alias{cusum}
\alias{cusum.default}
\alias{ccusum}
\alias{bcusum}
\alias{plot.CCUSUM}
\alias{plot.BCUSUM}
\alias{print.CCUSUM}
\alias{print.BCUSUM}
\alias{CCUSUM}
\alias{BCUSUM}
\alias{Restart}
\alias{Restart.default}
\alias{Restart.CCUSUM}
\title{
  Cumulative Sum Control Chart.
}
\encoding{UTF-8}
\description{
  To construct the CUSUM plot.}
\usage{
  cusum(x, ...)
  ccusum(x,mu=NULL, k=0.5,h=4, ...)
  bcusum(x,prob=NULL, R0=1.0,Ra=2, ...)
}
\arguments{
  \item{x}{A sequnence of results.}
  \item{mu}{mean values if 'x' is continous.}
  \item{h}{A parameter to control recalculating the new mean to restart the CUSUM 
  	     chart.  Default: 4.}
  \item{k}{Number of shift in sd unit, which controls the reference value
           K=k*s.  Default: 1/2.}
  \item{prob}{ predicted probabilities for risk-adjusted CUSUM.}
  \item{R0, Ra}{Odds ratios under H0 and Ha, respectively.Default: \code{R0=1}, \code{Ra=2}}
  \item{...}{Controls}
}

\details{
  To perform either continuous or binary CUSUM analysis based on the input "x".

  In \code{bcusum}, the probability c0 under H0 is computed based on
current data, which can also be given as a parameter \code{prob}.

  Choose lower and upper limits h0 and h1: h1=1.5, 2.0, 2.5, 3, 10; h0=-2.5

  Wald: h0 = -ln(1-alpha)/beta, h1 = ln(1-beta)/alpha.  Alert lines:
  alpha=beta=0.1.  Alarm/action lines: alpha=beta=0.01.

  Rogers et al defined a scale adjustment ln(OR), say log(1.5)

  de Leval et al: alert line with alpha=0.05, beta=0.2; alarm line
  with alpha=0.01, beta=0.1


  }

\value{ If the binary CUSUM is used, a class "BCUSUM" will be
  produced, otherwise, a class "CCUSUM" will be produced.  

  To detect increases we set \code{Ra > R0}, otherwise \code{Ra < R0}.
}

\references{
de Leval, M. R., Francois, K., Bull, C., Brawn, W., Spiegelhalter, D., Mar 1994. Analysis of a cluster of surgical failures. Application to a series of neonatal arterial switch operations. J. Thorac. Cardiovasc. Surg. 107, 914–923.

Page, E. S., 1954. Continuous inspection schemes. Biometrika 41 (1/2), pp. 100–115. URL http://www.jstor.org/stable/2333009

Parsonnet, V., Dean, D., Bernstein, A. D., Jun 1989. A method of uniform stratification of risk for evaluating the results of surgery in acquired adult heart disease. Circulation 79, 3–12. 

Steiner, S. H., Cook, R. J., Farewell, V. T., Jan 1999. Monitoring paired binary surgical outcomes using cumulative sum charts. Stat Med 18, 69–86. 

Steiner, S. H., Cook, R. J., Farewell, V. T., Treasure, T., Dec 2000. Monitoring surgical performance using risk-adjusted cumulative sum charts. Biostatistics 1, 441–452.
}

\note{
  Risk-adjusted CUSUM through the second parameter.
}

\seealso{
 \code{\link{robot}}
}

\examples{
  data(robot)
  hst = robot$TMHSTMIN
  out = cusum(hst,mu=45)
  plot(out)
}
\keyword{htest}
