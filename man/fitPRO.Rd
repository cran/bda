% This is a template R document file
% Distributed under GPL 3 or later

\name{fit.PRO}
\alias{fit.PRO}

\title{Fitting distributions to Patient-reported Outcome Data}
\description{
 To fit PRO data to distributions including GLD, PD, GPD, Weibull, EWD
 and other families.
}
\usage{
fit.PRO(x,dist,x.range,nclass)
}
\arguments{
  \item{x}{'x' can be a vector or a matrix, a 'histogram', or binned data.}
  \item{dist}{distribution type for 'x' and/or 'y'. Options include \code{Weibull}, \code{gpd} -- generalized Pareto distribution, \code{pd} or
 \code{Pareto}, \code{EWD}-- exponentiated Weibull distribtion.}
  \item{x.range}{Specifies the range of data. Used only for interval PRO data.}
  \item{nclass}{Number of classes/bins. Used only for interval PRO data.}
}

\value{
 \item{x.fit, y.fit}{Fitted marginal distribution for row- and column-data when 'x' is a matrix.}
 \item{Psi}{Plackett estimate of the Psi. ONLY for 'fit.FSD2'}

 If each of \code{x.fit} and \code{y.fit}, the following values are available:

 \item{xhist}{histogram.}
 \item{dist}{distribution type. Options include \code{Weibull},
   \code{gpd} -- generalized Pareto distribution, \code{pd} or
   \code{Pareto}, \code{EWD}-- exponentiated Weibull distribtion.}
 \item{size}{Total number of observations (sample size)}
 \item{pars}{Estimates of parameters.}
 \item{y,y2,x}{the pdf (\code{y}) and cdf (\code{y2}) values evaluated on a grid \code{x}}
}
\examples{
x <- rweibull(1000,2,1)
(out <- fit.PRO(x))

}

\keyword{distribution}

