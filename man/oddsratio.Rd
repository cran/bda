\name{oddsratio}
\alias{oddsratio}
\alias{print.odds}
\title{Odds Ratio and Relative Risk}

\description{To compute the odds ratio and relative risk based on 
a 2 X 2 table.}

\usage{
  oddsratio(x, alpha=0.05,...)
}
\arguments{
  \item{x}{a 2x2 table.}
  \item{alpha}{The significance level. Default: 0.05.}
  \item{...}{Controls}
}
\details{

\code{x} can be a matrix or a data.frame: the first columns showing
the number of events and the second column showing the sample sizes.

Exact confidence limits for the odds ratio by using an algorithm based
on Thomas (1971). See also Gart (1971).  If the sample sizes are too
large, the exact confidence interval may not work due to overflow
problem.

Asymptotic confidence limits are computed according to SAS/STAT(R) 9.2
User's Guide, Second Edition.

Score method: code has been published for generating confidence
intervals by inverting a score test. It is available from
http://web.stat.ufl.edu/~aa/cda/R/two_sample/R2/

See also "riskratio"  and "oddsratio" in R package \emph{epitools}.

}

\value{

\item{OR}{an estimate of odds ratio;}
\item{RR}{an estimate of realtive risk;}
\item{ORCI}{A table showing various (1-alpha)\% confidence limits for OR;}
\item{RRCI}{A table showing various (1-alpha)\% confidence limits for RR;}

}

\references{


    Agresti, A. (1990) _Categorical data analysis_.  New York: Wiley.
     Pages 59-66.

Agresti, A. (1992), A Survey of Exact Inference for Contingency Tables
Statistical Science, Vol. 7, No. 1. (Feb., 1992), pp. 131-153.

Agresti, A. (2002), Categorical Data Analysis, Second Edition, New
York: John Wiley \& Sons.

     Fisher, R. A. (1935) The logic of inductive inference.  _Journal
     of the Royal Statistical Society Series A_ *98*, 39-54.

     Fisher, R. A. (1962) Confidence limits for a cross-product ratio.
     _Australian Journal of Statistics_ *4*, 41.

     Fisher, R. A. (1970) _Statistical Methods for Research Workers._
     Oliver & Boyd.

     Mehta, C. R. and Patel, N. R. (1986) Algorithm 643. FEXACT: A
     Fortran subroutine for Fisher's exact test on unordered r*c
     contingency tables.  _ACM Transactions on Mathematical Software_,
     *12*, 154-161.

     Clarkson, D. B., Fan, Y. and Joe, H. (1993) A Remark on Algorithm
     643: FEXACT: An Algorithm for Performing Fisher's Exact Test in r
     x c Contingency Tables.  _ACM Transactions on Mathematical
     Software_, *19*, 484-488.

     Patefield, W. M. (1981) Algorithm AS159.  An efficient method of
     generating r x c tables with given row and column totals.
     _Applied Statistics_ *30*, 91-97.

Stokes, M. E., Davis, C. S., and Koch, G. G. (2000), Categorical Data
Analysis Using the SAS System, Second Edition, Cary, NC: SAS Institute
Inc.

}

\seealso{
\code{\link{fisher.test}}, \code{\link{chisq.test}}
}

\keyword{htest}
