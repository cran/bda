% This is a template R document file
% Distributed under GPL 3 or later

\name{perm}
\alias{perm}
\alias{perm.test}
\alias{perm.default}
\title{To perform a permutation test to compare two samples/populations.}
\description{
  To perform a permutation test to compare two samples/populations.
}
\usage{
  perm(x,y,fun,alternative = "two.sided", trials = 1000,...)
}
\arguments{
  \item{x}{Numerical or logical vector.}
  \item{y}{Numerical or logical vector.}
  \item{fun}{User-defined function to compare the difference of a function for two samples/data sets. If missing, compare the edf's.}
  \item{alternative}{Type of test.  Either 'one-sided' or 'two-sided'.}
  \item{trials}{Iteration times, default: 10,000.}
  \item{...}{Controls.}
}
\details{
To test a null hypothesis "H0: the two distributions are the same." against "H1: the two distributions are different". The test is performed by (1) computing the difference of the sample statistic; (2) pooling the two samples together to form a permutated population; (3) repeatedly bootstrapping/resampling from the permutated population without replacement, and computing the sample difference. (4) finding the p-value by check the percentage of times that the absoluted sample statistics differences are greater than the observed difference from (1).  
}
\value{

Missing values are not allowed.
}

\author{
  B. Wang  \email{bwang@jaguar1.usouthal.edu}
}

\references{
R. A. Fisher, The Design of Experiment, New York: Hafner, 1935.

Pitman, E. J. G., "Significance tests which may be applied to samples from any population", Royal Statistical Society Supplement, 1937; 4: 119-130 and 225-32 (parts I and II).

Pitman, E. J. G., "Significance tests which may be applied to samples from any population. Part III. The analysis of variance test", Biometrika, 1938; 29: 322-335.

E. S. Edgington, Randomization tests, 3rd ed. New York: Marcel-Dekker, 1995.
Phillip I. Good, Permutation, Parametric and Bootstrap Tests of Hypotheses, 3rd ed., Springer, 2005. ISBN 0-387-98898-X

Good, P. (2002) Extensions of the concept of exchangeability and their applications, J. Modern Appl. Statist. Methods, 1:243-247.

Lunneborg, Cliff. Data Analysis by Resampling, Duxbury Press, 1999. ISBN 0-534-22110-6.

Pesarin, F. 2001. Multivariate Permutation Tests, John Wiley & Sons.

Welch, W. J., Construction of permutation tests, Journal of American Statistical Association, 85:693-698, 1990.

  Wang, B. and Wertelecki, W. (2012) 
  Density Estimation for Data With Rounding Errors.
  \emph{Computational Statistics and Data Analysis}, 
  (in press), doi: 10.1016/j.csda.2012.02.016. 

}

\seealso{
 \code{\link{gof},\link{mediation.test}}.
}

\examples{

x1 = rnorm(100,34.5,1.5)
x2 = round(x1,0)
perm.test(x1,x2,trials=100)

## use user-defined function
foo = function(x,y) mean(x)-mean(y)
perm.test(x1,x2,fun=foo,trials=100)
## another user-defined function example
x3 = round(rnorm(10000,34.2,1.4),0)
foo = function(x,y){
 if(length(x)>length(y))
 cutoff = mean(x)-3*sd(x)
 else
 cutoff = mean(y)-3*sd(y)
 mean(x<cutoff)-mean(y<cutoff)
}
perm.test(x3,x2,fun=foo,trials=100)

}
\keyword{stats}

