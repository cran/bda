\name{robot}
\alias{robot}
\title{
  Robot surgery.
}
\encoding{UTF-8}
\description{Robot surgery results.}
\usage{
  data(robot)
}

\format{
A data frame with 388 rows and 111 columns (variables).
}

\references{

de Leval, M. R., Francois, K., Bull, C., Brawn, W., Spiegelhalter, D., Mar 1994. Analysis of a cluster of surgical failures. Application to a series of neonatal arterial switch operations. J. Thorac. Cardiovasc. Surg. 107, 914–923.

Page, E. S., 1954. Continuous inspection schemes. Biometrika 41 (1/2), pp. 100–115. URL http://www.jstor.org/stable/2333009

Parsonnet, V., Dean, D., Bernstein, A. D., Jun 1989. A method of uniform stratification of risk for evaluating the results of surgery in acquired adult heart disease. Circulation 79, 3–12. 

Steiner, S. H., Cook, R. J., Farewell, V. T., Jan 1999. Monitoring paired binary surgical outcomes using cumulative sum charts. Stat Med 18, 69–86. 

Steiner, S. H., Cook, R. J., Farewell, V. T., Treasure, T., Dec 2000. Monitoring surgical performance using risk-adjusted cumulative sum charts. Biostatistics 1, 441–452.

}

\seealso{
 \code{\link{cusum}}
}

\examples{
  data(robot)
  dim(robot)
}
\keyword{datasets}
