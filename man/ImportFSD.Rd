% This is a template R document file
% Distributed under GPL 3 or later

\name{ImportFSD}
\alias{ImportFSD}

\title{Import Firm Size and Firm Age Data}
\description{
 To read firm size and/or firm age data from built-in datasets.
}
\usage{
 ImportFSD(x,type,year)
}
\arguments{
  \item{x}{A built-in firm size and/or firm age dataset.}
  \item{type}{type of data: "size" or "age". If missing, read the age data from rows and the size data from columns.}
  \item{year}{The number year the firm size/age data to be read. If missing, assume the year is 2014.}
}

\value{
 \item{xy}{a matrix of joint frequency distribution table.}
 \item{breaks}{the class boundaries for firm age as a component \code{age}, and for firm size as a component \code{size}.}
 \item{size}{a 'bdata' subject of the firm size data}
 \item{age}{a 'bdata' object of the firm age data}
}

\examples{

 data(FSD)
 ## bivariate data
 xy = ImportFSD(Firm2)
 ## firm age of 2013
 x = ImportFSD(FirmAge, type="age", year=2013)
 ## firm size of 2013
 y = ImportFSD(FirmSize, type="size", year=2013)

}

\keyword{distribution}

