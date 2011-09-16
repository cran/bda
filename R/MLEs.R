mle.gamma <- function(x){
  x = x[!is.na(x)]
  .Fortran(.F_FitGamma, as.double(x),
           as.integer(length(x)), as.double(rep(0,2)))[[3]]
}

mle.weibull <- function(x){
  x = x[!is.na(x)]
  .Fortran(.F_FitWeibull, as.double(x),
           as.integer(length(x)), as.double(rep(0,2)))[[3]]
}
