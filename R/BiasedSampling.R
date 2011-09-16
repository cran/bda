## Define some typical biased sampling functions
## 2012/07/10 by Bin Wang

## length-biasing
bs.w.length <- function(x){
  if(any(x<=0)) stop("'x' must be positive")
  x
}
bs.W.length <- function(x){
  if(any(x<=0)) stop("'x' must be positive")
  x^2/2
}

## area-biasing
bs.w.area <- function(x){
  if(any(x<=0)) stop("'x' must be positive")
  x^2
}
bs.W.area <- function(x){
  if(any(x<=0)) stop("'x' must be positive")
  x^3/3
}

## inverse length-biasing
bs.w.invlength <- function(x){
  if(any(x<=0)) stop("'x' must be positive")
  1/x
}
bs.W.invlength <- function(x){
  if(any(x<=0)) stop("'x' must be positive")
  log(x)
}
