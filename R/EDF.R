
### The folloowing functions are to compute the empirical distribution
### function and construct a confidence band.

edf <- function(x, weights, xgrid, na.rm = TRUE){
  if (!is.numeric(x)) 
    stop("'x' must be numeric")
  name <- deparse(substitute(x))
  if(missing(weights)) weights = rep(1, length(x))
  isna = is.na(x)
  if(any(isna)){
    if(na.rm){
      x = x[!isna]; weights = weights[!is.na]
    }else stop("'x' contains missing value(s)")
  }
  n=length(x); weights = weights/sum(weights)
  o = order(x); x = x[o]; weights = weights[o];

  if(missing(xgrid))
    xgrid = as.numeric(names(table(x)))
  else xgrid = sort(xgrid)  # in case xgrid is not sorted
  M = length(xgrid)
  
  y = .Fortran(.F_wedf, as.double(x), as.double(weights),
    as.integer(n), as.double(xgrid), F = double(M),
    as.integer(M))$F

  structure(list(
                 y = y, x = xgrid, n = n, 
                 data.name = name), class = "edf")
}

.cumprop <- function(x,y){
  mean(y<=x)
}

.recumprop <- function(x,y){
  x0 = round(x,0)
  sgnx = -1
  if(x0<x) sgnx = 1
  sele = x==x0
  mean(y<=x) - mean(sele)*(sgnx * 0.5 + x0 - x)
}


print.edf <- function (x, digits = NULL, ...) 
{
  cat("\nCall:\n\t", deparse(x$call), "\n\nData: ", x$data.name, 
      " (", x$n, " obs.);", "\n", sep = "")
  print(summary(as.data.frame(x[c("x", "y")])), digits = digits, 
        ...)
  invisible(x)
}

plot.edf  <- 
function (x, main = NULL, xlab = NULL, ylab = "EDF", type = "l", 
    zero.line = TRUE, ...) 
{
  if (is.null(xlab)) 
    xlab <- paste("N =", x$n)
  if (is.null(main)) 
    main <- deparse(x$data.name)
  plot.default(x, type = 'n', 
               main = main, xlab = xlab, ylab = ylab, ...)
  x0 = x$x; y0 = x$y; m = length(x0)
  x1 = c(x0[-1], x0[m] * 1.2);
  y1 = c(y0)
  segments(x0,y0,x1,y1,...)
  
  x1 = x1[-m]; y1 = y1[-m]
  x3 = x0[-1]; y3 = y0[-1]
  segments(x1,y1,x3,y3, col='gray',lty=3)
  if (zero.line) 
    abline(h = 0, lwd = 0.1, col = "gray")
  invisible(NULL)
}

lines.edf  <- function (x, ...) 
{
  x0 = x$x; y0 = x$y; m = length(x0)
  x1 = c(x0[-1], x0[m] * 1.2);
  y1 = c(y0)
  segments(x0,y0,x1,y1,...)
  
  x1 = x1[-m]; y1 = y1[-m]
  x3 = x0[-1]; y3 = y0[-1]
  segments(x1,y1,x3,y3, col='gray',lty=3)
  invisible(NULL)
}
