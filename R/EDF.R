
###  Empirical distribution function

### The folloowing functions are to compute the empirical distribution
### function and construct a confidence band.

### Last updated on 04/18/2013

edf <- function(x, weights, freq, from, to, gridsize, digits=0,
                rounding) UseMethod("edf")


edf.default <- function(x, weights, freq, from, to, gridsize,
                        digits=0, rounding){
  if(missing(rounding)) rounding <- "none"
  xr <- .data.rounding(x=x,freq=freq, weights=weights, digits=digits,
                  type=rounding)
  edf.bdata(xr,gridsize=gridsize,from=from,to=to)
}

edf.bdata <- function(x, weights, freq, from, to, gridsize,
                      digits=0, rounding){
  xr <- x;
  if(missing(from)) from <- min(xr$x+xr$a)
  if(missing(to)) to <- max(xr$x+xr$b)
  stopifnot(to > from)
  if(missing(gridsize)){
    xgrid <- xr$x
    gridsize <- length(xgrid)
    Fx <- cumsum(xr$f*xr$w)
    y <- Fx/Fx[gridsize];
  }else{
    if(gridsize < 10) gridsize <- 512L
    xgrid <- seq(from, to, length=gridsize)
    gcounts <- .data.binning(xr, type='edf',from=from,to=to,gridsize=gridsize)
    y <- cumsum(gcounts);
#    y <- Fx/Fx[gridsize];
  }
  
  structure(list(y = y, x = xgrid, n = gridsize, 
                 data.name = xr$data.name),
            class = "edf")
}

print.edf <- function (x, digits = NULL, ...) 
{
  cat("\nCall:\n\t", deparse(x$call), "\n\nData: ", x$data.name, 
      " (", x$n, " obs.);", "\n", sep = "")
  print(summary(as.data.frame(x[c("x", "y")])), digits = digits, 
        ...)
  invisible(x)
}

plot.edf <- function (x, main = NULL, xlab = NULL, ylab = "EDF", type
    = "l", zero.line = TRUE, ...)

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
