##  07/24/2012

##  To define a method to construct pointwise confidence bands.

pcb <- function(x, level=0.95, ...) UseMethod("pcb")

pcb.default <- function(x, level=0.95,...)
{
  stopifnot(level>0 & level <1)
  xbar = mean(x, na.rm=TRUE)
  s = sd(x, na.rm=TRUE)
  n = length(x)
  cv = abs(qt(0.5-0.5*level, n-1))
  B = cv * s/sqrt(n)
  out = c(xbar-B, xbar+B)
  cat(level*100, "% confidence interval (using t-distribution):\n\n")
  cat("\t", out, "\n")
  invisible(out)
}



pcb.edf <- function(x, level = 0.95,...){
  stopifnot(level>0 & level <1)
  n = x$n; x0 = x$x; y0 = x$y; alpha = 1 - level
  
  epsn = sqrt(.5/n * log(2./alpha))
  LFn = y0 - epsn; LFn[LFn<0] = 0
  UFn = y0 + epsn; UFn[UFn>1] = 1

  m = length(x0)
  x1 = c(x0[-1], x0[m] * 1.2);
  y1 = c(UFn)
  y2 = c(LFn)

  out = structure(list(
    x0 = x0, x1 = x1, y1 = y1, y2 = y2,
    data = x), class = "pcb")
}


pcb.hist <- function(x, level = 0.95,...){
  stopifnot(level>0 & level <1)
  n = x$n; x0 = x$x; y0 = x$y; alpha = 1 - level
  m = x$nclass
  
  epsn = 0.5 * qnorm(1-alpha/2/m) * sqrt(m/n)
  LFn = sqrt(y0) - epsn;
  LFn[LFn<0] = 0; LFn=LFn^2
  UFn = (sqrt(y0) + epsn)^2

  breaks = x$plot$breaks
  m = length(breaks)
  x0 = breaks[-m];
  x1 = breaks[-1];

  y1 = c(UFn)
  y2 = c(LFn)

  out = structure(list(
    x0 = x0, x1 = x1, y1 = y1, y2 = y2,
    data = x), class = "pcb")
}

print.pcb <- function (x, digits = NULL, ...) print(x$data)

plot.pcb  <- function (x, ...) 
{
  plot(x$data,...)
  lines(x,...)
  invisible(NULL)
}

lines.pcb  <- function (x, col='gray',...) 
{
##  segments(x$x0, x$y1, x$x1, x$y1, col=col,...)
##  segments(x$x0, x$y2, x$x1, x$y2, col=col,...)

  x0 = as.vector(rbind(x$x0, x$x1))
  y1 = as.vector(rbind(x$y1, x$y1))
  y2 = as.vector(rbind(x$y2, x$y2))

  x3 = c(x0, rev(x0)); y3 = c(y1, rev(y2))
  polygon(x3, y3, border='aliceblue')
  invisible(NULL)
}
