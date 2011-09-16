##  The following programs help to prepare data to be analyzed with
##  the FMMBD functions.  In summary, data should have three columns:
##  x = center of the bins, widths = width of the bins, and counts =
##  counts (frequencies) of the bins.


###  function binning(...) 11/25/2011
##  Data could be rounded to the nearest integers (method='nearest'),
##  or rounded up (method='up'), or rounded down (method = 'down').
##  Some data may be recorded and rounded to the nearest 100mg.  In
##  that case, we can specify scale=100 (default scale=1).

binning <- function(x, scale=1, method='nearest'){
  name <- deparse(substitute(x))
  if(scale<=0) stop("Wrong value for 'scale'")
  method = match.arg(tolower(method), c('nearest', 'up','down'))
  x0 = switch(method,
    nearest = round(x/scale),
    up      = ceiling(x/scale) - 0.5,
    down    = floor(x/scale) + 0.5
    )
  tmp = table(x0);
  y = as.numeric(names(tmp))
  f = as.numeric(tmp)
  w = rep(1,length(f))
  mu = .wmean(y,f);  s = .wsd(y, f);
  structure(list(x = y, widths = w, counts = f,mu=mu,s=s,
                 scale = scale, data.name = name),
            class='bdata')
}

print.bdata <- function (x, digits = NULL, ...) 
{
  cat("\nData:\t", deparse(x$data.name), sep = "")
  cat("Mean: ", format(x$mu,6),
      ";\t\tStd.Dev: ", format(x$s,6), "\n\n",sep = "")
  print(data.frame(Center=x$x, Counts=x$counts, BinWidth=x$widths))
  cat("\n\n")
  invisible(x)
}

hist.bdata <- function(x,...)
  {
    x0 = rep(x$x, x$counts)
    breaks = x$x + 0.5 * x$widths
    breaks = c(x$x[1] - 0.5*x$widths[1],breaks)
    hist(x0, breaks=breaks,probability=TRUE,...)
  }

