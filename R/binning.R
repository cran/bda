binning <- function(x, weights, breaks, binned = FALSE,
                    range.x, na.rm=TRUE) {
  
  name <- deparse(substitute(x))
  if (!is.numeric(x)) 
    stop("'x' must be numeric")
  if(!missing(weights)){
    if (!is.numeric(weights)) 
      stop("'weights' must be numeric")
    stopifnot(all(weights>0))
  }

  if(missing(weights)){
    x.na = is.na(x);
    if(any(x.na)){
      if(na.rm){
        x = x[!x.na]; 
      }else stop("Data contain missing value(s)")
    }
  }else{
    x.na = is.na(x); w.na = is.na(weights)
    if(any(x.na|w.na)){
      if(na.rm){
        sele = !x.na & !w.na
        x = x[sele]; weights = weights[sele]
      }else stop("Data contain missing value(s)")
    }
  }
  
  if (missing(range.x))
    range.x <- c(min(x), max(x))
  a = range.x[1L]; b = range.x[2L]

  if(binned){
    if(missing(weights)){  # 'x' is prebinned/rounded, but not tabled
      tmp = table(x)
      x0 = as.numeric(names(tmp));
      y0 = as.numeric(tmp);
      N = length(x)
    }else{ # frequencies are given
      x0 = x; y0 = weights;
      N = sum(weights)
    }

    widths = diff(x0); width = widths[1L]
    if(any(widths != width)){# not equally-spaced
      if (missing(breaks)) breaks = 'lscv'
      if (is.character(breaks)) {
        breaks <- match.arg(tolower(breaks),
                            c("sturges", "fd", "freedman-diaconis", "scott", "lscv"))
        nclass <- switch(breaks,
                         sturges = nclass.Sturges(x), 
                         `freedman-diaconis` = ,
                         fd = nclass.FD(x),
                         scott = nclass.scott(x),
                         lscv = .bwlscv(x, w = weights),
                         stop("unknown 'breaks' algorithm"))
        h = (b-a)/nclass
        breaks = seq(a-h/2, b+h/2, length = nclass + 1)
      }else if(length(breaks)==1){
        nclass = breaks
        h = (b-a)/nclass
        breaks = seq(a-h/2, b+h/2, length = nclass + 1)
      }
      
      y = .wbin(x, weights, breaks)
      x = breaks[-1] - h*0.5
    }else{
      y = y0; x = x0
    }
  }else{
    if(missing(weights)){  # raw data unweighted
      weights = rep(1, length(x))
    }else{
      weights = weights / sum(weights) * length(x)
    }
      
    if (missing(breaks)) breaks = 'lscv'
    if (is.character(breaks)) {
      breaks <- match.arg(tolower(breaks),
                          c("sturges", "fd", "freedman-diaconis", "scott", "lscv"))
      nclass <- switch(breaks,
                       sturges = nclass.Sturges(x), 
                       `freedman-diaconis` = ,
                       fd = nclass.FD(x),
                       scott = nclass.scott(x),
                       lscv = .bwlscv(x, w = weights),
                       stop("unknown 'breaks' algorithm"))
      h = (b-a)/nclass
      breaks = seq(a, b, length = nclass + 1)
    }else if(length(breaks)==1){
      nclass = breaks
      h = (b-a)/nclass
      breaks = seq(a, b, length = nclass + 1)
    }

    y = .wbin(x, weights, breaks)
    x = breaks[-1-nclass]
  }
  
  list(y = y, x = x)
}
         
.wbin <- function(X, W, gpoints, truncate = TRUE)
{
  n <- length(X)
  if(missing(W)) W = rep(1, n)
  M <- length(gpoints)
  if (truncate) trun <- 1L else 0L
  a <- gpoints[1L]
  b <- gpoints[M]
  out = .Fortran(.F_wbin, as.double(X), as.double(W), as.integer(n),
    as.double(a), as.double(b), as.integer(M),
    as.integer(trun), wts = double(M))$wts
  out[-M]
}

## For application of linear binning to a regression
## data set.

##  07/30/2012: rlbin from KernSmooth 
.rlbin <- function(X, Y, gpoints, truncate = TRUE)
{
    n <- length(X)
    M <- length(gpoints)
    trun <- if (truncate) 1L else 0L
    a <- gpoints[1L]
    b <- gpoints[M]
    out <- .Fortran(.F_rlbin, as.double(X), as.double(Y), as.integer(n),
                    as.double(a), as.double(b), as.integer(M), as.integer(trun),
                    double(M), double(M))
    list(xcounts = out[[8L]],  ycounts = out[[9L]])
}
