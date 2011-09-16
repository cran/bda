##  Part of R package BDA
##  Copyright (C) 2009-2010 Bin Wang
##
##  Unlimited use and distribution (see LICENCE).

           
lprde <- function(x, weights, bandwidth, gridsize=256L, degree=1L,
                  kernel="gaussian", error = NULL, range.x,
                  na.rm=TRUE, binned = FALSE) {
  ## 07/30/2012: we develop algorithms for local polynomial with
  ## degree = 1.  The program can be developed for data with
  ## measurement errors.
  
  ## find the regression line and prepare information for
  ## simultanous confidence bands/pointwise confidence interval
  
  name <- deparse(substitute(x))
  
  ##  prepare data for nplpr.  Yi = sqrt(Qi+1/4), Qi is the count in
  ##  the i-th cell (bin)
  out = binning(x=x, weights=weights, binned=binned,
    range.x=range.x,na.rm=na.rm)
  stopifnot(all(out$y>=0))
  y0 = sqrt(out$y+0.25); x0 = out$x
  out = lpreg(x=x0, y=y0, bandwidth = bandwidth, gridsize = gridsize,
    degree=degree, kernel = kernel, error = error, binned=FALSE)

  x1 = out$x; y1 = out$y^2
  kappa = sum(y1)*diff(x1)[1L]
  out = structure(list(x = x1, y = y1/kappa,
    bw = out$bw, sp = NULL, data = x,
    call = match.call(), data.name = name),
    class='kde')
  invisible(out)
}

lpreg <- function(x, y, bandwidth, gridsize = 256L, degree = 1L,
                  kernel="gaussian", error = NULL, range.x, na.rm =
                  TRUE, binned = FALSE)
  
{
  name <- deparse(substitute(x))
  ## Install safeguard against non-positive bandwidths:
  if (!missing(bandwidth) && bandwidth <= 0)
    stop("'bandwidth' must be strictly positive")
  
  x.na = is.na(x); y.na = is.na(y)
  if(any(x.na | y.na)){
    if(na.rm){
      sele = !x.na & !y.na
      x = x[sele]; y = y[sele]
    }else stop("Data contain missing value(s)")
  }
  
  if (missing(range.x))
    range.x <- c(min(x), max(x))
  a = range.x[1L]; b = range.x[2L]

  kernel <- match.arg(tolower(kernel),
                      c("normal", "gaussian", "support"))
  if(is.null(error)) error = 'none'
  error <- match.arg(tolower(error),
                     c("none", "normal", "gaussian", "laplacian"))

  ##  in the current version, we support lp1 nop with Gaussian kernel
  ##  and no measurement errors (08/30/2012).
  if(kernel != "normal") kernel = "normal"
  if(error != "none") error = "none"
  if (!is.numeric(degree)) 
    stop("'degree' must be numeric")
  if(degree != 1L)     degree = 1L; # lp1 regression

  ## data prebinned over equally-spaced grid need to be re-binned
  ## differently.
  
  bwdisc = 25
    
  ## Rename common variables
  M <- gridsize
  Q <- as.integer(bwdisc)
  pp <- degree + 1L
  ppp <- 2L*degree + 1L
  tau <- 4

  truncate = TRUE  # pesudo r.v.
  
  if(binned){
    xcounts <- x
    ycounts <- y
    M <- length(xcounts)
    gpoints <- seq(a, b, length = M)
  }else{
    gpoints <- seq(a, b, length = M)
    out <- .rlbin(x, y, gpoints, truncate)
    xcounts <- out$xcounts
    ycounts <- out$ycounts
  }
  
  n = length(x)
  
  ## Set the bin width
  delta <- (b-a)/(M-1L)
  if(missing(bandwidth)){
    ## tmp = histo(x, weights=weights,binned=binned)
    ## bandwidth = dpill(tmp$x, tmp$y)
    bandwidth = bw.wnrd(x)
  }else{
    if (!is.numeric(bandwidth)) 
      stop("'bandwidth' must be numeric")
    stopifnot(length(bandwidth) == 1L)
  }

  ##  following we ingnored the case that bandwidth is a vector
  indic <- rep(1, M)
  Q <- 1L
  Lvec <- rep(floor(tau*bandwidth/delta), Q)
  hdisc <- rep(bandwidth, Q)
  
  stopifnot(min(Lvec) > 0)
  
  ## Allocate space for the kernel vector and final estimate
  
  dimfkap <- 2L * sum(Lvec) + Q
  fkap <- rep(0, dimfkap)
  curvest <- rep(0, M)
  midpts <- rep(0, Q)
  ss <- matrix(0, M, ppp)
  tt <- matrix(0, M, pp)
  Smat <- matrix(0, pp, pp)
  Tvec <- rep(0, pp)
  ipvt <- rep(0, pp)
  drv = 0L
  
  ## Call FORTRAN routine "locpol"
  
  out <- .Fortran(.F_locpol, as.double(xcounts), as.double(ycounts),
                  as.integer(drv), as.double(delta), as.double(hdisc),
                  as.integer(Lvec), as.integer(indic), as.integer(midpts),
                  as.integer(M), as.integer(Q), as.double(fkap), as.integer(pp),
                  as.integer(ppp), as.double(ss), as.double(tt),
                  as.double(Smat), as.double(Tvec), as.integer(ipvt),
                  as.double(curvest))
  
  curvest <- gamma(drv+1) * out[[19L]]
  curvest[curvest<0] = 0
  curvest = curvest/sum(curvest*delta)
  
  out = structure(list(x = gpoints, y = curvest,
    bw = bandwidth, sp = NULL, data = x,
    call = match.call(), data.name = name),
    class='kde')
  invisible(out)
}

