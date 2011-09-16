##  Part of R package BDA
##  Copyright (C) 2009-2010 Bin Wang
##
##  Unlimited use and distribution (see LICENCE).

histolpr <- function(x, weights, bandwidth,  gridsize = 512L, range.x,
                     binned = FALSE, truncate = TRUE)
  
{
    name <- deparse(substitute(x))
    ## Install safeguard against non-positive bandwidths:
    if (!missing(bandwidth) && bandwidth <= 0)
        stop("'bandwidth' must be strictly positive")
    
    kernel = "normal"
    drv <- 0L;
    degree = 1L; # lp1 regression
    bwdisc = 25
    
    ##    if (missing(range.x) && !binned)
    if (missing(range.x))
      range.x <- c(min(x), max(x))

    ## Rename common variables
    M <- gridsize
    Q <- as.integer(bwdisc)
    a <- range.x[1L]
    b <- range.x[2L]
    pp <- degree + 1L
    ppp <- 2L*degree + 1L
    tau <- 4

    if(binned){
      if(missing(weights)){  # 'x' is prebinned/rounded, but not tabled
        tmp = table(x)
        x0 = as.numeric(names(tmp));
        y0 = as.numeric(tmp);
        N = length(x)
      }else{
        x0 = x; y0 = weights;
        N = sum(weights)
      }
      bw = diff(x0); h = bw[1]
      stopifnot(all(bw==h))  # stop if not equally binned
      ycounts = y0/(N*h)
      M = length(ycounts)
      xcounts <- rep(1, M)
      gpoints = x0
      width = h
    }else{
      if(missing(weights)) weights = rep(1, length(x))
      N <- length(x)
      gpoints <- seq(a, b, length = M)
      xcounts <- .wlinbin(x, weights, gpoints, truncate)
      width = (b-a)/(M-1)
      ycounts <- xcounts/(N*width)
      xcounts <- rep(1, M)
    }
    
    n = N
    ## Set the bin width
    delta <- (b-a)/(M-1L)
    if(missing(bandwidth)){
      ## tmp = histo(x, weights=weights,binned=binned)
      ## bandwidth = dpill(tmp$x, tmp$y)
      bandwidth = bw.wnrd(x, weights)
    }else{
      if (!is.numeric(bandwidth)) 
        stop("'bandwidth' must be numeric")
      stopifnot(length(bandwidth) == 1L)
    }

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
    curvest = curvest/sum(curvest*width)
    
    out = structure(list(x = gpoints, y = curvest,
      bw = bandwidth, sp = NULL, data = x,
      call = match.call(), data.name = name),
      class='kde')
    invisible(out)
}

