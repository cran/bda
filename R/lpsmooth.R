
npr <- function(y,x,bw,from,to,gridsize)
  {
    if(missing(from)) from <- min(x)
    if(missing(to)) to <- max(x)
    stopifnot(to > from)
    
    if(missing(bw)){
      lscv <- 1
      bw <- bw.nrd(x)
      adaptive = 1  # use adaptive bandwidth selector
    }else{
      stopifnot(bw>0)
      lscv <- 0
      adaptive = 0  # don't use adaptive bandwidth selector
    }

    ##  Compute the variance based on the raw data
    ox = order(x)
    oy = y[ox]; ox = sort(x)
    sy = sqrt(0.5* mean((diff(oy))^2))

    if(missing(gridsize)) gridsize <- 512L
    stopifnot(gridsize>10)
    gpoints <- seq(from, to, length=gridsize)
    n <- length(x)
    stopifnot(length(y) == n)
    if(any(is.na(x)|is.na(y)))
      stop("Missing value(s) in 'x' and/or 'y'")
    if(any(!is.finite(x)|!is.finite(y)))
      stop("Inifite value(s) in 'x' and/or 'y'")
    out <- .Fortran(.F_lpsmooth,
                    fx = as.double(gpoints), as.integer(gridsize),
                    as.double(x), as.double(y), as.integer(n),
                    bw = as.double(bw), as.integer(lscv),
                    as.double(c(from, to)), as.integer(adaptive),
                    ellx = double(gridsize), kappa=double(1))

    ## print(out$kappa)
    cv <-  .Fortran(.F_tubecv,
                    cv=as.double(out$kappa), as.double(0.95))$cv
    ## print(cv)
    
    y = out$fx
    sele1 = is.na(y) | !is.finite(out$fx)
    if(any(sele1)) y[sele1] = 0.0
    MOE <- cv * out$ellx * sy
    ll = y - MOE; ll[ll<0] <- 0
    ul = y + MOE; 
    
    structure(list(y = y, x = gpoints,
                   bw=out$bw, ucb=ul, lcb=ll,
                   call = match.call()
                   ),
              class = 'smooth')
  }

