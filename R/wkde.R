
######   functions for birth defect data analysis

##  Objects: BD -- Sex, Head, Region
print.BD  <- function (x, digits = NULL, ...) 
{
  cat("Birth Defect Data:\n\n", sep = "")
  out1 <- summary(x$Head)
  out2 <- tapply(x$Head, x$Region, summary)
  out3 <- tapply(x$Head, x$Sex, summary)
  out <- rbind(out1,out2[[1]],out2[[2]],out3[[1]],out3[[2]])
  row.names(out) <- c("All", names(out2), names(out3))
  print(out)
  invisible(x)
}

wkde <- function(x, weights, freq, bw, from, to, gridsize, digits=0,
                 rounding) UseMethod("wkde")

wkde.default <- function(x, weights, freq, bw, from, to, gridsize, digits=0,
                         rounding){
  xr <- .data.rounding(x=x, freq=freq, weights=weights, digits=digits,
                       type=rounding)
  out <- wkde(xr, gridsize=gridsize, bw=bw)
  invisible(out)
}

wkde.bdata <- function(x, weights, freq, bw, from, to, gridsize, digits=0,
                       rounding){
#  xhist <- histogram(x, from=from,to=to)
    name <- x$data.name
  y <- x
  x <- x$x;
  n <- length(x);
  size <- sum(y$f)
  w <- y$w * y$f
  
  ## Set kernel support values    
  tau <-  4;
  h0 <- .bw.wtdnrd0(y);
  if(missing(gridsize)) gridsize <- 512L
  gridsize = max(512L, gridsize)  # make sure the grid is fine enough
  M <- 2^(ceiling(log2(gridsize)))
  if(missing(from)) from <- min(y$x+y$a)
  if(missing(to)) to <- max(y$x+y$b)
  stopifnot(to > from)
  ## Set up grid points and bin the data
  gpoints <- seq(from, to, length = M)
  
  ## Set default bandwidth and compute the estimate
  ## RR = 0; ##  no reflection used in this version 06/25/2012.  Features
  ##  can be added to the subroutines later.
  pars = c(h0,0)  ## for adaptive bandwidth selector: awkde
  
  ## Install safeguard against non-positive bandwidths:
  if(missing(bw)) bw <- "nrd0"
  if(is.numeric(bw)){
    stopifnot(bw > 0);
    h <- bw;
  } else if(is.character(bw)){
    bw <- match.arg(tolower(bw),
                    c('nrd', 'nrd0','lscv','mise','amise','erd'))
    h <- switch(bw,
                "nrd"   = .bw.wtdnrd(y),
                "nrd0"  = .bw.wtdnrd0(y),
                "erd"   = .bw.wtderd(y),
                "lscv"  = .bw.wtdlscv(y,h0,gridsize=M),
                "mise"  = .Fortran(.F_wkde,as.double(x),as.double(w),
                  as.integer(n),x=as.double(gpoints),y=as.double(gpoints),
                  Fy=as.double(gpoints),as.integer(M), bw=as.double(pars)),
                "amise" = .Fortran(.F_awkde,as.double(x),as.double(w),
                  as.integer(n),x=as.double(gpoints),y=as.double(gpoints),
                  Fy=as.double(gpoints),as.integer(M), bw=as.double(pars)),          
                stop(paste("Bandwidth selector '", bw, "'is not supported"))
                )
  }else stop("Invalid 'bandwidth'")

  if(bw=='mise'||bw=='amise'){
    fx <- h$y
    pars <- h$bw
  }else{
    pars <- c(h,0)
    fx <- .wkdefft(x=y,h=h,a=from,b=to,M=M)
  }
  
  totMass <- sum(fx) *(to-from)/(M-1)

  return(structure(list(y=fx/totMass, x=gpoints,
                        bw = pars,
                        ucb=NULL, lcb=NULL,
                        call = match.call(),
                        data.name = name,
#                        hist = xhist,
                        type=bw),
                   class='smooth'))
}

.wtdmean <- function(x){
  sum(x$x * x$w * x$f)
}

.wtdvar <- function(x){
  sum((x$x-.wtdmean(x))^2 * x$w * x$f);
}

.wtdsd <- function(x){
  sqrt(.wtdvar(x))
}

.wtdiqr <- function(x){
  x0 <- rep(x$x, x$f)
  as.numeric(diff(quantile(x0, c(0.25,0.75))))
}

.bw.wtderd <- function(x)
{
  size <- sum(x$f)
  if (size < 2L) 
    stop("need at least 2 data points")
  s <- .wtdmean(x)
  iqr <- .wtdiqr(x);
  0.9*min(s,iqr/1.34)*size^-.2
}


.bw.wtdnrd <- function(x)
{
  size <- sum(x$f)
  if (size < 2L) 
    stop("need at least 2 data points")
  s = .wtdsd(x); iqr = .wtdiqr(x);
  1.06*min(s,iqr/1.34)*size^-.2
}

.bw.wtdnrd0 <- function(x)
{
  size <- sum(x$f)
  if (size < 2L) 
    stop("need at least 2 data points")
  s = .wtdsd(x); iqr = .wtdiqr(x);
  0.90*min(s,iqr/1.34)*size^-.2
}

.bw.wtdlscv <- function(x,h0,gridsize=256)
{
  x0 <- x; # save a copy of the binned data
  x <- x0$x; n <- length(x); size <- sum(x0$f)
  w <- x0$f * x0$w 
  hl <- h0 * 0.5; hh <- 2.5*h0
  ## not needed if use optimize
  iter <- 50;   y <- rep(0,iter) # y stores the lscv's
  hs <- seq(hl,hh,length=iter)
  
  ## define grid and bin the weighted data
  a <- min(x); b <- max(x);
  M <- 2^(ceiling(log(gridsize)/log(2)))
  gpoints <- seq(a, b, length = M)
  binwidth <- (b-a)/(M-1)
##  gcounts <- .wbin(x, w, a, binwidth, M,linbin=TRUE,
##                   truncate=TRUE)$gcounts
##  gcounts <- gcounts * M/(b-a) 
  gcounts <- .data.binning(x0, type='kde', from=a, to=b,gridsize=M)
  gcounts <- gcounts * M/(b-a) #* size
  ##  compute the FT(gpoints) and |y_l|^2
  Myl <- .Fortran(.F_yldist, as.double(gcounts), as.integer(M),
                  y = double(M/2))$y
  ell <- 1:(M/2)
  sl <-  2*pi*ell/(b-a)
  lscv.temp <- function(h){
    tmp2 <- n*h*sqrt(2*pi)
    hs2 <- h^2*sl^2
    tmp <- exp(-hs2)-2.*exp(-0.5*hs2)
    sum(tmp*Myl)*(b-a)+1./tmp2
  }
  opt <- optimise(f = lscv.temp,
                  interval = c(0.25 * h0, 5 * h0,
                    tol = .Machine$double.eps))
  return(opt$minimum)
}

.wkdefft <- function(x, h, a, b, M)
{
  n <- length(x$x)
  ## Compute kernel weights
  delta  <- (b - a)/(h * (M-1L))
  if (M <5) warning("Binning grid too coarse!")
  gpoints <- seq(a, b, length = M)

  gcounts <- .data.binning(x, type='kde', from=a, to=b,gridsize=M)
  gcounts <- gcounts * sum(x$f)
  
  lvec <- 0L:M
  kappa <- dnorm(lvec*delta)/(n*h)
  
  ## Now combine weight and counts to obtain estimate 
  ## we need P >= 2L+1L, M: L <= M.
  P <- 2^(ceiling(log(2L*M+1L)/log(2)))
  kappa <- c(kappa, rep(0, P-2L*M-1L), rev(kappa[-1L]))
  tot <- sum(kappa) * (b-a)/(M-1L) * n # should have total weight one
  gcounts <- c(gcounts, rep(0L, P-M))
  kappa <- fft(kappa/tot)
  gcounts <- fft(gcounts)
  (Re(fft(kappa*gcounts, TRUE))/P)[1L:M]
}


##  The following programs help to prepare data to be analyzed with
##  the FMMBD functions.  In summary, data should have three columns:
##  x = center of the bins, widths = width of the bins, and counts =
##  counts (frequencies) of the bins.


###  function binning(...) 11/25/2011
##  Data could be rounded to the nearest integers (method='nearest'),
##  or rounded up (method='up'), or rounded down (method = 'down').
##  Some data may be recorded and rounded to the nearest 100mg.  In
##  that case, we can specify scale=100 (default scale=1).


##   digits: integer indicating the number of decimal places (‘round’)
##   or significant digits (‘signif’) to be used.  Negative values are
##   allowed.  Rounding to a negative number of digits means rounding
##   to a power of ten, so for example ‘round(x, digits = -2)’ rounds
##   to the nearest hundred.

.data.binning <- function(x, type='kde', from, to, gridsize=512L){
  ## x is an object produced by .rounding
  y <- x$x; f <- x$f; w <- x$w; b <- x$b; a <- x$a; n <- length(y);
  if(missing(from)) from <- min(y+a);
  if(missing(to)) to <- max(y+b);
  type <- match.arg(tolower(type), c('edf','hist','kde'))
  itype <- switch(type, edf=1, hist=-1, kde=0)

  .Fortran(.F_binning, as.double(y), as.double(f), as.double(w),
           as.double(a), as.double(b), as.integer(n), as.double(from),
           as.double(to), as.integer(gridsize), y=double(gridsize),
           as.integer(itype))$y

}

## this should be an internal function.  Users can use ceiling, round,
## floor or other functions to manupulate the data.
.data.rounding <- function(x, weights, freq, digits=0, type='none'){
  name <- deparse(substitute(x))
  stopifnot(is.numeric(digits))
  stopifnot(length(digits)==1)
  if(missing(type)) type <- 'none'
  method <- match.arg(tolower(type),
                      c('nearest', 'up','down','none'))

  stopifnot(is.numeric(x))
  n <- length(x)
  
  if(missing(freq)) freq <- rep(1, n)
  else{
    stopifnot(length(freq)==n)
    stopifnot(is.numeric(freq))
    if(any(is.na(freq)))
      stop("Missing value(s) in 'freq'")
    if(any(freq<=0))
      stop("Non-positive value(s) in 'freq'")
    if(any(freq != round(freq)))
      stop("Non-integer value(s) in 'freq'")
    if(any(!is.finite(freq)))
      stop("Infinite value(s) in 'freq'")
  }
  size <- sum(freq)
  
  if(missing(weights)) weights <- rep(1/size, n)
  else{
    stopifnot(length(weights)==n)
    stopifnot(is.numeric(weights))
    if(any(is.na(weights)))
      stop("Missing value(s) in 'weights'")
    if(any(weights<0.0)) stop("Invalid value(s) in 'weights'")
    if(sum(weights*freq) != 1)
      weights <- (weights * freq)/sum(weights * freq)
  }

  if(any(is.na(x))) stop("Missing value(s) in 'x'")

  x0 <- x*10^digits
  x0 <- switch(method,nearest=round(x0),up=ceiling(x0),
               down=floor(x0), x0)
  a <- switch(method,nearest=-0.5,up=-1, down=0, 0)
  b <- switch(method,nearest=0.5, up=0, down=1, 0)
  
  x0 <- x0*10^(-digits)
  a <- a*10^(-digits)
  b <- b*10^(-digits)

  x.finite <- is.finite(x0)
  x0 <- x0[x.finite]; totMass <- mean(x.finite)
  freq <- freq[x.finite]; weights <- weights[x.finite]
  w <- as.numeric(tapply(weights*freq, x0, sum))
  out <- tapply(freq, x0, sum)
  f <- as.numeric(out)
  x <- as.numeric(names(out))
  w <- w/f
  
  structure(list(x = x, f=f, w=w,
                 a=a, b=b,
                 type = method,
                 data.name = name),
            class='bdata')
}

print.bdata  <- function (x, ...) 
{
  cat("\nData: ", x$data.name, "\tRounding type: ", x$type, "\n", sep = "")
  out <- cbind(x=x$x,f=x$f,w=x$w)
  print(out)
  invisible(x)
}



