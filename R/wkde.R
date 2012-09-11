### Functions:
## (a) wkde (b) lprde

wkde <- function(x, w, bandwidth, freq=FALSE, gridsize = 512L,
                 range.x, truncate = TRUE, na.rm=TRUE,...)
  UseMethod("wkde")

wkde.default <- function(x, w, bandwidth, freq=FALSE, gridsize = 512L,
                 range.x, truncate = TRUE, na.rm=TRUE,...)
{
  if(!missing(bandwidth)&&is.character(bandwidth)){
    bw <- match.arg(tolower(bandwidth),
                    c('nrd', 'nrd0','lscv','mise','amise','erd','amae','mae'))
    if(bw=='amae'||bw=='mae'){
      out <- .wkde.mae(x,w,gridsize=gridsize)
    }else{
      x.wt <- .weighting(x,w,freq=freq,na.rm=na.rm)
      out <-   .wkde(x.wt=x.wt, bandwidth=bandwidth,gridsize=gridsize,
                     range.x=range.x, truncate=truncate,...)
    }
  }else{
    x.wt <- .weighting(x,w,freq=freq,na.rm=na.rm)
    out <-   .wkde(x.wt=x.wt, bandwidth=bandwidth,gridsize=gridsize,
                   range.x=range.x, truncate=truncate,...)
  }
  invisible(out)
}

wkde.wtdata <- function(x, w, bandwidth, freq=FALSE, gridsize = 512L,
                 range.x, truncate = TRUE, na.rm=TRUE,...)
{
  if(class(x)!= 'wtdata') x <- .weighting(x,w,freq=freq,na.rm=na.rm)
  .wkde(x.wt=x, bandwidth=bandwidth,gridsize=gridsize,
        range.x=range.x, truncate=truncate,...)
}

.wkde.mae <- function(x,w,gridsize=512L)
  {
    if(any(x<0)) stop("Invalid lifetime data in 'x'")
    n <- length(x)
    if(missing(w)) w <- rep(1,n)
    if(any(w!=0&&w!=1)) stop("Invalid censoring status in 'w'")
    if(length(w)!=n) stop("'x' and 'w' have different lengths")
    gridsize = max(512L, gridsize)  # make sure the grid is fine enough
    M <- 2^(ceiling(log2(gridsize)))
    range.x <- c(min(x), max(x))
    a <- range.x[1L];  b <- range.x[2L]    
    ## Set up grid points and bin the data
    gpoints <- seq(a, b, length = M)
    ## sort the data
    ox <- order(x); x <- x[ox]; w <- w[ox]

    y <- .Fortran(.F_wkdemae,as.double(x),as.double(w),
                    as.integer(n), y=as.double(gpoints),
                    as.integer(M))$y
    totMass <- sum(y) *(b-a)/(M-1)
    
    x.wt <- .weighting(x,w,freq=FALSE,na.rm=TRUE)
    out = structure(list(y=y/totMass, x=gpoints,
      bw= NULL, sp = NULL,
      call = match.call(),
      data=x.wt, data.name = x.wt$name),
      class = "kde")
    invisible(out)

  }
                      
.wkde <- function(x.wt, bandwidth, gridsize = 512L,
                  range.x, truncate = TRUE,...)
{
  x <- x.wt$x; w <- x.wt$w; n <- x.wt$n; size <- x.wt$size

  ## Set kernel support values    
  tau <-  4;
  h0 <- .bw.wnrd0(x.wt);
  gridsize = max(512L, gridsize)  # make sure the grid is fine enough
  M <- 2^(ceiling(log2(gridsize)))
  if (missing(range.x)) range.x <- c(min(x)-tau*h0, max(x)+tau*h0)
  a <- range.x[1L];  b <- range.x[2L]    
  ## Set up grid points and bin the data
  gpoints <- seq(a, b, length = M)
  
  ## Set default bandwidth and compute the estimate
  ## RR = 0; ##  no reflection used in this version 06/25/2012.  Features
  ##  can be added to the subroutines later.
  pars = c(h0,0)  ## for adaptive bandwidth selector: awkde

  ## Install safeguard against non-positive bandwidths:
  if(missing(bandwidth)) bandwidth <- "MISE"
  else if(is.numeric(bandwidth)){
    stopifnot(bandwidth>0);
    h <- bandwidth;
  }
  else if(is.character(bandwidth)){
    bw <- match.arg(tolower(bandwidth),
                    c('nrd', 'nrd0','lscv','mise','amise','erd'))
    h <- switch(bw,
                "nrd"   = .bw.wnrd(x.wt),
                "nrd0"  = .bw.wnrd0(x.wt),
                "erd"   = .bw.erd(x.wt),
                "lscv"  = .bw.blscv(x.wt,h0,gridsize=M),
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
    bw <- h$bw[1]; sp <- h$bw[2]
    fx <- h$y
  }else if(bw=='amae'){
    bw <- NULL; sp <- NULL
    fx <- h$y
  }else{
    bw <- h; sp <- NULL;
    out <- .fftwkde(x,w,h,n,a,b,M,tau=tau,truncate=truncate)
    fx <- out$y
  }

  totMass <- sum(fx) *(b-a)/(M-1)

  out = structure(list(y=fx/totMass, x=gpoints, bw= bw, sp = sp,
    call = match.call(),
    data=x.wt, data.name = x.wt$name),
    class = "kde")
  invisible(out)
}


 

##  Oct 27, 2012: de via lpr

lprde <- function(x, w, bandwidth, nclass,  binwidth, lb, gridsize = 512L,
                  range.x, freq=FALSE, truncate = TRUE, na.rm=TRUE)
{
  x.bin <- .binning(x=x, w=w, nclass=nclass, binwidth=binwidth, lb=lb,
                    range.x=range.x, freq=freq, truncate = truncate)

  x.wt <- .weighting(x,w,freq=freq,na.rm=na.rm)

  h0 <- .bw.wnrd0(x.wt)
  if (missing(range.x))
    range.x <- c(min(x.wt$x) - 4. * h0, max(x.wt$x) + 4. * h0)
  a <- range.x[1L];  b <- range.x[2L]    
  gridsize = max(512L, gridsize)  # make sure the grid is fine enough
  M <- 2^(ceiling(log2(gridsize)))
  xgrid <- seq(a, b, length = M)
  
  ## Install safeguard against non-positive bandwidths:
  lscv <- FALSE
  if(missing(bandwidth)){
    bandwidth <- h0; lscv <- TRUE
  }else if(!is.numeric(bandwidth)){
    bandwidth <- h0; lscv <- TRUE
  }else if(bandwidth <= 0){
    stop("'bandwidth' must be strictly positive")
  }

  y1 <- x.bin$counts; x1 <- x.bin$mids; nbin <- x.bin$nclass
  y1 <- sqrt(nbin/x.wt$size*(y1 + 0.25)) # transform data
  
  tmp <- .lpreg(x=x1, y=y1, x0=xgrid, bandwidth = bandwidth,lscv)
  
  h <- tmp$bw; fy <- tmp$y;
  fy[fy<0] <- 0; fy <- fy^2
  totMass <- sum(fy) * (b-a)/M 
  structure(list(y=fy/totMass,x = xgrid,bw = h, sp=NULL,nclass=nbin,
                 hist=x.bin, call = match.call(),
                 data = x.wt, data.name = x.wt$name),
            class='nprde')
}

weighting <- function(x,w,freq=FALSE,na.rm=TRUE, type, method,...)
  UseMethod("weighting")

weighting.default <- function(x,w,freq=FALSE,na.rm=TRUE,type,
                              method, ...)
{
  ##  support plain weighted data and random censoring data
  ##  (2012/11/04)
  if(missing(type)) out <- .weighting(x=x,w=w,freq=freq,na.rm=na.rm)
  else{
    stopifnot(is.character(type))
    type = match.arg(tolower(type), c('rc'))
    out <- switch(type,
                  rc = .rcweighting(x=x,w=w,freq=freq,na.rm=na.rm,
                    method=method),
                  stop(paste("Data type '",type,"'not supported"))
                  )
  }
}


## 2012/10/30: to prepare data for analysis.  Check for NA and
## infinite values.  If weights are missing, tabulate the data in case
## there are duplicates. Even if the weights are given, it is still
## necessary to tabulate the data to reduce the computational burdens.

.weighting <- function(x,w,freq=FALSE,na.rm=TRUE){
  name <- deparse(substitute(x))
  if(missing(w)){
    xt <- table(x); x <- as.numeric(names(xt));
    w <- as.numeric(xt); freq <- TRUE
  }
  if(any(w<0)) stop("Negative weights not allowed")
  if(any(!is.finite(w))) stop("Infinite weights not allowed")
  x.na <- is.na(x) | is.na(w)
  if(any(x.na)){
    if(na.rm){
      x <- x[!x.na]; w <- w[!x.na]
    }
    else stop("'x' and/or 'w' contains missing value(s)")
  }
  ##  sample size should consider the Inf's (i.e. survival data
  ##  analysis)
  size <- ifelse(freq, sum(w), length(x))

  x.finite <- is.finite(x)
  totMass <- mean(x.finite)
  x <- x[x.finite]; w <- w[x.finite];
  out <- tapply(w,x,sum)
  x <- as.numeric(names(out));
  w <- as.numeric(out)
  n <- length(x);  
  
  structure(list(x = x, w = w/sum(w), n = n, size = size,
                 totMass = totMass, type=NULL, method=NULL,
                 pars=NA, data.name = name),
            class='wtdata')
}

## last updated: 2012/11/04

.rcweighting <- function(x,w,freq=FALSE,na.rm=TRUE,method="Nelson"){
  name <- deparse(substitute(x))
  if(missing(w)) stop("Censoring information 'w' is missing")
  if(any(!(w==1|w==0)))
    stop("Invalid censoring status in 'w'")
  x.na <- is.na(x) | is.na(w)
  if(any(x.na)){
    if(na.rm){
      x <- x[!x.na]; w <- w[!x.na]
    }
    else stop("'x' contains missing value(s)")
  }
  ##  sample size should consider the Inf's (i.e. survival data
  ##  analysis)
  size <- length(x)
  x.finite <- is.finite(x)
  totMass <- mean(x.finite)
  x <- x[x.finite]; w <- w[x.finite];
  lambda.hat <- sum(x)/sum(w)
  if(missing(method)) method <- "Nelson"
  tmp <- .Sx(x,w,method=method)
  w <- diff(rev(c(1,tmp$y)))
  x <- tmp$x
  out <- tapply(w,x,sum)
  x <- as.numeric(names(out));
  w <- as.numeric(out)
  n <- length(x);  
  totMass <- tmp$totMass
  structure(list(x = x, w = w/sum(w), n = n, size = size,
                 totMass = totMass, type="Random Right Censored",
                 method=tmp$method, pars=lambda.hat,
                 data.name = name),
            class='wtdata')
}

.Sx <- function(x,w,method="Nelson"){
  method <- match.arg(tolower(method), c('km', 'pl','product','kaplan',
                                         'na','aalen','nelson'))
  switch(method,
         km=, pl=, product=,
         kaplan = .Sx.KM(x,w),
         na=,aalen=, nelson=.Sx.NA(x,w),
         stop(paste("Method '", method, "'is not supported"))
         )
}

## x <- c(10,7,32,23,22,6,16,34,32,25,11,20,19,6,17,35,6,13,9,6,10)
## w <- c(1,1,0,1,1,1,1,0,0,0,0,0,0,1,0,0,1,1,0,0,0)

.Sx.KM <- function(x,w){
  x.censor <- w==0 # 0 censored, 1 events
  x.events <- x[!x.censor]
  n.events <- tapply(w[!x.censor], x.events, length)
  t.events <- as.numeric(names(n.events))
  n.events <- as.numeric(n.events)
  n.atrisk <-   apply(as.matrix(t.events,ncol=1),1,.KM.count,y=x)
  y <- cumprod(1-n.events/n.atrisk)
  totMass <- 1-((rev(n.atrisk - n.events))[1])/length(x)
  list(x=t.events,y=y,totMass=totMass,method="Kaplan-Meier")
}

.Sx.NA <- function(x,w){
  x.censor <- w==0 # 0 censored, 1 events
  x.events <- x[!x.censor]
  n.events <- tapply(w[!x.censor], x.events, length)
  t.events <- as.numeric(names(n.events))
  n.events <- as.numeric(n.events)
  n.atrisk <-   apply(as.matrix(t.events,ncol=1),1,.KM.count,y=x)
  y <- cumsum(n.events/n.atrisk)
  totMass <- 1-((rev(n.atrisk - n.events))[1])/length(x)
  list(x=t.events,y=exp(-y),totMass=totMass,method="Nelson-Aalen")
}

.KM.count <- function(x,y) sum(y>=x)

## internal functions for weighted data.  Will be called only by wkde
## functions.  2012/10/31

.wmean <- function(x,w){
  if(class(x) != 'wtdata')# stop("'x' not weighted")
    x <- .weighting(x,w)
  sum(x$x * x$w)
}

.wvar <- function(x,w){
#  if(class(x) != 'wtdata') stop("'x' not weighted")
  if(class(x) != 'wtdata')# stop("'x' not weighted")
    x <- .weighting(x,w)
  sum((x$x-.wmean(x))^2 * x$w);
  ##  sum((x$x-.wmean(x))^2 * x$w) * x$size/(x$size-1)
}

.wsd <- function(x,w){
#  if(class(x) != 'wtdata') stop("'x' not weighted")
  if(class(x) != 'wtdata')# stop("'x' not weighted")
    x <- .weighting(x,w)
  sqrt(.wvar(x))
}

## .wquantile should not be used to compute the quantiles with (very)
## low/high quantile levels, or for data with small sizes.

.wquantile <- function(x,w,level){
#  if(class(x) != 'wtdata') stop("'x' not weighted")
  if(class(x) != 'wtdata')# stop("'x' not weighted")
    x <- .weighting(x,w)
  w <- x$w; x <- x$x; Fw <- cumsum(w)
  if(any(level<0|level>1)) stop("Invalid quantile level(s)")
  out = rep(0, length(level))
  level.na <- is.na(level)
  out[level.na] <- NA
  sele = level > min(Fw) & level < max(Fw)
  if(sum(sele)>0)
    out[sele] = approx(Fw,x,level[sele])$y
  if(sum(!sele)>0)
    out[!sele] = as.numeric(quantile(x,level[!sele], na.rm=TRUE))
  out  
}

.wiqr <- function(x,w){
  if(class(x) == 'wtdata')
    out <- IQR(x$x)
  else out <- IQR(x)
  out
#  if(class(x) != 'wtdata')# stop("'x' not weighted")
#    x <- .weighting(x,w)
#  w <- x$w; x <- x$x; Fw <- cumsum(w)
#  level <- c(.25,.75)
#  sele = level > min(Fw) & level < max(Fw)
#  if(sum(sele)>0)
#    out[sele] = approx(Fw,x,level[sele])$y
#  if(sum(!sele)>0)
#    out[!sele] = as.numeric(quantile(x,level[!sele], na.rm=TRUE))
#  abs(diff(out))  
}

  
.lpreg <- function(x,y, x0, bandwidth, lscv=FALSE){
  lscv <- ifelse(lscv, 1, 0)
  n <- length(x)
  stopifnot(length(y) == n)
  m <- length(x0)
  if(missing(bandwidth)) stop("'bandwidth' missing")
  ## fitting llr with lp1: no measurement error, adaptive bandwidth
  ## selection using LSCV
  fhat  = .Fortran(.F_llrGauss, as.double(x), as.double(y),
    as.integer(n), y=as.double(x0), as.integer(m),
    bw=as.double(bandwidth), as.integer(lscv))
  list(x = x0,y = fhat$y, bw=fhat$bw)
} 


print.nprde  <- function (x, digits = NULL, ...) 
{
  cat("\nCall:\n\t", deparse(x$call), "\n\nData: ", x$data.name, 
      " (", x$n, " obs.); nclass = ", x$nclass,
      "\n\n", "Bandwidth: \t'bw' = ",
      formatC(x$bw,digits = digits), "\n\n", sep = "")
  print(summary(as.data.frame(x[c("x", "y")])), digits = digits, 
        ...)
  invisible(x)
}

plot.nprde  <- function (x, hist=TRUE,...) 
{
  if(hist){
    plot(x$hist,...)
    lines(x,...)
  }
  else plot(x,...)
  invisible(x)
}

print.kde  <- function (x, digits = NULL, ...) 
{
  cat("\nCall:\n\t", deparse(x$call), "\n\nData: ", x$data.name, 
      " (", x$n, " obs.);", "\n\n", "Bandwidth: \t'bw' = ",
      formatC(x$bw,digits = digits), "\n\n", sep = "")
  print(summary(as.data.frame(x[c("x", "y")])), digits = digits, 
        ...)
  invisible(x)
}

.bw.erd <- function(x)
{
  if(class(x) != 'wtdata')
    x <- .weighting(x)
  size <- x$size
  if (size < 2L) 
    stop("need at least 2 data points")
  s <- ifelse(is.na(x$pars), .wmean(x), x$pars);
  stopifnot(s>0)
  iqr <- .wiqr(x);
  0.9*min(s,iqr/1.34)*size^-.2
}


.bw.wnrd <- function(x)
{
  if(class(x) != 'wtdata')
    x <- .weighting(x)
  size <- x$size
  if (size < 2L) 
    stop("need at least 2 data points")
  s = .wsd(x); iqr = .wiqr(x);
  1.06*min(s,iqr/1.34)*size^-.2
}

.bw.wnrd0 <- function(x)
{
  if(class(x) != 'wtdata')
    x <- .weighting(x)
  size <- x$size
  if (size < 2L) 
    stop("need at least 2 data points")
  s = .wsd(x); iqr = .wiqr(x);
  0.9*min(s,iqr/1.34)*size^-.2
}

.bw.wmise <- function(x){
  if(class(x) != 'wtdata')
    x <- .weighting(x)
  h0 <- .bw.wnrd(x)
  
  w <- x$w; n <- x$n; x <- x$x;
  m = 50; h = seq(h0/10, by = h0/10, length=m)
  out = .Fortran(.F_wmise, as.double(x), as.double(w), as.integer(n),
    h = as.double(h), mise = as.double(h), as.integer(m))
  sele = which(out$mise == min(out$mise))
  bw = out$h[sele[1]]
  list(bw = bw, x = out$h, y = out$mise)
}

.bw.blscv <- function(x,h,gridsize=256)
{
  if(class(x) != "wtdata") x.wt <- .weighting(x)
  else x.wt <- x

  n <- x.wt$n; x <- x.wt$x; w <- x.wt$w
  h0 <- ifelse(missing(h), .bw.wnrd0(x.wt), h)
  hl <- h0 * 0.5; hh <- 2.5*h0
  ## not needed if use optimize
  iter <- 50;   y <- rep(0,iter) # y stores the lscv's
  hs <- seq(hl,hh,length=iter)
  
  ## define grid and bin the weighted data
  a <- min(x); b <- max(x);
  M <- 2^(ceiling(log(gridsize)/log(2)))
  gpoints <- seq(a, b, length = M)
  binwidth <- (b-a)/(M-1)
  gcounts <- .wbin(x, w, a, binwidth, M,linbin=TRUE,
                   truncate=TRUE)$gcounts
  gcounts <- gcounts * M/(b-a)
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





.wedf <- function(y){
  stopifnot(class(y)=="wdata")
  out = structure(
    list(x = y$x, y = cumsum(y$w), data=y
         ), class = "edf")
}

.fftwkde <- function(x,w, h, n, a,b, M, tau=4, truncate=FALSE)
{
  ## Compute kernel weights
  delta  <- (b - a)/(h * (M-1L))
  if (M == 0) warning("Binning grid too coarse!")
  gpoints <- seq(a, b, length = M)
  if(any(w>1)) w <- w/sum(w)
  binwidth <- (b-a)/(M-1)
  gcounts <- .wbin(x, w*n, a, binwidth, M, linbin = TRUE,
                   truncate=truncate)$gcounts
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
  list(x = gpoints, y = (Re(fft(kappa*gcounts, TRUE))/P)[1L:M], bw = h)
}

##  Biasing functions for simulation

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

