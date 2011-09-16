## modified based on bkde(...) in KernSmooth 2.23-8 on 07/02/2012 by
## Bin Wang.  We use gaussian kernels and provide three types of
## bandwidth selectors, namely bw.wnrd, bw.wnrd0 and bw.wmise.

## file KernSmooth/R/all.R
## original file Copyright (C) M. P. Wand
## modifications for use with R copyright (C) B. D. Ripley
## Unlimited use and distribution (see LICENCE).

wkde <- function(x, w, bandwidth, freq=FALSE, gridsize = 401L,
                 range.x, truncate = TRUE, na.rm=TRUE)
{
  name <- deparse(substitute(x))
  ## Install safeguard against non-positive bandwidths:
  if (!missing(bandwidth) && bandwidth <= 0)
    stop("'bandwidth' must be strictly positive")
  if(freq){
    stopifnot(all(w>0))
    n = sum(w)
  }
  ## Check data, and rename common variables
  y = .weighting(x,w,n,na.rm)
  n <- y$n;  x <- y$x; w <- y$w
  ngrid = y$nx
  ## Set kernel support values    
  tau <-  4; h0 <- .whrot0(y);
  gridsize = max(512L, gridsize)  # make sure the grid is fine enough
  M <- 2^(ceiling(log2(gridsize)))
  if (missing(range.x)) range.x <- c(min(x)-tau*h0, max(x)+tau*h0)
  a <- range.x[1L];  b <- range.x[2L]    
  ## Set up grid points and bin the data
  gpoints <- seq(a, b, length = M)
  
  ## Set default bandwidth and compute the estimate
  RR = 0; ##  no reflection used in this version 06/25/2012.  Features
  ##  can be added to the subroutines later.
  pars = c(h0,0)  ## for adaptive bandwidth selector: awkde
  out <- if(missing(bandwidth))
    .nwkde(x,w,h0,ngrid,a,b,M,tau=tau,truncate=truncate)
  else if(is.numeric(bandwidth))
    .nwkde(x,w,bandwidth,ngrid,a,b,M,tau=tau,truncate=truncate)
  else if(is.character(bandwidth))
    switch(tolower(bandwidth),
           "wnrd" = .nwkde(x,w,.whrot(y),ngrid,a,b,M,tau=tau,truncate=truncate),
           "wnrd0" =  .nwkde(x,w,h0,ngrid,a,b,M,tau=tau,truncate=truncate),
           "wmise" = .Fortran(.F_wkde,as.double(x),as.double(w),as.integer(ngrid),
             x=as.double(gpoints),y=as.double(gpoints),Fy=as.double(gpoints),
             as.integer(M), bw=as.double(pars),as.integer(RR)),
           "awmise" = .Fortran(.F_awkde,as.double(x),as.double(w),as.integer(ngrid),
             x=as.double(gpoints),y=as.double(gpoints),Fy=as.double(gpoints),
             as.integer(M), bw=as.double(pars),as.integer(RR)),          
           stop("Bandwidth selector not supported."))
  else stop("Bandwidth selector not supported.")
  if(length(out$bw)==1){
    bw = out$bw; sp=NULL;
  }else{
    bw = out$bw[1]; sp = out$bw[2]
  }
  out = structure(list(y=out$y,x=out$x, bw= bw, sp = sp,
    data = y, call = match.call(), data.name = name),
    class = "kde")
    
  invisible(out)
}

bw.wmise <- function(x, w, n=length(x), na.rm=TRUE)
{
  ## Check data, and rename common variables
  y = .weighting(x,w,n,na.rm)
  .wmise(y) # a class with "h", y=mise, x = a grid of bandwidths
}

bw.wnrd <- function(x, w, n=length(x), na.rm=TRUE)
{
  ## Check data, and rename common variables
  y = .weighting(x,w,n=n,na.rm)
  .whrot(y)
}

bw.wnrd0 <- function(x, w, n=length(x), na.rm=TRUE)
{
  ## Check data, and rename common variables
  y = .weighting(x,w,n=n,na.rm)
  .whrot0(y)
}

##  Function: .weighting   06/14/2012
### This function is used to create an R object "wdata".  Subroutines
### will be developed to summarize, visualize the "wdata" objects
### later. The weights "w" could be a vector of non-negative values.
### An observation with weight zero will be ignored automatically in
### the analysis. If 'w' is missing, create one using relative
### frequencies.  The number of missing values, total size, name
### should be kept.

### 07/03/2012: add wtdata1() and wtdata2().

.weighting <- function(x,w,n=length(x),na.rm=TRUE)
  {
    out = NULL
    if(missing(w)) out = .wtdata1(x,na.rm=na.rm)
    else out =  .wtdata2(x,w,n,na.rm=na.rm)
    out
  }

.wtdata1 <- function(x, na.rm=TRUE){
  name <- deparse(substitute(x))
  N = length(x)  ## size of the raw data
  x.na = is.na(x); n.na = sum(x.na)
  if(any(x.na)){
    if(na.rm) x = x[!x.na]
    else stop("Data contain missing value(s)")
  }
  x.finite <- is.finite(x)
  totMass = mean(x.finite)
  x = x[x.finite]
  n = sum(x.finite)
  tmp = table(x)
  x = as.numeric(names(tmp))
  w = as.numeric(tmp)
  w = w/sum(w)
  
  out = structure(
    list(x = x, w = w, n = n, nx = length(x),
         N = N, n.rm = 0, n.na = n.na, 
         totMass = totMass, data.name = name),
    class = "wdata")
  out
}

.wtdata2 <- function(x,w,n=length(x),na.rm=TRUE){
  name <- deparse(substitute(x))
  stopifnot(n>=length(x))
  stopifnot(length(w)==length(x))  # 'x' and 'w' have different lengths
  stopifnot(all(is.finite(w)))  # invalid infinite "w" 
  N = n; w = w/sum(w)
  x.na = is.na(x)|is.na(w); n.na=0
  if(any(x.na)){
    if(na.rm){
      x = x[!x.na];
      w = w[!x.na]; w=w/sum(w)
      n.na = sum(w[x.na])*N
    } else stop("Data contain missing value(s)")
  }
  x.rm = w <= 0  # observations with zero weight will be removed
                 # (ignored)
  n.rm = 0
  if(any(x.rm)){
    n.rm = sum(w[x.rm])*(N-n.na);
    x = x[!x.rm]; w = w[!x.rm]; w=w/sum(w)
  }
  
  x.finite <- is.finite(x)
  totMass = 1.0;  n.infinite = 0;
  if(any(!x.finite)){
    totMass = sum(w[x.finite])/sum(w)
    x = x[x.finite]; w = w[x.finite]; w = w/sum(w)
    n.infinite = sum(w[!x.finite])*(N-n.na-n.rm)
  }
  n = N - n.na - n.rm - n.infinite
  ox = order(x); x = x[ox]; w = w[ox]; 
  
  out = structure(
    list(x = x, w = w, n = n, nx = length(x),
         N = N, n.rm = n.rm, n.na = n.na, 
         totMass = totMass, data.name = name),
    class = "wdata")
  out
}

.wmean <- function(y){
  stopifnot(class(y)=="wdata")
  sum(y$x * y$w)
}

.wvar <- function(y){
  stopifnot(class(y)=="wdata")
  mu = sum(y$x * y$w)
  sum((y$x-mu)^2 * y$w)*(y$n)/(y$n-1)
}

.wsd <- function(y){
  sqrt(.wvar(y))
}


.wquantile <- function(y, level){
  stopifnot(class(y)=="wdata")
  Fw = cumsum(y$w)
  Fx = y$x
  stopifnot(all(!is.na(level)))
  stopifnot(all(level <= 1 & level >= 0))
  
  sele = level > min(Fw) & level < max(Fw)
  out = rep(0, length(level))
  if(sum(sele)>0){
    out[sele] = approx(Fw,Fx,level[sele])$y
  }
  if(sum(!sele)>0)
    out[!sele] = as.numeric(quantile(y$x,level[!sele], na.rm=TRUE))
  out  
}

.wiqr <- function(y){
  stopifnot(class(y)=="wdata")
  diff(.wquantile(y, level=c(.25,.75)))
}

## rule of thumb bandwidth selector

.whrot0 <- function(y, size)
{
  stopifnot(class(y)=="wdata")
  s = .wsd(y); iqr = .wiqr(y);
  n = ifelse(missing(size), y$n, size)
  0.9*min(s,iqr/1.34)*n^-.2
}

.whrot <- function(y, size)
{
  stopifnot(class(y)=="wdata")
  s = .wsd(y); iqr = .wiqr(y);
  n = ifelse(missing(size), y$n, size)
  1.06*min(s,iqr/1.34)*n^-.2
}

##  file: .wmise(y)
###  "y" is an R object of type wdata
.wmise <- function(y)
{
  stopifnot(class(y)=="wdata")
  h0 <- .whrot0(y)  # initialize the bandwidth for grid searching
  X <- y$x; W <- y$w; n <- y$nx; N = y$n
  m = 50; h = seq(h0/10, by = h0/10, length=m)
  out = .Fortran(.F_wmise, as.double(X), as.double(W), as.integer(n),
    h = as.double(h), mise = as.double(h), as.integer(m))
  sele = which(out$mise == min(out$mise))
  bw = out$h[sele[1]]
  list(bw = bw, x = out$h, y = out$mise)
}

## For application of linear binning to a univariate data set.
.wlinbin <- function(X, W, gpoints, truncate = TRUE)
{
    n <- length(X)
    M <- length(gpoints)
    if (truncate) trun <- 1L else 0L
    a <- gpoints[1L]
    b <- gpoints[M]
    .Fortran(.F_wlinbin, as.double(X), as.double(W), as.integer(n),
             as.double(a), as.double(b), as.integer(M),
             as.integer(trun), wts = double(M))$wts
}

.wedf <- function(y){
  stopifnot(class(y)=="wdata")
  out = structure(
    list(x = y$x, y = cumsum(y$w), data=y
         ), class = "edf")
}

.plot.edf  <- 
function (x, main = NULL, xlab = NULL, ylab = "Probability", 
          zero.line = TRUE, ...) 
{
  if (is.null(xlab)) 
    xlab <- paste("N =", x$data$n.total, "  NA =", x$data$n.rm)
  if (is.null(main)) 
    main <- ifelse(is.null(x$data$data.name),
                   "Empirical Distribution Function",
                   deparse(x$data$data.name))
  plot.default(x, main = main, xlab = xlab, ylab = ylab,
               ylim=c(0,1),pch=20,cex=.5,...)
  n = length(x$x)
  for(i in 1: (n-1))
    segments(x$x[i], x$y[i], x$x[i+1], x$y[i],...)
  if (zero.line) 
    abline(h = 0, lwd = 0.1, col = "gray")
  invisible(NULL)
}


.lines.edf  <- function (x, ...) 
{
  points(x,pch=20,cex=.5,...)
  n = length(x$x)
  for(i in 1: (n-1))
    segments(x$x[i], x$y[i], x$x[i+1], x$y[i],...)
  segments(x$x[n], x$y[n], x$x[n]*10, x$y[n],...)
  invisible(NULL)
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

.nwkde <- function(x,w, h, n, a,b, M,tau=4,truncate=FALSE)
{
    gpoints <- seq(a, b, length = M)
    gcounts <- .wlinbin(x, w*n, gpoints, truncate)
    ## Compute kernel weights
    delta  <- (b - a)/(h * (M-1L))
    L <- min(floor(tau/delta), M)
    if (L == 0)
      warning("Binning grid too coarse for current (small) bandwidth: consider increasing 'gridsize'")

    lvec <- 0L:L
    kappa <- dnorm(lvec*delta)/(n*h)
    
    ## Now combine weight and counts to obtain estimate 
    ## we need P >= 2L+1L, M: L <= M.
    P <- 2^(ceiling(log(M+L+1L)/log(2)))
    kappa <- c(kappa, rep(0, P-2L*L-1L), rev(kappa[-1L]))
    tot <- sum(kappa) * (b-a)/(M-1L) * n # should have total weight one
    gcounts <- c(gcounts, rep(0L, P-M))
    kappa <- fft(kappa/tot)
    gcounts <- fft(gcounts)
    list(x = gpoints, y = (Re(fft(kappa*gcounts, TRUE))/P)[1L:M], bw = h)
}
