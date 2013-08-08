### Functions:
## (a) wkde (b) lprde

.wkde.mae <- function(x,w,range.x,gridsize=512L)
  {
    if(any(x<0)) stop("Invalid lifetime data in 'x'")
    n <- length(x)
    if(missing(w)) w <- rep(1,n)
    if(any(w!=0&&w!=1)) stop("Invalid censoring status in 'w'")
    if(length(w)!=n) stop("'x' and 'w' have different lengths")
    gridsize = max(512L, gridsize)  # make sure the grid is fine enough
    M <- 2^(ceiling(log2(gridsize)))
    if (missing(range.x)) range.x <- c(min(x), max(x))
    if(range.x[1L]<0) range.x[1L] <- 0.0  # life-time data analysis only
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

.rmwquantile <- function(x,w,level){
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

.loqntl <- function(x,w){
  ## x has been sorted
  wcum <- cumsum(w)
  n <- length(x)
  stopifnot(length(w) == n)
  l1 <- rev(which(wcum<=0.25))[1]
  l2 <- which(wcum > 0.25)[1]
  p1 <- 0.25 - wcum[l1]
  x[l1] + p1 *(x[l1+1]-x[l1])
}

.wiqr <- function(x,w){
  if(class(x) != 'wtdata') xwt <- .weighting(x,w)
  else xwt <- x
  q1 <- .loqntl(xwt$x,xwt$w)
  q3 <- .loqntl(rev(xwt$x),rev(xwt$w))
  q3 - q1
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
         ), class = "edf2")
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



bootkde <- function(x, method='z.score',scale=1, rounding = 'nearest',
                    from, to, alpha=0.05, gridsize=512L,na.rm=TRUE,iter=100)
{
  method <- match.arg(tolower(method),c("z.score","quantile","cdf"))
  name <- deparse(substitute(x))
  sele = is.na(x)
  if(any(sele)){
    if(na.rm) x = x[!sele]
    else stop("'x' contains missing values.")
  }
  n = length(x)
  out = .rounding(x, scale=scale, method=rounding)

  X = out$x; F=out$counts; B=out$width/2; A=-B; N=length(X)

  if(missing(alpha))alpha=0.05
  if(alpha<=0|alpha>=1)stop("Invalid confidence level.")
  alpha = min(alpha,1-alpha)

  ucb=NULL; lcb=NULL;
  ## Set default bandwidth
  h0 = bw.nrd(rep(X,F)+runif(sum(F),-B,B))
  h = .Fortran(.F_hbmise,as.double(X), as.double(F),as.double(2*B),
    as.integer(N), hopt=as.double(h0))$hopt

  if(missing(from)) from = min(X) - 3*h
  if(missing(to)) to = max(X) + 3*h
  a <- from;  b <- to
  ## set grid to evaluate the pdf/cdf
  M <- 2^(ceiling(log2(gridsize)))
  gpoints <- seq(a, b, length = M)

  
  y = .Fortran(.F_ofcpdf,
    as.double(X), as.double(F),as.double(-B),as.double(B),
    as.integer(N), y=as.double(gpoints), as.integer(M),
    para=as.double(h))$y
  y=cumsum(y)
  
  tmp = as.matrix(rep(sum(F),iter), ncol=1);
  ## simulate the rounding errors by a pilot estimate of F(x) based on a BME
  out = apply(tmp,1, .bootemp,x=gpoints,y=X,f=F,lb=A,ub=B,
    Fx=y,from=from,to=to)

  if(method=="z.score"){
    sigs = apply(out,1,sd);
    z0 = qnorm(1-alpha/2.);
    ym = apply(out,1,median);
    lcb = ym-sigs*z0; lcb[lcb<0]=0;
    ucb = ym+sigs*z0; ucb[ucb>1]=1
    y=out[,1];
    type="Density"
  }else if(method=="quantile"){
    ym = apply(out,1,median);
    lcb = as.numeric(apply(out,1,quantile, alpha/2)); lcb[lcb<0]=0;
    ucb = as.numeric(apply(out,1,quantile, 1-alpha/2)); ucb[ucb>1]=1
    y=out[,1];
    type="Density"
  }else{
    out = apply(out,2,cumsum)*(b-a)/M
    ym = apply(out,1,median);
    lcb = as.numeric(apply(out,1,quantile, alpha/2)); lcb[lcb<0]=0;
    ucb = as.numeric(apply(out,1,quantile, 1-alpha/2)); ucb[ucb>1]=1
    y=out[,1];
    type="Probability"
  }
  
  return(structure(list(y=y,x=gpoints,ym=ym,n = sum(F),
                        ucb=ucb,lcb=lcb,alpha=alpha,
                        call = match.call(), data.name = name
                        ), class = "bcb"))
}

.bootemp <- function(n,x,y,f,lb,ub,Fx,from,to){
  tmp = cbind(y,f,lb,ub);M=length(x);u=runif(n);N=length(f);

  smpl = .Fortran(.F_remp, as.integer(N), as.double(y), as.double(f),
    as.double(lb), as.double(ub), as.integer(M), as.double(Fx),
    as.double(x), smpl=as.double(u), as.integer(n))$smpl

  density(smpl,from=from,to=to)$y
}

plot.bcb <-
  function (x, main = NULL, xlab = "x", ylab = "Probability",
            lwd=1, col=1, lty=1, zero.line = TRUE,
            bgcol='gray',scb=FALSE, ...) 
{
  if (is.null(main)) main <- deparse(x$data.name)
  plot.default(x$x,x$ucb, main = main,
               xlab = xlab, ylab = ylab, type='n',...)
  if(scb){
    cord.x = c(x$x,rev(x$x))
    cord.y = c(x$lcb, rev(x$ucb))
    polygon(cord.x,cord.y,border=col,col=bgcol)
  }
  lines(x$x,x$y, lwd=lwd, col=col,lty=lty,...)
  
  if (zero.line) 
    abline(h = 0, lwd = 0.1, col = "gray")
  invisible(NULL)
}

lines.bcb <-
  function (x, lwd=1, col=1, lty=1, bgcol='gray', scb=FALSE, ...) 
{
  if(scb){
    cord.x = c(x$x,rev(x$x))
    cord.y = c(x$lcb, rev(x$ucb))
    polygon(cord.x,cord.y,border=col,col=bgcol)
  }
  lines(x$x,x$y, lwd=lwd, col=col,lty=lty,...)
  
  invisible(NULL)
}

print.bcb <- function (x, digits = NULL, ...) 
{
  cat("\nData: ", x$data.name, 
      " (", x$n, " obs.);\n", sep = "")
  print(summary(as.data.frame(x[c("x", "y","lcb","ucb")])),
        digits = digits, 
        ...)
  invisible(x)
}


##  Part of R package BDA
##  Copyright (C) 2009-2010 Bin Wang
##
##  Unlimited use and distribution (see LICENCE).

##  Part of R package BDA
##  Copyright (C) 2009-2010 Bin Wang
##
##  Unlimited use and distribution (see LICENCE).

.nkde <- function(x,x0,h)
{
    gpoints <- x0
    n = length(x)
    M <- length(gpoints)
    a <- gpoints[1L]
    b <- gpoints[M]
    if(missing(h))
      h <-   (1/(4*pi))^(1/10)*(243/(35*n))^(1/5)*sqrt(var(x))*15^(1/5)
    gcounts <- .Fortran(.F_linbin, as.double(x), as.integer(n),
                        as.double(a), as.double(b), as.integer(M), double(M))[[6]]
    ## Compute kernel weights
    delta  <- (b - a)/(h * (M-1L))
    L <- min(floor(4./delta), M)
    if (L == 0)
      warning("Binning grid too coarse for current (small) bandwidth: consider increasing 'gridsize'")

    lvec <- 0L:L
    kappa <-  dnorm(lvec*delta)/(n*h)

    ## Now combine weight and counts to obtain estimate
    ## we need P >= 2L+1L, M: L <= M.
    P <- 2^(ceiling(log(M+L+1L)/log(2)))
    kappa <- c(kappa, rep(0, P-2L*L-1L), rev(kappa[-1L]))
    tot <- sum(kappa) * (b-a)/(M-1L) * n # should have total weight one
    gcounts <- c(gcounts, rep(0L, P-M))
    kappa <- fft(kappa/tot)
    gcounts <- fft(gcounts)
    list(x = gpoints, y = (Re(fft(kappa*gcounts, TRUE))/P)[1L:M])
}


smkde <- function(x,f,bw,gridsize=512L,na.rm=TRUE,
                  just="center",binned=FALSE, scale=1.0){
  name <- deparse(substitute(x))
  if(scale <= 0) stop("Wrong 'scale' level.")
  x = x/scale
  ##  bin data if x is not binned.
  if(binned){
    if(missing(f))stop("Frequencies missing!")
    ijust = match.arg(tolower(just),
      c("center","left","right"))
    a = switch(ijust,center=0,left=0.5,right=-0.5)
    x1 = data.frame(x=x+a,f=f,b=0.5)
    X = x1$x; F = f;
  }else{
    if(any(is.na(x))) x = x[!is.na(x)]
    x0 = .discretize(x,na.rm=na.rm,just=just)
    X=x0$x; F=x0$f;
  }
  m = length(X);B=0.5;n=sum(F)
  ## define the grid where the densities will be evaluated
  range.x <- c(X[1] - 0.5, X[m] + 0.5)
  a <- range.x[1L];  b <- range.x[2L]
  ## set grid to evaluate the pdf/cdf
  M <- 2^(ceiling(log2(gridsize)))
  gpoints <- seq(a, b, length = M)
  ## to approximate the initial estimate of f(x) by simulation
  x0 = rep(X,F) + runif(n,-B,B)
  if(missing(bw))
    h <-  bw.SJ(x0)# (1/(4*pi))^(1/10)*(243/(35*n))^(1/5)*sqrt(var(x0))*15^(1/5)
  else h = bw
  
  f0 = .nkde(x0,gpoints,h)$y;

##  cv0 = -9999.0
##  hs = seq(h/10,3*h,length=100)
##  for(h0 in hs){
##    fx = .nkde(x0,gpoints,h0)$y;
##    cv = sum(fx^2)*(b-a)/M-2/(n-1)*(sum(fx)-1/sqrt(2*pi)/h0)
##    if(cv0<0 || cv < cv0){
##      f0 = fx; # initial estimate
##      cv0 = cv
##      h = h0
##    }
##  }
  iter=100;
  out <- .Fortran(.F_iterfx, fx=as.double(f0), as.double(gpoints), as.integer(M),
                  as.double(X), as.double(F), as.integer(m), as.double(B),
                  h = as.double(h),iter=as.integer(iter))
  if(out$iter>=iter) warning("Fixed point not found!")
  y = out$fx
  sele1 = is.na(y) | !is.finite(out$fx)
  y[sele1] = 0.0
  ##  list(out,jout,jfail)
  return(structure(list(y=y/scale,x=gpoints*scale,bw=out$h,scale=scale,
                        call = match.call(), data.name = name),
                   class='bde'))
}

print.bde <- function (x, digits = NULL, ...) 
{
  cat("\nCall:\t", deparse(x$call), "\n\n",sep = "")
  print(summary(as.data.frame(x[c("x", "y")])), digits = digits, 
        ...)
  invisible(x)
}

plot.bde  <- 
function (x, main = NULL, xlab = NULL, ylab = "Density", type = "l", 
    zero.line = TRUE, ...) 
{
  if (is.null(xlab)){
    if(!is.null(x$bw))
      xlab <- paste(deparse(x$call), "bw =",round(x$bw,3))
    else xlab <- 'X'
  }
  if (is.null(main)) 
    main <- deparse(x$data.name)
  plot.default(x, type = type, 
               main = main, xlab = xlab, ylab = ylab, ...)
  if (zero.line) 
    abline(h = 0, lwd = 0.1, col = "gray")
  invisible(NULL)
}


### The folloowing functions are to compute the empirical distribution
### function and construct a confidence band.  Last updated on
### 11/04/2012



.edf <- function(x, weights, xgrid, na.rm = TRUE){
  if (!is.numeric(x)) 
    stop("'x' must be numeric")
  name <- deparse(substitute(x))
  if(missing(weights)) weights <- rep(1, length(x))
  isna <- is.na(x)
  if(any(isna)){
    if(na.rm){
      x <- x[!isna]; weights <- weights[!is.na]
    }else stop("'x' contains missing value(s)")
  }
  n <- length(x); weights <- weights/sum(weights)
  o <- order(x); x <- x[o]; weights <- weights[o];

  Fwx <- cumsum(weights) 
  
  if(missing(xgrid))
    xgrid <- seq(min(x), max(x), length=256)
    
  structure(list(y = approx(x, Fwx, xgrid)$y,
                 x = xgrid, n = n, 
                 data.name = name), class = "myedf")
}

.cumprop <- function(x,y){
  mean(y<=x)
}

.recumprop <- function(x,y){
  x0 = round(x,0)
  sgnx = -1
  if(x0<x) sgnx = 1
  sele = x==x0
  mean(y<=x) - mean(sele)*(sgnx * 0.5 + x0 - x)
}



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

##  function histo(...).  Return an R object 'hist'.  Update the
##  object returned by histo(...).

.histogram <- function(x, w, nclass, binwidth, lb,
                      range.x, freq=FALSE, truncate = TRUE)
{
  x.bin <- .binning(x=x, w=w, nclass=nclass, binwidth=binwidth,
                    lb=lb, range.x=range.x, freq=freq,
                    truncate = truncate)
}



histospline <-
  function(x,f,gridsize=512L,na.rm=TRUE,just="center",
           binned=FALSE, scale=1.0){
    name <- deparse(substitute(x))
    if(scale <= 0) stop("Wrong 'scale' level.")
    x = x/scale
    ##  bin data if x is not binned.
    if(binned){
      if(missing(f))stop("Frequencies missing!")
      ijust = match.arg(tolower(just),
        c("center","left","right"))
      a = switch(ijust,center=0,left=0.5,right=-0.5)
      x1 = data.frame(x=x+a,f=f,b=0.5)
      X = x1$x; F = f;
    }else{
      if(any(is.na(x))) x = x[!is.na(x)]
      x0 = .discretize(x,na.rm=na.rm,just=just)
      X=x0$x; F=x0$f;
    }
    rf = F/sum(F)
    out = spline(X,rf,n=gridsize)
    f0 = out$y
    f0[f0<0]=0
    gpoints = out$x
    return(structure(list(y=f0/scale,x=gpoints*scale,bw=NULL,scale=scale,
                          call = match.call(), data.name = name),
                     class='bde'))
}


##  the following codes are to compute the optimal bin numbers




#####################################################################
## Construct histogram  ('histo')

.hist <- function(x,m){
  x=x[!is.na(x)]; n=length(x);  
  h = diff(range(x))/m
  x0 = seq(min(x),max(x),by=h)
  fhat = diff(c(.edf(x,xgrid=x0)$y,1))
  Jh = 2/h/(n-1)-(n+1)/h/(n-1)*sum(fhat^2)
  fhat = fhat/h
  list(Jh=Jh,y=fhat,x=x0,m=m)
}

.binhist <- function(x,just="center",scale=1.0){
  x=x[!is.na(x)]; x = round(x,0);
  if(scale <= 0) stop("Invalid scale.")
  if(scale != 1.0) x = x/scale
  n=length(x); h = 1
  ijust = match.arg(tolower(just),
    c("center","left","right"))
  a = switch(ijust,center=-0.5,left=0,right=-1)
  x0 = seq(min(x)+a, max(x)+1.0+a, by=h)
  m = length(x0)-1
  fhat = diff(c(.edf(x,xgrid=x0)$y,1))
  Jh = 2/h/(n-1)-(n+1)/h/(n-1)*sum(fhat^2)
  fhat = fhat/h
  list(Jh=Jh,y=fhat,x=x0,m=m)
}

.binwidth <- function(x){
  n = length(x)/5
  Jh=NULL; m=NULL
  for(i in 5:n){
    out = .hist(x,i)
    Jh = c(Jh,out$Jh)
    m = c(m,i)
  }
  list(Jh=Jh,m=m)
}

.histo <- function(x, m=NULL, alpha=0.05, binned=FALSE,just="center",scale=1.0){
  name <- deparse(substitute(x))
  if(!binned){
    out = .binhist(x,just,scale)
    Jh = out$Jh
    m = out$m
    ms = m
    Jhs=Jh;
  }else{
    if(is.null(m)){
      out = .binwidth(x)
      Jhs = out$Jh
      ms = out$m
      sele = which(Jhs==min(Jhs))[1]
      m = ms[sele]
      Jh = Jhs[sele]
      out = .hist(x,m)
    }else{
      out = .hist(x,m)
      Jh = out$Jh
      ms = m
      Jhs=Jh;
    }
  }
  n = length(x)
  fhat = out$y
  if(alpha>1|alpha<0)stop("Invalid confidence level.")
  if(alpha>0.5) alpha=1-alpha
  epsn = 0.5 * qnorm(1 - 0.5 * alpha/m)*sqrt(m/n)
  LFn = sqrt(fhat) - epsn;
  LFn[LFn<0] = 0; LFn=LFn^2
  UFn = (sqrt(fhat) + epsn)^2

 out = structure(list(y=fhat,x=out$x, l=LFn,u=UFn,
   Jhs = Jhs, Jh=Jh, ms=ms, m=m, n = n,
   data = x, plot=NULL,
   data.name = name), class = "hist")
  out
}

## pre-bin the (weighted) data to a sequence of equally-spaced bins

.binning <- function(x, w, nclass, binwidth, lb, range.x,
                     freq=FALSE, truncate = TRUE){
  x.wt <- .weighting(x,w,freq=freq,na.rm=TRUE)
  if (missing(range.x))
    range.x <- c(min(x.wt$x), max(x.wt$x))
  a <- range.x[1L];  b <- range.x[2L]    
  if(!missing(binwidth)&&missing(lb)) lb <- a-binwidth/2
  else if(missing(lb)) lb <- a;
  if(missing(binwidth)){
    if(missing(nclass)) nclass <- "Sturges"
    if (is.character(nclass)){ 
      nclass <- match.arg(tolower(nclass),
                          c("sturges","fd", "freedman-diaconis",
                            "scott", "lscv"))
      nclass <- switch(nclass,
                       sturges = .nclass.Sturges(x.wt$size), 
                       `freedman-diaconis` = ,
                       fd = .nclass.FD(x.wt),
                       scott = .nclass.scott(x.wt),
                       lscv = .nclass.lscv(x.wt),
                       stop("unknown 'breaks' algorithm"))
    } else if (is.function(nclass)) {
      nclass <- nclass(x.wt)
    }
    if (!is.numeric(nclass) || !is.finite(nclass) || nclass < 1) 
      stop("invalid number of 'nclass'")
    
    binwidth <- (b-a)/nclass
  }else{ # if binwidth is specified, nclass need to be recalculated
    nclass <- ceiling((b-a)/binwidth)
  }

  out <- .wbin(X=x.wt$x, W=x.wt$w, a=lb, bw=binwidth, ngrid=nclass+1,
               linbin = FALSE, truncate = truncate)

  x.hist <- structure(list(breaks = out$breaks,
                           counts = out$gcounts * x.wt$size,
                           intensities = out$gcounts/binwidth,
                           density = out$gcounts/binwidth,
                           mids = out$gpoints,
                           bw=binwidth, nclass=nclass,
                           xname = x.wt$name, equidist=FALSE),
                      class="histogram")  
}

.nclass.scott <- function (x) 
{
  if(class(x) != "wtdata") x.wt <- .weighting(x)
  else x.wt <- x
  size <- x$size
  h <- 3.5 * .wsd(x) * size^(-1/3)
  if (h > 0) 
    ceiling(diff(range(x$x))/h)
  else 1L
}

.nclass.FD <- function (x)  
{
  x <- ifelse(class(x) != 'wtdata',.weighting(x),x)
  size <- x$size
  h <- .wiqr(x)
  if (h == 0){
    cf <- abs(cumsum(x$w)-0.5)
    h <- x$x[which(cf==min(cf))[1]]
  }
  if (h > 0) 
    ceiling(diff(range(x$x))/(2 * h * size^(-1/3)))
  else 1L
}

.nclass.Sturges <- function (n) 
  ceiling(log2(n) + 1)

.nclass.lscv <- function(x){
  x.wt <- ifelse(class(x) != 'wtdata',.weighting(x),x)
  n <- x$size; x <- x.wt$x; w <- x.wt$w
  xrange = diff(range(x))
  if(n<=5) nclass = 2
  else{
    n1 = min(round(n * 0.2),5)
    n2 = min(round(n * 0.5),30)
    Jh=NULL; m=NULL
    for(i in n1:n2){
      h = xrange/i
      Jh = c(Jh,.Jh(x,w,h));
      m = c(m,i)
    }
    nclass = m[which(Jh==min(Jh))[1]]
  }
  nclass
}

.Jh <- function(x,w,h){
  n = length(x)
  x0 = seq(min(x),max(x),by=h)
  fhat = diff(c(.edf(x,weights=w, xgrid=x0)$y,1))
  Jh = 2/h/(n-1)-(n+1)/h/(n-1)*sum(fhat^2)
}

## For application of linear binning to a univariate data set.
.wbin <- function(X, W, a, bw, ngrid, linbin=TRUE, truncate = TRUE)
{
  n <- length(X)
  if(missing(W)) W <- rep(1,n)
  stopifnot(length(W)==n)
  if(missing(a)) a <- min(X)
  stopifnot(bw>0)
  if(truncate) trun <- 1L else trun <- 0L
  if(linbin) lbin <- 1L else lbin <- 0L
  
  gcounts <- .Fortran(.F_GridBinning, as.double(X), as.double(W),
                      as.integer(n), as.double(a), as.double(bw),
                      as.integer(ngrid), as.integer(trun),
                      as.integer(lbin), wts = double(ngrid))$wts
  xgrid <- seq(from=a, by=bw, length=ngrid)
  gpoints <- xgrid[-1] - bw*0.5
  if(!linbin) gcounts <- gcounts[-ngrid];
  
  list(breaks=xgrid, gcounts=gcounts,gpoints=gpoints)
}


##  2012/11/01

bfmm <- function(x,m=2,mu,type='gaussian',method='nelder',range.x, ...)
UseMethod("bfmm")

bfmm.default <- function(x,m=2,mu,type='gaussian',method='nelder',range.x, ...)
{  
  stopifnot(class(x) == 'bdata')
  method = match.arg(tolower(method), c('nelder', 'newton'))
  switch(method,
         nelder = .bfmm(x,m=m,mu=mu,type=type,range.x=range.x,...),
         newton = .binem(x,m=m,mu=mu,...)
         )
}

bfmm.data.frame <- function(x,m=2,mu,type='gaussian',method='nelder',range.x,...){
  name <- deparse(substitute(x))
  if(is.null(x$x)||is.null(x$widths)||is.null(x$counts))
    stop("Information missing in data 'x'. Try 'help(rounding)'.")
  
  out = structure(list(x = x$x, widths = x$widths, counts = x$counts,
    scale = 1.0, data.name = name), class='bdata')
  bfmm(out,m=m,mu=mu,type=type,method=method,...)
}

bfmm.numeric <- function(x,m=2,mu,type='gaussian',method='nelder',
                         range.x,scale=1, rounding = 'nearest',...){
  out = .rounding(x, scale=scale, method=rounding)
  bfmm(out,m=m,mu=mu,type=type,method=method,...)  
}

print.mm <- function(x,...)
  {
    if(!is.null(x$range.x))
      cat("\nData were rescaled from [",
          x$range.x[1], ",", x$range.x[2], "] to [0,1].",sep = "")
    
    cat("\nFitted a mixed model of ( ", x$m, " ) ", x$type,
        " component(s).\n", sep = "")
    print(x$llk)
    cat("\n\tParamters:\n")
    print(x$para)
    cat("\n\nUse density(x) to compute the fitted density function.\n\n")
    invisible(x)

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

.rounding <- function(x, scale=1, method='nearest'){
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
  mu = .wtmean(y,f);  s = .wtsd(y, f);
  structure(list(x = y, widths = w, counts = f,mu=mu,s=s,
                 scale = scale, data.name = name),
            class='bdata')
}

##########################################################################
.bfmm <- function(x,m=2,mu,type='gaussian',range.x, ...)
{  
  type = match.arg(tolower(type),
    c('gaussian', 'normal','beta','gamma','weibull'))
  idist = switch(type,
    gaussian = 0,
    normal = 0,
    beta = 1,
    gamma = 2,
    weibull = 3
    )
  
  n = length(x$x)
  if (missing(range.x))
    range.x <- c(x$x[1] - 0.5 * x$widths[1], x$x[n] + 0.5 * x$widths[n])
  xmin <- range.x[1L]
  xmax <- range.x[2L]
  xspan = xmax - xmin;

  x0 = (x$x - xmin)/xspan;
  x0wd = x$widths/xspan

  llk = 0 # compute llk.
  if(missing(mu)){
    m = round(m)
    if(m<1) stop("Invalid number of components.")
    if(n < m)
      warning("Number of mixing components is larger than the distinct bin centers!")
    m = min(m, n - 2)
    type = ifelse(m>1,'normal',type)
    idist = ifelse(m>1,0,idist)
    par = rep(0,3*m)  ## works only for two parameter distribution families.
    Iter = min(100, 2*choose(n-2,m))
    x0s = x0[-c(1,n)]
    par[3*(1:m)-1] = sample(x0s, m)
    out = .Fortran(.F_fitmm,
      as.double(x0), as.double(x$counts), as.double(x0wd),
      as.integer(n), as.integer(idist), 
      as.integer(m), par = as.double(par), llk=as.double(llk)
      )
    for(i in 1:Iter){
      par[3*(1:m)-1] = sample(x0s, m)
      out2 = .Fortran(.F_fitmm,
        as.double(x0), as.double(x$counts), as.double(x0wd),
        as.integer(n), as.integer(idist), 
        as.integer(m), par = as.double(par), llk=as.double(llk)
        )
      if(out2$llk[1] > out$llk[1]) out = out2
    }
  }else{
    m = length(mu);
    type = ifelse(m>1,'normal',type)
    idist = ifelse(m>1,0,idist)
    par = rep(0,3*m)  ## works only for two parameter distribution families.
    mu = (mu-xmin)/xspan
    par[3*(1:m)-1] = mu

    out =  .Fortran(.F_fitmm,
      as.double(x0), as.double(x$counts), as.double(x0wd),
      as.integer(n), as.integer(idist), 
      as.integer(m), par = as.double(par), llk=as.double(llk))
  }

  pars =  matrix(out$par, nrow=3, ncol=m)
  sele = pars[1,] > 0.00001
  pars = pars[,sele]
  pars = as.data.frame(pars)
  row.names(pars) <- c("Proportion","Para.1","Para.2")
  m = sum(sele)
  
  llk0 = out$llk
  K = ifelse(m == 1, 2., 3. * m - 1.);
  N = sum(x$counts);
  llktransform = 2.0 * sum(x$counts) * log(xspan)
  AIC = -2.0 * llk0 + 2.0 * K + llktransform;
  AICc = AIC + 2.0 * K * (K+1.) / (N - K - 1.0);
  BIC = -2.0 * llk0 + log(N) * K + llktransform;
  structure(list(para = pars, type=type, m=m,
                 range.x=c(xmin,xmax), data = x,
                 llk = c("AIC"=AIC, "BIC"=BIC,"AICc"=AICc)),
            class='mm')
}


density.mm <- function(x,x0,gridsize=500,...)
  {
    if(is.null(x$range.x)){
      xmin = min(x$data$x); xmax = max(x$data$x);
      xspan = xmax - xmin;
      if(missing(x0)){
        x.ext = 1.5 * x$data$s; 
        x0 = seq(xmin-x.ext,xmax+x.ext, length=gridsize)
      }else{
        gridsize = length(x0)
      }
      xpts = x0
      fx = rep(0,gridsize)
      for(i in 1:x$m){
        fx0 = dnorm(x0,x$para[2,i],x$para[3,i])
        fx = fx + x$para[1,i] * fx0;
      }
    }else{
      xmin = x$range.x[1]; xmax = x$range.x[2];
      xspan = xmax - xmin;
      if(missing(x0)){
        s = x$data$s; z = x$data$x
        x.ext = 1.5*s/(max(z)-min(z)+1) 
        x0 = switch(x$type,
          gaussian = seq(0-x.ext,1+x.ext, length=gridsize),
          normal = seq(0,1+x.ext, length=gridsize),
          beta = seq(0,1, length=gridsize),
          gamma = seq(0,1+x.ext, length=gridsize),
          weibull = seq(0,1+x.ext, length=gridsize)
          )
        xpts = x0*xspan + xmin
      }else{
        gridsize = length(x0)
        xpts = x0 #save raw x0
        x0 = (x0 - xmin)/xspan
      }
      fx = rep(0,gridsize)
      for(i in 1:x$m){
        fx0 = switch(x$type,
          gaussian = dnorm(x0,x$para[2,i],x$para[3,i]),
          normal = dnorm(x0,x$para[2,i],x$para[3,i]),
          beta = dbeta(x0,x$para[2,i],x$para[3,i]),
          gamma = dgamma(x0,shape=x$para[2,i],scale=x$para[3,i]),
          weibull = dweibull(x0,shape=x$para[2,i],scale=x$para[3,i])
          )
        fx = fx + x$para[1,i] * fx0;
      }
      fx = fx/xspan
    }
    list(x = xpts, y = fx)
  }

##  Part of R package BDA
##  Copyright (C) 2009-2010 Bin Wang
##
##  Unlimited use and distribution (see LICENCE).

.discretize <- function(x, na.rm=TRUE,just="center"){
  if(na.rm){
    x=x[!is.na(x)];
    if(length(x)<1)stop("Invalid data!");
  }else{
    if(any(is.na(x)))stop("Data contain missing value(s)!");
  }
  ijust = match.arg(tolower(just),
    c("center","left","right"))
  a = switch(ijust,center=0,left=0.5,right=-0.5)
  x = round(x,0)
  y = table(x);
  x = as.numeric(names(y));
  N = length(x);
  ul = rep(.5,N);
  n=as.numeric(y);
  data.frame(x=x+a,f=n,b=ul);
}


.wtmean <- function(x,f){
  if(length(x)!=length(f))stop("'x' and 'f' have different lengths.")
  sele1 = !is.na(x)
  sele2 = !is.na(f)
  sele = sele1&sele2
  x=x[sele]
  f=f[sele]
  sum(x*f)/sum(f)
}

.wtsd <- function(x,f){
  if(length(x)!=length(f))stop("'x' and 'f' have different lengths.")
  sele1 = !is.na(x)
  sele2 = !is.na(f)
  sele = sele1&sele2
  x=x[sele]
  f=f[sele]
  mu = sum(x*f)/sum(f)
  sqrt(sum(f*(x-mu)^2)/(sum(f)-1))
}

.renorm <- function(X,F,B){
  theta = c(.wtmean(X,F),.wtsd(X,F));
  N=length(X)
  if(length(B)==1) B = rep(B,N)
  .Fortran(.F_remlenorm,
           as.double(X), as.double(F),as.double(B),
           iter=as.integer(N),theta=as.double(theta))$theta
}


.emnorm <- function(X,F,B,N,m,mu,s){
  p = rep(1/m,m);  sigs = rep(s,m);  N = length(X)
  .Fortran(.F_reemnorm,
           as.double(X), as.double(F),as.double(B/2),as.integer(N),
           iter=as.integer(m), p=as.double(p),mu=as.double(mu),
           sig=as.double(sigs),llk=as.double(0))
}


## Input: x which are the OFC data. m=2, number of components.
  
## Output: AIC/AICc/BIC, and the corresponding parameter settings.

##  There are only a few distinct OFC values after rounding.  Here
##  we first round x to the nearest centimeters, and get a list of
##  distinct values: x1.  Pick the one with the largest frequency
##  and m-1 from the rest.  Add U~runif(0,1)-.5 to the selected
##  value.

.binem <- function(x,m=2,mu,...){

  n = length(x$x); sx = x$s
  
  if(missing(mu)){
    m = round(m)
    Iter = min(100, 2*choose(n-2,m))
  }else{
    m = length(mu)
    Iter = 1
  }
  if(m<1) stop("Invalid number of components.")
  if(n < m)
    warning("Number of mixing components is larger than the distinct bin centers!")
  m = min(m,n - 2)
  llk0 = 0;  
  if(m == 1){
    theta = .renorm(x$x,x$counts,0.5)  # two-parameter estimate
    pars = matrix(c(1,theta), ncol=1,nrow=3)
    fx = dnorm(x$x, theta[1], theta[2])
    sele = fx > 0
    llk0 = sum(x$counts[sele] * log(fx[sele]))
  }else{
    if(Iter==1){
      out = .emnorm(x$x,x$counts,x$widths,n,m,mu,sx)
      pars =  rbind(out$p,out$mu,out$sig)
      sele = pars[1,] > 0.00001
      pars = pars[,sele]
      m = sum(sele)
      llk0 = out$llk
    }else{
      x0s = x$x[-c(1,n)]
      mu = sample(x0s, m)
      out = .emnorm(x$x,x$counts,x$widths,n,m,mu,sx)
      for(i in 1:Iter){
        mu = sample(x0s, m)
        out2 = .emnorm(x$x,x$counts,x$widths,n,m,mu,sx)
        if(out2$llk > out$llk) out = out2
      }
      pars =  rbind(out$p,out$mu,out$sig)
      sele = pars[1,] > 0.00001
      pars = pars[,sele]
      m = sum(sele)
      llk0 = out$llk
    }
  }
  
  pars = as.data.frame(pars)
  row.names(pars) <- c("Proportion","Para.1","Para.2")
  
  K = ifelse(m == 1, 2., 3. * m - 1.);
  N = sum(x$counts);
  AIC = -2.0 * llk0 + 2.0 * K;
  AICc = AIC + 2.0 * K * (K+1.) / (N - K - 1.0);
  BIC = -2.0 * llk0 + log(N) * K;
  structure(list(para = pars, type = 'normal', m=m,
                 range.x = NULL, data = x,
                 llk = c("AIC"=AIC, "BIC"=BIC,"AICc"=AICc)),
            class='mm')
}


##  Part of R package meada
##  Copyright (C) 2009-2010 Bin Wang
##
##  Unlimited use and distribution (see LICENCE).

##  2012/10/31: define four functions for a mixture of k normal
##  components f(x) = \sum p_i * dnorm(x,mu_i, sigma_i)

.dmnorm <- function(x,p,mu,s){
  k <- length(p)
  res <- 0
  for(i in 1:k){
    res <- res + p[i] * dnorm(x,mu[i],s[i])
  }
  res
}

dmixnorm <- function(x,p,mu,s){
  if(missing(p)) p <- 1
  if(missing(mu)) mu <- 0
  if(missing(s)) s <- 1
  ndim <- length(p)
  if(length(mu) != ndim | length(s) != ndim)
    stop("Parameters have different lengths")
  if(any(p>1|p<0) | sum(p) != 1) stop("Wrong mixing coefficients")
  if(any(s<=0)) stop("Invalid standard deviation(s)")
  
  sapply(x,.dmnorm,p=p,mu=mu,s=s)
}

.pmnorm <- function(x,p,mu,s){
  k <- length(p)
  res <- 0
  for(i in 1:k){
    res <- res + p[i] * pnorm(x,mu[i],s[i])
  }
  res
}

pmixnorm <- function(q,p,mu,s){
  if(missing(p)) p <- 1
  if(missing(mu)) mu <- 0
  if(missing(s)) s <- 1
  ndim <- length(p)
  if(length(mu) != ndim | length(s) != ndim)
    stop("Parameters have different lengths")
  if(any(p>1|p<0) | sum(p) != 1) stop("Wrong mixing coefficients")
  if(any(s<=0)) stop("Invalid standard deviation(s)")
  
  sapply(q,.pmnorm,p=p,mu=mu,s=s)
}

.rmnorm <- function(x,p){
  k <- length(p)
  cump <- cumsum(p)
  x[which(runif(1)-cump<=0)[1]]
}

rmixnorm <- function(n,p,mu,s){
  if(missing(p)) p <- 1
  if(missing(mu)) mu <- 0
  if(missing(s)) s <- 1
  ndim <- length(p)
  if(length(mu) != ndim | length(s) != ndim)
    stop("Parameters have different lengths")
  if(any(p>1|p<0) | sum(p) != 1) stop("Wrong mixing coefficients")
  if(any(s<=0)) stop("Invalid standard deviation(s)")
  n <- ceiling(n)
  stopifnot(n>0)
  
  tmp <- NULL
  k <- length(p)
  for(i in 1:k){
    tmp <- cbind(tmp, rnorm(n,mu[i], s[i]))
  }
  res <- apply(tmp,1,.rmnorm,p=p)
  as.numeric(res)
}

qmixnorm <- function(prob,p,mu,s){
  if(missing(p)) p <- 1
  if(missing(mu)) mu <- 0
  if(missing(s)) s <- 1
  sele <- !is.na(prob)
  if(any(prob[sele]>1|prob[sele]<0))
    stop("Invalid 'prob' value(s)")

  ndim <- length(p)
  if(length(mu) != ndim | length(s) != ndim)
    stop("Parameters have different lengths")
  if(any(p>1|p<0) | sum(p) != 1) stop("Wrong mixing coefficients")
  if(any(s<=0)) stop("Invalid standard deviation(s)")

  mu.pool <- sum(p*mu)
  s.pool <- sqrt(sum(p^2*s^2))  
  x <- seq(mu.pool-4*s.pool,mu.pool+4*s.pool,length=401L)
  Fx <- sapply(x,.pmnorm,p=p,mu=mu,s=s)
  approx(Fx,x,prob)$y
}
## This function is developed to compute the MLEs for different
## families of distributions based on data that are weighted, rounded
## or censored.

##  2012/11/06.
##  
mle <- function(x,w, type="Weighted", family="Gaussian"){
  if(!is.character(family)) family <- 'normal'
  if(!is.character(type)) type='Weighted'
  if(missing(w)) type='Weighted'
  type <- match.arg(tolower(type),
                    c('rc','right.censoring','re','rounding',
                      'wt','weighting',"weighted"))
  switch(type,
         wt = ,  weighted  = ,
         weighting = .wtmle(x=x,w=w,family=family),
         re = ,
         rounding  = .remle(x=x,w=w,family=family),
         rc = ,
         right.censoring  = .rcmle(x=x,w=w,family=family),
         stop("Data type  not supported")
         )
}


.remle <- function(x,w,family='normal'){
  family <- match.arg(tolower(family),
                      c('normal','gaussian', 'weibull',
                        'gamma'))
  stop("MLE for data with rounding error: not supported yet")
}

.wtmle <- function(x,w,family='normal'){
  family <- match.arg(tolower(family),
                      c('normal','gaussian', 'weibull',
                        'gamma'))
  if(missing(w)) w <- rep(1,length(x))
  x.wt <-  .weighting(x,w)
  switch(family,
         gaussian = ,
         normal = c(.wmean(x.wt),.wsd(x.wt)),
         stop("Family of distribution not supported")
         )
}

.rcmle <- function(x,w,family){
  if(any(x<=0)) stop("Invalid lifetime data in 'x'")
  n <- length(x)
  if(missing(w)) w <- rep(1,n)
  if(any(w!=0&&w!=1)) stop("Invalid censoring status in 'w'")
  if(length(w)!=n) stop("'x' and 'w' have different lengths")
  family <- match.arg(tolower(family),
                      c('normal','gaussian', 'weibull',
                        'gamma',"exponential"))
  switch(family,
         exponential = sum(x)/sum(w),
         weibull  =  .rcmle.weibull(x,w),
         stop("Family of distribution not supported")
         )
}

.rcmle.weibull <- function(x,w){
  n <- length(x);
  sele <- w==1
  stopifnot(sum(sele) > 5) # two few data points
  x0 <- x[sele]; rx <- rank(x0);
  Fhat <- (rx-0.3)/(n+0.4)
  lx <- log(x0); ly <- log(-log(1-Fhat))
  out <- lm(ly~lx)
  kappa <- out$coef[[2]];
  lambda <- exp(-out$coef[[1]]/kappa);
  pars <- c(kappa, lambda);
  .Fortran(.F_RcMleWeibull, as.double(x),
           as.double(w),as.integer(n),
           pars = as.double(pars))$pars
}
gof.test <- function(object, x, ...)
UseMethod("gof")

gof.default <- function(object, x, type='chisq',...)
{
  nam <- deparse(substitute(x));
  ##  check data
   if (!is.numeric(x) && !is.complex(x) && !is.logical(x)) {
    warning("'x' is not numeric or logical: returning NA")
    return(NA_real_)
  }
  x = x[!is.na(x)]; n=length(x)
  if (!is.character(object)) {
    warning("'object' is not a distribution name: returning NA")
    return(NA_real_)
  }
  dist = match.arg(tolower(object),
    c("normal","gamma","weibull"))

  anam = paste("Data does not come from a ", dist, "distribution")

  out = .bincount(x)
  nbins = length(out$y)
  np = switch(dist,
    normal = 2,
    gamma  = 2,
    weibull  = 2,
    stop("Not supported yet.")
    )
  DF = nbins - np -1
  if(DF>0 && type == 'chisq'){
    Fx = switch(dist,
      normal = pnorm(out$x[-1],mean(x),sd(x)),
      gamma  = pgamma(out$x[-1], .mle2.gamma(x)),
      weibull  = pweibull(out$x[-1], .mle2.weibull(x))
      )
    EXs = c(Fx[1],diff(Fx))
    es = n*c(EXs,1-sum(EXs))
    D2 = sum((out$y-es)^2/es)
    p2 = ifelse(DF>0,pchisq(D2,DF,lower.tail=FALSE),NA)
    mname = "Chi-Square Test for goodness-of-fit.";
  }else{
    mname = "Kolmogorov-Smirnov Test for goodness-of-test."
    x0 = seq(min(x),max(x),length=100)
    Fx = switch(dist,
      normal = pnorm(x0,mean(x),sd(x)),
      gamma  = pgamma(x0, .mle2.gamma(x)),
      weibull  = pweibull(x0, .mle2.weibull(x))
      )
    Fn = .edf(x,xgrid=x0)
    D2 = max(abs(Fx-Fn$y))
    ft = D2*sqrt(n);
    p2 = .Fortran(.F_kspvalue, p=as.double(ft))$p
  }
  ##  If paraetric, using chisq-test instead of KS-test
  RVAL <- list(statistic = c(ChiSq=D2),#c(D = D,Chisq=D2),
               p.value = p2,#c(p1,p2),
               parameter = c(df=DF),
               method = mname, 
               alternative = anam,
               data.name = nam)
  class(RVAL) <- "htest"
  return(RVAL)
}

gof.em <- function(object, x,type='chisq',...) {
  nam <- object$data.name
  if(missing(x)) x = rep(object$X, object$F)
  ##  check data
   if (!is.numeric(x) && !is.complex(x) && !is.logical(x)) {
    warning("'x' is not numeric or logical: returning NA")
    return(NA_real_)
  }
  x = x[!is.na(x)]; n=length(x)
  type = match.arg(tolower(type), c("ks","chisq"))

  anam = paste("The EM doesn't well fit the data")
  out = .bincount(x)
  nbins = length(out$y)
  DF = nbins-prod(dim(object$Para))
  if(DF>0 && type == 'chisq'){
    Fx = pmixnorm(out$x[-1], object$Para[,1], object$Para[,2],object$Para[,3])
    EXs = c(Fx[1],diff(Fx))
    es = n*c(EXs,1-sum(EXs))
    D2 = sum((out$y-es)^2/es)
    p2 = ifelse(DF>0,pchisq(D2,DF,lower.tail=FALSE),NA)
    mname = "Chi-Square Test for goodness-of-fit.";
  }else{
    mname = "Kolmogorov-Smirnov Test for goodness-of-test."
    x0 = seq(min(x),max(x),length=100)
    Fx = pmixnorm(x0, object$Para[,1], object$Para[,2],object$Para[,3])
    Fn = .edf(x,xgrid=x0)
    D2 = max(abs(Fx-Fn$y))
    ft = D2*sqrt(n);
    p2 = .Fortran(.F_kspvalue, p=as.double(ft))$p
  }
  ##  If paraetric, using chisq-test instead of KS-test
  RVAL <- list(statistic = c(ChiSq=D2),#c(D = D,Chisq=D2),
               p.value = p2,#c(p1,p2),
               parameter = c(df=DF),
               method = mname, 
               alternative = anam,
               data.name = nam)
  class(RVAL) <- "htest"
  return(RVAL)
}

gof.bde <- function(object, x,...) {
  nam <- object$data.name
  if(missing(x)) x = rep(object$X, object$F)
  ##  check data
   if (!is.numeric(x) && !is.complex(x) && !is.logical(x)) {
    warning("'x' is not numeric or logical: returning NA")
    return(NA_real_)
  }
  x = x[!is.na(x)]; n=length(x)
  
  x0 = object$x; bws = diff(x0); bws = c(bws[1],bws)
  F0 = cumsum(object$y*bws)
  Fn = .edf(x,xgrid=x0)$y;
  D = max(abs(F0-Fn))
  ft = D*sqrt(n);
  out = .Fortran(.F_kspvalue, p=as.double(ft))
  anam = paste("The binned density estimate doesn't well fit the data")
  ##conf.int; null.value; parameter

  RVAL <- list(statistic = c(D = D), p.value = out$p,
               method = "Kolmogorov-Smirnov Test for goodness-of-fit.", 
               alternative = anam,
               data.name = nam)
  class(RVAL) <- "htest"
  return(RVAL)
}


###  subroutines to be called


.lecount <- function(x,y) sum(y<=x)

.bincount <- function(x){
  x=x[!is.na(x)]; s = sd(x); mu = mean(x);
  n=length(x)
  x0 = seq(mu-5.7*s,mu+6*s,by=0.3*s)
  ps = apply(matrix(x0,ncol=1),1,FUN=.lecount,y=x)
  cts = diff(ps)
  bounds = NULL; counts=NULL;
  k = length(cts)
  for(i in 1:(k-1)){
    if(cts[i]>=5){
      bounds = c(bounds,x0[i])
      counts = c(counts,cts[i])
      ul = x0[i+1]
    }else{
      cts[i+1] = cts[i+1]+cts[i]
    }
  }
  m = n - sum(counts)
  if(m>5){
    bounds = c(bounds,ul)
    counts = c(counts,m)
  }else{
    counts[length(counts)] = counts[length(counts)]+m
  }
  list(x=bounds,y=counts);
}


.mle2.gamma <- function(x){
  x = x[!is.na(x)]
  .Fortran(.F_FitGamma, as.double(x),
           as.integer(length(x)), as.double(rep(0,2)))[[3]]
}

.mle2.weibull <- function(x){
  x = x[!is.na(x)]
  .Fortran(.F_FitWeibull, as.double(x),
           as.integer(length(x)), as.double(rep(0,2)))[[3]]
}


