

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



edf <- function(x, weights, xgrid, na.rm = TRUE){
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
                 data.name = name), class = "edf")
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


print.edf <- function (x, digits = NULL, ...) 
{
  cat("\nCall:\n\t", deparse(x$call), "\n\nData: ", x$data.name, 
      " (", x$n, " obs.);", "\n", sep = "")
  print(summary(as.data.frame(x[c("x", "y")])), digits = digits, 
        ...)
  invisible(x)
}

plot.edf  <- function (x, main = NULL, xlab = NULL, ylab = "EDF", type = "l", 
    zero.line = TRUE, ...) 

{
  if (is.null(xlab)) 
    xlab <- paste("N =", x$n)
  if (is.null(main)) 
    main <- deparse(x$data.name)
  plot.default(x, type = 'n', 
               main = main, xlab = xlab, ylab = ylab, ...)
  x0 = x$x; y0 = x$y; m = length(x0)
  x1 = c(x0[-1], x0[m] * 1.2);
  y1 = c(y0)
  segments(x0,y0,x1,y1,...)
  
  x1 = x1[-m]; y1 = y1[-m]
  x3 = x0[-1]; y3 = y0[-1]
  segments(x1,y1,x3,y3, col='gray',lty=3)
  if (zero.line) 
    abline(h = 0, lwd = 0.1, col = "gray")
  invisible(NULL)
}

lines.edf  <- function (x, ...) 
{
  x0 = x$x; y0 = x$y; m = length(x0)
  x1 = c(x0[-1], x0[m] * 1.2);
  y1 = c(y0)
  segments(x0,y0,x1,y1,...)
  
  x1 = x1[-m]; y1 = y1[-m]
  x3 = x0[-1]; y3 = y0[-1]
  segments(x1,y1,x3,y3, col='gray',lty=3)
  invisible(NULL)
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
