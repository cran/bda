.onUnload <- function(libpath)
  library.dynam.unload("bda",  libpath)

.bdaConnect <- NULL


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
  xhist <- histogram(x, from=from,to=to)

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
##    out <- .fftwkde(x,w,h,n,from,to,M,tau=tau,truncate=truncate)
##    fx <- out$y
  }
  
  totMass <- sum(fx) *(to-from)/(M-1)

  return(structure(list(y=fx/totMass, x=gpoints,
                        bw = pars,
                        ucb=NULL, lcb=NULL,
                        call = match.call(),
                        data.name = xhist$xname,
                        hist = xhist,
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


###  Permutation test

##perm.test <-
##  function(x,y,fun,alternative = "two.sided", trials = 1000,...)
##  UseMethod("perm")

##  if fun is missing, define it as comparing distributional
##  difference; otherwise, user-defined.
perm.test <- function(x,y,fun,alternative = "two.sided", 
                      trials = 1000,...)
{
  xnam <- deparse(substitute(x));
  ynam <- deparse(substitute(y));
  nam = paste(xnam, '(',length(x), ') vs ', ynam, '(',length(y),")")
  if (!is.numeric(x) && !is.complex(x) && !is.logical(x)) {
    warning("'x' is not numeric or logical: returning NA")
    return(NA_real_)
  }
  x <- x[!is.na(x)]
  if (!is.numeric(y) && !is.complex(y) && !is.logical(y)) {
    warning("'y' is not numeric or logical: returning NA")
    return(NA_real_)
  }
  y <- y[!is.na(y)]

  alternative = match.arg(tolower(alternative),
      c("one.sided","two.sided"))

    
  ## use of replicate() with parameters:
  if(missing(fun)) fun = .edfperm

  D = fun(x,y)
  rfun <- function(x,y){
    nx = length(x); ny = length(y); n = nx + ny;
    n0 = sample(n)
    xy = c(x,y)
    fun(xy[n0[1:nx]], xy[n0[(nx+1):n]])
  }
  bar <- function(n, x,y) replicate(n, rfun(x=x,y=y))
  z = bar(trials, x=x,y=y)

  ##  z = replicate(trials, fun(x=x,y=y))
  ##  z = apply(as.matrix(rep(1,trials),ncol=1),1,fun,x=x,y=y)
  pv = switch(alternative,
    two.sided = mean(abs(D) <= abs(z)),
    one.sided = mean(D <= z),
    stop("Wrong test type!")
    )

  RVAL <- list(statistic = c(D = D), p.value = pv,
               method = "Permutation test to compare two samples/populations.", 
               data.name = nam)
  class(RVAL) <- "htest"
  return(RVAL)

}

.edfperm <- function(x,y){
  ## no randomization within this function.  It could be defined by
  ## users. The radomization should be a step in the main program. 
  rngx <- range(c(x,y))
  from <- rngx[1]; to <- rngx[2]
  Fx1 <- edf(x,from=from,to=to,gridsize=512L);
  Fx2 <- edf(y,from=from,to=to,gridsize=512L);
  max(abs((Fx1$y-Fx2$y)))
}


###  Mediation test:


mediation.test <- function(mv,iv,dv){
  ## mx: mediation variable
  ## iv:  indep. variable
  ## dv: dep. var.
  if(any(is.na(mv)))stop("Mediator contains missing value(s)")
  if(any(is.na(iv)))stop("Mediator contains missing value(s)")
  if(any(is.na(dv)))stop("Mediator contains missing value(s)")
  nm = length(mv); ni = length(iv); nd = length(dv);
  if(nm!=ni | nm!=nd | ni!=nd) stop("Variables have different lengths.")
  tmp = summary(lm(mv~iv));
  a = tmp$coef[2,1];sa=tmp$coef[2,2];
  tmp = summary(lm(dv~mv+iv));
  b = tmp$coef[2,1];sb=tmp$coef[2,2];
  tmp1 = b^2*sa^2+a^2*sb^2
  tmp2 = sa^2*sb^2
  zsob = a*b/sqrt(tmp1);
  psob = pnorm(-abs(zsob))*2;
  zaro = a*b/sqrt(tmp1+tmp2);
  paro = pnorm(-abs(zaro))*2;
  if(tmp1>tmp2){
    zgm = a*b/sqrt(tmp1-tmp2)
    pgm = pnorm(-abs(zgm))*2;
  }else{
    zgm=NA
    pgm=NA;
  }
  p.value=c(psob,paro,pgm)
  z.value = c(zsob,zaro,zgm)
  out = data.frame(rbind(z.value,p.value));
  names(out)=c("Sobel","Aroian","Goodman")
  out
}

###  Empirical distribution function

### The folloowing functions are to compute the empirical distribution
### function and construct a confidence band.

### Last updated on 04/18/2013

edf <- function(x, weights, freq, from, to, gridsize, digits=0,
                rounding) UseMethod("edf")


edf.default <- function(x, weights, freq, from, to, gridsize,
                        digits=0, rounding){
  if(missing(rounding)) rounding <- "none"
  xr <- .data.rounding(x=x,freq=freq, weights=weights, digits=digits,
                  type=rounding)
  edf.bdata(xr,gridsize=gridsize,from=from,to=to)
}

edf.bdata <- function(x, weights, freq, from, to, gridsize,
                      digits=0, rounding){
  xr <- x;
  if(missing(from)) from <- min(xr$x+xr$a)
  if(missing(to)) to <- max(xr$x+xr$b)
  stopifnot(to > from)
  if(missing(gridsize)){
    xgrid <- xr$x
    gridsize <- length(xgrid)
    Fx <- cumsum(xr$f*xr$w)
    y <- Fx/Fx[gridsize];
  }else{
    if(gridsize < 10) gridsize <- 512L
    xgrid <- seq(from, to, length=gridsize)
    gcounts <- .data.binning(xr, type='edf',from=from,to=to,gridsize=gridsize)
    y <- cumsum(gcounts);
#    y <- Fx/Fx[gridsize];
  }
  
  structure(list(y = y, x = xgrid, n = gridsize, 
                 data.name = xr$data.name),
            class = "edf")
}

print.edf <- function (x, digits = NULL, ...) 
{
  cat("\nCall:\n\t", deparse(x$call), "\n\nData: ", x$data.name, 
      " (", x$n, " obs.);", "\n", sep = "")
  print(summary(as.data.frame(x[c("x", "y")])), digits = digits, 
        ...)
  invisible(x)
}

plot.edf <- function (x, main = NULL, xlab = NULL, ylab = "EDF", type
    = "l", zero.line = TRUE, ...)

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


### histogram

## 2013/04/19
histogram <- function(x, weights, freq, from, to, nclass,
                      digits=0, rounding)
  UseMethod("histogram")


histogram.default <- function(x, weights, freq, from, to, nclass,
                              digits=0, rounding){
  if(missing(rounding)) rounding <- "none"
  xr <- .data.rounding(x=x,freq=freq, weights=weights, digits=digits,
                  type=rounding)
  histogram(xr, from=from,to=to,nclass=nclass)
}

histogram.bdata <- function(x, weights, freq, from, to, nclass,
                            digits=0, rounding){
  xr <- x
  if(missing(from)) from <- min(xr$x+xr$a)
  if(missing(to)) to <- max(xr$x+xr$b)
  stopifnot(to > from)

  if(missing(nclass)){
    bw <- xr$b - xr$a
    if(bw == 0){# not binned/compute nclass using lscv
      nclass <- .histclass.lscv(xr, xrange=to-from)
      bw <- (to-from)/nclass
    }else{
      nclass <- round((to-from)/bw)
    }
  }else{
    bw <- (to-from)/nclass
  }
  gcounts <- .data.binning(xr, type='hist',from=from,to=to,gridsize=nclass+1)
  gcounts <- gcounts[-1]; 
  y <- gcounts * sum(xr$f)
  breaks <- seq(from,to,length=nclass+1)
  structure(list(breaks = breaks,
                 counts = y,
                 intensities = gcounts/bw,
                 density = gcounts/bw,
                 mids = breaks[-1]-0.5*bw,
                 bw = bw, nclass=nclass,
                 xname = xr$data.name,
                 equidist=FALSE),
            class="histogram")  
}


.histclass.lscv <- function(y,xrange){
  n <- sum(y$f); 
  if(n<=5) nclass <- 2
  else{
    n1 <- min(round(n * 0.2),5)
    n2 <- min(round(n * 0.5),30)
    Jh <- NULL; m <- NULL
    for(i in n1:n2){
      h <- xrange/i
      Jh <- c(Jh,.Jhscore(y,h,i));
      m <- c(m,i)
    }
    nclass <- m[which(Jh==min(Jh))[1]]
  }
  nclass
}

.Jhscore <- function(x,h,nclass){
  n <- sum(x$f)
  out <- edf(x, gridsize=nclass+1)
  fhat <- diff(c(out$y,1))
  Jh <- 2/h/(n-1)-(n+1)/h/(n-1)*sum(fhat^2)
}


histosmooth <- function(x, weights, freq, type='spline', bw,
                    from, to, gridsize, digits=0, rounding)
  UseMethod("histosmooth")

histosmooth.default <- function(x, weights, freq, type='spline', bw,
                                from, to, gridsize, digits=0, rounding)
{

  if(missing(rounding)) rounding <- "none"
  xhist <- histogram(x=x,freq=freq, weights=weights, from=from,
                  to=to, digits=digits,
                  rounding=rounding)
  histosmooth(x=xhist,freq=freq, weights=weights, from=from,
              to=to, gridsize=gridsize, digits=digits,
              rounding=rounding, type=type,bw=bw)
}

histosmooth.histogram <- function(x, weights, freq, type='spline', bw,
                                  from, to, gridsize, digits=0, rounding)
{
  xhist <- x
  bws <- diff(xhist$breaks)
##  if(any(bws != bws[1]))
##    stop("Not applicable to histograms with unequal bin widths")

  if(missing(gridsize)) gridsize <- 512L
  
  type <- match.arg(tolower(type),
                    c('spline',"lp",'localpolynomial',
                      "kde"))
  res <- switch(type,
                spline = .histospline(xhist, gridsize),
                kde = .histokde(xhist,bw,gridsize),
                lp = ,
                localpolynomial = .histonpr(xhist,bw,gridsize)
                )

  return(structure(list(y = res$y, x = res$x,
                        ucb=res$ub, lcb=res$lb,
                        call = match.call(),
                        bw = res$bw,
                        data.name = xhist$xname,
                        hist = xhist,
                        type=type),
                   class='smooth'))
}


print.smooth <- function (x, digits = NULL, ...)
{
  cat("Call:\n  ", deparse(x$call),
      "\n  Data: ", x$data.name, sep = "")
  if(!is.null(x$bw)){
    cat(";\tbw =[", round(x$bw[1],3),",",
        round(x$bw[2],3),"];\tselector := '",
        x$type, "'\n", sep = "")
  }else if(!is.null(x$pars)){
    tmp <- rbind(x$f,x$a,x$b)
    f <- x$f; m <- x$m
    nbin <- length(f)
    if(f[nbin]==0) m[nbin] <- m[nbin-1]
    n <- sum(x$f)
    mu <- ifelse(is.null(x$mean), sum(f*m)/n, x$mean)
    s <- ifelse(is.null(x$sigma), sqrt(sum((m-mu)^2*f)/(n-1)),x$sigma);
    cat("\n  Dist := '",
        x$type, "',\tPara := [", sep = "")
    k <- length(x$pars)
    if(k==1){
      cat(x$pars, "]\n",sep='')
    }else{
      for(i in 1:(k-1)){
        cat(x$pars[i], ", ",sep='')
      }
      cat(x$pars[k], "]\n", sep='')
    }
    cat("  mean = ", round(mu,3),",\ts.d. = ",
        round(s,3),".\n", sep='')
  }else{
    cat("\n", sep="")
  }
  print(summary(as.data.frame(x[c("x", "y")])),
        digits = digits, ...)
  invisible(x)
}

plot.smooth <- function (x, hist=TRUE, ...)
{
  if(!is.null(x$hist)&&hist){
    plot(x$hist)
    lines(x,...)
  }else{
    plot(x$y~x$x,...)
  }
  invisible(x)
}

hist.smooth <- function(x,plot=TRUE,main,...){
  if(missing(main)) main <- x$dist
  if(plot)
    plot(x$hist,main=main,...)
  invisible(x$hist)
}

.histonpr <- function(x,bw,gridsize)
  {
    F <- x$counts; X <- x$mids;
    nbin <- x$nclass; n <- sum(F)
    a <- x$breaks[1]; b <- x$breaks[nbin+1];
    if(missing(bw)){
      lscv <- 1
      A <- x$breaks[1:nbin]; B <- x$breaks[-1]
      x0 <- runif(n,rep(A,F),rep(B,F))
      bw <- bw.nrd(x0)
      adaptive = 1  # use adaptive bandwidth selector
    }else{
      lscv <- 0
      adaptive = 0  # use adaptive bandwidth selector
    }
    M <- 2^(ceiling(log2(gridsize)))
    gpoints <- seq(a, b, length = M)
    ##    kernel <- 1  # gaussian kernel to ensure "absolute continuous"
    Y = sqrt(nbin/n*(F+0.25))
    out <- .Fortran(.F_lpsmooth,
                    fx = as.double(gpoints), as.integer(M),
                    as.double(X), as.double(Y), as.integer(nbin),
                    bw = as.double(bw), as.integer(lscv),
                    as.double(c(a,b)), as.integer(adaptive),
                    ellx = double(M), kappa=double(1))

    ## print(out$kappa)
    cv <-  .Fortran(.F_tubecv,
                    cv=as.double(out$kappa), as.double(0.95))$cv
    ## print(cv)
    
    y = out$fx
    sele1 = is.na(y) | !is.finite(out$fx)
    if(any(sele1)) y[sele1] = 0.0
    MOE <- 0.5*cv * out$ellx * sqrt(nbin/n)
    ll = y - MOE; ll[ll<0] <- 0
    ul = y + MOE; 
    
    y <- y*y; tot.mass <- sum(y)*(b-a)/M; 
    ll <- ll*ll/tot.mass
    ul <- ul*ul/tot.mass
    
    list(y = y/tot.mass, x = gpoints, bw=out$bw, ub=ul, lb=ll)
  }

.histospline <- function(x,gridsize)
  {
    out <- spline(x$mids, x$density, n=gridsize,
                  xmin=x$breaks[1], xmax=x$breaks[x$nclass+1])
    f0 <- out$y
    f0[f0<0] <- 0
    gpoints <- out$x
    list(y = f0, x = gpoints, bw=NULL)
  }

.histokde <- function(x,bw,gridsize)
  {
    ## set grid to evaluate the pdf/cdf
    nbin <- x$nclass
    a <- x$breaks[1]; b <- x$breaks[nbin+1];
    M <- 2^(ceiling(log2(gridsize)))
    gpoints <- seq(a, b, length = M)
    ## to approximate the initial estimate of f(x) by simulation
    F <- x$counts; X <- x$mids; n <- sum(F)
    A <- x$breaks[1:nbin]; B <- x$breaks[-1]
    x0 <- runif(n,rep(A,F),rep(B,F))
    f0 <- density(x0, n=M)$y
    
    if(missing(bw))
      h <-  (1/(4*pi))^(1/10)*(243/(35*n))^(1/5)*sqrt(var(x0))*15^(1/5)
    else h = bw
    ## iterate to find a smoothed KDE
    iter <- 100;
    out <- .Fortran(.F_smoothkde, fx=as.double(f0), as.double(gpoints),
                    as.integer(M), as.double(X), as.double(F),
                    as.integer(nbin), as.double(x$bw/2),
                    bw = as.double(h),iter=as.integer(iter))
    y = out$fx
    if(out$iter>=iter) y = rep(0, M)
    sele1 = is.na(y) | !is.finite(y)
    if(any(sele1))
      y[sele1] = 0.0
    tot.mass <- sum(y)*(b-a)/(M-1)
    
    list(y = y/tot.mass, x = gpoints, bw=out$bw)
  }


##################################################################
### Functions:
##  bdde: binned data density estimation
## create an internal R object "bindata" with
## f -> frequencies, a -> lower limits, b -> upper limits

### dist: weibull=0,

bde <- function(f, breaks, dist="ewd", gridsize=512L,iter=100,...)
{
  name <- deparse(substitute(f))

  ## NA values are not allowed
  if(any(is.na(f)))
    stop("'f' contains missing value(s).")
  if(any(is.na(breaks)))
    stop("'f' contains missing value(s).")
  if(any(!is.finite(f)))
    stop("Infinite value of 'f' not allowed.")
  ## frequencies should to be positive integers or zero's.
  if(any(f != round(f))) stop("Non-integer frequency found.")
  if(any(f <0)) stop("Negative frequency found.")

  if(any(diff(breaks) <= 0))
    stop("'breaks' is not strictly increasing.")
  nbin <- length(f) # bin number
  if(length(breaks) != nbin + 1)
    stop("Lengths not match for 'f' and 'breaks'")
  if(min(breaks)<0) stop("Negative not allowed")
  if(!is.finite(breaks[nbin])) stop("Invalid value in 'breaks'")
  dist <- match.arg(tolower(dist),
                    c('ewd','normal', 'gaussian','dagum',
                      'weibull','gb','gld'))

  ## if top coded, compute an upper bound.
  xmax <- max(breaks)
  if(!is.finite(xmax)){
    top.coded <- TRUE
  }else{
    top.coded <- FALSE
  }
  a <- breaks[1:nbin]; b <- breaks[-1]
  b0 <- b; # copy data for display
  if(top.coded){
    b[nbin] <- 2 * max(a)
    tmp <- .bweibullplot(f,a,b,top.coded)
    lambda <- tmp$pars[[1]]; kappa <- tmp$pars[[2]]
    mu <- lambda * gamma(1+1/kappa)
    s2 <- lambda^2*gamma(1+2/kappa)-mu^2
    xmax <- mu + 3.0 * sqrt(s2)
    xmax <- max(xmax, 2*max(a))
    b[nbin] <- xmax; breaks[nbin+1] <- xmax
  }
  
  ##  draft a histogram
  n0 <- which(match(rev(f)==0,FALSE)==1)[1]-1
  if(n0>0){
    nbin <- nbin - n0
    f <- f[1:nbin]
    a <- a[1:nbin]
    b <- b[1:nbin]; b0 <- b  # update date for display
    breaks <- breaks[1:(nbin+1)]
    xmax <- max(breaks)
  }
  m <- 0.5 * (b + a); 
  binwidth <- diff(breaks)
  den <- f/(binwidth*sum(f))

  x.hist <- structure(list(breaks = breaks,
                           counts = f,
                           intensities = den,
                           density = den,
                           mids = m, nclass=nbin,
                           xname = name,
                           equidist=FALSE),
                      class="histogram")  

  #####################################################
  if(gridsize<10) stop("Gridsize might be too small")
  x <- seq(breaks[1L],breaks[nbin+1], length=gridsize)

  if(dist=="weibull"){
    fixed <- TRUE
    pars <- .ewdbde(f,a,b,iter,top.coded,fixed)
    lambda <- pars[[2]]; kappa <- pars[[1]]
    pars <- c(kappa, lambda)
    y <- dweibull(x,kappa,lambda)
    sele <- is.finite(y)
    if(any(!sele)){
      y <- y[sele]; x <- x[sele]
    }
    sele <- y>0
    if(any(sele)){
      y <- y[sele]; x <- x[sele]
    }
    mu <- lambda * gamma(1+1/kappa)
    s2 <- lambda^2*gamma(1+2/kappa)-mu^2
  }else if(dist=='ewd'){
    fixed <- FALSE
    pars <- .ewdbde(f,a,b,iter,top.coded,fixed)
    alpha <- pars[[3]];
    lambda <- pars[[2]];
    kappa <- pars[[1]]
    pars <- c(kappa, lambda,alpha)
    y  <- .dewd(x, kappa,lambda,alpha)
    sele <- is.finite(y)
    if(any(!sele)){
      y <- y[sele]; x <- x[sele]
    }
    sele <- y>0
    if(any(sele)){
      y <- y[sele]; x <- x[sele]
    }
    Fx <- .pewd(x, kappa,lambda,alpha)
    Px <- diff(Fx); Px <- c(Px,1-sum(Px))
    mu <- sum(x*Px)
    s2 <- sum((x-mu)^2*Px)
  }else if(dist=='dagum'){
    pars <- .dagumbde(f,a,b,iter,top.coded)
    alpha <- pars[[3]];
    lambda <- pars[[2]];
    kappa <- pars[[1]]
    y <- .ddagum(x,kappa,lambda,alpha)
    sele <- y>0
    if(any(sele)){
      y <- y[sele]; x <- x[sele]
    }
    Fx <- .pdagum(x,kappa,lambda,alpha)
    Px <- diff(Fx); Px <- c(Px,1-sum(Px))
    mu <- sum(x*Px)
    s2 <- sum((x-mu)^2*Px)
  }else if(dist=='gb'){
    sele <- f==0  # bins with zero counts
    f1 <- f[!sele]
    a1 <- a[!sele]
    b1 <- b[!sele]
    xmin <- min(a1);
    x <- seq(xmin,xmax, length=gridsize)
    a2 <- (a1-xmin)/(xmax-xmin)
    b2 <- (b1-xmin)/(xmax-xmin)
    pars <- .betabde(f1,a2,b2,iter,top.coded)
    pa <- pars[[1]];
    pb <- pars[[2]];

    y <- dbeta((x-xmin)/(xmax-xmin),pa,pb)/(xmax-xmin)
    sele <- is.finite(y)
    if(any(!sele)){
      y <- y[sele]; x <- x[sele]
    }
    sele <- y>0
    if(any(sele)){
      y <- y[sele]; x <- x[sele]
    }
    mu <- pa/(pa+pb)*(xmax-xmin)
    s2 <- pa*pb/(pa+pb)^2/(pa+pb+1)*(xmax-xmin)^2
  }else stop("Distribution not supported")

  out <- structure(list(y=y, x=x,
                        pars=pars,
                        mean=mu, sigma=sqrt(s2),
                        ucb=NULL, lcb=NULL,
                        f=f, a=a, b=b0,
                        call = match.call(),
                        data.name = name,
                        hist=x.hist,
                        type = dist
                        ),
                   class = "smooth")
  invisible(out)
}

######################################################
######################################################
.betabde <- function(f,a,b,iter,top.coded)
{
  alphas <- rep(1,iter)
  tmp <- apply(as.matrix(alphas,ncol=1),1,.betaplot,
               f=f,a=a,b=b,top.coded=top.coded)
  tmp[-1,which(tmp[1,]==max(tmp[1,]))[[1]]]
}

.betaplot <- function(x,f,a,b,top.coded)
  {
    n <- sum(f); nbin <- length(f)
    x <- runif(n,rep(a,f),rep(b,f))
    mu <- mean(x); s2 <- var(x)
    alpha <- mu*(mu*(1-mu)/s2-1)
    pbeta <- alpha/mu-alpha
    Fa <- pbeta(a,alpha,pbeta)
    if(top.coded) Fa <- c(Fa, 1)
    else Fa <- c(Fa, pbeta(b[nbin],alpha,pbeta))
    llk <- sum(log(diff(Fa)))
    c(llk,alpha, pbeta)
  }

######################################################
######################################################

.dagumbde <- function(f,a,b,iter,top.coded){
  alphas <- c(runif(iter), runif(iter,1,5))
  tmp <- apply(as.matrix(alphas,ncol=1),1,.dagumplot,
               f=f,a=a,b=b)
  tmp[1:3,which(tmp[4,]==max(tmp[4,]))[[1]]]
}

.dagumplot <- function(alpha,f,a,b)
  {
    tol <- 0.000001
    n <- sum(f)
    x <- runif(n,rep(a,f),rep(b,f))
    Fhat <- (rank(x)-0.3)/(n+0.4)
    ly <- log(-1-log(Fhat)/alpha+tol)
    lx <- log(x+tol)
    out <- lm(ly~lx)
    kappa <- -out$coef[[2]]
    lambda <- exp(out$coef[[1]]/kappa)
    llk <- .dagumllk(f,a,b,alpha,kappa,lambda)
    pars <- c(kappa, lambda, alpha,llk)
  }

.dagumllk <- function(f,a,b,alpha,kappa,lambda){
  if(is.na(lambda)||is.na(kappa)){
    llk <- -9999999999999.99
  }else if(lambda<=0||kappa<=0){
    llk <- -9999999999999.99
  }else{
    tol <- 0.0000001
    Fa <- .pdagum(a, kappa,lambda,alpha)
    Fb <- .pdagum(b, kappa,lambda,alpha)
    llk <- sum(f*log(Fb-Fa+tol))
  }
}

.ewdbde <- function(f,a,b,iter,top.coded,fixed){
  if(fixed)
    alphas <- rep(1,iter)
  else
    alphas <- c(runif(iter), runif(iter,1,5))
  tmp <- apply(as.matrix(alphas,ncol=1),1,.ewdplot,
               f=f,a=a,b=b)
  tmp[1:3,which(tmp[4,]==max(tmp[4,]))[[1]]]
}

.ewdplot <- function(alpha,f,a,b)
  {
    tol <- 0.000001
    n <- sum(f)
    x <- runif(n,rep(a,f),rep(b,f))
    Fhat <- (rank(x)-0.3)/(n+0.4)
    ly <- log(-log(1-Fhat^(1/alpha))+tol)
    lx <- log(x+tol)
    out <- lm(ly~lx)
    kappa <- out$coef[[2]]
    lambda <- exp(-out$coef[[1]]/kappa)
    llk <- .ewdllk(f,a,b,alpha,kappa,lambda)
    pars <- c(kappa, lambda, alpha,llk)
  }

.ewdllk <- function(f,a,b,alpha,kappa,lambda){
  if(is.na(lambda)||is.na(kappa)){
    llk <- -9999999999999.99
  }else if(kappa<=0||lambda<=0){
    llk <- -999999999999
  }else{
    tol <- 0.00000000001
    Fa <- .pewd(a, kappa,lambda,alpha)
    Fb <- .pewd(a, kappa,lambda,alpha)
    llk <- sum(f*log(Fb-Fa+tol))
  }
  llk
}


##################################################
##  Weibull(x,lambda,kappa)
##  f(x) = k/l*(x/l)^(k-1)*exp(-(x/l)^k)
##  in R Weibill(x,kappa,lambda)
.bweibullplot <- function(f,a,b,top.coded){
  nbin <- length(f)
  n <- sum(f)
  y <- runif(n,rep(a,f),rep(b,f))
  pars <- .weibullplot(y)
  Fa <- pweibull(a,pars[[2]],pars[[1]])
  if(top.coded) Fa <- c(Fa, 1)
  else Fa <- c(Fa, pweibull(b[nbin],pars[[2]],pars[[1]]))
  llk <- sum(log(diff(Fa)))
  list(pars=pars,llk=llk)
}

.weibullplot <- function(x){
  tol <- 0.0000001
  if(any(x<0))
    stop("Negative value not allowed in Weibull distribution")
  n <- length(x)
  Fhat <- (rank(x)-0.3)/(n+0.4)
  ly <- log(-log(1-Fhat))
  lx <- log(x+tol)
  out <- lm(ly~lx)
  kappa <- out$coef[[2]]
  lambda <- exp(-out$coef[[1]]/kappa)
  c(lambda, kappa)
}

.weibullbde <- function(f,a,b){
  nbin <- length(f)
  if(!is.finite(b[nbin])) b[nbin] <- a[nbin]*2-a[nbin-1]
  n <- sum(f)
  y <- runif(n,rep(a,f),rep(b,f))
  pars <- .weibullplot(y)
  ##.bdmle(f,a,b,dist=0,pars=pars)
}

.bdmle <- function(f,a,b,dist,pars){
  nbin <- length(f); npar <- length(pars)
  .Fortran(.F_BDMLE, as.double(f), as.double(a),
           as.double(b),as.integer(nbin),
           pars = as.double(pars), as.integer(npar),
           as.integer(dist))$pars
}

#############################################################


.pewd <- function(x,kappa,lambda,alpha)
  {
    n <- length(x)
    res <- rep(0,n)
    sele1 <- x > 0
    sele2 <- is.finite(x); 
    sele3 <- is.na(x);
    if(any(sele3)) res[sele3] <- NA
    if(any(!sele3&!sele2)) res[!sele3&!sele2] <- 1.0
    x <- x[sele1&sele2]
    stopifnot(kappa>0&&lambda>0&&alpha>0)
    tmp <- (1-exp(-(x/lambda)^kappa))^alpha
    res[sele1&sele2] <- tmp
    res
  }

.dewd <- function(x,kappa,lambda,alpha)
  {
    n <- length(x)
    res <- rep(0,n)
    sele1 <- x > 0
    sele2 <- is.finite(x); 
    sele3 <- is.na(x);
    if(any(sele3)) res[sele3] <- NA
    x <- x[sele1&sele2]
    stopifnot(kappa>0&&lambda>0&&alpha>0)
    tmp <- alpha * (1-exp(-(x/lambda)^kappa))^(alpha-1) *
      exp(-(x/lambda)^kappa) * kappa * (x/lambda)^(kappa-1)/lambda
    res[sele1&sele2] <- tmp
    res
  }

.pdagum <- function(x,a,b,p)
  {
    n <- length(x)
    res <- rep(0,n)
    sele1 <- x > 0
    sele2 <- is.finite(x); 
    sele3 <- is.na(x);
    if(any(sele3)) res[sele3] <- NA
    if(any(!sele3&!sele2)) res[!sele3&!sele2] <- 1.0
    x <- x[sele1&sele2]
    stopifnot(a>0&&b>0&&p>0)
    tmp <- (1+(x/b)^(-a))^(-p)
    res[sele1&sele2] <- tmp
    res
  }

.ddagum <- function(x,a,b,p)
  {
    n <- length(x)
    res <- rep(0,n)
    sele1 <- x > 0
    sele2 <- is.finite(x); 
    sele3 <- is.na(x);
    if(any(sele3)) res[sele3] <- NA
    x <- x[sele1&sele2]
    stopifnot(a>0&&b>0&&p>0)
    tmp <- (a*p)/x*(x/b)^(a*p)*((x/b)^a+1)^(-p-1)
    res[sele1&sele2] <- tmp
    res
  }


