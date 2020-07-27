wdekde <- function(x, w, bandwidth,s.x)
{
    if (!missing(bandwidth) && bandwidth <= 0)
        stop("'bandwidth' must be strictly positive")
    if(missing(w)){
        w <- rep(1,length(x))
    }else{
        stopifnot(length(w)==length(x))
        stopifnot(all(w>0))
    }
    sele <- is.na(x)|is.na(w)
    if(any(sele)){
        x <- x[!sele]
        w <- w[!sele]
    }
    if(missing(bandwidth)){
        h0 <- bw.nrd(x);
    }else{
        stopifnot(is.numeric(bandwidth))
        stopifnot(length(bandwidth)==1)
        stopifnot(bandwidth > 0)
        h0 <- bandwidth
    }
    
    ## Set kernel support values    
    tau <-  4;
    gridsize = 512L  # make sure the grid is fine enough
    M <- 2^(ceiling(log2(gridsize)))
    range.x <- c(min(x)-tau*h0, max(x)+tau*h0)
    a <- range.x[1L];  b <- range.x[2L]    
    ## Set up grid points and bin the data
    gpoints <- seq(a, b, length = M)
  
    out <- .Fortran(.F_wdekde,
                    as.double(x),
                    as.double(w/sum(w)),
                    as.integer(length(x)),
                    y=as.double(gpoints),
                    as.integer(M),
                    as.double(h0),
                    as.double(s.x))
           
    list(y=out$y,x=gpoints)
}


wkde <- function(x, w, bandwidth, freq=FALSE, gridsize = 401L,
                 range.x, truncate = TRUE, na.rm=TRUE)
{
  name <- deparse(substitute(x))
  ## Install safeguard against non-positive bandwidths:
  if (!missing(bandwidth) && bandwidth <= 0)
    stop("'bandwidth' must be strictly positive")
  if(missing(w)){
    ##      w <- rep(1,n)
    
    if(any(is.na(x))){
      if(na.rm) x <- x[!is.na(x)]
      else stop("'x' contains missing value(s)")
    }
    xt <- table(x); 
    x <- as.numeric(names(xt))
    w <- as.numeric(xt)
    n <- length(x)
    freq <- TRUE
    totMass <- 1.0
  }else{
    tmp = .wtdata(x,w,na.rm=na.rm)
    x <- tmp$x; w <- tmp$w; n <- tmp$n; totMass <- tmp$totMass
  }
  ##  ngrid = n  ## ???? can be removed later

  ## Set kernel support values    
  tau <-  4;
  h0 <- bw.wnrd0(x,w, freq=freq);
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
  out <- if(missing(bandwidth))
    .nwkde(x,w,h0,n,a,b,M,tau=tau,truncate=truncate)
  else if(is.numeric(bandwidth))
    .nwkde(x,w,bandwidth,n,a,b,M,tau=tau,truncate=truncate)
  else if(is.character(bandwidth))
    switch(tolower(bandwidth),
           "wnrd" = .nwkde(x,w,bw.wnrd(x,w,freq=freq),n,a,b,M,
             tau=tau,truncate=truncate),
           "wnrd0" = .nwkde(x,w,h0,n,a,b,M,
             tau=tau,truncate=truncate),
           "blscv" = .nwkde(x,w,
             bw.blscv(x,w,h0,gridsize=M),
             n,a,b,M,tau=tau,truncate=truncate),
           "wmise" = .Fortran(.F_wkde,as.double(x),as.double(w/sum(w)),
             as.integer(n),x=as.double(gpoints),y=as.double(gpoints),
             Fy=as.double(gpoints),as.integer(M), bw=as.double(pars)),
           "awmise" = .Fortran(.F_awkde,as.double(x),as.double(w/sum(w)),
             as.integer(n),x=as.double(gpoints),y=as.double(gpoints),
             Fy=as.double(gpoints),as.integer(M), bw=as.double(pars)),          
           stop("Bandwidth selector not supported."))
  else stop("Bandwidth selector not supported.")

  if(length(out$bw)==1){
    bw = out$bw; sp=NULL;
  }else{
    bw = out$bw[1]; sp = out$bw[2]
  }

  list(y=out$y,x=out$x, bw= bw, sp = sp)
}

.bw.wmise <- function(x, w, na.rm=TRUE)
{
  tmp = .wtdata(x,w,na.rm=na.rm)
  x <- tmp$x; w <- tmp$w; n <- tmp$n;
  h0 <- bw.wnrd(x,w)
  m = 50; h = seq(h0/10, by = h0/10, length=m)
  out = .Fortran(.F_wmise, as.double(x), as.double(w), as.integer(n),
    h = as.double(h), mise = as.double(h), as.integer(m))
  sele = which(out$mise == min(out$mise))
  bw = out$h[sele[1]]
  list(bw = bw, x = out$h, y = out$mise)
}

 
bw.blscv <- function(x,w,h,gridsize=256)
{
  n <- length(x); 
  if(missing(w)) w = rep(1/n,n)
  h0 <- ifelse(missing(h), bw.wnrd0(x,w), h)
  hl <- h0 * 0.5; hh <- 2.5*h0
  ## not needed if use optimize
  iter <- 50;   y <- rep(0,iter) # y stores the lscv's
  hs <- seq(hl,hh,length=iter)

  ## define grid and bin the weighted data
  a <- min(x); b <- max(x);
  M <- 2^(ceiling(log(gridsize)/log(2)))
  gpoints <- seq(a, b, length = M)
  gcounts <- .wlinbin(x, w*n, gpoints, truncate=TRUE)
  gcounts <- gcounts * M/n/(b-a)
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
  
#  list(bw = h0, x=hs, y = Myl)
}

bw.wnrd <- function(x, w, freq=FALSE, na.rm=TRUE)
{
  size <- ifelse(freq, sum(w), length(x))
  tmp = .wtdata(x,w,na.rm=na.rm)
  x <- tmp$x; w <- tmp$w; n <- tmp$n; 
  s = .wsd(x,w); iqr = .wiqr(x,w);
  1.06*min(s,iqr/1.34)*size^-.2
}

bw.wnrd0 <- function(x, w, freq=FALSE, na.rm=TRUE)
{
  size <- ifelse(freq, sum(w), length(x))
  tmp = .wtdata(x,w,na.rm=na.rm)
  x <- tmp$x; w <- tmp$w; n <- tmp$n;
  s = .wsd(x,w); iqr = .wiqr(x,w);
  0.9*min(s,iqr/1.34)*size^-.2
}

##  Internal subroutines to be called

.wtdata <- function(x,w,na.rm=TRUE){
  n <- length(x)
  if(missing(w)) w <- rep(1/n,n)
  stopifnot(all(w>0))
  if(any(is.na(x))){
    if(na.rm){
      sele <- !is.na(x) 
      x <- x[sele]; w <- w[sele]
      n <- length(x)
    }else stop("Missing value(s) in 'x'")
    stopifnot(any(!is.na(w)))  # no 'na' allowed in 'w'
  }    
  x.finite <- is.finite(x)
  totMass <- mean(x.finite)
  x <- x[x.finite]; w <- w[x.finite];
  n <- length(x); w <- w/sum(w)
  list(x=x,w=w,n=n,totMass=totMass)
}

.wmean <- function(x,w) 
  ifelse(missing(w), mean(x), sum(x*w)/sum(w))

.wvar <- function(x,w)
  sum((x-.wmean(x,w))^2 * w/sum(w))

.wsd <- function(x,w)
  sqrt(.wvar(x,w))

## .wquantile should not be used to compute the quantiles with (very)
## low/high quantile levels, or for data with small sizes.

.wquantile <- function(x,w,level){
  if(missing(w)) w = rep(1,length(x))/length(x)
  else stopifnot(length(x)==length(w))
  w <- w/sum(w)
  ox <- order(x); w <- w[ox]; x <- x[ox]
  Fw <- cumsum(w)/sum(w)
  stopifnot(all(!is.na(level)))
  stopifnot(all(level <= 1 & level >= 0))
  
  sele = level > min(Fw) & level < max(Fw)
  out = rep(0, length(level))
  if(sum(sele)>0){
    out[sele] = approx(Fw,x,level[sele])$y
  }
  if(sum(!sele)>0)
    out[!sele] = as.numeric(quantile(x,level[!sele], na.rm=TRUE))
  out  
}

.wiqr <- function(x,w){
  diff(.wquantile(x,w,level=c(.25,.75)))
}

##  class number selectors

.nclass.scott <- function (x,w) 
{
  size <- length(x)
  h <- 3.5 * .wsd(x,w) * size^(-1/3)
  if (h > 0) 
    ceiling(diff(range(x))/h)
  else 1L
}

.nclass.FD <- function (x,w)  
{
  size <- length(x)
  w <- w/sum(w)
  h <- .wiqr(x,w)
  if (h == 0){
    cf <- abs(cumsum(w)-0.5)
    h <- x[which(cf==min(cf))[1]]
  }
  if (h > 0) 
    ceiling(diff(range(x))/(2 * h * size^(-1/3)))
  else 1L
}

.nclass.Sturges <- function (n) 
  ceiling(log2(n) + 1)


## For application of linear binning to a univariate data set.
.wlinbin <- function(X, W, gpoints, truncate = TRUE)
{
    n <- length(X)
    M <- length(gpoints)
    trun <- ifelse(truncate, 1L, 0L)
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


.nwkde <- function(x,w, h, n, a,b, M,tau=4,truncate=FALSE)
{
    gpoints <- seq(a, b, length = M)
    if(any(w>1)) w <- w/sum(w)
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

