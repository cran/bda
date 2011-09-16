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
