
##  function histo(...).  Return an R object 'hist'.  Update the
##  object returned by histo(...).

histogram <- function(x, w, nclass, binwidth, lb,
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
  fhat = diff(c(edf(x,xgrid=x0)$y,1))
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
  fhat = diff(c(edf(x,xgrid=x0)$y,1))
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
  
  if(missing(lb)) lb <- a;
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
  fhat = diff(c(edf(x,weights=w, xgrid=x0)$y,1))
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


