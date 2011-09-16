
##  function histo(...).  Return an R object 'hist'.  Update the
##  object returned by histo(...).

histo <- function(x, weights, nclass, binned = FALSE, range.x)
{
  name <- deparse(substitute(x))
  if (!is.numeric(x)) 
    stop("'x' must be numeric")
  if(any(is.na(x))) stop("'x' contains missing value(s)")
  if(!missing(weights)){
    if (!is.numeric(weights)) 
      stop("'weights' must be numeric")
    stopifnot(all(weights>0))
  }
 
  if(missing(range.x)){
    a = min(x); b = max(x);
  }else{
    a = range.x[1L]; b = range.x[2L]
  }
  if(binned){
    if(missing(weights)){  # 'x' is prebinned/rounded, but not tabled
      tmp = table(x)
      x0 = as.numeric(names(tmp));
      y0 = as.numeric(tmp);
      N = length(x)
    }else{
      x0 = x; y0 = weights;
      N = sum(weights)
    }
    nclass = length(x0)
    bw = diff(x0); h = bw[1]
    stopifnot(all(bw==h))  # stop if not equally binned
    y0 = y0/(N*h)
    breaks = c(x0 - h/2, x0[nclass] + h/2)
  }else{
    if (missing(nclass)) nclass = 'lscv'
    if (is.character(nclass)) {
      nclass <- match.arg(tolower(nclass),
                          c("sturges", "fd", "freedman-diaconis", "scott", "lscv"))
      nclass <- switch(nclass,
                       sturges = nclass.Sturges(x), 
                       `freedman-diaconis` = ,
                       fd = nclass.FD(x),
                       scott = nclass.scott(x),
                       lscv = .bwlscv(x, w = weights),
                       stop("unknown 'breaks' algorithm"))
    }
    h = (b-a)/nclass
    breaks = seq(a, b, length = nclass + 1)
    y0 = .wbin(x, weights, breaks)
    x0 = breaks[-1] - h*0.5
    y0 = y0/(sum(y0)*h)
    N = length(x)
  }

  tmp = structure(
    list(breaks = breaks, counts = y0 * N,
         intensities = y0, density = y0, mids = x0,
         xname = name, equidist=FALSE), class="histogram")

  out = structure(
    list(y = y0, x = x0, n = N, nclass = nclass, 
         data.name = name, plot=tmp), class = "hist")
  return(out)
}

print.hist <- function (x, digits = NULL, ...) 
{
  cat("\nData: ", x$data.name, 
      " (", x$n, " obs.)\n", sep = "")
  print(summary(as.data.frame(x[c("x", "y")])), digits = digits, 
        ...)
  invisible(x)
}

plot.hist  <- function (x,  ...)   plot(x$plot,...)

lines.hist  <- function (x,  ...)   lines(x$plot,...)



.bwlscv <- function(x, w){
  n = length(x)
  xrange = diff(range(x))
  if(missing(w)) w = rep(1, n)
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

