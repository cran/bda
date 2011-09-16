## Construct histogram

.hist <- function(x,m){
  x=x[!is.na(x)]; n=length(x);  
  h = diff(range(x))/m
  x0 = seq(min(x),max(x),by=h)
  fhat = diff(c(edf(x,x0)$y,1))
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
  fhat = diff(c(edf(x,x0)$y,1))
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

histo <- function(x, m=NULL, alpha=0.05, binned=FALSE,just="center",scale=1.0){
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
  epsn = qnorm(1-alpha/2/m)/2*sqrt(m/n)
  LFn = sqrt(fhat) - epsn;
  LFn[LFn<0] = 0; LFn=LFn^2
  UFn = (sqrt(fhat) + epsn)^2

  return(structure(list(y=fhat,x=out$x, l=LFn,u=UFn,
                        Jhs = Jhs, Jh=Jh, ms=ms, m=m, n = n,
                        data = x, data.name = name
                        ), class = "hist"))
}

print.hist <- function (x, digits = NULL, ...) 
{
  cat("\nData: ", x$data.name, 
      " (", x$n, " obs.); J(h)=",x$Jh, ", m=",x$m,"\n", sep = "")
  print(summary(as.data.frame(x[c("x", "y","l","u")])), digits = digits, 
        ...)
  invisible(x)
}

plot.hist  <- 
function (x, main = NULL, xlab = NULL, ylab = "Probability", lwd=1,col=1,
    zero.line = TRUE, scb=FALSE, ...) 
{
  if (is.null(xlab)) 
    xlab <- paste("N =", x$n, "J(h)=",formatC(x$Jh), "m=",x$m)
  if (is.null(main)) 
    main <- deparse(x$data.name)
  if(scb){
    plot.default(x$x,x$u, ylim=c(0,max(x$u)),
                 main = main, xlab = xlab, ylab = ylab, type='n',...)
    n = length(x$y)
    segments(x$x[1],0,x$x[1],x$y[1], lwd=lwd,col=col)
    segments(x$x[1],x$y[1],x$x[2],x$y[1], lwd=lwd,col=col)
    ul = x$u[1];ll=x$l[1];x0=x$x[1]
    for(i in 2:(n-1)){
      ul = c(ul,x$u[i-1],x$u[i])
      ll = c(ll,x$l[i-1],x$l[i])
      x0 = c(x0,x$x[i],x$x[i])
      segments(x$x[i],x$y[i-1],x$x[i],x$y[i], lwd=lwd,col=col)
      segments(x$x[i],x$y[i],x$x[i+1],x$y[i], lwd=lwd,col=col)
    }
    ul = c(ul,x$u[n])
    ll = c(ll,x$l[i])
    x0 = c(x0,2*x$x[n]-x$x[n-1])
    segments(x$x[n],x$y[n],x$x[n],0, lwd=lwd,col=col)
  
    cord.x = c(x0,rev(x0));
    cord.y = c(ll,rev(ul))
    polygon(cord.x,cord.y,border=col,lty=2)
  }else{
    hist(x$data,breaks=x$x,main=main, probability=TRUE,xlab=xlab,ylab=ylab,...)
  }
  
  if (zero.line) 
    abline(h = 0, lwd = 0.1, col = "gray")
  invisible(NULL)
}

lines.hist  <- 
function (x, lwd=1,col=1, zero.line = TRUE, ...) 
{
  n = length(x$y)
  n = length(x$y)
  segments(x$x[1],0,x$x[1],x$y[1], lwd=lwd,col=col)
  segments(x$x[1],x$y[1],x$x[2],x$y[1], lwd=lwd,col=col)
  ul = x$u[1];ll=x$l[1];x0=x$x[1]
  for(i in 2:(n-1)){
    ul = c(ul,x$u[i-1],x$u[i])
    ll = c(ll,x$l[i-1],x$l[i])
    x0 = c(x0,x$x[i],x$x[i])
    segments(x$x[i],x$y[i-1],x$x[i],x$y[i], lwd=lwd,col=col)
    segments(x$x[i],x$y[i],x$x[i+1],x$y[i], lwd=lwd,col=col)
  }
  ul = c(ul,x$u[n])
  ll = c(ll,x$l[i])
  x0 = c(x0,2*x$x[n]-x$x[n-1])
  segments(x$x[n],x$y[n],x$x[n],0, lwd=lwd,col=col)
  cord.x = c(x0,rev(x0));
  cord.y = c(ll,rev(ul))
  polygon(cord.x,cord.y,border=col,lty=2)

  invisible(NULL)
}
