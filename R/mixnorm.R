##  Part of R package meada
##  Copyright (C) 2009-2010 Bin Wang
##
##  Unlimited use and distribution (see LICENCE).

##  Apply an EM algorithm to fit a model of normal mixture.


dmixnorm <- function(x,p,mu,s,na.rm=TRUE,pobj=NULL){
  if(!is.null(pobj)){
    p = pobj$Prop
    mu = pobj$Mean
    s = pobj$Std.Dev
  }
  sele = is.na(x)
  if(any(sele)){
    if(na.rm) x = x[!sele]
    else stop("Missing value(s) in 'x'!")
  }
  if(missing(p)) p=1
  nc = length(p)
  if(any(is.na(p))) warning("Missing values in 'p'.")
  if(any(p<0 || p>1))
    stop("Invalid component proportions!")
  if(missing(mu)) mu=rep(0,nc)
  if(missing(s)) s=rep(1,nc)
  if(any(is.na(s))) warning("Missing values in 's'.")
  if(any(is.na(mu))) warning("Missing values in 'mu'.")
  if(length(mu)!=nc||length(s)!=nc)
    stop("Parameters of different sizes")
  if(any(s<0)) stop("Invalid standard deviation(s)!")

  res = rep(0,length(x));
  for(i in 1:nc)
    res = res + p[i]*dnorm(x,mu[i],s[i]);
  res;
}

pmixnorm <- function(x,p,mu,s,na.rm=TRUE,pobj=NULL){
  if(!is.null(pobj)){
    p = pobj$Prop
    mu = pobj$Mean
    s = pobj$Std.Dev
  }
  sele = is.na(x)
  if(any(sele)){
    if(na.rm) x = x[!sele]
    else stop("Missing value(s) in 'x'!")
  }
  if(missing(p)) p=1
  nc = length(p)
  if(any(is.na(p))) warning("Missing values in 'p'.")
  if(any(p<0 || p>1))
    stop("Invalid component proportions!")
  if(missing(mu)) mu=rep(0,nc)
  if(missing(s)) s=rep(1,nc)
  if(any(is.na(s))) warning("Missing values in 's'.")
  if(any(is.na(mu))) warning("Missing values in 'mu'.")
  if(length(mu)!=nc||length(s)!=nc)
    stop("Parameters of different sizes")
  if(any(s<0)) stop("Invalid standard deviation(s)!")

  res = rep(0,length(x));
  for(i in 1:nc)
    res = res + p[i]*pnorm(x,mu[i],s[i]);
  res;
}

rmixnorm <- function(n,p,mu,s,pobj=NULL){
  if(!is.null(pobj)){
    print(pobj)
    p = pobj$Prop
    mu = pobj$Mean
    s = pobj$Std.Dev
  }
  if(n<1) stop("Invalid sample size 'n'!")
  if(missing(p)) p=1
  nc = length(p)
  if(any(is.na(p))) warning("Missing values in 'p'.")
  if(any(p<0 || p>1))
    stop("Invalid component proportions!")
  if(missing(mu)) mu=rep(0,nc)
  if(missing(s)) s=rep(1,nc)
  if(any(is.na(s))) warning("Missing values in 's'.")
  if(any(is.na(mu))) warning("Missing values in 'mu'.")
  if(length(mu)!=nc||length(s)!=nc)
    stop("Parameters of different sizes")
  if(any(s<0)) stop("Invalid standard deviation(s)!")

  res = NULL;
  tmp = runif(n)
  sump = cumsum(p)

  for(i in 1:(nc-1)){
    m = sum(tmp<sump[nc+1-i] & tmp>=sump[nc-i])
    if(m>0) res = c(res, rnorm(m,mu[nc+1-i],s[nc+1-i]))
  }
  m = sum(tmp < sump[1])
  if(m>0) res = c(res, rnorm(m,mu[1],s[1]))
  res;
}

qmixnorm <- function(x,p,mu,s,na.rm=TRUE,pobj=NULL){
  if(!is.null(pobj)){
    print(pobj)
    p = pobj$Prop
    mu = pobj$Mean
    s = pobj$Std.Dev
  }
  sele = is.na(x)
  if(any(sele)){
    if(na.rm) x = x[!sele]
    else stop("Missing value(s) in 'x'!")
  }
  
  if(missing(p)) p=1
  nc = length(p)
  if(any(is.na(p))) warning("Missing values in 'p'.")
  if(any(p<0 || p>1) || sum(p)!=1)
    stop("Invalid component proportions!")
  if(missing(mu)) mu=rep(0,nc)
  if(missing(s)) s=rep(1,nc)
  if(any(is.na(s))) warning("Missing values in 's'.")
  if(any(is.na(mu))) warning("Missing values in 'mu'.")
  if(length(mu)!=nc||length(s)!=nc)
    stop("Parameters of different sizes")
  if(any(s<0)) stop("Invalid standard deviation(s)!")

  n = length(x); tmp = c(1:n)
  ox = order(x); x0 = x[ox]; xo = match(tmp,ox)
  
  res = rep(0,n); #initialize
  sele0 = x0==0; res[sele0] = -Inf;
  sele1 = x0==1; res[sele1] = Inf
  sele2 = x0<0|x0>1; res[sele2] = NA;
  sele = !(sele0|sele1|sele2)
  y1 = x0[sele]; # sorted, no zeros and ones
  m = sum(sele);
  xmin = NULL; xmax = NULL
  for(i in 1:nc){
    xmin = c(xmin, mu[i]-4*s[i])
    xmax = c(xmax, mu[i]+4*s[i])
  }
  l = 10000
  x0 = seq(min(xmin),max(xmax), length=l)
  y0 = pmixnorm(x0,p,mu,s)
  out = .Fortran(.F_findxyz,
    as.double(y0), as.double(x0), as.integer(l),
    x= as.double(y1),as.integer(m))
  res[sele] = out$x
  res[xo];
}
