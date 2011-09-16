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


.wmean <- function(x,f){
  if(length(x)!=length(f))stop("'x' and 'f' have different lengths.")
  sele1 = !is.na(x)
  sele2 = !is.na(f)
  sele = sele1&sele2
  x=x[sele]
  f=f[sele]
  sum(x*f)/sum(f)
}

.wsd <- function(x,f){
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
  theta = c(.wmean(X,F),.wsd(X,F));
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
