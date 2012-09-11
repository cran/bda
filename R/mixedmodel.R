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
