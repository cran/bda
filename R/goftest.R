gof.test <- function(object, x, ...)
UseMethod("gof")

gof <- function(object, x, ...)
UseMethod("gof")

gof.default <- function(object, x, type='chisq',...)
{
  nam <- deparse(substitute(x));
  ##  check data
   if (!is.numeric(x) && !is.complex(x) && !is.logical(x)) {
    warning("'x' is not numeric or logical: returning NA")
    return(NA_real_)
  }
  x = x[!is.na(x)]; n=length(x)
  if (!is.character(object)) {
    warning("'object' is not a distribution name: returning NA")
    return(NA_real_)
  }
  dist = match.arg(tolower(object),
    c("normal","gamma","weibull"))

  anam = paste("Data does not come from a ", dist, "distribution")

  out = .bincount(x)
  nbins = length(out$y)
  np = switch(dist,
    normal = 2,
    gamma  = 2,
    weibull  = 2,
    stop("Not supported yet.")
    )
  DF = nbins - np -1
  if(DF>0 && type == 'chisq'){
    Fx = switch(dist,
      normal = pnorm(out$x[-1],mean(x),sd(x)),
      gamma  = pgamma(out$x[-1], mle.gamma(x)),
      weibull  = pweibull(out$x[-1], mle.weibull(x))
      )
    EXs = c(Fx[1],diff(Fx))
    es = n*c(EXs,1-sum(EXs))
    D2 = sum((out$y-es)^2/es)
    p2 = ifelse(DF>0,pchisq(D2,DF,lower.tail=FALSE),NA)
    mname = "Chi-Square Test for goodness-of-fit.";
  }else{
    mname = "Kolmogorov-Smirnov Test for goodness-of-test."
    x0 = seq(min(x),max(x),length=100)
    Fx = switch(dist,
      normal = pnorm(x0,mean(x),sd(x)),
      gamma  = pgamma(x0, mle.gamma(x)),
      weibull  = pweibull(x0, mle.weibull(x))
      )
    Fn = edf(x,xgrid=x0)
    D2 = max(abs(Fx-Fn$y))
    ft = D2*sqrt(n);
    p2 = .Fortran(.F_kspvalue, p=as.double(ft))$p
  }
  ##  If paraetric, using chisq-test instead of KS-test
  RVAL <- list(statistic = c(ChiSq=D2),#c(D = D,Chisq=D2),
               p.value = p2,#c(p1,p2),
               parameter = c(df=DF),
               method = mname, 
               alternative = anam,
               data.name = nam)
  class(RVAL) <- "htest"
  return(RVAL)
}

gof.em <- function(object, x,type='chisq',...) {
  nam <- object$data.name
  if(missing(x)) x = rep(object$X, object$F)
  ##  check data
   if (!is.numeric(x) && !is.complex(x) && !is.logical(x)) {
    warning("'x' is not numeric or logical: returning NA")
    return(NA_real_)
  }
  x = x[!is.na(x)]; n=length(x)
  type = match.arg(tolower(type), c("ks","chisq"))

  anam = paste("The EM doesn't well fit the data")
  out = .bincount(x)
  nbins = length(out$y)
  DF = nbins-prod(dim(object$Para))
  if(DF>0 && type == 'chisq'){
    Fx = pmixnorm(out$x[-1], object$Para[,1], object$Para[,2],object$Para[,3])
    EXs = c(Fx[1],diff(Fx))
    es = n*c(EXs,1-sum(EXs))
    D2 = sum((out$y-es)^2/es)
    p2 = ifelse(DF>0,pchisq(D2,DF,lower.tail=FALSE),NA)
    mname = "Chi-Square Test for goodness-of-fit.";
  }else{
    mname = "Kolmogorov-Smirnov Test for goodness-of-test."
    x0 = seq(min(x),max(x),length=100)
    Fx = pmixnorm(x0, object$Para[,1], object$Para[,2],object$Para[,3])
    Fn = edf(x,xgrid=x0)
    D2 = max(abs(Fx-Fn$y))
    ft = D2*sqrt(n);
    p2 = .Fortran(.F_kspvalue, p=as.double(ft))$p
  }
  ##  If paraetric, using chisq-test instead of KS-test
  RVAL <- list(statistic = c(ChiSq=D2),#c(D = D,Chisq=D2),
               p.value = p2,#c(p1,p2),
               parameter = c(df=DF),
               method = mname, 
               alternative = anam,
               data.name = nam)
  class(RVAL) <- "htest"
  return(RVAL)
}

gof.bde <- function(object, x,...) {
  nam <- object$data.name
  if(missing(x)) x = rep(object$X, object$F)
  ##  check data
   if (!is.numeric(x) && !is.complex(x) && !is.logical(x)) {
    warning("'x' is not numeric or logical: returning NA")
    return(NA_real_)
  }
  x = x[!is.na(x)]; n=length(x)
  
  x0 = object$x; bws = diff(x0); bws = c(bws[1],bws)
  F0 = cumsum(object$y*bws)
  Fn = edf(x,xgrid=x0)$y;
  D = max(abs(F0-Fn))
  ft = D*sqrt(n);
  out = .Fortran(.F_kspvalue, p=as.double(ft))
  anam = paste("The binned density estimate doesn't well fit the data")
  ##conf.int; null.value; parameter

  RVAL <- list(statistic = c(D = D), p.value = out$p,
               method = "Kolmogorov-Smirnov Test for goodness-of-fit.", 
               alternative = anam,
               data.name = nam)
  class(RVAL) <- "htest"
  return(RVAL)
}


###  subroutines to be called


.lecount <- function(x,y) sum(y<=x)

.bincount <- function(x){
  x=x[!is.na(x)]; s = sd(x); mu = mean(x);
  n=length(x)
  x0 = seq(mu-5.7*s,mu+6*s,by=0.3*s)
  ps = apply(matrix(x0,ncol=1),1,FUN=.lecount,y=x)
  cts = diff(ps)
  bounds = NULL; counts=NULL;
  k = length(cts)
  for(i in 1:(k-1)){
    if(cts[i]>=5){
      bounds = c(bounds,x0[i])
      counts = c(counts,cts[i])
      ul = x0[i+1]
    }else{
      cts[i+1] = cts[i+1]+cts[i]
    }
  }
  m = n - sum(counts)
  if(m>5){
    bounds = c(bounds,ul)
    counts = c(counts,m)
  }else{
    counts[length(counts)] = counts[length(counts)]+m
  }
  list(x=bounds,y=counts);
}


