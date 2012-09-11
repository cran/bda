gof.test <- function(object, x, ...)
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
      gamma  = pgamma(out$x[-1], .mle2.gamma(x)),
      weibull  = pweibull(out$x[-1], .mle2.weibull(x))
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
      gamma  = pgamma(x0, .mle2.gamma(x)),
      weibull  = pweibull(x0, .mle2.weibull(x))
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


.mle2.gamma <- function(x){
  x = x[!is.na(x)]
  .Fortran(.F_FitGamma, as.double(x),
           as.integer(length(x)), as.double(rep(0,2)))[[3]]
}

.mle2.weibull <- function(x){
  x = x[!is.na(x)]
  .Fortran(.F_FitWeibull, as.double(x),
           as.integer(length(x)), as.double(rep(0,2)))[[3]]
}


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
  ##  if(trials<500) warning("Trial number might be too small.")
  ##  check data (x,y)
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
  x0 <- seq(rngx[1],rngx[2],length=256);
  Fx1 <- edf(x,xgrid=x0);
  Fx2 <- edf(y,xgrid=x0);
  sele <- is.na(Fx1$y) | is.na(Fx2$y)
  max(abs((Fx1$y-Fx2$y)[!sele]))
}

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
