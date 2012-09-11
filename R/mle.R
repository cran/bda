## This function is developed to compute the MLEs for different
## families of distributions based on data that are weighted, rounded
## or censored.

##  2012/11/06.
##  
mle <- function(x,w, type="Weighted", family="Gaussian"){
  if(!is.character(family)) family <- 'normal'
  if(!is.character(type)) type='Weighted'
  if(missing(w)) type='Weighted'
  type <- match.arg(tolower(type),
                    c('rc','right.censoring','re','rounding',
                      'wt','weighting',"weighted"))
  switch(type,
         wt = ,  weighted  = ,
         weighting = .wtmle(x=x,w=w,family=family),
         re = ,
         rounding  = .remle(x=x,w=w,family=family),
         rc = ,
         right.censoring  = .rcmle(x=x,w=w,family=family),
         stop("Data type  not supported")
         )
}


.remle <- function(x,w,family='normal'){
  family <- match.arg(tolower(family),
                      c('normal','gaussian', 'weibull',
                        'gamma'))
  stop("MLE for data with rounding error: not supported yet")
}

.wtmle <- function(x,w,family='normal'){
  family <- match.arg(tolower(family),
                      c('normal','gaussian', 'weibull',
                        'gamma'))
  if(missing(w)) w <- rep(1,length(x))
  x.wt <-  .weighting(x,w)
  switch(family,
         gaussian = ,
         normal = c(.wmean(x.wt),.wsd(x.wt)),
         stop("Family of distribution not supported")
         )
}

.rcmle <- function(x,w,family){
  if(any(x<=0)) stop("Invalid lifetime data in 'x'")
  n <- length(x)
  if(missing(w)) w <- rep(1,n)
  if(any(w!=0&&w!=1)) stop("Invalid censoring status in 'w'")
  if(length(w)!=n) stop("'x' and 'w' have different lengths")
  family <- match.arg(tolower(family),
                      c('normal','gaussian', 'weibull',
                        'gamma',"exponential"))
  switch(family,
         exponential = sum(x)/sum(w),
         weibull  =  .rcmle.weibull(x,w),
         stop("Family of distribution not supported")
         )
}

.rcmle.weibull <- function(x,w){
  n <- length(x);
  sele <- w==1
  stopifnot(sum(sele) > 5) # two few data points
  x0 <- x[sele]; rx <- rank(x0);
  Fhat <- (rx-0.3)/(n+0.4)
  lx <- log(x0); ly <- log(-log(1-Fhat))
  out <- lm(ly~lx)
  kappa <- out$coef[[2]];
  lambda <- exp(-out$coef[[1]]/kappa);
  pars <- c(kappa, lambda);
  .Fortran(.F_RcMleWeibull, as.double(x),
           as.double(w),as.integer(n),
           pars = as.double(pars))$pars
}
