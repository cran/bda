### Functions:
##  bdde: binned data density estimation
## create an internal R object "bindata" with
## f -> frequencies, a -> lower limits, b -> upper limits

### dist: weibull=0,

bdde <- function(f, breaks, dist="normal", gridsize=512L,...)
{

  ## NA values are not allowed
  if(any(is.na(f)))
    stop("'f' contains missing value(s).")
  if(any(is.na(breaks)))
    stop("'f' contains missing value(s).")
  if(any(!is.finite(f)))
    stop("Infinite value of 'f' not allowed.")
  ## frequencies should to be positive integers or zero's.
  if(any(f != round(f))) stop("Non-integer frequency found.")
  if(any(f <0)) stop("Negative frequency found.")

  if(any(diff(breaks) <= 0))
    stop("'breaks' is not strictly increasing.")
  nbin <- length(f) # bin number
  if(length(breaks) != nbin + 1)
    stop("Lengths not match for 'f' and 'breaks'")
  if(min(breaks)<0) stop("Negative not allowed")
  if(!is.finite(breaks[nbin])) stop("Invalid value in 'breaks'")
  ## define lower and upper limits
  a <- breaks[1:nbin]; b <- breaks[-1]
  m <- 0.5 * (b + a); m0 <- m # save a copy 

  dist <- match.arg(tolower(dist),
                    c('normal', 'gaussian','dagum','weibull','gb',
                      'gld'))


  ##  draft a histogram
  if(!is.finite(breaks[nbin+1])){
    ##    hmax <- breaks[nbin]-breaks[nbin-1]
    hmax <- breaks[nbin]-breaks[1]
    breaks[nbin+1] <- breaks[nbin] + hmax
    m[nbin] <- hmax
  }
  binwidth <- diff(breaks)
  den <- f/(binwidth*sum(f))

  if(gridsize<10) stop("Gridsize might be too small")
  x <- seq(breaks[1L],breaks[nbin+1], length=gridsize)

  ##pars <- switch(dist,
  ##               nornal =,
  ##               gaussian =,
  ##               .weibullbde(f,a,b,x)
  ##               )

  pars <- .weibullbde(f,a,b)
  lambda <- pars[[1]]; kappa <- pars[[2]]
  
  y <- dweibull(x,kappa,lambda)
  mu <- lambda * gamma(1+1/kappa)
  s2 <- lambda^2*gamma(1+2/kappa)-mu^2
  
  x.hist <- structure(list(breaks = breaks,
                           counts = f,
                           intensities = den,
                           density = den,
                           mids = m, nclass=nbin,
                           equidist=FALSE),
                      class="histogram")  
  
  out <- structure(list(x=x,y=y, dist = dist,
                        hist=x.hist, mean=mu, sigma=sqrt(s2),
                        f=f, a=a, b=b, m=m0,
                        call = match.call()),
                   class = "bindata")
  invisible(out)
}

.weibullbde <- function(f,a,b){
  nbin <- length(f)
  if(!is.finite(b[nbin])) b[nbin] <- a[nbin]*2-a[nbin-1]
  n <- sum(f)
  y <- runif(n,rep(a,f),rep(b,f))
  pars <- .weibullplot(y)
  ##  pars <- .bdmle(f,a,b,dist=0,pars=pars0)
}

.bdmle <- function(f,a,b,dist,pars){
  nbin <- length(f); npar <- length(pars)
  .Fortran(.F_BDMLE, as.double(f), as.double(a),
           as.double(b),as.integer(nbin),
           pars = as.double(pars), as.integer(npar),
           as.integer(dist))$pars
}

.weibullplot <- function(x){
  tol <- 0.0000001
  if(any(x<0))
    stop("Negative value not allowed in Weibull distribution")
  n <- length(x)
  Fhat <- (rank(x)-0.3)/(n+0.4)
  ly <- log(-log(1-Fhat))
  lx <- log(x+tol)
  out <- lm(ly~lx)
  kappa <- out$coef[[2]]
  lambda <- exp(-out$coef[[1]]/kappa)
  c(lambda, kappa)
}

print.bindata  <- function(x, digits = NULL,...) 
{
  cat("\nCall:\n\t", deparse(x$call), sep = "")
  cat("\n\nTotal frequency:", round(sum(x$f)), sep="")
  tmp <- rbind(x$f,x$a,x$b,x$m)
  f <- x$f; m <- x$m
  nbin <- length(f)
  if(f[nbin]==0) m[nbin] <- m[nbin-1]
  n <- sum(f)
  mu <- ifelse(is.null(x$mean), sum(f*m)/n, x$mean)
  
  cat("\n\tMean:", round(mu,3), sep="")
  s <- ifelse(is.null(x$sigma), sqrt(sum((m-mu)^2*f)/(n-1)),x$sigma);
  cat("\n\tStdDev:", round(s,3), "\n\n", sep="")
  row.names(tmp) <- c("Freq","Lower.limit","Upper.limit","Midpoint")
  print(t(tmp))
  cat("\n\nDensity estimates\n")
  print(summary(as.data.frame(x[c("x", "y")])), digits = digits, 
        ...)
  invisible(x)
}

hist.bindata <- function(x,plot=TRUE,...){
  if(plot)
    plot(x$hist,...)
  invisible(x$hist)
}
