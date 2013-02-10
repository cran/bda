### Functions:
##  bdde: binned data density estimation
## create an internal R object "bindata" with
## f -> frequencies, a -> lower limits, b -> upper limits

### dist: weibull=0,

bdde <- function(f, breaks, dist="ewd", gridsize=512L,iter=100,...)
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
  dist <- match.arg(tolower(dist),
                    c('ewd','normal', 'gaussian','dagum',
                      'weibull','gb','gld'))

  ## if top coded, compute an upper bound.
  xmax <- max(breaks)
  if(!is.finite(xmax)){
    top.coded <- TRUE
  }else{
    top.coded <- FALSE
  }
  a <- breaks[1:nbin]; b <- breaks[-1]
  b0 <- b; # copy data for display
  if(top.coded){
    b[nbin] <- 2 * max(a)
    tmp <- .bweibullplot(f,a,b,top.coded)
    lambda <- tmp$pars[[1]]; kappa <- tmp$pars[[2]]
    mu <- lambda * gamma(1+1/kappa)
    s2 <- lambda^2*gamma(1+2/kappa)-mu^2
    xmax <- mu + 3.0 * sqrt(s2)
    xmax <- max(xmax, 2*max(a))
    b[nbin] <- xmax; breaks[nbin+1] <- xmax
  }
  
  ##  draft a histogram
  n0 <- which(match(rev(f)==0,FALSE)==1)[1]-1
  if(n0>0){
    nbin <- nbin - n0
    f <- f[1:nbin]
    a <- a[1:nbin]
    b <- b[1:nbin]; b0 <- b  # update date for display
    breaks <- breaks[1:(nbin+1)]
    xmax <- max(breaks)
  }
  m <- 0.5 * (b + a); 
  binwidth <- diff(breaks)
  den <- f/(binwidth*sum(f))

  x.hist <- structure(list(breaks = breaks,
                           counts = f,
                           intensities = den,
                           density = den,
                           mids = m, nclass=nbin,
                           equidist=FALSE),
                      class="histogram")  

  #####################################################
  if(gridsize<10) stop("Gridsize might be too small")
  x <- seq(breaks[1L],breaks[nbin+1], length=gridsize)

  if(dist=="weibull"){
    fixed <- TRUE
    pars <- .ewdbde(f,a,b,iter,top.coded,fixed)
    lambda <- pars[[2]]; kappa <- pars[[1]]
    pars <- c(kappa, lambda)
    y <- dweibull(x,kappa,lambda)
    sele <- is.finite(y)
    if(any(!sele)){
      y <- y[sele]; x <- x[sele]
    }
    sele <- y>0
    if(any(sele)){
      y <- y[sele]; x <- x[sele]
    }
    mu <- lambda * gamma(1+1/kappa)
    s2 <- lambda^2*gamma(1+2/kappa)-mu^2
  }else if(dist=='ewd'){
    fixed <- FALSE
    pars <- .ewdbde(f,a,b,iter,top.coded,fixed)
    alpha <- pars[[3]];
    lambda <- pars[[2]];
    kappa <- pars[[1]]
    pars <- c(kappa, lambda,alpha)
    y  <- .dewd(x, kappa,lambda,alpha)
    sele <- is.finite(y)
    if(any(!sele)){
      y <- y[sele]; x <- x[sele]
    }
    sele <- y>0
    if(any(sele)){
      y <- y[sele]; x <- x[sele]
    }
    Fx <- .pewd(x, kappa,lambda,alpha)
    Px <- diff(Fx); Px <- c(Px,1-sum(Px))
    mu <- sum(x*Px)
    s2 <- sum((x-mu)^2*Px)
  }else if(dist=='dagum'){
    pars <- .dagumbde(f,a,b,iter,top.coded)
    alpha <- pars[[3]];
    lambda <- pars[[2]];
    kappa <- pars[[1]]
    y <- .ddagum(x,kappa,lambda,alpha)
    sele <- y>0
    if(any(sele)){
      y <- y[sele]; x <- x[sele]
    }
    Fx <- .pdagum(x,kappa,lambda,alpha)
    Px <- diff(Fx); Px <- c(Px,1-sum(Px))
    mu <- sum(x*Px)
    s2 <- sum((x-mu)^2*Px)
  }else if(dist=='gb'){
    sele <- f==0  # bins with zero counts
    f1 <- f[!sele]
    a1 <- a[!sele]
    b1 <- b[!sele]
    xmin <- min(a1);
    x <- seq(xmin,xmax, length=gridsize)
    a2 <- (a1-xmin)/(xmax-xmin)
    b2 <- (b1-xmin)/(xmax-xmin)
    pars <- .betabde(f1,a2,b2,iter,top.coded)
    pa <- pars[[1]];
    pb <- pars[[2]];

    y <- dbeta((x-xmin)/(xmax-xmin),pa,pb)/(xmax-xmin)
    sele <- is.finite(y)
    if(any(!sele)){
      y <- y[sele]; x <- x[sele]
    }
    sele <- y>0
    if(any(sele)){
      y <- y[sele]; x <- x[sele]
    }
    mu <- pa/(pa+pb)*(xmax-xmin)
    s2 <- pa*pb/(pa+pb)^2/(pa+pb+1)*(xmax-xmin)^2
  }else stop("Distribution not supported")


  
  out <- structure(list(x=x,y=y, dist = dist,
                        hist=x.hist,
                        dist=dist,pars=pars,
                        mean=mu, sigma=sqrt(s2),
                        f=f, a=a, b=b0,
                        call = match.call()),
                   class = "bindata")
  invisible(out)
}

######################################################
######################################################
.betabde <- function(f,a,b,iter,top.coded)
{
  alphas <- rep(1,iter)
  tmp <- apply(as.matrix(alphas,ncol=1),1,.betaplot,
               f=f,a=a,b=b,top.coded=top.coded)
  tmp[-1,which(tmp[1,]==max(tmp[1,]))[[1]]]
}

.betaplot <- function(x,f,a,b,top.coded)
  {
    n <- sum(f); nbin <- length(f)
    x <- runif(n,rep(a,f),rep(b,f))
    mu <- mean(x); s2 <- var(x)
    alpha <- mu*(mu*(1-mu)/s2-1)
    pbeta <- alpha/mu-alpha
    Fa <- pbeta(a,alpha,pbeta)
    if(top.coded) Fa <- c(Fa, 1)
    else Fa <- c(Fa, pbeta(b[nbin],alpha,pbeta))
    llk <- sum(log(diff(Fa)))
    c(llk,alpha, pbeta)
  }

######################################################
######################################################

.dagumbde <- function(f,a,b,iter,top.coded){
  alphas <- c(runif(iter), runif(iter,1,5))
  tmp <- apply(as.matrix(alphas,ncol=1),1,.dagumplot,
               f=f,a=a,b=b)
  tmp[1:3,which(tmp[4,]==max(tmp[4,]))[[1]]]
}

.dagumplot <- function(alpha,f,a,b)
  {
    tol <- 0.000001
    n <- sum(f)
    x <- runif(n,rep(a,f),rep(b,f))
    Fhat <- (rank(x)-0.3)/(n+0.4)
    ly <- log(-1-log(Fhat)/alpha+tol)
    lx <- log(x+tol)
    out <- lm(ly~lx)
    kappa <- -out$coef[[2]]
    lambda <- exp(out$coef[[1]]/kappa)
    llk <- .dagumllk(f,a,b,alpha,kappa,lambda)
    pars <- c(kappa, lambda, alpha,llk)
  }

.dagumllk <- function(f,a,b,alpha,kappa,lambda){
  if(is.na(lambda)||is.na(kappa)){
    llk <- -9999999999999.99
  }else if(lambda<=0||kappa<=0){
    llk <- -9999999999999.99
  }else{
    tol <- 0.0000001
    Fa <- .pdagum(a, kappa,lambda,alpha)
    Fb <- .pdagum(b, kappa,lambda,alpha)
    llk <- sum(f*log(Fb-Fa+tol))
  }
}

.ewdbde <- function(f,a,b,iter,top.coded,fixed){
  if(fixed)
    alphas <- rep(1,iter)
  else
    alphas <- c(runif(iter), runif(iter,1,5))
  tmp <- apply(as.matrix(alphas,ncol=1),1,.ewdplot,
               f=f,a=a,b=b)
  tmp[1:3,which(tmp[4,]==max(tmp[4,]))[[1]]]
}

.ewdplot <- function(alpha,f,a,b)
  {
    tol <- 0.000001
    n <- sum(f)
    x <- runif(n,rep(a,f),rep(b,f))
    Fhat <- (rank(x)-0.3)/(n+0.4)
    ly <- log(-log(1-Fhat^(1/alpha))+tol)
    lx <- log(x+tol)
    out <- lm(ly~lx)
    kappa <- out$coef[[2]]
    lambda <- exp(-out$coef[[1]]/kappa)
    llk <- .ewdllk(f,a,b,alpha,kappa,lambda)
    pars <- c(kappa, lambda, alpha,llk)
  }

.ewdllk <- function(f,a,b,alpha,kappa,lambda){
  if(is.na(lambda)||is.na(kappa)){
    llk <- -9999999999999.99
  }else if(kappa<=0||lambda<=0){
    llk <- -999999999999
  }else{
    tol <- 0.00000000001
    Fa <- .pewd(a, kappa,lambda,alpha)
    Fb <- .pewd(a, kappa,lambda,alpha)
    llk <- sum(f*log(Fb-Fa+tol))
  }
  llk
}


##################################################
##  Weibull(x,lambda,kappa)
##  f(x) = k/l*(x/l)^(k-1)*exp(-(x/l)^k)
##  in R Weibill(x,kappa,lambda)
.bweibullplot <- function(f,a,b,top.coded){
  nbin <- length(f)
  n <- sum(f)
  y <- runif(n,rep(a,f),rep(b,f))
  pars <- .weibullplot(y)
  Fa <- pweibull(a,pars[[2]],pars[[1]])
  if(top.coded) Fa <- c(Fa, 1)
  else Fa <- c(Fa, pweibull(b[nbin],pars[[2]],pars[[1]]))
  llk <- sum(log(diff(Fa)))
  list(pars=pars,llk=llk)
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

.weibullbde <- function(f,a,b){
  nbin <- length(f)
  if(!is.finite(b[nbin])) b[nbin] <- a[nbin]*2-a[nbin-1]
  n <- sum(f)
  y <- runif(n,rep(a,f),rep(b,f))
  pars <- .weibullplot(y)
  ##.bdmle(f,a,b,dist=0,pars=pars)
}

.bdmle <- function(f,a,b,dist,pars){
  nbin <- length(f); npar <- length(pars)
  .Fortran(.F_BDMLE, as.double(f), as.double(a),
           as.double(b),as.integer(nbin),
           pars = as.double(pars), as.integer(npar),
           as.integer(dist))$pars
}

#############################################################
print.bindata  <- function(x, digits = NULL,...) 
{
  cat("\nCall:\n\t", deparse(x$call), sep = "")
  cat("\n\nDistribution:\t", x$dist, "(",  sep="")
  cat(round(x$pars,3), ")",  sep=" ")

  tmp <- rbind(x$f,x$a,x$b)
  f <- x$f; m <- x$m
  nbin <- length(f)
  if(f[nbin]==0) m[nbin] <- m[nbin-1]
  n <- sum(f)
  mu <- ifelse(is.null(x$mean), sum(f*m)/n, x$mean)
  
  cat("\n Size: \t", round(sum(x$f)), sep="")
  cat("\n Mean: \t", round(mu,3), sep="")
  s <- ifelse(is.null(x$sigma), sqrt(sum((m-mu)^2*f)/(n-1)),x$sigma);
  cat("\n StDev: ", round(s,3), "\n\n", sep="")
  row.names(tmp) <- c("Freq","Lower.limit","Upper.limit")
  print(t(tmp))
  cat("\n\nDensity estimates\n")
  print(summary(as.data.frame(x[c("x", "y")])), digits = digits, 
        ...)
  invisible(x)
}

hist.bindata <- function(x,plot=TRUE,main,...){
  if(missing(main)) main <- x$dist
  if(plot)
    plot(x$hist,main=main,...)
  invisible(x$hist)
}

.pewd <- function(x,kappa,lambda,alpha)
  {
    n <- length(x)
    res <- rep(0,n)
    sele1 <- x > 0
    sele2 <- is.finite(x); 
    sele3 <- is.na(x);
    if(any(sele3)) res[sele3] <- NA
    if(any(!sele3&!sele2)) res[!sele3&!sele2] <- 1.0
    x <- x[sele1&sele2]
    stopifnot(kappa>0&&lambda>0&&alpha>0)
    tmp <- (1-exp(-(x/lambda)^kappa))^alpha
    res[sele1&sele2] <- tmp
    res
  }

.dewd <- function(x,kappa,lambda,alpha)
  {
    n <- length(x)
    res <- rep(0,n)
    sele1 <- x > 0
    sele2 <- is.finite(x); 
    sele3 <- is.na(x);
    if(any(sele3)) res[sele3] <- NA
    x <- x[sele1&sele2]
    stopifnot(kappa>0&&lambda>0&&alpha>0)
    tmp <- alpha * (1-exp(-(x/lambda)^kappa))^(alpha-1) *
      exp(-(x/lambda)^kappa) * kappa * (x/lambda)^(kappa-1)/lambda
    res[sele1&sele2] <- tmp
    res
  }

.pdagum <- function(x,a,b,p)
  {
    n <- length(x)
    res <- rep(0,n)
    sele1 <- x > 0
    sele2 <- is.finite(x); 
    sele3 <- is.na(x);
    if(any(sele3)) res[sele3] <- NA
    if(any(!sele3&!sele2)) res[!sele3&!sele2] <- 1.0
    x <- x[sele1&sele2]
    stopifnot(a>0&&b>0&&p>0)
    tmp <- (1+(x/b)^(-a))^(-p)
    res[sele1&sele2] <- tmp
    res
  }

.ddagum <- function(x,a,b,p)
  {
    n <- length(x)
    res <- rep(0,n)
    sele1 <- x > 0
    sele2 <- is.finite(x); 
    sele3 <- is.na(x);
    if(any(sele3)) res[sele3] <- NA
    x <- x[sele1&sele2]
    stopifnot(a>0&&b>0&&p>0)
    tmp <- (a*p)/x*(x/b)^(a*p)*((x/b)^a+1)^(-p-1)
    res[sele1&sele2] <- tmp
    res
  }
