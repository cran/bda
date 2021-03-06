## To do: for y = r(x) + e when the data is very large, we can group
## the data and then fit lpr based on the grouped data.  Use a
## masurement error model instead.  That says, we fit y.bar = r(x.bar)
## + e, and x = x.bar + e2.  Laplace error can be assumed for e2 to
## speed up the algorithm.

## we use four method to compute the variance of r(x).

## method 1) Larry Wasserman--nearly unbiased.  This method based on
## an lps object;

## method 2) Rice 1984

## method 3) Gasser et al (1986) -- a variation of method 3.

## method 4) For heteroscedastic errors. Need to estimate based on an
## lpr object. Yu and Jones (2004).

lps.variance <- function(y,x, bw, method="Rice"){
    method <- match.arg(tolower(method),
                        c("wasserman","rice","gasser",
                          "heteroscedastic"))
    sele <- is.na(y)|is.na(x)
    y <- y[!sele]; x <- x[!sele]
    
    n <- length(y)
    ox <- order(x)
    oy <- y[ox];
    ox <- sort(x)
    s2 <- 0.5* mean((diff(oy))^2)

    if(method == "wasserman"){
        ## we fit lpr
        if(missing(bw)){
            fhat <- lpsmooth(y=y,x=x,lscv=TRUE,sd.y=sqrt(s2))
            bw <- fhat$pars$bw
        }else{
            fhat <- lpsmooth(y=y,x=x,bw=bw,sd.y=sqrt(s2))
        }
        tmp <- .DesignMatrix(x,bw)
        nu1 <- sum(diag(tmp))
        dl2 <- apply(tmp^2,2,sum)
        nu2 <- sum(dl2)
        res <- sum((y-fhat$y0)^2)/(n-2*nu1+nu2)
    }else if(method == "rice"){
        res <- s2
    }else if(method == "gasser"){
        ox <- order(x)
        oy <- y[ox];
        ox <- sort(x)
        ai <- (ox[3:n] - ox[2:(n-1)])/(ox[3:n] - ox[1:(n-2)])
        bi <- (ox[2:(n-1)] - ox[1:(n-2)])/(ox[3:n] - ox[1:(n-2)])
        ci2 <- 1/(ai^2+bi^2+1)
        di <- ai * oy[1:(n-2)] + bi * oy[3:n] - oy[2:(n-1)]
        res <- mean(ci2*di^2)
    }else{
        ## we fit an lpr to get temporary results.  sd.y does not matter.
        fhat <- lpsmooth(y=y,x=x,lscv=TRUE,sd.y=sqrt(s2))
        bw1 <- fhat$pars$bw
        from <- min(x); to <- max(x)
        fhat <- .Fortran(.F_lpsmooth,
                         fx = as.double(x), as.integer(n),
                         as.double(x),y0=as.double(y), as.integer(n),
                         bw = as.double(bw1), as.integer(0),
                         as.double(c(from, to)), as.integer(0),
                         ellx = double(n), kappa=double(1))
        z <- log((y - fhat$y0)^2)
        if(missing(bw)){
            qhat <- lpsmooth(y=z,x=x,lscv=TRUE)
        }else{
            qhat <- lpsmooth(y=z,x=x,bw=bw,lscv=FALSE)
        }
        res <- list(x=qhat$x, y=exp(qhat$y),bw=qhat$pars$bw)
    }
    res
}

## this function can be used to compute the design (smoothing) matrix
.DesignMatrix <- function(x,bw){
    stopifnot(bw > 0)
    x <- x[!is.na(x)]; # remove missing values
    ndim = length(x)
    
    L = rep(0,ndim*ndim)
    DM = .Fortran(.F_DesignMatrix, as.double(x), as.integer(ndim),
        as.double(bw), dm=as.double(L))$dm
    matrix(DM, ncol = ndim)
}



## 2015/08/3: add two options lscv=T/F and adaptive=T/F.

## If lscv = FALSE, use the given bandwidth to fit lpr directly.

## If lscv = TRUE and adaptive = FALSE, compute lscv bandwidth and fit
## lpr.  Initial bandwidth should be given.

## If lscv = TRUE and adaptive = TURE, compute lscv bandwidth and then
## compute varying smoothing parameter, then fit lpr.  This algorithm
## could be extremeely slow when the sample size is very large.

lpsmooth <-
    function(y,x,bw,sd.y,lscv=FALSE,adaptive=FALSE,
             from,to,gridsize,conf.level=0.95)
{
    out <- .lps(y=y,x=x,bw=bw,sd.y=sd.y,sd.x=0,
                lscv=lscv,adaptive=adaptive,
                from=from,to=to,gridsize=gridsize,
                conf.level=conf.level)
}


.lps <- function(y,x,bw,sd.y,sd.x,lscv=FALSE, adaptive=FALSE,
                 from,to,gridsize,conf.level=0.95)
{
    stopifnot(conf.level<1 & conf.level >0)
    if(length(y) != length(x))
        stop("'x' and 'y' have different lengths")
    sele <- is.na(y)|is.na(x)
    y <- y[!sele]; x <- x[!sele]
    
    size.limit <- 1000
    n <- length(y)
    if(missing(from)) from <- min(x)
    if(missing(to)) to <- max(x)
    stopifnot(to > from)
    
    if(missing(bw)){
        if(n < size.limit){
            lscv <- 1
            bw <- bw.nrd(x)
            ##adaptive = 1  # use adaptive bandwidth selector
            lscv <- 0
        }else{
            stop("'bw' is missing.")
        }
    }else{
        stopifnot(bw>0)
    }
    
    if(n >= size.limit){
        if(lscv || adaptive)
            stop("sample size is too large.")
        lscv <- 0
        adaptive <- 0
    }else{
        lscv <- ifelse(lscv, 1, 0)
        adaptive <- ifelse(adaptive, 1, 0)  
    }
    
    ##  Compute the variance based on the raw data
    if(missing(sd.y)){
        sd.y <- sqrt(lps.variance(y=y,x=x,bw=bw))
    }else{
        stopifnot(is.numeric(sd.y))
        stopifnot(!any(sd.y < 0))
    }

    sd.x <- 0
    ##if(missing(sd.x)){
    ##    sd.x <- 0
    ##}else{
    ##    stopifnot(is.numeric(sd.x))
    ##    stopifnot(sd.x >= 0)
    ##}
    
    if(missing(gridsize)) gridsize <- 512L
    stopifnot(gridsize > 10)
    gpoints <- seq(from, to, length=gridsize)

    n <- length(x)
    stopifnot(length(y) == n)
    if(any(is.na(x)|is.na(y)))
        stop("Missing value(s) in 'x' and/or 'y'")
    if(any(!is.finite(x)|!is.finite(y)))
        stop("Inifite value(s) in 'x' and/or 'y'")

    out <- .Fortran(.F_lpsmooth,
                    fx = as.double(gpoints),
                    as.integer(gridsize),
                    as.double(x),
                    y=as.double(y),
                    as.integer(n),
                    bw = as.double(bw),
                    as.integer(lscv),
                    as.double(c(from, to)),
                    as.integer(adaptive),
                    ellx = double(gridsize),
                    kappa=double(1))
    y0 <- out$y
    y1 <- out$fx
    
    if(is.null(y1)){
        if(adaptive==0)
            cat("bw=",bw,"\n")
        stop("fail to converge")
    }
    if(length(y1)!=gridsize){
        if(adaptive==0)
            cat("bw=",bw,"\n")
        stop("fail to converge")
    }

    
    ellx <- out$ellx
    kappa <- out$kappa
    
    cv <-  .Fortran(.F_tubecv,
                    cv=as.double(out$kappa),
                    as.double(conf.level))$cv

    bw <- out$bw
    pars <- list(cv=cv, kappa=kappa, bw=bw,sigma=sd.y)

    y <- y1
    sele1 = is.na(y) | !is.finite(y)
    if(any(sele1)) y[sele1] = 0.0
    
    MOE <- cv * ellx * sd.y
    ll = y - MOE; ##ll[ll<0] <- 0
    ul = y + MOE;
    
    
    structure(list(y = y, x = gpoints,
                   x0 = x, y0 = y0,
                   conf.level = conf.level,
                   pars = pars,
                   ucb=ul, lcb=ll,
                   call = match.call()
                   ),
              class = 'scb')
}

print.scb <- function(x,...){
    print(x$pars,...)
}

.histonpr <- function(x,from,to,gridsize,conf.level=0.95)
    {
        if(missing(conf.level)) conf.level <- 0.95
                
        F <- x$counts;
        X <- x$mids;
        nbin <- length(X);
        n <- sum(F)
        a <- x$breaks[1]; b <- x$breaks[nbin+1];

        ##        if(missing(bw)){
        lscv <- 1
        A <- x$breaks[1:nbin]; B <- x$breaks[-1]
        x0 <- runif(n,rep(A,F),rep(B,F))
        bw <- bw.nrd(x0)
        adaptive = 1  # use adaptive bandwidth selector
        ##        }else{
        ##            lscv <- 0
        ##            adaptive = 0  # don't use adaptive bandwidth selector
        ##        }

        if(missing(from)) from <- a#min(X)#a - 3*bw
        if(missing(to)) to <- b#max(X)#b + 3*bw
        stopifnot(to > from)

        if(missing(gridsize)) gridsize <- 512L
        stopifnot(gridsize > 10)
        ##  root-transformation
        if(!x$equalwidth){
            ##stop("Not applicable for unequal width bins")
            wd <- diff(x$breaks)
            wdmean <- mean(wd)
            Y <- sqrt(nbin/n*(F*wdmean/wd+0.25))
        }else{
            Y <- sqrt(nbin/n*(F+0.25))
        }
        ## add two more points on the two sides for the boundary
        ## effects
        tmp <- sqrt(nbin/n*0.5)
        Y <- c(tmp, Y, tmp)
        tmp <- diff(X)
        X <- c(X[1]-tmp[1], X, max(X)+rev(tmp)[1])
        
        ##print(x$counts)
        
        ## call npr(...) to estimate the nonparametric regression
        ## function
        out <- lpsmooth(y=Y, x=X, bw=bw, sd.y= 0.5*sqrt(nbin/n),
                   from=from, to=to, gridsize=gridsize,
                   conf.level=conf.level)
        y <- out$y; x <- out$x
        ll <- out$lcb; ul <- out$ucb;
        y <- y*y; tot.mass <- sum(y)*(to-from)/gridsize;
        ll[ll<0] <- 0 # density can not be negative
        ll <- ll*ll/tot.mass
        ul <- ul*ul/tot.mass
        
        f0 <- y/tot.mass
        x0 <- x;

        pars <- list(cv=out$cv, kappa=out$kappa, bw=out$bw, npar=NULL)

        structure(list(y = y/tot.mass, x = x,
                       conf.level = out$conf.level,
                       pars = pars,
                       ucb=ul, lcb=ll,
                       call = match.call()
                       ),
                  class = 'scb')
    }

wlpsmooth <- function(y,x,w,s.x,bw,from,to,gridsize,conf.level=0.95)
{

    if(missing(s.x)){
        s.x <- 0
    }else{
        stopifnot(is.numeric(s.x))
        stopifnot(length(s.x)==1)
        stopifnot(s.x>=0)
    }
    
    n <- length(y)
    if(n<5)
        stop("too few data points to fit a regression model")

    stopifnot(length(x)==n)
    if(missing(w)){
        w <- rep(1/n,n)
    }else{
        stopifnot(length(w)==n)
    }
    sele <- is.na(x) | is.na(y) | is.na(w)
    if(any(sele)){
        if(sum(!sele) < 5)
            stop("too many missing values")
        x <- x[!sele]
        y <- y[!sele]
        w <- w[!sele]
        n <- sum(!sele)
    }
    
    if(missing(bw)){
        bw <- bw.nrd(x)
    }else{
        stopifnot(is.numeric(bw))
        stopifnot(length(bw)==1)
        stopifnot(bw>0)
    }
    
    if(missing(from)) from <- min(x)
    if(missing(to)) to <- max(x)
    stopifnot(to > from)
    if(missing(gridsize)) gridsize <- 512L
    stopifnot(gridsize > 10)
    gpoints <- seq(from, to, length=gridsize)

    out <- .Fortran(.F_wnpr,
                    rx = as.double(gpoints),
                    as.integer(gridsize),
                    as.double(x),
                    as.double(y),
                    as.double(w),
                    as.integer(n),
                    bw = as.double(bw),
                    as.double(s.x))
    list(x=gpoints,y=out$rx,bw=out$bw)    
}

.bwdmise <- function(y,sig)
  {
      hnorm2 <- function(h1,sig,Rfx,n,grid=100,ub=2){
          out <- .Fortran(.F_NormLap2,
                          as.integer(n),
                          as.double(Rfx),
                          as.double(sig),
                          h = as.double(h1),
                          as.double(grid),
                          as.double(ub))
          out$h
      }
      
      sig <- sig^2;
      s2bar <- mean(sig);
      sbar <- sqrt(s2bar);
      s2y <- var(y);
      Rfx <- 0.09375*(abs(s2y-s2bar))^(-2.5)*pi^(-.5); # R(f'')/4
      ##cat("\nRfx=", Rfx,"\n")
      n <- length(y);
      h1 <- bw.nrd(y)
      ##print(h1)
      grid <- 100;
      ub <- 2
      
      result <- hnorm2(h1,sig,Rfx,n,grid=grid,ub=ub)
      ##print(result)
      return(result);
  }
