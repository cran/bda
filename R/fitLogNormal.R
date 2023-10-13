###  Lognormal distribution
#####################################################################
## Created on Feb 8, 2019 by Bin Wang
## Laste updated on Feb 8, 2019

fit.lognormal <- function(x, k=1, normal=FALSE){
    if(inherits(x,'numeric')){
        out <- .fitlnorm.raw(x,k=1,normal)
    }else{
        stop("data type not supported. See 'fit.FDS' for more detail")
    }
    out
}

.fitlnorm.raw <- function(x, k=1,normal=FALSE){
    ## only for k=1. Using other method to fit finite (log-)normal mixture
    if(normal) stop("using mean and SD for univariate normal")
    k <- round(k)
    if(k != 1) stop("'k<>1' not supported")
    if(any(x<0)) stop("invalid value(s) in 'x'")
    n <- length(x) 
    n0 <- sum(x==0)
    
    if(n0 > 0){
        out <- .mlnormEM0(x=x)
        npar <- 3*k
        p0 <- out$p0
        p <- out$p
        mu <- out$mu
        sig <- out$sig
        llk <- out$llk
    }else{
        out <- .mlnormEM(x=x)
        npar <- 3*k-1
        p0 <- 0
        p <- out$p
        mu <- out$mu
        sig <- out$sig
        llk <- out$llk
    }
    structure(list(p0=p0,p=p, mean=mu, sigma=sig,
                   n=n, npar=npar, llk=llk),
              class="mixlognormal")
}

.mlnormEM <- function(x){
    mu <- mean(log(x))
    sig <- sd(log(x))
    p <- 1
    llk <- .mnormllk(x,p0=0,p,mu,sig)
    list(p=p,mu=mu,sig=sig,llk=llk)
}

.mnormllk <- function(x,p0,p,mu,s){
    x1 <- x[x>0]
    p1 <- p/sum(p)
    f1 <- dmnorm(x1,p1,mu,s)*p
    res <- sum(log(f1))
    if(p0>0){
        res <- res+ log(p0)*sum(x==0)
    }
    res
}

.mlnormEM0 <- function(x){
    ## x = c(0's, x>0)
    n0 <- sum(x==0)
    x <- sort(log(x[x>0]))
    n <- length(x)
    xt <- unique(x)
    nt <- length(xt)
    ft <- as.numeric(table(x))
    r2a <- 0
    for(i in 1:n0){
        ## fit SLR to find rough estimates
        Fn <- cumsum(c(i,ft))/(n+i)
        Zn <- qnorm(Fn[-(nt+1)])
        lm0 <- lm(xt~Zn)
        mu <- lm0$coef[[1]]
        sig <- lm0$coef[[2]]
        p0 <- (n0-i)/(n+n0)
        p <- 1-p0
        llk <- .mnormllk(x,p0,p,mu,sig)
        r2b <- summary(lm0)$r.square
        if(r2b > r2a){
            r2a <- r2b
            lmout <- lm0
        }
    }
    mu <- lmout$coef[[1]]
    sig <- lmout$coef[[2]]

    n1 <- pnorm(xt[1], mu,sig)*n
    if(n1>n0){
        p0 <- 0
        p <- 1
    }else{
        p0 <- (n0 - n1)/(n+n0)
        p <- 1-p0
    }
    llk <- 1000
    list(p0=p0,p=p,mu=mu,sig=sig,llk=llk)
}

print.mixlognormal <- function(x,...){
    tmp <- data.frame(Prop=signif(c(x$p0,x$p),3),
                      Mean=signif(c(0,x$mean),3),
                      SD=signif(c(NA,x$sigma),3))
    print(tmp)
    cat("\n")
}

.fitLogNormal <- function(x, x.limit, k=1){
    k <- round(k)
    if(k<0) stop("invalid 'k' value")
    
    if(k==1){
        out <- .fitnormk1(x=x, x.limit=x.limit)
#    }else if(k<=3){
#        out <- .fitnormk2(x=x,x.limit=x.limit, k=k)        
    }else
        stop("'k>1' not supported in this version")
    out
}

.fitnormk1 <- function(x, x.limit){
    mu <- NA
    sig <- NA
    
    ## initial estimate using LS method
    f <- x
    x0 <- x.limit
    k <- length(x0)
    if(any(x0 < 0))
        stop("negative class limit(s) not allowed")
    if(any(f < 0))
        stop("negative counts not allowed")
    if(sum(f) == 0)
        stop("no observation found")
    ## zero frequencies are allowed. But we need to handle the
    ## classes with zero counts specially: keep more details.
    sele <- f == 0
    if(any(sele)){
        f <- f[!sele]
        x0 <- x0[-which(sele)+1]
    }
    Fn <- cumsum(f)/sum(f)
    qi <- qnorm(Fn[-length(Fn)])
    mi <- log(x0[-c(1,length(x0))])
    xbar <- mean(qi)
    ybar <- mean(mi)
    ssxy <- sum((qi-xbar)*(mi-ybar))
    ssxx <- sum((qi-xbar)^2)
    sig <- ssxy/ssxx
    mu <- ybar - xbar * sig
    ## find MLE
    w <- x
    nclass <- length(w)
    a <- x0[-(nclass+1)]
    b <- x0[-1]
    if(!is.finite(b[nclass]))
        b[nclass] <- 999.999
    sele <- w == 0
    if(any(sele)){
        a <- a[!sele]
        b <- b[!sele]
        w <- w[!sele]
        nclass <- sum(!sele)
    }
    res <- .Fortran(.F_lnormBinMLE3,
                    as.double(a),
                    as.double(b),
                    as.double(w),
                    as.integer(nclass),
                    mu=as.double(mu),
                    s=as.double(sig))
    mu <- res$mu
    sig <- res$s
    
    structure(list(p0=0,p=1, mean=mu, sigma=sig,
                   n=sum(w), npar=2, llk=NA),
              class="mixlognormal")
}

.fitnormk2 <- function(x, mle=TRUE,k=2, x.limit){
    out <- .fitnormk2mle(x=x,k=k,x.limit=x.limit)
    ##if(mle){
    ##out <- .fitnormk2mle(x=x,k=k,x.limit=x.limit)
    ##}else{
    ##out <- .fitnormk2lse(x=x,k=k)
    ##}
    out
}

.fitnormk2mle <- function(x, k=2, x.limit){
    mu <- NA
    sig <- NA
    ps <- NA
    ncomp <- round(k)
    
    if(inherits(x,"numeric")){
        if(any(is.na(x))){
            x <- x[!is.na(x)]
            warning("missing value(s) removed")
        }
        if(any(!is.finite(x)))
            stop("'x' values must be finite")
        if(any(x < 0))
            stop("'x' cannot be negative")

        xmin <- min(x[x>0])
        if(missing(x.limit)){
            xlmt <- xmin
        }else{
            stopifnot(x.limit > 0)
            stopifnot(is.numeric(x.limit))
            xlmt <- min(x.limit, xmin)
        }

        ## least square estimates by default; or to estimate the
        ## initial values for the maximum likelihood estimates.
        xt <- table(x)
        Fn <- cumsum(xt)/sum(xt)
        mi <- as.numeric(names(xt))
        k <- length(mi)
        if(mi[1]==0){
            mi <- log(mi[-c(1,k)])
            qi <- qnorm(Fn[-c(1,k)])
        }else{
            mi <- log(mi[-k])
            qi <- qnorm(Fn[-k])
        }
        xbar <- mean(qi)
        ybar <- mean(mi)
        ssxy <- sum((qi-xbar)*(mi-ybar))
        ssxx <- sum((qi-xbar)^2)
        sig <- ssxy/ssxx
        mu <- ybar - xbar * sig
        
        n0 <- sum(x==0)
        if(any(x == 0)){
            x1 <- x[x>0]
        }else{
            x1 <- x
        }
        xtmp <- table(x1)
        xt <- as.numeric(names(xtmp))
        xn <- as.numeric(xtmp) # can be weights
        n <- length(xn)
        ps <- rep(1/ncomp, ncomp)
        mus <- seq(.8,1.2,length=ncomp)*mu
        sigs <- rep(sig, ncomp)

        res <- .Fortran(.F_lnormMixK, 
                        as.double(log(xt)),
                        as.double(xn),
                        as.integer(c(n,ncomp)),
                        as.double(c(n0,log(xlmt))),
                        p=as.double(ps),
                        mu=as.double(mus),
                        sig=as.double(sigs))
        ## update the estimates using MLE
        ps <- res$p
        mu <- res$mu
        sig <- res$sig
    }else
        stop("data type not supported")
    list(p=ps, mu=mu, sigma=sig)
}


.fitnormk2lse <- function(x, k=2){
    mu <- NA
    sig <- NA
    ps <- NA
    ncomp <- round(k)
    
    if(inherits(x,"numeric")){
        if(any(is.na(x))){
            x <- x[!is.na(x)]
            warning("missing value(s) removed")
        }
        if(any(!is.finite(x)))
            stop("'x' values must be finite")
        if(any(x < 0))
            stop("'x' cannot be negative")

        ## least square estimates by default; or to estimate the
        ## initial values for the maximum likelihood estimates.
        xt <- table(x)
        Fn <- cumsum(xt)/sum(xt)
        mi <- as.numeric(names(xt))
        x0 <- mi #lognormal
        y0 <- xt
        k <- length(mi)
        if(mi[1]==0){
            mi <- log(mi[-c(1,k)])
            qi <- qnorm(Fn[-c(1,k)])
            x0 <- x0[-1]
            y0[2] <- y0[2] + y0[1]
            y0 <- y0[-1]
        }else{
            mi <- log(mi[-k])
            qi <- qnorm(Fn[-k])
        }
    }else if(inherits(x,"histogram")){
        xc <- x$counts
        y0 <- xc
        x0 <- x$breaks[-1]
        k <- length(xc)
        Fn <- cumsum(xc)/sum(xc)
        qi <- qnorm(Fn[-k])
        mi <- log(x$breaks[-c(1,k+1)])
    }else
        stop("data type not supported")

    xbar <- mean(qi)
    ybar <- mean(mi)
    ssxy <- sum((qi-xbar)*(mi-ybar))
    ssxx <- sum((qi-xbar)^2)
    sig <- ssxy/ssxx
    mu <- ybar - xbar * sig

    ps <- rep(1/ncomp, ncomp)
    mus <- seq(.8,1.2,length=ncomp)*mu
    sigs <- rep(sig, ncomp)
    n <- length(x0)

    stopifnot(n>=3)
    
    res <- .Fortran(.F_lnormLSEK, 
                    as.double(log(x0)),
                    as.double(y0),
                    as.integer(c(n,ncomp)),
                    p=as.double(ps),
                    mu=as.double(mus),
                    sig=as.double(sigs))
        ## update the estimates using MLE
    ps <- res$p
    mu <- res$mu
    sig <- res$sig
    
    list(p=ps, mu=mu, sigma=sig)
}

