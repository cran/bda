##  Part of R package BDA
##  Copyright (C) 2009-2010 Bin Wang
##
##  Unlimited use and distribution (see LICENCE).

## created on 2018/12/08.
## updated on 2018/12/09
## updated on 2019/11/02
## updated on 2019/11/11 ## Cleaned

fit.Pareto <- function(x, xm, method='mle'){
    method <- match.arg(tolower(method),
                        c("mle","ls"))
    if(class(x) == 'histogram')
        x <- binning(x)
    
    if(class(x) == 'bdata'){
        x <- list(counts=x$freq, breaks=c(x$ll[1],x$ul))
        if(method=='mle'){
            if(missing(xm)){
                f <- x$counts
                nclass <- length(f)
                xbrks <- x$breaks[-c(1,nclass+1)]
                delta <- (x$breaks[2]-x$breaks[1])*0.01
                res <- .Fortran(.F_mle2Pareto,
                                as.double(f),
                                llk=as.double(xbrks),
                                as.integer(nclass),
                                xm=as.double(delta),
                                alpha=as.double(0.0))
                alpha <- res$alpha
                xm <- res$xm
                out <- list(xm=xm, alpha=alpha)
            }else{
                out <- .fitPareto1(x, xm)
            }
        }else{
            out <- .fitPareto2(x)
        }
    }else{
        if(any(x <= 0)){
            out <- NULL
            warning("negative value(s) not allowed")
        }
        
        if(method == "mle"){
            xm <- min(x)
            n <- length(x)
            tmp <- sum(log(x/xm))
            alpha <- (n-2)/tmp
            out <- list(xm=xm, alpha=alpha)
        }else{
            out <- .lsePareto(x)
        }
    }
    out
}

.fitPareto1 <- function(x, lbound){
    if(lbound <= 0){
        warning("invalid value of 'lower.bound'")
        out <- NULL
    }else{
        bk <- x$breaks[-1]
        sele <- bk > lbound
        if(any(sele)){
            nclass <- sum(sele)
            if(nclass<3){
                warning("Fitting failed: too few classes")
                out <- NULL
            }else{
                ##xhist2 <- binning(counts=x$counts[which(sele)],
                ##                  breaks=c(lbound, bk[which(sele)]))
                ##out <- .funcPareto1(x=xhist2,lbound=lbound)
                f <- x$counts[which(sele)]
                nclass <- length(f)
                xbrks <- bk[which(sele)]
                xbrks <- xbrks[-nclass]
                ##x$breaks[-c(1,nclass+1)]
                res <- .Fortran(.F_mle1Pareto,
                                as.double(f),
                                as.double(xbrks),
                                as.integer(nclass),
                                llk=as.double(lbound),
                                alpha=as.double(0.0))
                alpha <- res$alpha
                out <- list(xm=lbound, alpha=alpha)
            }
        }else{
            warning("Fitting failed: 'lbound' is too large")
            out <- NULL
        }
    }
    out
}

.funcPareto1 <- function(x, lbound){
    nclass <- length(x$counts)
    xm <- lbound
    Fn <- cumsum(x$counts)/sum(x$counts)
    dFn <- abs(Fn-0.5)
    isele <- which(dFn==min(dFn))[1] #in case multiple minimums
    p <- Fn[isele]
    q <- x$breaks[isele+1]
    alpha <- log(1-p)/log(xm/q) # initial point estimate near median
    
    options(warn=-1)
    res <- .mlePareto1(x, xm, alpha)
    options(warn=0)
    
    alpha <- res$par[1]
    out <- list(xm=xm, alpha=alpha)
}

.mlePareto1 <- function(x, xm, alpha){
    f <- function(pars,x,xm){
        n <- x$counts
        q <- x$breaks
        
        k <- length(q)
        Fx <- pPareto(q[-c(1,k)], xm, pars)            
        tmp <- diff(c(0,Fx,1))
        llk <- -sum(log(tmp)*n)
        if(!is.finite(llk)||is.nan(llk))
            llk <- 1.e9+runif(1)
        llk
    }
    optim(alpha,f, x=x, xm=xm,
          lower=alpha*0.1, upper=alpha*10)
}

.fitPareto2 <- function(x){
    ## lm.est based on the Zipf plot
    nclass <- length(x$counts)
    if(nclass<3){
        warning("Fitting failed: too few classes")
        out <- NULL
    }else{
        p <- cumsum(x$counts)/sum(x$counts)
        p <- p[-nclass]
        q <- x$breaks[-c(1,nclass+1)]
        sele <- p < 1
        if(sum(sele) > 3){
            p <- p[sele]
            q <- q[sele]
            lF <- log(1-p)
            lx <- log(q)
            lmout <- lm(lF~lx)
            pidx <- -coef(lmout)[[2]]
            xm <- exp(coef(lmout)[[1]]/pidx)
            out <- list(xm=xm, alpha=pidx)
        }else{
            out <- NULL
            warning("not enough data points")
        }
    }
    out
}

.lsePareto <- function(x){
    out <- NULL
    ## LSE for raw data
    x <- x[!is.na(x)]
    xt <- table(x)
    f <- as.numeric(xt)
    nclass <- length(f)
    q <- as.numeric(names(xt))
    if(nclass<3){
        warning("Fitting failed: too few classes")
    }else{
        p <- cumsum(f)/(1+sum(f))
        lF <- log(1-p)
        lx <- log(q)
        mu.y <- mean(lF)
        mu.x <- mean(lx)
        ssxy <- sum((lF-mu.y)*(lx-mu.x))
        ssxx <- sum((lx-mu.x)^2)
        pidx <- - ssxy/ssxx
        yint <- mu.y - mu.x * pidx
        xm <- exp(yint/pidx)
        ##n <- length(lx)
        ##res <- .Fortran(.F_lsePareto,
        ##as.double(lF),
        ##as.double(lx),
        ##as.integer(n),
        ##a=as.double(0),
        ##b=as.double(0))
        ##lmout <- lm(lF~lx)
        ##pidx <- -coef(lmout)[[2]]
        ##xm <- exp(coef(lmout)[[1]]/pidx)
        ##pidx <- -res$b
        ##xm <- exp(res$a/pidx)
        out <- list(xm=xm, alpha=pidx)
    }
    out
}

dPareto <- function(x, xm, alpha){
    if(any(is.na(x)))
        stop("missing value not allowed")
    if(xm<=0) stop("invalid parameter 'xm'")
    if(alpha<=0) stop("invalid parameter 'alpha'")
    n <- length(x)
    res <- rep(0,n)
    sele <- x >= xm
    if(any(sele)){
        y <- x[sele]
        tmp <- alpha * (xm/y)^alpha /y
        res[sele] <- tmp
    }
    res
}

pPareto <- function(q, xm, alpha){
    if(any(is.na(q)))
        stop("missing value not allowed")
    if(xm<=0) stop("invalid parameter 'xm'")
    if(alpha<=0) stop("invalid parameter 'alpha'")
    n <- length(q)
    res <- rep(0,n)
    sele <- q >= xm
    if(any(sele)){
        y <- q[sele]
        tmp <- 1 - (xm/y)^alpha
        res[sele] <- tmp
    }
    res
}

qPareto <- function(p, xm, alpha){
    if(any(is.na(p)))
        stop("missing value not allowed")
    if(xm<=0) stop("invalid parameter 'xm'")
    if(alpha<=0) stop("invalid parameter 'alpha'")
    n <- length(p)
    res <- rep(NA,n)
    res[p==0] <- xm
    res[p==1] <- Inf
    sele <- p>0 & p<1
    if(any(sele)){
        y <- p[sele]
        tmp <- xm * (1-y)^(-1/alpha)
        res[sele] <- tmp
    }
    res
}

rPareto <- function(n, xm, alpha){
    if(xm<=0) stop("invalid parameter 'xm'")
    if(alpha<=0) stop("invalid parameter 'alpha'")
    u <- runif(n)
    qPareto(u,xm=xm, alpha=alpha)
}

pmixPU <- function(x, xm, alpha, a, b, p){
    stopifnot(xm > 0)
    stopifnot(alpha > 0)
    stopifnot(b>a)
    stopifnot(a>0)
    stopifnot(p>=0 && p<=1)
    if(x<=xm){
        F1 <- 0
    }else{
        F1 <- 1-(xm/x)^alpha
    }

    if(x<=a){
        F2 <- 0
    }else if(x>=b){
        F2 <- 1
    }else{
        F2 <- (x-a)/(b-a)
    }
    p*F1 + (1-p)*F2
}

qmixPU <- function(u, xm, alpha, a, b, p){
    stopifnot(xm > 0)
    stopifnot(alpha > 0)
    stopifnot(b>a)
    stopifnot(a>0)
    stopifnot(p>=0 && p<=1)
    stopifnot(u > 0 && u < 1)
    f <- function(x, u, xm, alpha, a, b, p)
        pmixPU(x, xm, alpha, a, b, p) - u
    out <- uniroot(f,c(min(a,xm),max(2*b,10000)),
                   extendInt="yes",
                   u=u,xm=xm,alpha=alpha,a=a,b=b,p=p)
    out$root
}
