
##  Part of R package BDA
##  Copyright (C) 2009-2010 Bin Wang
##
##  Unlimited use and distribution (see LICENCE).

##  2018/12/08:

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

