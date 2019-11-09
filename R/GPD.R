
##  Part of R package BDA
##  Copyright (C) 2009-2010 Bin Wang
##
##  Unlimited use and distribution (see LICENCE).

##  2018/12/08:

dGPD <- function(x, mu,sig,eta){
    if(any(is.na(x)))
        stop("missing value not allowed")
    if(sig<=0) stop("invalid parameter 'sig'")
    z <- (x-mu)/sig
    if(eta==0){
        res <- exp(-z)
    }else{
        tmp <- eta*z+1
        sele <- tmp <= 0
        res <- tmp^(-1/eta-1)
        res[sele] <- 0
    }
    res[res<0] <- 0
    res/sig
}

pGPD <- function(q,mu,sig,eta){
    if(any(is.na(q)))
        stop("missing value not allowed")
    if(sig<=0) stop("invalid parameter 'sig'")

    z <- (q-mu)/sig
    if(eta==0){
        res <- 1-exp(-z)
    }else{
        tmp <- eta*z+1
        sele <- tmp <= 0
        res <- 1-tmp^(-1/eta)
        res[sele] <- 1
    }
    res[res<0] <- 0
    res[res>1] <- 1
    res
}

qGPD <- function(p, mu,sig,eta){
    if(any(is.na(p)))
        stop("missing value not allowed")
    if(sig<=0) stop("invalid parameter 'sig'")
    if(any(p<0 | p>1)) stop("invalid value(s) in 'p'")
    
    if(eta==0){
        res <- mu - log(1-p)*sig
    }else{
        res <- ((1-p)^(-1/eta)-1)/eta*sig + mu
    }

    res
}

rGPD <- function(n, mu,sig,eta){
    n <- round(n)
    stopifnot(n>0)
    if(sig<=0) stop("invalid parameter 'sig'")

    u <- runif(n)
    if(eta==0){
        res <- mu - log(u)*sig
    }else{
        res <- (u^(-eta) - 1)/eta*sig + mu
    }

    res
    
}

