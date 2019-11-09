##  Part of R package BDA
##  Copyright (C) 2009-2010 Bin Wang
##
##  Unlimited use and distribution (see LICENCE).

##  2018/12/08:

fit.GPD <- function(x){
    if(class(x) != 'bdata')
        x <- binning(x)

    f2 <- function(pars,x,mu){
        n <- x$counts
        q <- x$breaks
        tmp <- diff(pGPD(q, mu, pars[1], pars[2]))
        if(any(is.na(tmp))){
            llk <- 999999999999.+runif(1)
        }else if(any(tmp==0)){
            llk <- 999999999999.+runif(1)
        }else{
            llk <- -sum(log(tmp)*n)
        }
        llk
    }

    mu <- 0
    nclass <- length(x$counts)
   
    options(warn=-1)
    icount <- 0
    s <- 1;
    eta <- 1
    while(icount < 10){
        out <- optim(c(s,eta),f2, x=x, mu=mu,
                     lower=c(0.1*s, -2*abs(eta)),
                     upper=c(10*s,2*abs(eta)))
        s <- out$par[1]
        eta <- out$par[2]
        if(out$converge==1){
            break;
        }else{
            icount <- icount + 1
        }
    }
    options(warn=0)

    if(icount > 10)
        cat("\nwarning:", out$message,"\n")

    ModSel <- .aic(-out$value,2,sum(x$counts))
   
    list(pars=c(mu, s, eta), ModSel=ModSel)
}


.fit.GPD <- function(x){
    if(class(x) != 'bdata')
        x <- binning(x)
    cat("\nMethod of maximum likelihood:\n")

    f <- function(pars,x,mu){
        n <- x$counts
        q <- x$breaks
        tmp <- diff(pGPD(q, mu, pars[1], pars[2]))
        if(any(is.na(tmp))){
            llk <- 999999999999.+runif(1)
        }else if(any(tmp==0)){
            llk <- 999999999999.+runif(1)
        }else{
            llk <- -sum(log(tmp)*n)
        }
        llk
    }

#    options(warn=-1)
    mu <- 0
    nclass <- length(x$counts)
    s <- x$breaks[nclass]*.25
    eta <- 1.0
    out <- optim(c(s,eta),f, x=x, mu=mu,
                 lower=c(0.0001, 0.0001),
                 upper=c(10*s,Inf))
#    options(warn=0)
    
    list(mu=mu, sig=out$par[1], eta=out$par[2])
}

