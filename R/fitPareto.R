
##  Part of R package BDA
##  Copyright (C) 2009-2010 Bin Wang
##
##  Unlimited use and distribution (see LICENCE).

## created on 2018/12/08.
## updated on 2018/12/09
## updated on 2019/11/02
fit.Pareto <- function(x, lbound, method='mle'){
    method <- match.arg(tolower(method),
                        c("mle","percentile"))
    if(class(x) != 'bdata')
        x <- binning(x)
    
    if(method=='percentile'){
        out <- .fitPareto2(x)
    }else{
        if(missing(lbound)){
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
            ModSel <- .aic(-res$llk[1],1,sum(x$counts))
            out <- list(pars=c(xm, alpha), ModSel=ModSel)
        }else{
            out <- .fitPareto1(x, lbound)
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
                ModSel <- .aic(-res$llk,1,sum(x$counts))
                out <- list(pars=c(lbound, alpha), ModSel=ModSel)
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
    ModSel <- .aic(-res$value,2,sum(x$counts))
    out <- list(pars=c(xm, alpha), ModSel=ModSel)
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
        lF <- log(1-p)
        lx <- log(q)
        lmout <- lm(lF~lx)
        pidx <- -coef(lmout)[[2]]
        xm <- exp(coef(lmout)[[1]]/pidx)
        tmp <- pPareto(x$breaks, xm, pidx)
        llk <- sum(x$counts * log(diff(tmp)))
        ModSel <- .aic(llk,2,sum(x$counts))
        out <- list(pars=c(xm, pidx), ModSel=ModSel)
    }
    out
}




