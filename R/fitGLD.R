###  To fit a smooth curve to histogram-type data using FMKL GLD

## We use percentile method (or quantile matching method to fit FMKL
## GLD to find initial values, then find MLE numerically.

## Created on 2014/06/12 updated 2019/11/11 ## cleaned modified on
## 03/13/2020: use exact percentiles if possible. Otherwise, stop. MLE
## will also be provided based on the initial MOP estimates. This
## function will be called by fit.FSD.
## Input 'x' must be ".bdata".

fit.GLD <- function(x,qtl,qtl.levels,lbound,ubound){
    
    if(inherits(x,'numeric')){
        xh <- hist(x, plot=FALSE)
        x <- binning(xh)
    }else if(inherits(x,'histogram')){
        x <- binning(x)
    }

    if(missing(lbound)) lbound <- NULL
    if(missing(ubound)) ubound <- NULL
    
    if(!inherits(x,'bdata'))
        stop("data type not supported")

    xhist <- x
    xbrks <- x$breaks
    x <- x$freq

    if(length(x)<6)
        stop("too few classes to fit FMKL-GLD")

    a <- xbrks[1]; b <- rev(xbrks)[1]

    if(missing(qtl)||missing(qtl.levels)){
        ## find different combinations
        k <- length(xbrks)-2; 
        pmat <- NULL
        for(i1 in 1:(k-4)){
            for(i2 in (i1+1):(k-3)){
                for(i3 in (i2+1):(k-2)){
                    for(i4 in (i3+1):(k-1)){
                        for(i5 in (i4+1):k){
                            y2 <- c(i1,i2,i3,i4,i5)
                            pmat <- rbind(pmat, y2)
                        }
                    }
                }
            }
        }
        
        Fn <- cumsum(x)/sum(x)
        k <- length(Fn)
        x0 <- xbrks[-c(1,length(xbrks))];
        lx0 <- log(x0)
        
        par.mop <- NULL
        par.mle <- NULL
        Dmax <- 1e9
        ibest <- NULL
        for(i in 1:nrow(pmat)){
            pts <- pmat[i,]
            lmd.mop <- .gld.mop(Fn[pts],xbrks[pts+1],a,b,lbound,ubound)
            lmd.mle <- .gld.mle(lmd.mop,xhist,a,b,lbound,ubound)
            F2 <- .pgld(x0, lmd.mop)
            E <- diff(c(0,F2,1)) * sum(x)
            mse <- sum((x-E)^2/x)
            if(mse < Dmax){
                ibest <- i
                Dmax <- mse
                par.mle <- lmd.mle
                par.mop <- lmd.mop
            }
        }
        cat("\nQuantiles for best fit:")
        pts <- pmat[ibest,]
        tmp <- data.frame(levels=Fn[pts],quantile=xbrks[pts+1])
        print(tmp)
    }else{
        if(length(qtl) != length(qtl.levels))
            stop("quantiles and quantile levels have different lengths")
        par.mop <- .gld.mop(qtl.levels,qtl,a,b,lbound,ubound)
        par.mle <- .gld.mle(par.mop,xhist,a,b,lbound,ubound)
    }
    if(!is.finite(b)) xmax <- rev(xbrks)[2]
    else xmax <- rev(xbrks)[1]
    if(a==0) a <- 0.01
    lx0 <- seq(log(a),log(xmax),length=501)
    x0 <- exp(lx0)
    y <- .dgld(x0,par.mop)
    y2 <- .pgld(x0,par.mop)
    
    out <- list(pars=par.mop,aux=par.mle,
                x=x0,y=y,y2=y2)
}

.findQTL <- function(x)
{
    if(!inherits(x,'bdata'))
        stop("data type not supported")
    xbrks <- x$breaks
    Fn <- cumsum(x$freq)/sum(x$freq)
    M <- length(xbrks)
    a <- xbrks[1]
    b <- xbrks[M]
    if(M<7){
        stop("too few classes to fit FMKL-GLD")
    }else{
        q1 <- xbrks[2]; lvl1 <- Fn[1]
        q5 <- xbrks[M-1]; lvl5 <- Fn[M-2]
        out <- .midpoint(2,M-1,Fn,xbrks)
        q3 <- out$q;lvl3 <- out$lvl;k3 <- out$k
        out <- .midpoint(2,k3,Fn,xbrks)
        q2 <- out$q;lvl2 <- out$lvl;
        out <- .midpoint(k3,M-1,Fn,xbrks)
        q4 <- out$q;lvl4 <- out$lvl;
        qlevels <- c(lvl1,lvl2,lvl3,lvl4,lvl5)
        qtls <- c(q1,q2,q3,q4,q5)
    }
    list(qtl=qtls,qtl.levels=qlevels)
}

.midpoint <- function(n1,n2,Fn,xbrks){
    M <- length(xbrks)
    stopifnot(n1<n2&&n1>0&&n2<M)
    ## pts available
    n <- n2-n1-1
    if(n - (n%/%2)*2 ==0){# two in the middle
        k1 <- n/2 + n1;k2 <- k1+1
        if(abs(Fn[k1]-0.5)<abs(Fn[k2]-0.5)) k <- k1
        else k <- k2
    }else{
        k <- (n+1)/2 + n1
    }
    q <- xbrks[k]
    lvl <- Fn[k-1]
    list(k=k,q=q,lvl=lvl)
}


.gld.mle <- function(pars, x,a,b,lbound,ubound){
    .g0.mle <- function(lmd){
        breaks <- x$breaks
        counts <- x$freq
        Fx <- diff(.pgld(breaks, lmd))
        sele <- Fx <= 0
        Fx[sele] <- 1e-10
        -sum(counts * log(Fx))
    }
    .g1.mle <- function(lmd){
        l1 <- lbound + 1/(lmd[1]*lmd[2])
        lambdas <- c(l1, lmd)
        breaks <- x$breaks
        counts <- x$freq
        Fx <- diff(.pgld(breaks, lambdas))
        sele <- Fx <= 0
        Fx[sele] <- 1e-10
        -sum(counts * log(Fx))
    }
    .g2.mle <- function(lmd){
        l1 <- ubound - 1/(lmd[1]*lmd[3])
        lambdas <- c(l1, lmd)
        breaks <- x$breaks
        counts <- x$freq
        Fx <- diff(.pgld(breaks, lambdas))
        sele <- Fx <= 0
        Fx[sele] <- 1e-10
        -sum(counts * log(Fx))
    }
    .g3.mle <- function(lmd){
        l3 <- lmd[1]; l4 <- lmd[2]
        l2 <- (1/l3+1/l4)/(ubound-lbound)
        l1 <- ubound - 1/(l2*l4)
        lambdas <- c(l1, l2,l3,l4)
        breaks <- x$breaks
        counts <- x$freq
        Fx <- diff(.pgld(breaks, lambdas))
        sele <- Fx <= 0
        Fx[sele] <- 1e-10
        -sum(counts * log(Fx))
    }

    
    ##  constraints:
    if(is.null(lbound)){
        l3l <- -Inf; l3u <- Inf; lfinite <- FALSE
    }else if(is.finite(lbound)){
        l3l <- 1e-8; l3u <- Inf; lfinite <- TRUE
    }else{
        l3l <- -Inf; l3u <- -1e-8; lfinite <- FALSE
    }
    
    if(is.null(ubound)){
        l4l <- -Inf; l4u <- Inf; ufinite <- FALSE
    }else if(is.finite(ubound)){
        l4l <- 1e-8; l4u <- Inf; ufinite <- TRUE
    }else{
        l4l <- -Inf; l4u <- -1e-8; ufinite <- FALSE
    }
    ## find the MLEs
    if(lfinite&&ufinite){
        out <- optim(pars[-c(1,2)], .g3.mle, method="L-BFGS-B",
                     lower=c(l3l, l4l),
                     upper=c(l3u,l4u),
                     control = list(maxit=10000))
        pars <- as.numeric(out$par)
        l3 <- pars[1]; l4 <- pars[2]
        l2 <- (1/l3+1/l4)/(ubound-lbound)
        l1 <- ubound - 1/(l2*l4)
        pars <- c(l1, l2,l3,l4)
    }
    
    if(!lfinite&&ufinite){
        out <- optim(pars[-1], .g2.mle, method="L-BFGS-B",
                     lower=c(1e-10,l3l, l4l),
                     upper=c(Inf,l3u,l4u),
                     control = list(maxit=10000))
        pars <- as.numeric(out$par)
        l1 <- ubound - 1/(pars[1]*pars[3])
        pars <- c(l1, pars)
    }

    if(lfinite&&!ufinite){
        out <- optim(pars[-1], .g1.mle, method="L-BFGS-B",
                     lower=c(1e-10,l3l, l4l),
                     upper=c(Inf,l3u,l4u),
                     control = list(maxit=10000))
        pars <- as.numeric(out$par)
        l1 <- lbound + 1/(pars[1]*pars[2])
        pars <- c(l1, pars)
    }

    if(!lfinite&&!ufinite){
        out <- optim(pars, .g0.mle, method="L-BFGS-B",
                      lower=c(-Inf,1e-10,l3l, l4l),
                      upper=c(Inf,Inf,l3u,l4u),
                     control = list(maxit=10000))
        pars <- as.numeric(out$par)
    }



    if(pars[2] <= 0){
        warning("No valid parameters found")
        pars <- NULL
    }
    pars
}


## l1 can be computed based on the median, which is more stable.  So
## we can just adjust lambda2 (dispersion parameter to cover the whole
## range).

.fexp <- function(p,lmd){
    stopifnot(p>=0&&p<=1)
    if(lmd==0)
        res <- 1e10
    else{
        if(p==0){
            if(lmd>=0)
                res <- 0
            else
                res <- 1e10
        }else
            res <- (p^lmd -1)/lmd
    }
}

.gld.mop <- function(p,q,a,b,lbound,ubound){
    .g.mop <- function(lmd,q){
        rho3hat <- (q[3] - q[1])/(q[5] - q[3])
        rho4hat <- (q[4]-q[2])/(q[5]-q[1])
        Q1 <- .fexp(p[1],lmd[1]) - .fexp(1-p[1],lmd[2])
        Q2 <- .fexp(p[2],lmd[1]) - .fexp(1-p[2],lmd[2])
        Q3 <- .fexp(p[3],lmd[1]) - .fexp(1-p[3],lmd[2])
        Q4 <- .fexp(p[4],lmd[1]) - .fexp(1-p[4],lmd[2])
        Q5 <- .fexp(p[5],lmd[1]) - .fexp(1-p[5],lmd[2])

        rho3 <- ifelse(Q5-Q3==0, 1e10, (Q3-Q1)/(Q5-Q3))
        rho4 <- ifelse(Q5-Q1==0, 1e10, (Q4-Q2)/(Q5-Q1))
        out <- (rho3-rho3hat)^2+(rho4-rho4hat)^2
        if(!is.finite(out)) out <- 999
        if(is.na(out)) out <- 999
        if(is.nan(out)) out <- 999
        out
    }
    lambdas <- c(1,1)

    lb <- FALSE
    if(is.null(lbound)){
        l3l <- -Inf; l3u <- Inf; 
    }else if(is.finite(lbound)){
        l3l <- 1e-8; l3u <- Inf; lb <- TRUE
    }else{
        l3l <- -Inf; l3u <- -1e-8;
    }

    ub <- FALSE
    if(is.null(ubound)){
        l4l <- -Inf; l4u <- Inf
    }else if(is.finite(ubound)){
        l4l <- 1e-8; l4u <- Inf; ub <- TRUE
    }else{
        l4l <- -Inf; l4u <- -1e-8
    }

    pars <- optim(lambdas, .g.mop, q=q,
                  method="L-BFGS-B",
                  lower=c(l3l,l4l),upper=c(l3u,l4u),
                  control = list(maxit=10000))$par

    if(all(is.finite(pars))){
        l3 <- pars[1]; l4 <- pars[2]
        Q5 <- (p[5]^l3-1)/l3 - ((1-p[5])^l4-1)/l4
        Q3 <- (p[3]^l3-1)/l3 - ((1-p[3])^l4-1)/l4
        Q1 <- (p[1]^l3-1)/l3 - ((1-p[1])^l4-1)/l4
        l2 <- (Q5-Q1)/(q[5]-q[1])
        l1 <- q[3] - Q3/l2
        
        if(lb){
            if(ub){
                l2 <- (1/l3+1/l4)/(ubound-lbound)
                l1 <- ubound - 1/(l2*l4)
            }else{
                l1 <- lbound + 1/(l2*l3)
            }
        }else{
            if(ub)
                l1 <- ubound - 1/(l2*l4)
        }
        res <- c(l1,l2,l3,l4)    
    }else{
        res <- c(1,2,3,4)
    }
    res
}


.rgld <- function(n, lambdas){
    if(any(!is.finite(lambdas)) || lambdas[2]<=0)
        stop("Invalid FMKL-GLD parameters")
    l1 <- lambdas[1]
    l2 <- lambdas[2]
    l3 <- lambdas[3]
    l4 <- lambdas[4]
    u <- runif(n)
    l1+((u^l3-1)/l3-((1-u)^l4-1)/l4)/l2
}

.qgld <- function(p, lambdas){
    sele <- p < 0 | p > 1
    out <- rep(NA, length(p))
    if(any(!is.finite(lambdas)) || lambdas[2]<=0)
        stop("Invalid FMKL-GLD parameters")
    l1 <- lambdas[1]
    l2 <- lambdas[2]
    l3 <- lambdas[3]
    l4 <- lambdas[4]
    if(any(sele)){
        if(!all(sele)){
            u <- p[!sele]
            tmp <- l1+((u^l3-1)/l3-((1-u)^l4-1)/l4)/l2
            out[!sele] <- tmp
        }
    }else{
        u <- p
        out <- l1+((u^l3-1)/l3-((1-u)^l4-1)/l4)/l2
    }
    out
}

.fungld <- function(u,q,lambdas){
    .qgld(u,lambdas) - q
}

.gld.proot <- function(q, lambdas){
    if(is.null(lambdas))
        stop("Empty GLD parameters")
    if(is.null(q))
        stop("'q' value missing")
    if(lambdas[2]<=0){
        stop("Invalid GLD parameters")
    }
    xmax <- .qgld(1,lambdas)
    xmin <- .qgld(0,lambdas)
    if(q == Inf)
        res <- 1
    else if(q==-Inf)
        res <- 0
    else{
        if(q <= xmin)
            res <- 0
        else if(q >= xmax)
            res <- 1
        else{
            res <- .gld.root(q,lambdas=lambdas)
        }
    }
    res
}

.pgld <- function(q, lambdas){
    if(any(!is.finite(lambdas)) || lambdas[2]<=0)
        stop("Invalid FMKL-GLD parameters")
    as.numeric(lapply(q, .gld.proot, lambdas=lambdas))
}


.gld.droot <- function(q, lambdas){
    xmax <- .qgld(1,lambdas)
    xmin <- .qgld(0,lambdas)
    if(!is.finite(q)){
        res <- 0
    }else if(q < xmin)
        res <- 0
    else if(q > xmax)
        res <- 0
    else{
        u <- .gld.root(q,lambdas=lambdas)
        l1 <- lambdas[1]
        l2 <- lambdas[2]
        l3 <- lambdas[3]
        l4 <- lambdas[4]
        res <- u^(l3-1)+(1-u)^(l4-1)
        res <- l2/res
    }
    res
}

.dgld <- function(x, lambdas){
    if(any(!is.finite(lambdas)) || lambdas[2]<=0)
        stop("Invalid FMKL-GLD parameters")
    as.numeric(lapply(x, .gld.droot, lambdas=lambdas))
}

.gld.root <- function(q, lambdas){
    .Fortran(.F_rootGldFmklBisection,
             u=as.double(q), as.double(lambdas))$u
}


.lnorm.lm <- function(x){
    xc <- x$counts
    bk <- x$breaks
    stopifnot(all(bk>=0))
    k <- length(bk)
    Fn <- cumsum(xc)/sum(xc);
    bk <- bk[-1]
    sele <- which(diff(Fn)==0)
    if(any(sele)){
        Fn <- Fn[-c(sele+1)]
        bk <- bk[-c(sele+1)]
    }
    sele <- which(Fn==1)
    if(any(sele)){
        Fn <- Fn[-sele]
        bk <- bk[-sele]
    }
    z <- qnorm(Fn)
    lx <- log(bk)
    lmout <- lm(lx~z)
    coef(lmout)
}
#(lnout <- myest(xhist)$par)

.rmzero <- function(x){
    sele <- x$count==0
    xc <- x$count
    bk <- x$breaks
    xbrks <- bk
    if(any(sele)){
        sele2 <- which(sele)
        xc <- xc[-sele2]
        bk <- bk[-(sele2+1)]
        if(any(sele2==length(xbrks)-1)){
            bk[length(bk)] <- Inf
        }
    }
    binning(counts=xc, breaks=bk)
}
