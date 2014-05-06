## 2014-04-10: add an option that the number of components are unknown

## x of size nr+1, sorted, class limits; y: frequencies/counts, of
## size nr.  mu: initial mean values of the components.

## k is the number of parameters; L is the likelihood; n is the sample
## size

## AIC = 2*k - 2*ln(L)
## BIC = k*[ln(n)+ln(2*pi)] - 2*ln(L)
## BIC = k*ln(n) - 2*ln(L) => for large n
## AICc = AIC + 2*k*(k+1)/(n-k-1)

## 2014/04/19: if 'k' is not specified, we try k=1,2,...,10 and choose
## the best fit using AICc.

bnmm <- function(breaks,freq,mu,s,p,k, trunc,lognormal=FALSE,
                 from,to,gridsize=512L)
    UseMethod("bnmm")

bnmm.default <- function(breaks,freq,mu,s,p,k,trunc,lognormal=FALSE,
                         from,to,gridsize=512L){
    f.call <- match.call()
    if(missing(mu)){# && missing(k)){
        nr <- length(freq)
        if(missing(k)){
            ngroup <- max(1, round(0.5 * nr))
            ks <- c(1:ngroup)
        }else{
            stopifnot(is.numeric(k))
            k <- round(k)
            stopifnot(!(any(k<1 | k>nr)))
            ngroup <- length(k)
            ks <- sort(k)
        }
        res <- .cemmix(breaks=breaks,freq=freq,mu=mu,s=s,p=p,
                       k=ks[1],trunc=trunc,lognormal=lognormal)
        aic0 <- res$AICc
        AIC <- res$AIC;
        BIC <- res$BIC;
        AICc <- res$AICc

        if(ngroup>1){
            for(j in 2:ngroup){
                res1 <- .cemmix(breaks=breaks,freq=freq,mu=mu,s=s,p=p,
                                k=ks[j],trunc=trunc,lognormal=lognormal)
                if(res1$ifault==0){
                    if(aic0<0 || aic0 > res1$AICc){
                        res <- res1
                        aic0 <- res1$AICc
                    }
                    AIC <- c(AIC, res1$AIC)
                    BIC <- c(BIC, res1$BIC)
                    AICc <- c(AICc, res1$AICc)
                }else{
                    AIC <- c(AIC, NA)
                    BIC <- c(BIC, NA)
                    AICc <- c(AICc, NA)
                }
            }
            model.sel <- data.frame(g=ks,AIC=AIC,BIC=BIC, AICc=AICc)
            res$model.sel = model.sel
        }
    }else{
        res <- .cemmix(breaks=breaks,freq=freq,mu=mu,s=s,p=p,k=k,
                       trunc=trunc, lognormal=lognormal)
    }

    ## define fine grid point
    stopifnot(gridsize > 10)
    nx <- length(breaks)
    if(missing(from)){
        if(is.finite(breaks[1]))
            from <- breaks[1]
        else
            from <- 2 * breaks[2] - breaks[3] 
    }
    if(missing(to)){
        if(is.finite(breaks[nx]))
            to <- max(breaks)
        else
            to <- 2* breaks[nx-1] - breaks[nx-2]
    }
    x0 <- seq(from, to, length=gridsize)

    if(lognormal){
        lx0 <- log(x0 - min(x0)) 
        if(length(res$p)==1)
            ly0 <- dnorm(lx0, mean=res$mu, sd=res$s)
        else
            ly0 <- dmixnorm(lx0,p=res$p, mean=res$mu, sd=res$s)
        ly0[1] <- 0;
        tot.mass <- sum(ly0)*(to-from)/gridsize
        y0 <- ly0/tot.mass
    }else{
        if(length(res$p)==1)
            y0 <- dnorm(x0, mean=res$mu, sd=res$s)
        else
            y0 <- dmixnorm(x0,p=res$p, mean=res$mu, sd=res$s)
    }
    res$x <- x0;
    res$y <- y0;
    res$mean <- sum(x0*y0)*diff(x0)[1]
    Fx <- cumsum(y0)*diff(x0)[1]
    res$median <- approx(Fx, x0, 0.5)$y
    res$call <- f.call
    res
}

bnmm.histogram <- function(breaks,freq,mu,s,p,k,
                           trunc, lognormal=FALSE,
                           from,to,gridsize=512L){
    f.call <- match.call()
   if(missing(trunc)){
        trunc <- "none"
        if(breaks$top.coded == -1) trunc <- "left"
        if(breaks$top.coded == 1) trunc <- "right"
        if(breaks$top.coded == 2) trunc <- "both"
    }
    
    res <- bnmm(breaks=breaks$breaks,freq=breaks$counts,
                mu=mu,s=s,p=p,k=k,trunc=trunc,
                lognormal=lognormal,
                from=from,to=to,gridsize=gridsize)
    res$call <- f.call
    res
}


.cemmix <- function(breaks,freq,mu,s,p,k,trunc,lognormal=FALSE){

    stopifnot(is.numeric(breaks))
    stopifnot(is.numeric(freq))
    stopifnot(!any(is.na(breaks)))
    stopifnot(!any(is.na(freq)))

    freq <- round(freq) # in case not integer(s)
    if(any(freq < 0)) stop("Counts cannot be negative")

    ## breaks should be monotone increasing
    stopifnot(any(diff(breaks) > 0))

    nr <- length(freq)
    stopifnot(nr>3)
    stopifnot(length(breaks) == nr +1)

    ## McLachlan & Jones (1988): in the case of not truncated
    ## data, the equations require two extra classes (-\infty, x0)
    ## and (xr,+\infty) with corresponding frequencies n_{l} and
    ## n_{u}.
    if(missing(trunc)) trunc <- "none"
    ctrunc <- match.arg(tolower(trunc),
                       c("left", "right","both","none"))
    itrunc <- 0; trunc <- FALSE # output
    nl <- 0; nu <- 0;
    if(!is.finite(breaks[nr+1]) || ctrunc=="right"||ctrunc=="both"){
        itrunc <- 1;
        trunc <- TRUE;
        nu <- freq[nr];
        breaks <- breaks[-(nr+1)]; freq <- freq[-nr];
    }
    if(!is.finite(breaks[1]) || ctrunc=="left"||ctrunc=="both"){
        itrunc <- 1; nl <- freq[1];
        trunc <- TRUE;
        breaks <- breaks[-1]; freq <- freq[-1];
    }

##    cat("Breaks\n")
##    print(breaks)
##    cat("Frequencies\n")
##    print(freq)
##    cat("(nl,nu)=(", nl, ",",nu,")")
    
    ## renew the number of groups (bins)
    nr <- length(freq)

    ## whether to log-transform data
    c.glog <- 0;  # a constant for glog-transformation
    if(lognormal){
        c.glog <- ifelse(min(breaks)<0.25, 0.25-min(breaks), min(breaks))
        breaks <- log(breaks + c.glog)
    }
    
    ## Apply EM for model fitting
    x0 <- breaks[1]; x1 <- breaks[-1]
    xc <- (breaks[-1]+breaks[-(nr+1)])*0.5

    if(missing(k)){
        stopifnot(is.numeric(mu))
        ng <- length(mu)
        stopifnot(ng <= nr-2)
    }else{
        ng <- round(k)
        stopifnot(ng>0)
        if(missing(mu)){
            if(ng==1) mu <- mean(rep(xc, freq))
            else mu <- runif(ng, min(breaks), max(breaks));
        }
    }

    if(missing(p)) p <- rep(1/ng,ng)  # initial proportions
    else{
        stopifnot(!any(p<=0))
        p <- p/sum(p)
    }
    
    if(missing(s)) v <- rep(var(rep(xc, freq)), ng)
    else v <- s^2

    stopifnot(length(p) == ng)
    stopifnot(length(v) == ng)

    ## input as WK (workspace size) and output as error info
    wk <- max(500, 3*ng*(nr+2+3*ng)+(nr+2))

    out <- .Fortran(.F_emmix,
                    as.integer(freq), as.integer(nr), as.integer(ng),
                    as.double(x0), as.double(x1),
                    p=as.double(p), mu=as.double(mu), v=as.double(v),
                    xlogl=double(1), ifault=as.integer(wk),
                    iter = as.integer(itrunc), 
                    nl=as.integer(nl), nu=as.integer(nu))
    llk0 <- out$xlogl
    K <- ifelse(ng == 1, 2., 3. * ng - 1.);
    N <- sum(freq);
    xspan <- diff(range(breaks)) 
    AIC = -2.0 * llk0 + 2.0 * K;
    AICc = AIC + 2.0 * K * (K+1.) / (N - K - 1.0);
    BIC = -2.0 * llk0 + log(N*2.0*pi) * K;

    model.sel <- data.frame(g=ng, AIC=AIC,BIC=BIC, AICc=AICc)
    structure(list(x=breaks, y=freq,#faked.
                   ng=ng,p=out$p, mu=out$mu-c.glog, s=sqrt(out$v),
                   llk = out$xlogl, nl = out$nl, nu = out$nu,
                   iter = out$iter, ifault = out$ifault,
                   "AIC"=AIC, "BIC"=BIC,"AICc"=AICc,
                   lognormal = lognormal, c.glog=c.glog,
                   mean=NULL, median=NULL,
                   call = match.call(),
                   trunc = trunc, model.sel=model.sel),
              class='nmix')
}

print.nmix <- function(x,...)
    {
        cat("\nCall:\n\t", deparse(x$call), "\n")
        tmp.names <- paste("Comp.",c(1:length(x$mu)),sep='')
        tmp <- cbind(x$p, x$mu, x$s)
        row.names(tmp) <- tmp.names
        tmp <- t(tmp)
        row.names(tmp) <- c("Proportion","Mean", "Std.Dev")
        print(x$model.sel)
        cat("\nMean=",x$mean,"\tMedian=",x$median,"\n",sep='')
        ##        cat("\nAIC :=", x$AIC, "\nBIC :=", x$BIC, "\nAICc :=", x$AICc,"\n")
        cat("\nParameter estimates:\n")
        print(tmp, ...)
        if(!is.null(x$trunc)){
            if(!x$trunc){
                cat("\nEstimates of the lower and upper class frequencies:\n")
                cat("Lower class: nl=", x$nl, "\tUpper class: nu=", x$nu, "\n")
            }
        }
    }


