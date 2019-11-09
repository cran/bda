## two versions of lognormal distributions fittings.  One for NGS data
## and one for FSD.

## FSD main function: fit.lnmix
## Others: fit.mlnorm

fit.mlnorm <- function(x, method="mixed", k=1, gridsize=1000,optim=TRUE){
    method <- match.arg(tolower(method),
                        c("mixed","lognormal","mop","mom",
                          "lnorm","fnmm","percentile"))
    options(warn=-1)
    ## number of components
    k <- round(k)
    stopifnot(k>0)

    ## gridsize
    gridsize <- round(gridsize)
    stopifnot(gridsize > 10)
    xgrid <- seq(min(x)+0.0000001,max(x)+1, length=gridsize)

    xnam <- deparse(substitute(x))
    x <- as.numeric(x)
    xmax <- max(x)
    nx <- length(x)
    nx0 <- sum(x==0)
    p0 <- nx0/nx
    if(any(x<0)) stop("negative 'x' value(s) not allowed")


    ## create a histogram to visualize the data.  If the data is
    ## highly discretized, use single valued classes or treat them as
    ## truncated values.
    tmp <- table(x)
    xt <- as.numeric(names(tmp))
    fn <- as.numeric(tmp)
    if(length(xt) <= 10){
        xbrks <- xt
        xbrks[-1] <- xbrks[-1] - 0.000001
        xbrks <- c(xbrks, max(xbrks)+1); 
    }else{
        ## if there are a lot of distinct values we create a
        ## histogram by dealing the zero separately.
        xp <- log(x[x>0])
        xh <- hist(xp, plot=FALSE)
        xbrks <- exp(xh$breaks)
        ## check for outliers
        ##xbrks <- xh$breaks
        ##cutoff <- mean(xp)+5*sd(xp)
        ##isele <- xbrks > cutoff
        ##if(any(isele)){
        ##    xbrks <- c(xbrks[-which(isele)], max(xbrks))
        ##    xgrid <- seq(min(x)+0.0000001,exp(cutoff), length=gridsize)
        ##}
        ##xbrks <- exp(xbrks)
        ##xh <- hist(xp, breaks=xbrks, plot=FALSE)
        sele <- xh$counts == 0
        if(any(sele)){
            isele <- which(sele) + 1
            xbrks <- xbrks[-isele]
        }
        if(p0 > 0){
            xbrks <- c(0, xbrks)
        }
    }
    
    if(length(xt) <= 2){
        if(method != "mop" && method != "percentile"){
            method <- "mop"
            warning("sample too small. switched to 'method of percentile'")
        }
    }
    
    if(method=="mixed"){
        meth <- "Mixed lognormal mixture"
        out <- .em1mmmlnorm(x=x,k=k)
        if(k==1){
            pmix <- 1
            mu <- out$mu[2];
            s <- out$s[2];
            K <- 3
        }else{
            pmix <- out$p
            mu <- out$mu;
            s <- out$s;
            K <- 3*k
        }
    }else if(method=="lognormal" || method=="lnorm"){
        K <- 2
        if(optim){
            out <- .fitlnormOPT(x=x)
        }else{
            out <- .fitlnorm(x=x)
        }
        pmix <- out$p
        mu <- out$mu;
        s <- out$s;
        meth <- out$meth
    }else if(method=="fnmm"){
        if(k==1){
            K <- 2
            out <- .fitlnorm(x=x)
            pmix <- out$p
            mu <- out$mu;
            s <- out$s;
            meth <- out$meth
        }else{
            meth <- "Lognormal mixture"
            if(optim){
                out <- .em4mlnorm(x=x,k=k)
            }else{
                out <- .em5mlnorm(x=x,k=k)
            }
            pmix <- out$p
            mu <- out$mu;
            s <- out$s;
            K <- out$npar
        }
    }else if(method=="mop"||method=="percentile"){
        meth <- "lognormal (mop)"
        out <- .moplnorm(x)
        pmix <- 1.0;
        mu <- out$mu;
        s <- out$s;
        K <- 2
    }else{
        meth <- "lognormal (mop)"
        out <- .moplnorm(x)
        pmix <- 1.0;
        mu <- out$mu;
        s <- out$s;
        K <- 2
    }

    
    ## summary statistics ############################
    ## data mean and sd
    MU <- sum(pmix * exp(mu+0.5*s^2))
    SD <- sqrt(sum(pmix^2 * (exp(s^2)-1)*exp(2*mu+s^2)))

    ## initialize the output variables
    ygrid <- dmlnorm(xgrid, pmix, mu, s)

    ## compute log-likelihood
    llk <- .llkmlnorm(xt,fn,pmix,mu,s)
    options(warn=0)

    structure(list(
        x.name = xnam, n = nx, n0 = nx0,
        data = x, method = meth,
        xhist = NULL,
        breaks = xbrks,
        x = xgrid, y = ygrid, 
        p = pmix,
        meanlog = mu,
        sdlog = s,
        npara = K, llk = llk,
        Mean = MU, SD = SD,
        xbar = mean(x, na.rm=TRUE), ### additional
        stdev = sd(x, na.rm=TRUE) ### addintional
    ),
    class='NGS.lognormal')
}



.llklnorm <- function(x,f,mu,s){
    if(length(x) == 1){
        Px <- log(plnorm(x,mu,s))*f
    }else{
        if(x[1]==0){
            x <- x[-1]
            Fx <- plnorm(x,mu,s)
            Px <- sum(log(diff(c(0,Fx,1)))*f)
        }else{
            Fx <- plnorm(x,mu,s)
            Px <- sum(log(diff(c(Fx,1)))*f)
        }
    }
    Px
}

.llkmlnorm <- function(x,f,p,mu,s){
    if(length(x) == 1){
        Px <- log(pmlnorm(x+1,p,mu,s))*f
    }else{
        if(x[1]==0){
            x <- x[-1]
            Fx <- pmlnorm(x,p,mu,s)
            Px <- sum(log(diff(c(0,Fx,1)))*f)
        }else{
            Fx <- pmlnorm(x,p,mu,s)
            Px <- sum(log(diff(c(Fx,1)))*f)
        }
    }
    Px
}

plot.NGS.lognormal <- function(x,main=NULL,xlab=NULL,xlim=NULL,...){
    xhist <- x$xhist
    if(is.null(main)) main <- x$x.names
    if(is.null(xlab)) xlab <- 'expression'
    if(is.null(xlim)){
        xlim <- c(0, max(x$data)+1)
        x0 <- x$x; y0 <- x$y
    }else{
        x0 <- seq(xlim[1], xlim[2], length=400)
        y0 <- dmlnorm(x0, x$p, x$meanlog, x$sdlog)
    }
    
    if(is.null(xhist)){
        if(is.null(x$breaks))
            hist(x$data, prob=TRUE,main=main,
                 xlab=xlab,xlim=xlim,...)
        else
            hist(x$data, breaks=x$breaks, prob=TRUE,
                 main=main,xlab=xlab,xlim=xlim,...)
        lines(x0, y0, ...)
    }else{
        xbrks <- log(xhist$breaks)
        xbrks[1] <- 0
        k <- length(xbrks)
        if(!is.finite(xbrks[k]))
            xbrks[k] <- 2*xbrks[k-1] - xbrks[k-2]
        xhist$breaks <- xbrks
        plot(xhist,main=main,xlab=xlab,xlim=xlim,...)
        lines(log(x$x), x$y,...)
    }
}

lines.NGS.lognormal <- function(x,...){
    xhist <- x$xhist
    if(is.null(xhist)){
        lines(x$x, x$y, ...)
    }else{
        lines(log(x$x), x$y,...)
    }
}

print.NGS.lognormal <- function(x,...){
    ##cat("\n", x$method, " for '",x$x.name,"'",sep='')
    cat("\n  n.obs = ", x$n, ", n.zero = ", x$n0,
        " (",round(100*x$n0/x$n,2),"%)",sep='')        
    cat("\n  Sample mean (SD): ", round(x$xbar,3), "(",
        round(x$stdev,3), "). \n  Fitted mean (SD): ",
        round(x$Mean,3), "(", round(x$SD,3),")")

    cat("\nParameter estimates:\n")
    tmp <- data.frame("P(%)"=100*x$p, meanlog=x$meanlog, sdlog=x$sdlog)

    tmp.mu <- NULL
    tmp.s <- NULL
    for(i in 1:nrow(tmp)){
        mu <- tmp[i,2]
        s <- tmp[i,3]
        if(is.finite(mu)){
            MU <- exp(mu+0.5*s^2)
            SD <- sqrt((exp(s^2)-1)*exp(2*mu+s^2))
        }else{
            MU <- 0
            SD <- 0
        }
        tmp.mu <- c(tmp.mu, MU)
        tmp.s <- c(tmp.s, SD)
    }
    tmp <- data.frame(tmp, Mean=tmp.mu, SD=tmp.s)
    print(round(tmp,3))
}


.lnormllk <- function(parms,x,x.limits)
{    
    if(!is.finite(parms[1])) parms[1] <- 1
    if(!is.finite(parms[2])) parms[2] <- 1
    if(is.na(parms[1])) parms[1] <- 1
    if(is.na(parms[2])) parms[2] <- 1
    dF <- plnorm(x+x.limits, meanlog=parms[1],
                 sdlog=parms[2], log.p = FALSE) -
        plnorm(x, meanlog=parms[1],
               sdlog=parms[2], log.p = FALSE)
    sele1 <- dF > 0
    sele2 <- is.finite(dF)
    logl <- sum(log(dF[sele1&sele2]))
    dF <- dlnorm(x,meanlog=parms[1],
                 sdlog=parms[2], log=TRUE)
    sele2 <- is.finite(dF)
    logl <- logl + sum(dF[sele2])
    -logl                         # return negative log likelihood
}

    
.lnorm2 <- function(parms, x) {
    sele <- x == 0
    if(sum(sele) == length(x)){
        res <- NA
    }else{
        n1 <- sum(sele) # model zero-measures separately
        lx <- log(x[!sele])
        llk1 <- n1*dnorm(min(lx),parms[1], parms[2])
        ## check whether there are extreme outliers
        q3 <- as.numeric(quantile(lx, prob=0.75))
        sele <- lx > q3 + 2.0 * IQR(lx)
        if(any(sele)){
            llk3 <- sum(sele)*
                pnorm(max(lx[!sele]),parms[1], parms[2],
                      lower.tail=FALSE,log.p=TRUE)
            llk2 <- sum(dnorm(lx[!sele],parms[1],parms[2],log=TRUE))
            llk2 <- llk2 + llk3
        }else{
            llk2 <- sum(dnorm(lx,parms[1],parms[2],log=TRUE))
        }
        
        res <- -llk1-llk2
    }
    res
}

.lnorm <- function(parms, x) {
    if(any(x<0)) stop("negative value found in 'x's")
    if(!is.finite(parms[1])) parms[1] <- runif(1)
    if(!is.finite(parms[2])) parms[2] <- runif(1,1,2)
    if(is.na(parms[1])) parms[1] <- runif(1)
    if(is.na(parms[2])) parms[2] <- runif(1,1,2)
    if(parms[2]<=0) parms[2] <- runif(1,1,2)
    
    xt <- table(x)
    x0 <- as.numeric(names(xt))
    n0 <- as.numeric(xt)
    delta0 <- min(diff(x0))[1]; 

    llk <- sum(n0 * log(pnorm(log(x0+delta0),parms[1],parms[2]) -
                        pnorm(log(x0),parms[1],parms[2])))    
}


.gofbin <- function(x,K=5){
    ## group the data into four classes with reasonably large
    ## frequencies.
    ## K <- round(K)
    ## stopifnot(K>3)
    ## don't use 'K' as a parameter.  We fix K=5 if the number of
    ## classes are large enough.
    K <- 5
    
    nx <- length(x)
    x.table <- table(x)
    x.low <- as.numeric(names(x.table))
    fn <- as.numeric(x.table)
    ## initialize the groups for GOF test
    ng <- length(fn)
    xbreaks <- c(x.low, max(x.low)+1)-0.00001
    Fb <- fn

    if(all(fn >= K)){
        Large <- TRUE
    }else{
        Large <- FALSE
    }
    while(Large == FALSE && ng > 3){
        fmin <- min(Fb)
        if(fmin >= K){
            break
        }else{
            sele <- which(Fb == fmin)#[1] # will never fail
            ## if there are more than one classes with the smallest
            ## frequency, we start with the rightmost one and merge to
            ## the left.
            gid <- rev(sele)[1]
            Fb[gid-1] <- Fb[gid-1] + Fb[gid]
            Fb <- Fb[-gid]
            xbreaks <- xbreaks[-gid]
            ng <- ng - 1
        }
    }
    xbreaks <- xbreaks[-c(1,ng+1)]

    list(x=x.low,f=fn,breaks=xbreaks,counts=Fb,nclass=ng)
}

.gofbin1 <- function(x,K=5){
    ## group the data into four classes with reasonably large
    ## frequencies.
    K <- round(K)
    stopifnot(K>3)

    nx <- length(x)
    x.table <- table(x)
    x.low <- as.numeric(names(x.table))
    fn <- as.numeric(x.table)
    Fn <- cumsum(fn)
    ng <- length(fn)
    n1 <- ceiling(nx/K)
    xb <- 0
    xbreaks <- NULL
    gid1 <- 0
    for(i in 1:K){
        n2 <- min(i * n1, nx)
        gid <- which(Fn >= n2)[1] # will never fail
        if(gid!=gid1){
            if(gid==ng){
                xb <- c(xb, nx)
                break;
            }else{
                xbreaks <- c(xbreaks, x.low[gid+1])
                xb <- c(xb, Fn[gid])
            }
            gid1 <- gid
        }
    }
    Fb <- diff(xb)
    ng <- length(Fb)
    list(x=x.low,f=fn,breaks=xbreaks,counts=Fb,nclass=ng)
}

.gofbin2 <- function(x,K=5){
    ## Bin the data such that all classes have frequency 5+
    K <- round(K)
    stopifnot(K>4)
    
    x.table <- table(x)
    x.low <- as.numeric(names(x.table))
    fn <- as.numeric(x.table)
    nclass <- length(fn)
    if(nclass > 30){
        xhist <- hist(x, nclass=10,plot=FALSE)
        xbreaks <- xhist$breaks
        fn <- xhist$count
    }else{
        xbreaks <- c(x.low, max(x)+1)
    }
    
    for(i in 1:100){
        nclass <- length(fn)
        xmin <- min(fn)
        if(xmin < K){
            gid <- which(fn == xmin)[1]
            if(gid==1){
                fn[2] <- fn[2] + fn[1]
                fn <- fn[-1]
                xbreaks <- xbreaks[-2]
            }else if(gid==nclass){
                fn[nclass-1] <- fn[nclass-1] + fn[nclass]
                fn <- fn[-nclass]
                xbreaks <- xbreaks[-(nclass-1)]
            }else{
                if(fn[gid-1] <= fn[gid+1]){
                    fn[gid-1] <- fn[gid-1] + fn[gid]
                    fn <- fn[-gid]
                    xbreaks <- xbreaks[-gid]
                }else{
                    fn[gid+1] <- fn[gid+1] + fn[gid]
                    fn <- fn[-gid]
                    xbreaks <- xbreaks[-(gid+1)]
                }
            }
        }else{
            break;
        }
    }
    list(breaks=xbreaks, counts=fn,nclass=length(fn))
}

.goflnorm <- function(x,parms){
    nx <- length(x)
    xgrp <- .gofbin(x,5)
    Fx <- plnorm(xgrp$breaks,meanlog = parms[1], sdlog = parms[2])
    ##    print(Fx)
    Fx <- c(0,Fx,1.0)
    Fa <- diff(Fx) * nx
    ##    print(as.numeric(Fa))
    Fb <- xgrp$counts
    ##    print(as.numeric(Fb))
    out <- chisq.test(cbind(Fa,Fb))
    p.gof <- out$p.value
}

.gofmlnorm <- function(x,parms){
    nx <- length(x)
    xgrp <- .gofbin(x,5)
    Fx <- pmlnorm(xgrp$breaks,parms[,1],parms[,2],parms[,3])
    Fx <- c(0,Fx,1.0)
    Fa <- diff(Fx) * nx
    Fb <- xgrp$counts
    out <- chisq.test(cbind(Fa,Fb))
    p.gof <- out$p.value
}

.mmelnorm <- function(x){
    nx <- length(x)
    x.table <- table(x)
    x.low <- as.numeric(names(x.table))
    fn <- as.numeric(x.table)
    n <- length(fn)

    ## counts
    x0 <- round(x.low)
    if(any(x0 != x.low)){
        delta1 <- diff(x.low)
        if(x.low[1] == 0){
            delta <- c(x.low[1], rep(min(delta1), n-1))
        }else{
            delta <- rep(min(delta1),n)
        }
    }else{
        delta <- rep(1,n)
    }
    xbar0 <- mean(x)
    xbar1 <- sum(fn*(x.low+delta))/sum(fn)
    var1 <- var(x) + (max(delta))^2/12
    s <- sqrt(log(1+var1/xbar0^2))
    mu <- log(xbar1)-s^2/2
    list(mean=mu, sd=s,delta=delta)
}

.rrmlelnorm <- function(x){
    nx <- length(x)
    x.table <- table(x)
    x.low <- as.numeric(names(x.table))
    fn <- as.numeric(x.table)
    n <- length(fn)

    ## counts
    x0 <- round(x.low)
    if(any(x0 != x.low)){
        delta1 <- diff(x.low)
        if(x.low[1] == 0){
            delta <- c(x.low[1], rep(min(delta1), n-1))
        }else{
            delta <- rep(min(delta1),n)
        }
    }else{
        delta <- rep(1,n)
    }
    y <- NULL
    for(i in 1:n){
        y1 <- x.low[i] + runif(fn[i],0,delta[i])
        y <- c(y,y1)
    }
    if(any(y<=0))
        y[y<=0] <- 0.0001
    y <- log(y)
    s <- sd(y)
    mu <- mean(y)
    list(mean=mu, sd=s,delta=delta)
}


##################
pmlnorm <- function(q,p,mean,sd){
    sapply(q,FUN=.pmlnorm, p=p,mu=mean,s=sd)
}

.pmlnorm <- function(x,p,mu,s){
    k <- length(p)
    stopifnot(length(x) == 1)
    stopifnot(length(mu) == k)
    stopifnot(length(s) == k)
    stopifnot(all(p>0))
    p <- p/sum(p)
    .myfun <- function(x){
        if(is.finite(x[2])){
            out <- plnorm(x[1],meanlog = x[2], sdlog = x[3])
        }else{
            out <- 1.0
        }
    }
    x0 <- rep(x, k)
    tmp <- cbind(x0,mu,s)
    sum(p*apply(tmp, 1, .myfun))
}

dmlnorm <- function(x,p,mean,sd){
    sapply(x,FUN=.dmlnorm, p=p,mu=mean,s=sd)
}

.dmlnorm <- function(x,p,mu,s){
    k <- length(p)
    stopifnot(length(x) == 1)
    stopifnot(length(mu) == k)
    stopifnot(length(s) == k)
    stopifnot(all(p>0))
    .myfun <- function(x){
        out <- 0
        if(is.finite(x[2])){
            out <- dlnorm(x[1],meanlog = x[2], sdlog = x[3])
        }
        out
    }
    x0 <- rep(x, k)
    tmp <- cbind(x0,mu,s)
    sum(p*apply(tmp, 1, .myfun))
}

rmlnorm <- function(n,p,mean,sd){
    mu <- mean
    s <- sd
    k <- length(p)
    stopifnot(length(n) == 1)
    n <- round(n)
    stopifnot(n>0)
    stopifnot(length(mu) == k)
    stopifnot(length(s) == k)
    stopifnot(all(p>0))
    .myfun <- function(x){
        if(is.finite(x[2])){
            out <- rlnorm(x[1],meanlog = x[2], sdlog = x[3])
        }else{
            out <- rep(0, x[1])
        }
        out
    }
    x0 <- rep(n, k)
    tmp <- cbind(n,mu,s)
    out <-  apply(tmp, 1, .myfun)
    gid <- sample(1:k,size=n,replace=TRUE,prob=p)
    gid <- cbind(1:n,gid)
    out[gid]
}


qmlnorm <- function(prob,p,mean,sd){
    sapply(prob,FUN=.qmlnorm, p=p,mu=mean,s=sd)
}

.qmlnorm <- function(x,p,mu,s){
    k <- length(p)
    stopifnot(length(x) == 1)
    stopifnot(length(mu) == k)
    stopifnot(length(s) == k)
    stopifnot(all(p>0))
    p <- p/sum(p)
    if(x <= 0){
        res <- 0
    }else if(x==1){
        res <- Inf
    }else{
        ## find xa and xb such that F(xa) > x and F(xb) <= x
        xa <- 0.1
        xb <- 10
        Fxa <- pmlnorm(xa,p,mu,s)
        Fxb <- pmlnorm(xb,p,mu,s)
        if(Fxa > x){
            xb <- xa
            for(i in 1:10000){
                xa <- 0.5 * xb
                Fxa <- pmlnorm(xa,p,mu,s)
                if(Fxa > x){
                    xb <- xa
                    xa <- 0.5 * xb
                }else{
                    ##cat("\n iter1=", i, "xa=",xa, "xb=",xb)
                    break
                }
            }
        }else if(Fxb < x){
            xa <- xb
            for(i in 1:10000){
                xb <- 2 * xa
                Fxb <- pmlnorm(xb,p,mu,s)
                if(Fxb < x){
                    xa <- xb
                    xb <- 2.0 * xb
                }else{
                    ##cat("\n iter1=", i, "xa=",xa, "xb=",xb)
                    break
                }
            }
        }
        for(i in 1:10000){
            xc <- (xa+xb)*0.5
            Fxc <- pmlnorm(xc,p,mu,s)
            if(Fxc < x){
                xa <- xc
            }else{
                xb <- xc
            }
            if(is.na(xb - xa)){
                res <- xa
                break
            }else if(xb - xa < 0.000001){
                res <- xc
                break
            }
        }
    }
    return(res)
}

.EMflnmm <- function(n0,x,f,k,mu,sig2)
{
    ## renew the number of groups (bins)
    nr <- length(f)
    ng <- k
    x0 <- x[-nr]; x1 <- x[-1]
    ## McLachlan & Jones (1988): in the case of not truncated
    ## data, the equations require two extra classes (-\infty, x0)
    ## and (xr,+\infty) with corresponding frequencies n_{l} and
    ## n_{u}.
    if(n0 > 0){
        itrunc <- 1;
        nl <- n0;
    }else{
        itrunc <- 0;
        nl <- 0
    }
    nu <- 0;

    

    ## initialize the parameters: p, mu and var (variance)
    p <- rep(1/ng, ng)
    mu <- (1+c(0:(ng-1))*0.1) * mu
    v <- (1+c(0:(ng-1))*0.1) * sig2
    ## input as WK (workspace size) and output as error info
    wk <- max(500, 3*ng*(nr+2+3*ng)+(nr+2))

    out <- .Fortran(.F_emmix,
                    as.integer(f), as.integer(nr), as.integer(ng),
                    as.double(x0), as.double(x1),
                    p=as.double(p), mu=as.double(mu), v=as.double(v),
                    xlogl=double(1), ifault=as.integer(wk),
                    as.integer(itrunc), 
                    nl=as.integer(nl), nu=as.integer(nu))
    
    llk0 <- out$xlogl
    
    ## sort the fitted components
    x.ord <- order(out$mu)
    
    pmix <- out$p[x.ord]
    mu <- out$mu[x.ord]
    s <- sqrt(out$v)[x.ord]
    sele <- (pmix <= 0) | (s <= 0)
    if(mean(sele) == 1)
        stop("model fitting failed")
    if(any(sele)){
        s <- s[!sele]
        mu <- mu[!sele]
        pmix <- pmix[!sele]
    }
    ## re-calculate npar
    ng <- length(mu)
    K <- ifelse(ng == 1, 2., 3. * ng - 1.);

    list(pmix=pmix,mu=mu,s=s, npar=K, llk=llk0)    
}

.em2lnorm <- function(x=x,k=k){
    x <- as.numeric(x)
    xmax <- max(x)
    nx <- length(x)
    nx0 <- sum(x==0)
    p0 <- nx0/nx

    x.table <- table(x)
    x.low <- as.numeric(names(x.table))
    fn <- as.numeric(x.table)

    out <- .mmelnorm(x)
    delta <- out$delta  ## use for MLE
    s <- out$sd
    
    if(p0==0){
        out <- fit.NGS(x)
        pmix <- out$fitted$p
        mu <- out$fitted$mu
        s <- out$fitted$s
    }else{
        xlow2 <- log(x.low[-1])
        fn2 <- fn[-1]
        nclass <- length(fn2)
        xmin <- xlow2[1]
        
        qp <- qnorm(p0)
        pmix <- rep(1/k,k)
        s <- rep(s,k)
        out <- .Fortran(.F_mlemixTN,
                        as.double(xlow2),
                        as.double(delta[-1]),
                        as.double(fn2),
                        converge=as.integer(nclass),
                        as.double(xmin),
                        as.double(qp),
                        sd = as.double(s),
                        pmix = as.double(pmix),
                        as.integer(k))
        
        pmix <- out$pmix
        pmix <- pmix/sum(pmix)
        s <- out$sd
        sele <- (pmix <= 0) | (s <= 0)
        if(any(sele)){
            s <- s[!sele]
            pmix <- pmix[!sele]
        }
        pmix <- pmix/sum(pmix)
        mu <- xmin - qp * s
    }

    K <- 3*length(pmix) - 1
    list(p=pmix,mu=mu,s=s,npar=K)
}
                      
.em1mmmlnorm <- function(x=x,k=k)
{
    p0 <- mean(x==0)
    ## create a histogram
    x.table <- table(x)
    x.low <- as.numeric(names(x.table))
    fn <- as.numeric(x.table)

    ## initial estimate by MME
    out <- .mmelnorm(x)
    delta <- out$delta  ## use for MLE
    mu0 <- out$mean
    v0 <- (out$sd)^2

    ## McLachlan & Jones (1988): in the case of not truncated
    ## data, the equations require two extra classes (-\infty, x0)
    ## and (xr,+\infty) with corresponding frequencies n_{l} and
    ## n_{u}.
    if(x.low[1]==0){
        itrunc <- 0;
        nl <- fn[1];
        x0 <- x.low[-1]
        x1 <- x0 + delta[-1]
        freq <- fn[-1]
    }else{
        itrunc <- 0;
        nl <- 0
        x0 <- x.low
        x1 <- x0 + delta
        freq <- fn
    }
    nu <- 0;

    ## renew the number of groups (bins)
    nr <- length(freq)
    ng <- k

    ## initialize the parameters: p, mu and var (variance)
    p <- rep(1/k, k)
    mu <- (1+c(0:(k-1))*0.1) * mu0
    v <- (1+c(0:(k-1))*0.1) * v0
    ## input as WK (workspace size) and output as error info
    wk <- max(500, 3*ng*(nr+2+3*ng)+(nr+2))
    ## if the percentage of zeros is larger than 10%, we model the
    ## data using a mixed mixture model with a zero component.
    if(p0 <= 0.1){
        out <- .Fortran(.F_emmix,
                        as.integer(freq), as.integer(nr), as.integer(ng),
                        as.double(log(x0)), as.double(log(x1)),
                        p=as.double(p), mu=as.double(mu), v=as.double(v),
                        xlogl=double(1), ifault=as.integer(wk),
                        as.integer(itrunc), 
                        nl=as.integer(nl), nu=as.integer(nu))
        llk.max <- out$xlogl
        K <- ifelse(ng == 1, 2., 3. * ng - 1.);
        ## sort the fitted components
        x.ord <- order(out$mu)
        
        pmix <- out$p[x.ord]
        mu <- out$mu[x.ord]
        s <- sqrt(out$v)[x.ord]
        sele <- (pmix <= 0) | (s <= 0)
        if(mean(sele) == 1)
            stop("model fitting failed")
        if(any(sele)){
            s <- s[!sele]
            mu <- mu[!sele]
            pmix <- pmix[!sele]
        }

    }else{
        M <- 50
        pstep <- seq(0.001,0.999,length=M)
        nl1 <- pstep * nl
        pnl2 <- (1 - pstep) * p0
        nl2 <- nl - nl1

        ## initialize

        out <- .Fortran(.F_emmix,
                        as.integer(freq), as.integer(nr), as.integer(ng),
                        as.double(log(x0)), as.double(log(x1)),
                        p=as.double(p), mu=as.double(mu), v=as.double(v),
                        xlogl=double(1), ifault=as.integer(wk),
                        as.integer(itrunc), 
                        nl=as.integer(nl1[1]), nu=as.integer(nu))
        llk1 <- out$xlogl
        llk.max <- llk1 + nl2[1] * log(pnl2[1])
        llkM <- 0
        out.opt <- out
        pnl <- pnl2[1]
        
        for(i in 2:M){
            out <- .Fortran(.F_emmix,
                            as.integer(freq), as.integer(nr), as.integer(ng),
                            as.double(log(x0)), as.double(log(x1)),
                            p=as.double(p), mu=as.double(mu), v=as.double(v),
                            xlogl=double(1), ifault=as.integer(wk),
                            as.integer(itrunc), 
                            nl=as.integer(nl1[i]), nu=as.integer(nu))
            llk1 <- out$xlogl
            llk2 <-llk1
            llk3 <- llk1 + nl2[i] * log(pnl2[i])
            if(llk2 > llk.max){
                llk.max <- llk2
                llkM <- llk3
                out.opt <- out
                pnl <- pnl2[i]
            }
        }
        out <- out.opt
        K <- ifelse(ng == 1, 2., 3. * ng - 1.);
        K <- K + 1 ## p of 0-comp

        ## sort the fitted components
        x.ord <- order(out$mu)
    
        pmix <- out$p[x.ord]
        mu <- out$mu[x.ord]
        s <- sqrt(out$v)[x.ord]
        sele <- (pmix <= 0) | (s <= 0)
        if(mean(sele) == 1)
            stop("model fitting failed")
        if(any(sele)){
            s <- s[!sele]
            mu <- mu[!sele]
            pmix <- pmix[!sele]
        }
        pmix <- c(pnl, (1-pnl)*pmix)
        mu <- c(-Inf, mu)
        s <- c(0, s)
    }
    list(p=pmix,mu=mu,s=s,npar=K)
}

.fitlnormPTest <- function(x){
    x <- as.numeric(x)
    xmax <- max(x)
    nx <- length(x)
    nx0 <- sum(x==0)
    p0 <- nx0/nx
    if(p0 < 0.10)
        stop("not applicable when percentage of zero measures is small")

    ## create a histogram
    x.table <- table(x)
    x.low <- as.numeric(names(x.table))
    fn <- as.numeric(x.table)

    ## estimate the parameters ################
    mu1 <- mean(log(0.001+x), na.rm=TRUE)
    s1 <- sd(log(0.001+x), na.rm=TRUE)
    p.gof1 <- .goflnorm(x,c(mu1,s1))
    meth <- "normal approximation"
    mu <- mu1; s <- s1; p.gof <- p.gof1
    
    
    ## MME ####################################
    out <- .mmelnorm(x)
    delta <- out$delta  ## use for MLE
    mu2 <- out$mean
    s2 <- out$sd
    
    
    p.gof2 <- .goflnorm(x,c(mu2,s2))
    if(p.gof < p.gof2){
        mu <- mu2; s <- s2; p.gof <- p.gof2;
        meth <- "method of moments"
    }

    nclass <- length(fn)
    
    out <- .Fortran(.F_mclnorm2,
                    as.double(x.low),
                    F=as.double(fn),
                    as.double(delta),
                    as.integer(nclass),
                    mean = as.double(mu),
                    sd = as.double(s))
    mu2 <- out$mean
    s2 <- out$sd
    p0hat <- out$F[1]
    if(p0hat < p0){
        pmix <- c(p0-p0hat, 1+p0hat-p0)
        mu <- c(-Inf,mu2)
        s <- c(0,s2)
        K <- 3
    }else{
        pmix <- 1.0
        mu <- mu2
        s <- s2
        K <- 2
    }
    list(p=pmix,mu=mu,s=s,npar=K)    
}

.fitlnorm <- function(x){
    x <- as.numeric(x)
    xmax <- max(x)
    nx <- length(x)
    nx0 <- sum(x==0)
    p0 <- nx0/nx

    ## create a histogram
    x.table <- table(x)
    x.low <- as.numeric(names(x.table))
    fn <- as.numeric(x.table)

    ## estimate the parameters ################
    mu1 <- mean(log(0.001+x), na.rm=TRUE)
    s1 <- sd(log(0.001+x), na.rm=TRUE)
    p.gof1 <- .goflnorm(x,c(mu1,s1))
    meth <- "normal approximation"
    mu <- mu1; s <- s1; p.gof <- p.gof1

    
    ## MME ####################################
    out <- .mmelnorm(x)
    delta <- out$delta  ## use for MLE
    mu2 <- out$mean
    s2 <- out$sd
    

    p.gof2 <- .goflnorm(x,c(mu2,s2))
    if(p.gof < p.gof2){
        mu <- mu2; s <- s2; p.gof <- p.gof2;
        meth <- "method of moments"
    }
    ## MLE: MC approach #######################

    Fn <- cumsum(fn)/nx
    nclass <- length(Fn)

    out <- .Fortran(.F_mclnorm,
                    as.double(x.low[-1]), as.double(Fn[-nclass]),
                    as.integer(nclass-1),
                    mean = as.double(mu),
                    sd = as.double(s))
    mu2 <- out$mean
    s2 <- out$sd
    
    p.gof2 <- .goflnorm(x,c(mu2,s2))
    if(p.gof < p.gof2){
        mu <- mu2; s <- s2; p.gof <- p.gof2
        meth <- "Monte Carlo method"
    }

    if(p0 == 0){
        out <- .rrmlelnorm(x)
        s2 <- out$sd
        mu2 <- out$mean
        p.gof2 <- .goflnorm(x,c(mu2,s2))
        if(p.gof < p.gof2){
            mu <- mu2; s <- s2; p.gof <- p.gof2
            meth <- "MLE (replication and reflection)"
        }
    }else{
        xlow2 <- log(x.low[-1])
        fn2 <- fn[-1]
        nclass <- length(fn2)
        xmin <- xlow2[1]
        
        qp <- qnorm(p0)
        out <- .Fortran(.F_mleTN,
                        as.double(xlow2),
                        as.double(delta[-1]),
                        as.double(fn2),
                        converge=as.integer(nclass),
                        as.double(xmin),
                        as.double(qp),
                        sd = as.double(s))
        if(out$converge==0){
            s2 <- out$sd
            mu2 <- xmin - qp * s2
            p.gof2 <- .goflnorm(x,c(mu2,s2))
            if(p.gof < p.gof2){
                mu <- mu2; s <- s2; p.gof <- p.gof2
                meth <- "MLE for truncated normal"
            }
        }
    }
    llk0 <- .LLKfnmm(x,1.0,mu,s)
    
    list(p=1.0,mu=mu,s=s, meth=meth,npar=2,llk=llk0)    
}

.em3lnorm <- function(x=x,k=k){
    x <- as.numeric(x)
    xmin <- min(x[x>0])
    nx <- length(x)
    nx0 <- sum(x==0); p0 <- nx0/nx
    stopifnot(p0>0.05)
    tmp <- table(x)
    xt <- as.numeric(names(tmp))
    fn <- as.numeric(tmp)
    nclass <- length(fn)
    ## initial estimate by MME
    out <- .mmelnorm(x)
    delta <- out$delta  ## use for MLE

    out <- .fitlnorm(x=x)
    pmix <- out$p
    mu <- out$mu;
    s <- out$s;
    
    ##    xbar1 <- log(mean(x)); xbar2 <- log(xbar1 + 1.)
    z0 <- qnorm(p0);
    ## the mean expression mu <= 0.  When p0 is small, z0 is negative
    ## and mu could be large xbar = x(1) - z0*sigma. We define the range
    mu.min <- 0.5 * mu;  mu.max <- 1.5 * mu 
    mus <- seq(mu.min, mu.max, length=100)
    sig.min <- 0.1 * s;  sig.max <- 1.5 * s 
    sigs <- seq(sig.min, sig.max, length=30)
    Dmin <- 1000000
    MU0 <- NULL; S0 <- NULL; P0 <- NULL
    for(mu.hat in mus){
        for(s.hat in sigs){
            Fn <- plnorm(xt[-1], mu.hat, s.hat)
            Fn <- diff(c(0,Fn, 1.0)) * nx
            fn.hat <- round(fn - Fn)
            fn.hat <- fn.hat[-1] # exclude the class of [0,xmin)
            sele <- fn.hat > 0
            phat <- sum(fn.hat[sele])/nx
            if(sum(sele) >=3 &&  phat > 0.05){
                nu <- 0; nl <- 0
                itrunc <- 0;
                nr <- sum(sele)
                ng <- k - 1
                x0 <- (xt[-1])[sele]
                x1 <- ((xt+delta)[-1])[sele]
                xtemp <- rep(0.5*(x0+x1), fn.hat[sele])
                mu <- mean(log(xtemp))*runif(ng,0.5,0.99)
                v <- var(log(xtemp))*runif(ng,0.5,0.99)
                p <- rep(1/ng, ng)
                wk <- max(500, 3*ng*(nr+2+3*ng)+(nr+2))
                
                out <- .Fortran(.F_emmix,
                                as.integer(fn.hat[sele]), as.integer(nr),
                                as.integer(ng), as.double(log(x0)),
                                as.double(log(x1)), p=as.double(p),
                                mu=as.double(mu), v=as.double(v),
                                xlogl=double(1), ifault=as.integer(wk),
                                as.integer(itrunc), nl=as.integer(nl),
                                nu=as.integer(nu))
                x.ord <- order(out$mu)
                pmix <- c(1-phat, phat*out$p[x.ord])
                mu <- c(mu.hat, out$mu[x.ord])
                s <- c(s.hat, sqrt(out$v)[x.ord])
                ##            print(pmix)
                ## remove zero components
                sele <- pmix > 0
                pmix <- pmix[sele]
                pmix <- pmix/sum(pmix)
                mu <- mu[sele]
                s <- s[sele]
                
                D1 <- .L2mlnorm(x,pmix,mu,s)
                if(D1 < Dmin){
                    MU0 <- mu
                    S0 <- s
                    P0 <- pmix
                    Dmin <- D1
                }
            }
        }
    }
    
    K <- 3 * k - 1
    list(p=P0,mu=MU0,s=S0,npar=K)
}

.L2mlnorm <- function(x,p,mu,s){
    nx <- length(x)
    tmp <- table(x)
    xt <- as.numeric(names(tmp))
    fn <- as.numeric(tmp)
    nclass <- length(fn)
    Fx <- pmlnorm(xt[-1],p,mu,s)
    Fx <- c(0,Fx,1.0)
    Fa <- diff(Fx) * nx
    sum((Fa - fn)^2/fn)
}

.moplnorm <- function(x){
    tmp <- table(x)
    x <- as.numeric(names(tmp))
    f <- as.numeric(tmp)
    ## the first quntile should be at least 50%, then search for the
    ## second one as large as possible to ensure the variance is
    ## appropriately estimated.
    delta <- 0.00001
    n <- sum(f)
    k <- length(f)
    pf <- cumsum(f)/n
    if(k == 1){
        p1 <- 1.0 - delta*1.1/n
        x1 <- log(x[1]+1)
        p2 <- 1.0 - delta*.1/n
        x2 <- log(x[1]+2)
    }else if(k==2){
        p1 <- pf[1]
        x1 <- log(x[1]+1)
        p2 <- pf[2] - delta*1.1/n
        x2 <- log(x[2]+1)
    }else{
        p2 <- pf[k-1]
        x2 <- log(x[k])
        p0 <- 0.5 * p2
        isele <- min(k-2,which(pf>=p0)[1])
        p1 <- pf[isele]
        x1 <- log(x[isele]+1)
    }
    
    z1 <- qnorm(p1)
    z2 <- qnorm(p2)
    sig <- (x2-x1)/(z2-z1)
    mu <- x1 - sig * z1
    if(sig <= 0){
        x0 <- rep(x,f)+runif(sum(f))+0.0000001
        mu <- mean(log(x0)) - 1
        sig <- sd(log(x0))
    }
    llk0 <- .LLKfnmm(x,1.0,mu,sig)
    list(p=1,mu=mu,s=sig,llk=llk0)
}

.em4mlnorm <- function(x,k)
{
    stopifnot(k>1)
    n0 <- sum(x==0)
    p0 <- mean(x==0)
    n1 <- sum(x>0)
    nx <- length(x)
    
    ## create a histogram
    xp <- log(x[x>0])
    xmin <- min(xp)
    x.table <- table(xp)
    xbrks <- c(as.numeric(names(x.table)), max(x) + 1)
    fn <- as.numeric(x.table)
    ng <- length(fn)

    ## initial the estimates
    mu0 <- mean(xp)
    v0 <- var(xp)
    
    if(all(fn <= 2) && ng >= 10){
        xh <- hist(xp, plot=FALSE)
        xbrks <- xh$breaks
        fn <- xh$counts
        sele <- xh$counts == 0
        if(any(sele)){
            isele <- which(sele) + 1
            xbrks <- xbrks[-isele]
            xh <- hist(xp, breaks=xbrks, plot=FALSE)
            fn <- xh$counts
        }
        xbrks[1] <- xmin
    }

    ## we choose the best fit using BIC.  We first fit FNMM based on
    ## X>0.  Check if the zeros can be completely covered. if not, fit
    ## another model using the not covered data (MOP).
    ## fit unimodal lognormal as the baseline
    out <- .fitlnorm(x=x)
    pmix0 <- 1.0 
    mu0 <- out$mu;
    s0 <- out$s;
    npar0 <- out$npar
    llk0 <- out$llk

    ##print(c(pmix0, mu0, s0,npar0,llk0))

    
    ##    BIC.best <- npar0 * log(nx/2/pi) - 2.0 * llk0
    ##    K <- ifelse(ng == 1, 2., 3. * ng - 1.);
    ##    N <- nx;
    ##    AIC = -2.0 * llk0 + 2.0 * K;
    ##    AICc = AIC + 2.0 * K * (K+1.) / (N - K - 1.0);
    ##    BIC = -2.0 * llk0 + log(0.5*N/pi) * K; ## do we need a pi here?

    f.best <- out
    if(p0 < 0.01){
        for(i in 1:k){
            out <- .EMflnmm(n0=n0,x=xbrks,f=fn,k=i,mu=mu0,sig2=v0)
            pmix <- out$pmix 
            mu <- out$mu;
            s <- out$s;
            npar <- out$npar
            llk <- .LLKfnmm(x, pmix,mu,s)
            ##Fx <- pmlnorm(xbrks[-1],pmix,mu,s)
            ##llk <- sum(log(diff(c(0,Fx,1)))*fn)
            ## choose best model: if return TRUE, switch, if FALSE, don't
            ## do anything.
            swModel <- .chooseModel(llk0, npar0, llk, npar,n=nx)
            ##cat("p0=",p0,"ng=",i,"switch model?",swModel,"\n")
            if(swModel){
                f.best <- out
                pmix0 <- pmix
                mu0 <- mu
                s0 <- s
                llk0 <- llk
                npar0 <- npar
            }else{
                break;
            }            
        }
    }else{
        for(i in 1:(k-1)){
            ## i <- 1 # iterate for i = k-1 (one component is used for
            ## 0's)
            out <- .EMflnmm(n0=0,x=xbrks,f=fn,k=i,mu=mu0,sig2=v0)
            pmix <- out$pmix 
            mu <- out$mu;
            s <- out$s;
            npar <- out$npar
            
            Fx <- pmlnorm(xbrks[c(1,2)],pmix,mu,s)
            nhat <- diff(c(0,Fx)) * n1 # the obs covered
            if(nhat[1] > n0){
                if(nhat[2] > fn[1]){
                    ntot <- n0 + fn[1] - sum(nhat)
                    p1 <- n0/ntot
                    x1 <- xbrks[1]
                    p2 <- 1 - 0.000011/ntot
                    x2 <- xbrks[2]
                    pmix <- c(ntot/nx, (1-ntot/nx)*pmix)
                }else{
                    ntot <- n0 - nhat[1]
                    p1 <- 1 - 0.000011/ntot
                    x1 <- xbrks[1]
                    p2 <- 1 - 0.000001/ntot
                    x2 <- xbrks[2]
                    pmix <- c(ntot/nx, (1-ntot/nx)*pmix)
                }
                z1 <- qnorm(p1)
                z2 <- qnorm(p2)
                sig1 <- (x2-x1)/(z2-z1)
                mu1 <- x1 - sig1 * z1
                mu <- c(mu1, mu)
                s <- c(sig1, s)
                npar <- npar + 2
            }
            ##Fx <- pmlnorm(xbrks[-1],pmix,mu,s)
            ##llk <- sum(log(diff(c(0,Fx,1)))*fn)
            llk <- .LLKfnmm(x, pmix,mu,s)
            
            swModel <- .chooseModel(llk0, npar0, llk, npar,n=nx)
            ##cat("p0=",p0,"ng=",i,"switch model?",swModel,"\n")
            if(swModel){
                f.best <- out
                pmix0 <- pmix
                mu0 <- mu
                s0 <- s
                llk0 <- llk
                npar0 <- npar
            }else{
                break;
            }            
        }
    }
    
    
    list(p=pmix0,mu=mu0,s=s0,llk=llk0,npar=npar0)
}

.chooseModel <- function(llk0, npar0, llk, npar,n){
    ## if DF0 and DF are different, use LR-test, otherwise choose the
    ## one with smaller BIC    
    res <- FALSE
    DF <- npar - npar0
    D <- 2 * (llk - llk0)
        
    if(!is.na(D)){
        if(DF == 0){
            BIC0 <- npar0 * log(n/2/pi) - 2.0 * llk0
            BIC <- npar * log(n/2/pi) - 2.0 * llk
            if(BIC < BIC0){
                res <- TRUE
            }
        }else if(DF > 0){
            if(D > 0){
                p.v <- pchisq(D, DF, lower.tail=FALSE)
                if(p.v < 0.05){
                    res <- TRUE
                }
            }
        }else{
            if(D > 0){
                res <- TRUE
            }else{
                p.v <- pchisq(-D, abs(DF), lower.tail=FALSE)
                if(p.v >= 0.05){
                    res <- TRUE
                }
            }
        }
    }
    ##cat("LLk=(",llk0,",",llk,"), npar=c(",
    ##    npar0,",",npar,")\nD=", D, ", DF=", DF,
    ##    "Reject H0?:", res,"\n")

    res
}

.LLKfnmm <- function(x,pmix,mu,s){
    tmp <- c(pmix,mu,s)
    if(any(is.null(tmp))){
        llk0 <- -Inf
    }else if(any(is.na(tmp))){
        llk0 <- -Inf
    }else{        
        ## create a histogram
        x.table <- table(x)
        x.low <- as.numeric(names(x.table))
        fn <- as.numeric(x.table)
        if(length(x.low)==1){
            tmp <- x.low
            if(tmp==0) tmp <- 1
        }else{
            tmp <- x.low[-1]
        }
            
        Fx <- pmlnorm(tmp,pmix,mu,s)
        tmp <- log(diff(c(0,Fx,1)))
        isele <- which(!is.finite(tmp))
        if(any(isele)) tmp[isele] <- -309
        llk0 <- sum(tmp*fn)
    }
    llk0
}

.em5mlnorm <- function(x,k)
{
    stopifnot(k>1)
    n0 <- sum(x==0)
    p0 <- mean(x==0)
    n1 <- sum(x>0)
    nx <- length(x)
    
    ## create a histogram
    xp <- log(x[x>0])
    xmin <- min(xp)
    x.table <- table(xp)
    xbrks <- c(as.numeric(names(x.table)), max(x) + 1)
    fn <- as.numeric(x.table)
    ng <- length(fn)

    ## initial the estimates
    mu0 <- mean(xp)
    v0 <- var(xp)
    
    if(all(fn <= 2) && ng >= 10){
        xh <- hist(xp, plot=FALSE)
        xbrks <- xh$breaks
        fn <- xh$counts
        sele <- xh$counts == 0
        if(any(sele)){
            isele <- which(sele) + 1
            xbrks <- xbrks[-isele]
            xh <- hist(xp, breaks=xbrks, plot=FALSE)
            fn <- xh$counts
        }
        xbrks[1] <- xmin
    }

    ## we choose the best fit using BIC.  We first fit FNMM based on
    ## X>0.  Check if the zeros can be completely covered. if not, fit
    ## another model using the not covered data (MOP).
    ## fit unimodal lognormal as the baseline
    out <- .fitlnorm(x=x)
    pmix0 <- 1.0 
    mu0 <- out$mu;
    s0 <- out$s;
    npar0 <- out$npar
    llk0 <- out$llk

    ##print(c(pmix0, mu0, s0,npar0,llk0))

    
    ##    BIC.best <- npar0 * log(nx/2/pi) - 2.0 * llk0
    ##    K <- ifelse(ng == 1, 2., 3. * ng - 1.);
    ##    N <- nx;
    ##    AIC = -2.0 * llk0 + 2.0 * K;
    ##    AICc = AIC + 2.0 * K * (K+1.) / (N - K - 1.0);
    ##    BIC = -2.0 * llk0 + log(0.5*N/pi) * K; ## do we need a pi here?

    f.best <- out
    if(p0 < 0.01){
        i <- k
        out <- .EMflnmm(n0=n0,x=xbrks,f=fn,k=i,mu=mu0,sig2=v0)
        pmix <- out$pmix 
        mu <- out$mu;
        s <- out$s;
        npar <- out$npar
        llk <- .LLKfnmm(x, pmix,mu,s)
        swModel <- .chooseModel(llk0, npar0, llk, npar,n=nx)
#        if(swModel){
            f.best <- out
            pmix0 <- pmix
            mu0 <- mu
            s0 <- s
            llk0 <- llk
            npar0 <- npar
#        }
    }else{
        i <- k - 1
        out <- .EMflnmm(n0=0,x=xbrks,f=fn,k=i,mu=mu0,sig2=v0)
        pmix <- out$pmix 
        mu <- out$mu;
        s <- out$s;
        npar <- out$npar
        
        Fx <- pmlnorm(xbrks[c(1,2)],pmix,mu,s)
        nhat <- diff(c(0,Fx)) * n1 # the obs covered
        if(nhat[1] > n0){
            if(nhat[2] > fn[1]){
                ntot <- n0 + fn[1] - sum(nhat)
                p1 <- n0/ntot
                x1 <- xbrks[1]
                p2 <- 1 - 0.000011/ntot
                x2 <- xbrks[2]
                pmix <- c(ntot/nx, (1-ntot/nx)*pmix)
            }else{
                ntot <- n0 - nhat[1]
                p1 <- 1 - 0.000011/ntot
                x1 <- xbrks[1]
                p2 <- 1 - 0.000001/ntot
                x2 <- xbrks[2]
                pmix <- c(ntot/nx, (1-ntot/nx)*pmix)
            }
            z1 <- qnorm(p1)
            z2 <- qnorm(p2)
            sig1 <- (x2-x1)/(z2-z1)
            mu1 <- x1 - sig1 * z1
            mu <- c(mu1, mu)
            s <- c(sig1, s)
            npar <- npar + 2
        }
        llk <- .LLKfnmm(x, pmix,mu,s)
        swModel <- .chooseModel(llk0, npar0, llk, npar,n=nx)
#        if(swModel){
            f.best <- out
            pmix0 <- pmix
            mu0 <- mu
            s0 <- s
            llk0 <- llk
            npar0 <- npar
#        }            
    }
        
    list(p=pmix0,mu=mu0,s=s0,llk=llk0,npar=npar0)
}

.fitlnormOPT <- function(x){
    out <- .fitlnorm(x=x)
    out1 <- .moplnorm(x)
    if(out$llk < out1$llk){
        out <- out1
    }
    out
}


## Firm size data are usually of histogram-type. If not, we bin the
## data first then fit lognormal distributions.
fit.lnmix <- function(x, method='lognormal',
                    k=1, lbound, ubound){
    if(class(x)=="histogram")
        x <- binning(x)

    if(class(x) != 'bdata')
        stop("data type not supported")

    freq <- x$counts;
    breaks <- x$breaks
    ## data validation
    if(any(is.na(freq)))
        stop("missing value(s) found in 'freq'")
    if(any(freq<=0))
        stop("non-positive frequencies not allowed")
    tmp <- diff(breaks)
    if(any(is.na(tmp)))
        stop("missing value(s) found in class limits")
    if(any(tmp <= 0))
        stop("invalid values found in class limits")
    
    nclass <- length(freq)
    xhist <- binning(counts=freq, breaks=breaks)
    
    ##if(breaks[1] == 0) breaks[1] <- 0.001
    ##if(!is.finite(breaks[nclass+1]))
    ##breaks[nclass+1] <- 2*breaks[nclass]
    
    ll <- breaks[-(nclass+1)]; ul <- breaks[-1]
    nx <- sum(freq)
    rf <- freq/nx
    xhist1 <- data.frame(Lower.limit=ll, # show raw data
                         Upper.limit=ul,
                         Counts=round(freq),
                         Percent=round(rf*100,2)
                         )
    
    method <- match.arg(tolower(method),
                        c("lognormal", "normal", "gld","fmkl"))

    ntotal <- sum(freq)
    
    if(method=="lognormal"){
        k <- round(k)
        stopifnot(k>0)
        out <- .fitFSDbin(breaks, freq, k,1,
                          lbound=lbound,ubound=ubound)
        pars <- out$pars
        meth <- "lognormal"
        
        if(k == 1){
            K <- 2
        }else{
            K <- 3*k - 1
        }
        ## compute LLK
        tmp <- pmlnorm(breaks, pars$p, pars$mean,pars$sd);
        llk0 <- sum(log(abs(diff(tmp)))*freq)
        AIC = -2.0 * llk0 + 2.0 * K;
        AICc = AIC + 2.0 * K * (K+1.) / (ntotal - K - 1.0);
        BIC = -2.0 * llk0 + log(0.5*ntotal/pi) * K; 
        ModSel <- data.frame(AIC=AIC, BIC=BIC, AICc=AICc)
    }else if(method=="normal"){
        k <- round(k)
        stopifnot(k>0)
        out <- .fitFSDbin(breaks, freq, k,0,
                          lbound=lbound,ubound=ubound)
        pars <- out$pars
        meth <- "normal"

        if(k == 1){
            K <- 2
        }else{
            K <- 3*k - 1
        }
        ## compute LLK
        tmp <- pmixnorm(breaks, pars$p, pars$mean,pars$sd);
        llk0 <- sum(log(abs(diff(tmp)))*freq)
        AIC = -2.0 * llk0 + 2.0 * K;
        AICc = AIC + 2.0 * K * (K+1.) / (ntotal - K - 1.0);
        BIC = -2.0 * llk0 + log(0.5*ntotal/pi) * K; 
        ModSel <- data.frame(AIC=AIC, BIC=BIC, AICc=AICc)
    }else{
        k <- round(k)
        stopifnot(k>0)
        out <- .fitFSDbin(breaks, freq, k,1,
                          lbound=lbound,ubound=ubound)
        pars <- out$pars
        meth <- "FMKL-GLD"

        K <- 4
        ## compute LLK
        tmp <- pgld(breaks, pars);
        llk0 <- sum(log(abs(diff(tmp)))*freq)
        AIC = -2.0 * llk0 + 2.0 * K;
        AICc = AIC + 2.0 * K * (K+1.) / (ntotal - K - 1.0);
        BIC = -2.0 * llk0 + log(0.5*ntotal/pi) * K;
        ModSel <- data.frame(AIC=AIC, BIC=BIC, AICc=AICc)
    }

    structure(list(data = xhist1, ###
                   xhist = xhist,
                   pars = pars,
                   ModSel = ModSel,
                   method = meth
                   ),
              class='lnmix')
}



plot.lnmix <- function(x,xlim=NULL,ngrid=1000,...){
    if(is.null(xlim)){
        xbrks <- x$xhist$breaks
        k <- length(xbrks)
        if(is.finite(xbrks[1]))
            xmin <- xbrks[1]
        else
            xmin <- xbrks[2]
        if(is.finite(xbrks[k]))
            xmax <- xbrks[k]
        else
            xmax <- xbrks[k-1]
        
        xlim <- c(xmin, xmax)
    }
    plot(x$xhist,xlim=xlim,...)
    lines(x,xlim=xlim,ngrid=ngrid,...)
}

lines.lnmix <- function(x,xlim=NULL,ngrid=1000,...){
    ngrid <- round(ngrid)
    stopifnot(ngrid > 1)
    if(is.null(xlim)){
        xbrks <- x$xhist$breaks
        k <- length(xbrks)
        if(is.finite(xbrks[1]))
            xmin <- xbrks[1]
        else
            xmin <- xbrks[2]
        if(is.finite(xbrks[k]))
            xmax <- xbrks[k]
        else
            xmax <- xbrks[k-1]
        x0 <- seq(xmin, xmax, length=ngrid)
    }else{
        x0 <- seq(xlim[1], xlim[2], length=ngrid)
    }

    if(x$method=='normal'){
        f0 <- dmixnorm(x0,x$pars$p,x$pars$mean,x$pars$sd)
    }else if(x$method=='lognormal'){
        f0 <- dmlnorm(x0,x$pars$p,x$pars$mean,x$pars$sd)
    }else if(x$method=='FMKL-GLD'){
        f0 <- dgld(x0, x$pars)
    }else{
        f0 <- rep(1/diff(range(x0)),length(x0))
    }

    lines(x0,f0,...)
}

print.lnmix <- function(x,...){
    if(x$method=='normal'){
        cat("\nNormal Mixture:")
        k <- nrow(x$pars)
        cat(k, " components")
        cat("\nParameter estimates:\n")
        print(round(x$pars,3))
    }else if(x$method=='lognormal'){
        cat("\nLognormal Mixture:")
        k <- nrow(x$pars)
        cat(k, " components")
        cat("\nParameter estimates:\n")
        print(round(x$pars,3))
    }else{
        cat("\nFKML-GLD:")
        cat("\nParameter estimates:\n")
        print(round(x$pars,3))
    }
    cat("\nModel selection:\n")
    print(x$ModSel)
    cat("\n")
}

######################################################################

.gofmlnormBin <- function(xbrks, freq, pmix, mu, s){
    nx <- sum(freq, na.rm=TRUE)
    k <- length(xbrks)
    Fx <- pmlnorm(xbrks[-c(1,k)],pmix,mu,s)
    Fx <- c(0,Fx,1.0)
    Fa <- diff(Fx) * nx
    out <- chisq.test(cbind(Fa, freq))
    list(stat=out$stat, p.value=out$p.value)
}

.llkmlnormBin <- function(xbrks, freq, pmix, mu, s){
    nx <- sum(freq, na.rm=TRUE)
    k <- length(xbrks)
    Fx <- pmixnorm(log(xbrks[-c(1,k)]),pmix,mu,s)
    Fx <- c(0,Fx,1.0)
    sum(log(diff(Fx)) * freq)
}

.fitlnormBin <- function(x, method=1){
    nx <- sum(x$counts, na.rm=TRUE)
    xbrks <- x$breaks
    nclass <- length(xbrks) - 1
    xbrks <- log(xbrks)
    xbrks[1] <- -log(xbrks[nclass])
    ## range of mean
    xt0 <- xbrks[-(nclass+1)]
    mu0 <- sum(xt0 * x$counts) / nx
    xt1 <- xbrks[-1]; xt1[nclass] <- xt1[nclass-1] * 10
    mu1 <- sum(xt1 * x$counts) / nx
    ## variance
    xc <- 0.5 * (xt0 + xt1)
    mu <- sum(xc * x$counts) / nx
    sig2 <- sum(x$counts * xc^2)/nx - mu^2
    sig <- sqrt(sig2)

    ##cat("\nSample means", mu0, mu1, mu,"\n")
    if(method==1){
        out <- .Fortran(.F_lnormBinMLE,
                        as.integer(nclass-1),
                        as.double(xbrks[-c(1,nclass+1)]),
                        as.double(x$counts),
                        mean = as.double(c(mu0,mu1)),
                        sd = as.double(sig))
    }else{
        out <- .Fortran(.F_lnormBinChisq,
                        as.integer(nclass-1),
                        as.double(xbrks[-c(1,nclass+1)]),
                        as.double(x$counts),
                        mean = as.double(c(mu0,mu1)),
                        sd = as.double(sig))
    }
    
    list(p=1,mu=out$mean[1],s=out$sd)    
}

######################################################################

.fitFSD1 <- function(xbrks, freq, sig){
    nx <- sum(freq)
    nx0 <- freq[1]
    p0 <- nx0/nx
    fn <- freq
    Fn <- cumsum(fn)/nx

    ## create a histogram
    nclass <- length(fn)
    x.low <- xbrks[-(nclass + 1)]
    x.hi <- xbrks[-1]
    if(!is.finite(x.hi[nclass]))
        x.hi[nclass] <- 3*x.hi[nclass-1]-2*x.hi[nclass-2]

    xlow2 <- log(x.low[-1])
    fn2 <- fn[-1]
    nclass <- length(fn2)
    xmin <- xlow2[1]
        
    delta <- x.hi - x.low
    qp <- qnorm(p0)
    
    s <- log(sig)/5
    mu <- xmin - qp * s
    
    out <- .Fortran(.F_mclnorm2,
                    as.double(x.low),
                    F=as.double(fn),
                    as.double(delta),
                    as.integer(nclass),
                    mean = as.double(mu),
                    sd = as.double(s))
    mu <- out$mean
    s <- out$sd

    list(p=1.0,mu=mu,s=s)    
}


.fitFSDMLE <- function(xbrks, freq, mu,s2){
    nclass <- length(freq)
    ## infinite or mi values are not allowed to pass to C/Fortran. We
    ## need to remove the first and last breaks.
    out <- .Fortran(.F_lnormBinMLE,
                    as.integer(nclass),
                    as.double(xbrks[-c(1,nclass+1)]),
                    F=as.double(freq),
                    mean = as.double(log(mu)),
                    sd = as.double(log(sqrt(s2)))
                    )
    
    list(p=1.0,mu=out$mean,s=out$sd)    
}

.fitFSDk <- function(freq, xbrks, k){
    nx <- sum(freq)
    nclass <- length(freq)
    ## if there are five or less groups, we don't use this model
    ## fitting algorithm
    stopifnot(nclass >= 5)
    
    k <- round(k) - 1
    stopifnot(k >= 1)
    ## estimate parameters of the lognormal based on the first two
    ## classes.
    tmp <- cumsum(freq)/nx
    p1 <- tmp[1]; q1 <- qnorm(p1)
    p2 <- tmp[2]; q2 <- qnorm(p2)
    x1 <- log(xbrks[2])
    x2 <- log(xbrks[3])

    if(q2 - q1 == 0){
        s <- 0.3 * (xbrks[nclass] - xbrks[2]) # approx.
    }else{
        s <- (x2 - x1)/(q2 - q1)
    }
    mu <- x2 - q2 * s
    pmix <- 1

    ## W1 and W2
    W1 <- freq[1] + freq[2]
    icount <- 0
    ll <- NULL; ul <- NULL
    f2 <- NULL
    for(i in 1:(nclass-2)){
        a <- log(xbrks[i+2])
        b <- log(xbrks[i+3])
        n3 <- nx * (pnorm(b,mu,s) - pnorm(a,mu,s))
        if(n3 < freq[i+2]){
            icount <- icount + 1
            ll <- c(ll, a)
            ul <- c(ul, b)
            f2 <- c(f2,freq[i+2] - n3)
            W1 <- W1 + n3
        }else{
            W1 <- W1 + freq[i+2]
            ##cat("\n class=", i+2, "Freq=", freq[i+2], ", E(x)=", n3,"\n")
        }
    }
    if(!is.finite(ul[icount])){
        ul[icount] <- 5 * ul[icount-1]
    }

##    print(rbind(ll,ul,f2))
    
    if(icount >= 2){
        out <- .emmix(f2, ll, ul, k)
        P1 <- W1/nx
        pmix <- c(P1, (1-P1) * out$pmix)
        mu <- c(mu, out$mu)
        s <- c(s, out$s)
    }
    
    list(pmix=pmix, mu=mu, s=s)
}

## the following function is used to fit normal mixture models based
## on binned data. 'll' and 'ul' need to be log-transformed if it is
## used to fit lognormal mixtures.
.emmix <- function(freq, ll, ul, k){
    if(any(ll >= ul))
        stop("Wrong class limit(s)")
    if(any(freq <= 0))
        stop("Non-postive frequencies not allowed")
    ## renew the number of groups (bins)
    nr <- length(freq)
    nx <- sum(freq)
    ng <- round(k)
    if(ng<1)
        stop("Invalud number of components")
    if(length(ll) != nr)
        stop("wrong length of lower limits")
    if(length(ul) != nr)
        stop("wrong length of upper limits")
    
    ## initialize the parameters: p, mu and var (variance)
    p <- rep(1/k, k)
    ## sample mean using mid-points
    x0 <- 0.5 * (ll+ul)
    if(!is.finite(ll[1])){
        x0[1] <- ul[1]
        ll[1] = min(1, ul[1]*0.1)
    }
    if(!is.finite(ul[nr])){
        x0[nr] <- 2*ll[nr]
        ul[nr] <- 2*ll[nr]
    }

    xbar <- sum(freq*x0)/nx
    vbar <- sum(freq*x0^2)/nx - xbar^2

    if(k > 1){
        xbar <- (1+c(0:(k-1))*0.2) * xbar
        vbar <- rep(vbar, k)
    }
    
    ## input as WK (workspace size) and output as error info
    wk <- max(500, 3*ng*(nr+2+3*ng)+(nr+2))
    itrunc <- 0
    nl <- 0; nu <- 0

    out <- .Fortran(.F_emmix,
                    as.integer(freq), as.integer(nr), as.integer(ng),
                    as.double(ll), as.double(ul),
                    p=as.double(p), mu=as.double(xbar), v=as.double(vbar),
                    xlogl=double(1), ifault=as.integer(wk),
                    as.integer(itrunc), 
                    nl=as.integer(nl), nu=as.integer(nu))
    ## sort the fitted components
    x.ord <- order(out$mu)
    pmix <- out$p[x.ord]
    mu <- out$mu[x.ord]
    s <- sqrt(out$v)[x.ord]
    sele <- (pmix <= 0) | (s <= 0)
    if(any(sele)){
        s <- s[!sele]
        mu <- mu[!sele]
        pmix <- pmix[!sele]
        pmix <- pmix/sum(pmix)
    }
    list(pmix=pmix, mu=mu,s=s)
}


.musd.bdata <- function(breaks,f){
    ## compute sample mean and SD
    nclass <- length(f)
    stopifnot(length(breaks)==nclass+1)
    x0 <- breaks[-1] - 0.5 * diff(breaks)
    if(!is.finite(breaks[nclass+1])){
        x0 <- x0[-nclass]
        f <- f[-nclass]
    }
    if(!is.finite(breaks[1])){
        x0 <- x0[-1]
        f <- f[-1]
    }
    nx <- sum(f)
    xbar <- sum(f*x0)/nx
    s2 <- sum(f*x0^2)/nx - xbar^2
    list(mean=xbar,sd=sqrt(s2))
}


.fitFSDbin <- function(xbrks, freq, k,lognormal,lbound,ubound){
    stopifnot(is.numeric(lognormal))
    if(lognormal != 1) lognormal <- 0

    nr <- length(freq)
    ng <- round(k)
    nl <- 0; nu <- 0
    stopifnot(ng>0)
    
    if(!missing(lbound)){
        stopifnot(lbound < xbrks[2])
        xbrks[1] <- lbound
    }
    if(!missing(ubound)){
        stopifnot(ubound > xbrks[nr])
        xbrks[nr+1] <- ubound
    }
    if(lognormal==1){ # do logarithm transformation if lognormal
        if(any(xbrks < 0))
            stop("lognormal is not suitable for this data")
        xbrks <- log(xbrks)
    }

    ntotal <- sum(freq)
    n.5 <- round(ntotal * 0.05)
    rfreq <- cumsum(freq)/ntotal
    if(xbrks[1] == -Inf){
        if(rfreq[1] < 0.05){
            nl <- freq[1]
            xbrks <- xbrks[-1]
            freq <- freq[-1]
        }else{
            q1 <- xbrks[2]
            q2 <- xbrks[3]
            qn1 <- qnorm(rfreq[1])
            qn2 <- qnorm(rfreq[2])
            s <- (q1-q2)/(qn1-qn2)
            mu <- q1 - qn1 * s
            nl <- min(ceiling(freq[1]*0.01), n.5)
            freq[1] <- freq[1] - nl
            qhat <- qnorm(nl/ntotal,mu,s)
            xbrks[1] <- qhat
            ##cat("\nnl = ", nl, 'lower limit =',
            ##round(qhat,3), " (", round(exp(qhat),3),")")
        }
    }

    nr <- length(xbrks) - 1
    if(xbrks[nr+1] == Inf){
        if(rfreq[nr] < 0.05){
            nu <- freq[nr]
            xbrks <- xbrks[-(nr+1)]
            freq <- freq[-(nr+1)]
        }else{
            q1 <- xbrks[nr-1]
            q2 <- xbrks[nr]
            qn1 <- qnorm(rfreq[nr-2])
            qn2 <- qnorm(rfreq[nr-1])
            s <- (q1-q2)/(qn1-qn2)
            mu <- q1 - qn1 * s
            nu <- min(ceiling(freq[nr]*0.01), n.5)
            freq[nr] <- freq[nr] - nu
            qhat <- qnorm(1-nu/ntotal,mu,s)
            xbrks[nr+1] <- qhat
            ##cat("\nnu = ", nl, 'upper limit =',
            ##round(qhat,3), " (", round(exp(qhat),3),")\n")
        }
    }

    ll <- xbrks[-length(xbrks)]
    ul <- xbrks[-1]
    nr <- length(ll)

    res <- .musd.bdata(xbrks,freq)
    p <- rep(1/ng,ng)
    xbar <- runif(ng)*res$mean
    vbar <- rep((res$sd)^2,ng)

    nc <- (ll+ul)/2
    tmp <- rev(sort(freq))
    sele <- which(freq >= tmp[k])
    xbar <- nc[sele]
    p <- freq[sele]/sum(freq[sele])
    vbar <- rep(diff(range(xbrks))/4, ng)
    vbar <- vbar^2
    ## input as WK (workspace size) and output as error info
    wk <- max(500, 3*ng*(nr+2+3*ng)+(nr+2))
    if(nl > 0 || nu > 0){
        itrunc <- 1
    }else{
        itrunc <- 0
    }
    
    res <- .Fortran(.F_emmix,
                    as.integer(freq), as.integer(nr), as.integer(ng),
                    as.double(ll), as.double(ul),
                    p=as.double(p), mu=as.double(xbar), v=as.double(vbar),
                    xlogl=double(0),
                    ifault=as.integer(wk),
                    as.integer(itrunc), 
                    nl=as.integer(nl), nu=as.integer(nu))
    
    xord <- order(res$mu)
    out <- data.frame(p=res$p[xord],
                      mean=res$mu[xord],
                      sd=sqrt(res$v[xord]))
    k <- nrow(out)
    nams <- paste("comp.", 1:k, sep='')
    rownames(out) <- nams
    list(pars=out, breaks=xbrks, freq=freq)
}
