## Created in June 2016 for the UCAS talk Updated on 06/16/2017 --
## rewrite the density estimation function and define a new R object,
## add plot function.

## bin the NGS data.  The smallest (0 counts) are also included.  The
## last class is not top-coded here.

##2017/10/14: to guarantee the smoothness, we need to double check and
##make sure there is no gaps in the first n.first+1 classes.  This
##won't happen because if we use single-valued classes, not
##zero-frequency class will exist.  The only concern is that if we use
##sLog(x+c) and c is very small

bin.NGS <- function(x, sLog, n.first=10){
    x <- x[!is.na(x)]
    if(any(x<0)) stop("negative value not allowed")
    if(missing(sLog)){
        sLog <- 0
    }else{
        stopifnot(is.numeric(sLog))
        if(sLog < 0)
            stop("invalid constant for sLog-transormation")
        if(!is.finite(sLog))
            stop("constant not finite for sLog-transormation")
    }
    
    if(sLog>0){
        x <- log(x + sLog)
        p0 <- 0
    }else{
        p0 <- mean(x == 0)
        x <- log(x[x>0])
    }
    
    res <- .binNGS(x, n.first=n.first)
    if(p0>0)
        res$density <- res$density * (1-p0)
    structure(list(xhist=res, 
                   sLog = sLog,
                   p0 = p0),
              class='NGS.hist')
}

print.NGS.hist <- function(x,...){
    cat("Percent of zeros =", round(100*x$p0,2),"%\n")
    cat("Breaks:\n")
    print(x$xhist$breaks)
    cat("Bin counts:\n")
    print(x$xhist$count)
    cat("\n")
}

plot.NGS.hist <- function(x,ylim=NULL,...){
    if(is.null(ylim))
        ylim <- c(0, max(x$p0,x$xhist$density))
    plot(x$xhist,ylim=ylim,...)
    if(x$p0 > 0){
        x0 <- x$xhist$breaks[1]
        segments(x0,x$p0,x0,0,lwd=2,col='red')
    }
}

.binNGS <- function(x, n.first=10){
    ## to keep some details at the left side of the distribution, we
    ## define some classes using the distrinct log(x) values.
    xq <- sort(unique(x))
    nxq <- length(xq)
    if(nxq < 5) stop("Too few distinct data points")
    if(nxq < n.first){ #use 20% if too few unique values
        k <- round(nxq * 0.2)
        if(k<1) k <- 2
    }else{
        k <- n.first
    }
    y <- x[x>xq[k]]
    yhist <- hist(y, plot=FALSE)
    ##    print(yhist)
    ybreaks <- yhist$breaks
    tmp.breaks <- c(xq[1]-abs(xq[1]*.1),xq[2:k],ybreaks[-1])    
    ##tmp <- xq[2:(k+1)]*.9999
    ##tmp.breaks <- c(xq[1],tmp)
    ##xbreaks <- c(tmp.breaks, ybreaks[-1])
    xhist <- hist(x, breaks=tmp.breaks, plot=FALSE)
    xhist
}

## We assume all data are NGS data and without any transformation such
## as log- or glog-transformation.  We will transform the data within
## the fit.NGS function. 
fit.NGS <- function(x, iter.max=30,reflect=TRUE){
    if(class(x) == 'histogram'){
        xhist0 <- x;
        p0 <- 0;
        sLog <- -1
        x <- structure(list(xhist=xhist0, 
                            sLog = sLog,
                            p0 = p0),
                       class='NGS.hist')
        
    }else if(class(x) == 'NGS.hist'){
        xhist0 <- x$xhist
        p0 <- x$p0;
        sLog <- x$sLog
    }else{
        x <- bin.NGS(x, sLog=0)
        xhist0 <- x$xhist
        p0 <- x$p0;
        sLog <- x$sLog
    }
    
    ## To stablize the mixture model estimation, we repeat the
    ## procedue and carefully choose initial values.
    if(reflect)
        xhist0$breaks[1] <- 2*xhist0$breaks[1]-xhist0$breaks[2]

#    if(x$sLog==0 && x$p0 > 0 && !mixed)
#        xhist0$breaks[1] <- -Inf

    maxLLK <- NULL; fit0 <- NULL
    fnmm <- list() #to store different fitted results with k components
    LLK1 <- NULL
    for(k in 1:7){
        res <- .fitMM.NGS(xhist0,iter.max,k=k)
        llk0 <- res$llk
        if(p0>0){
            n0 <- round(sum(xhist0$counts)/(1-p0)*p0)
            llk0 <- llk0 + n0 * log(p0)
        }
        K <- res$K
        N <- res$N
        AIC = -2.0 * llk0 + 2.0 * K;
        AICc = AIC + 2.0 * K * (K+1.) / (N - K - 1.0);
        BIC = -2.0 * llk0 + log(N*2.0*pi) * K;
        tmp <- list(fitted=res$fit,AIC=AIC,BIC=BIC,AICc=AICc)
        fnmm[[k]] <- tmp
        ## we choose the model that minimize AIC/BIC/AICc
        ##print(AICc)
        ##print(maxLLK)
        if(length(AICc) == 0) AICc <- NA
        if(!is.na(AICc)){
            if(is.null(maxLLK)){
                np0 <- K
                LLK0 <- llk0;
                maxLLK <- AICc
                LLK1 <- llk0;
                np1 <- K;
                fit0 <- res$fit
            }else if(maxLLK > AICc){
                maxLLK <- AICc
                LLK1 <- llk0;
                np1 <- K
                fit0 <- res$fit
            }
        }
    }

    ## a likelihood-ratio test for homogeneity
    tmp <- -2 * (LLK0 - LLK1)
    ##    cat("\nLLK1=",LLK1, "LLK0=",maxLLK,"df=",np-np1)
    p.hetero <- pchisq(tmp, np1-np0, lower.tail=FALSE)
    if(np0==np1) p.hetero <- 1.0

    structure(list(xhist=x, 
                   fitted.model=fit0,
                   fnmm = fnmm,
                   p0 = p0, p.hetero = p.hetero,
                   sLog = sLog),
              class='NGS.nmix')
}

.fitMM.NGS <- function(xhist, iter.max,k=5){
    res <- NULL
    llk <- NULL;N <- NULL;K <- NULL
    if(missing(iter.max)) iter.max <- 10
    if(iter.max<10) iter.max <- 10
    for(i in 1:iter.max){
        tmp <- fit.mixnorm(xhist, k=k)
        tmp2 <- cdf.nmix(tmp, xhist$breaks)
        llk2 <- sum(log(diff(tmp2$y)) * xhist$counts)

        ## it doesn't matter we use AIC, BIC or ... we choose the
        ## model that maximize llk
        if(!is.na(llk2)){
            if(is.null(llk)){
                llk <- llk2
                res <- tmp
                N <- tmp$N
                K <- tmp$K
            }else if(llk < llk2){
                llk <- llk2
                res <- tmp
                N <- tmp$N
                K <- tmp$K
            }
        }
    }
    list(fit=res,llk=llk,N=N,K=K)
}

print.NGS.nmix <- function(x,...){
    cat("\nWeight of component zero: ", x$p0,"\n")
    print(x$fitted.model)
}

plot.NGS.nmix <- function(x,hist=TRUE,main=NULL,type=NULL,...){
    out <- pdf.NGS.nmix(x)
    
    if(hist){
        plot(x$xhist,main=main,...)
        lines(out)
    }else{
        plot(out, main=main,type='l',...)
    }
}

lines.NGS.nmix <- function(x,...){
    out <- pdf.NGS.nmix(x)
    lines(out,...)
}

summary.NGS.nmix <- function(object,...){
    x <- object
    res <- NULL
    n <- length(x$fnmm)
    for(k in 1:n){
        out <- x$fnmm[[k]]
        tmp <- c(out$AIC,out$BIC,out$AICc)
        res <- rbind(res,tmp)
    }
    res <- as.data.frame(res)
    names(res) <- c("AIC","BIC","AICc")
    row.names(res) <- paste("k=",1:n,sep='')
    res
}

pdf.NGS.nmix <- function(x, x0, from, to, gridsize=512){
    stopifnot(class(x)=='NGS.nmix')
    out <- pdf.nmix(x$fitted.model, x0=x0, from=from,
                    to=to,gridsize=gridsize)
    xmin <- x$xhist$xhist$breaks[1]
    x0 <- out$x
    y0 <- out$y

    if(x$sLog!=0 || x$p0==0){
        sele <- x0 < xmin
        if(any(sele)){
            n <- sum(sele)
            x0 <- x0[!sele]
            y1 <- y0[1:(n+1)]
            y0 <- y0[!sele]
            y0[1:(n+1)] <- y0[1:(n+1)] + rev(y1)
        }
    }
    if(x$p0 > 0){
        y0 <- (1-x$p0)*y0
    }
    
    list(x = x0, y = y0)
}

cdf.NGS.nmix <- function(x, x0, from, to, gridsize=512){
    stopifnot(class(x)=='NGS.nmix')
    out <- cdf.nmix(x$fitted.model, x0=x0, from=from,
                    to=to,gridsize=gridsize)
    xmin <- x$xhist$xhist$breaks[1]
    x0 <- out$x
    y0 <- out$y
    
    sele <- x0 < xmin
    if(any(sele)){
        x0 <- x0[!sele]
        y0 <- y0[!sele]
    }
    if(x$p0 > 0){
        y0 <- x$p0 + (1-x$p0)*y0 
    }
    list(x = x0, y = y0)
}

.rescale.mm <- function(fit1, fit2,plot=FALSE){
    ## we just need the 2nd and 3rd components
    ## fit2 is used as reference and we want to find a rescaling
    ## parameter mc such that the L1-distances between m+fit1 and fit2
    ## be the smallest.  m will change the two means, not the SD.  We
    ## use a numerical method to do this.
    stopifnot(class(fit1)=='nmix')
    stopifnot(class(fit2)=='nmix')
    ox <- order(fit1$mu)
    mu1 <- fit1$mu[ox]
    s1 <- fit1$s[ox]
    p1 <- fit1$p[ox]
    ox <- order(fit2$mu)
    mu2 <- fit2$mu[ox]
    s2 <- fit2$s[ox]
    p2 <- fit2$p[ox]
    ## determine the range
    ll <- min(mu1[2]-3*s1[2], mu2[2]-3*s2[2])
    ul <- max(mu1[3]+3*s1[3], mu2[3]-3*s2[3])
    x0 <- seq(ll, ul, length=401)
    xdiff <- diff(x0)[1]
    ## choose a range to search for m

    f2 <- p2[2]*dnorm(x0,mu2[2],s2[2]) +
        p2[3]*dnorm(x0,mu2[3],s2[3])
    if(plot){
        m <- 0
        f1 <- p1[2]*dnorm(x0,mu1[2]+m,s1[2]) +
            p1[3]*dnorm(x0,mu1[3]+m,s1[3])
        plot(f1~x0, type='l')
        lines(f2~x0, col='red')
    }

    ms <- range(mu2[-1]-mu1[-1])
    k <- 50
    ms <- seq(min(0.5,ms[1]*.85), max(ms[2]*2,1.5), length=k)
    l1s <- rep(0, k)
    for(i in 1:k){
        m <- ms[i]
        f1 <- p1[2]*dnorm(x0,mu1[2]+m,s1[2]) +
            p1[3]*dnorm(x0,mu1[3]+m,s1[3])
        l1s[i] <- sum(abs(f2-f1)) * xdiff
    }
    n.opt <- which(l1s == min(l1s))[1]
    mopt <- ms[n.opt]
    if(plot){
        f1 <- p1[2]*dnorm(x0,mu1[2]+mopt,s1[2]) +
            p1[3]*dnorm(x0,mu1[3]+mopt,s1[3])
        lines(f1~x0, col='blue',lwd=2)
        plot(l1s~ms, type='l', xlab="rescaling parameter",
             ylab="L1-distance")
    }
    list(opt=mopt,y=l1s,x=ms)
}

.GetID <- function(x, mmobj){
    if(missing(mmobj))
        mmobj <- fit.NGS(x)
    .findmax <- function(x)
        which(x==max(x))[1]
    out <- rep(0, length(x))
    sele <- x > 0 #cannot be negative
    lx <- log(x)
    mu <- mmobj$fitted.model$mu
    s <- mmobj$fitted.model$s
    xord <- order(mu)
    mu <- mu[xord]
    s <- s[xord]
    p1 <- dnorm(lx,mean=mu[1],sd=s[1])
    p2 <- dnorm(lx,mean=mu[2],sd=s[2])
    p3 <- dnorm(lx,mean=mu[3],sd=s[3])
    tmp1 <- cbind(p1,p2,p3)
    tmp2 <- apply(tmp1, 1, .findmax)
    out[sele] <- tmp2[sele]
    p3 <- pnorm(lx,mean=mu[3],sd=s[3])
    sele <- p3 > 0.95
    if(any(sele)) out[sele] <- 5
    out
}

.qnmatch <- function(x1, res1, res0){
    Fx0 <- cdf.NGS.nmix(res0, gridsize=2001)
    x0 <- Fx0$x
    y0 <- Fx0$y
    p0 <- res0$p0
    mu <- res0$fit$mu
    s <- res0$fit$s
    P <- res0$fit$p
    ox <- order(mu)
    mu <- rev(mu[ox])[1]
    s <- rev(s[ox])[1]
    P <- P[ox]
    
    sele <- x1>0
    Fx1 <- cdf.NGS.nmix(res1, log(x1)[sele])
    y1 <- Fx1$y
    p1 <- res1$p0
    .findx <- function(y,y0,x0,p0,p){
        if(y<=p0){
            res <- 0
        }else if(y>max(y0)){
            tmp <- (y+p[3]-1)/p[3]
            if(tmp>=1.0){
                tmp <- 0.9999999
            }
            res <- qnorm(tmp, mean=mu, sd=s)
        }else{
            tmp <- abs(y-y0)
            xid <- which(tmp == min(tmp))[1]
            res <- as.numeric(x0[xid])
        }
        res
    }

    tmp <- matrix(y1, ncol=1)
    out <- apply(tmp,1,.findx,y0=y0,x0=x0,p0=p0,p=P)
    x1[sele] <- exp(out)
    x1
}

normalize.NGS <- function(x,y, method="mixture"){
    if(class(x) != 'data.frame')
        stop("'x' must be a data.frame")
    nx <- ncol(x)
    if(missing(y)){
        xy.split <- FALSE
        xy <- x
    }else{
        if(class(y) != 'data.frame')
            stop("'y' must be a data.frame")
        if(nrow(x) != nrow(y))
            stop("'x' and 'y' have different number of rows")
        if(any(row.names(x) != row.names(y)))
            stop("'x' and 'y' have different genes")
        xy.split <- TRUE
        xy <- cbind(x,y) # merge the data
    }
    
    method = match.arg(tolower(method),
                       c("mixture", "quantile","fnmm",'nmix','qn'))

    xy.rank <- apply(xy,2,rank,ties.method="min")
    
    .ind2mean <- function(my_index, my_mean){
        return(my_mean[my_index])
    }
    
    if(method=='quantile' ||method=='qn'){
        xmin <- min(c(as.matrix(xy[xy>0])))
        xy.t <- log(xy+xmin*0.25)
        xy.sorted <- data.frame(apply(xy.t, 2, sort))
        xy.mean <- apply(xy.sorted, 1, mean)
        xy.mean <- exp(xy.mean)
    }else{
        xbar <- NULL;
        lmts <- as.numeric(apply(xy,2,.find.array.limit))
        for(i in 1:nrow(xy)){
            ##            cat("\n record #:", i)
            xtmp <- as.numeric(xy[i,]);
            fx <- fit.lognormal(xtmp,lmts)
            xbar <- c(xbar, fx$Mean)
        }
        xy.mean <- xbar
    }

    xy.final <- apply(xy.rank, 2, .ind2mean, my_mean=sort(xy.mean))
    xy.final <- as.data.frame(xy.final)
    row.names(xy.final) <- row.names(x)
    if(xy.split){
        xt <- xy.final[,c(1:nx)]
        yt <- xy.final[,-c(1:nx)]
        xt[x==0] <- 0
        yt[y==0] <- 0
    }else{
        xt <- xy.final
        xt[x==0] <- 0
        yt <- NULL
    }
    out <- list(x=xt,y=yt,mean.profile=xy.mean)
    invisible(out)
}

.kdediff <- function(x,y){
    px <- mean(x==0)
    py <- mean(y==0)
    lx <- log(x[x>0])
    ly <- log(y[y>0])
    xy <- c(lx,ly); s <- sd(xy)
    xmin <- min(xy)-s; xmax <- max(xy)+s
    fx <- density(lx, from=xmin,to=xmax)
    fy <- density(ly, from=xmin,to=xmax)
    xdiff <- diff(fx$x)[1]
    sum(abs((1-py)*fy$y-(1-px)*fx$y))*xdiff+abs(px-py)
}

.ecdf <- function(x,y){
    mean(y <= x, na.rm=TRUE)
}

.ecdfdiff <- function(x,y){
    xy <- c(x,y)
    x0 <- seq(min(xy), max(xy), length=401)
    Fx <- lapply(x0, FUN=.ecdf,y=x)
    Fy <- lapply(x0, FUN=.ecdf,y=y)
    D <- abs(max(as.numeric(Fy)-as.numeric(Fx)))
}

    
perm.test.NGS <- function(x,y, alternative = "two.sided", iter = 1001){
    xnam <- deparse(substitute(x))
    ynam <- deparse(substitute(y))
    nam = paste(xnam, "(", length(x), ") vs ",
                ynam, "(", length(y), ")")
    if (!is.numeric(x) && !is.complex(x) && !is.logical(x)) {
        warning("'x' is not numeric or logical: returning NA")
        return(NA_real_)
    }
    x <- x[!is.na(x)]
    if (!is.numeric(y) && !is.complex(y) && !is.logical(y)) {
        warning("'y' is not numeric or logical: returning NA")
        return(NA_real_)
    }
    y <- y[!is.na(y)]
    
    alternative = match.arg(tolower(alternative),
                            c("one.sided", "two.sided"))
    fun = .ecdfdiff
    D = fun(x, y)
    rfun <- function(x, y) {
        nx = length(x)
        ny = length(y)
        n = nx + ny
        n0 = sample(n)
        xy = c(x, y)
        fun(xy[n0[1:nx]], xy[n0[(nx + 1):n]])
    }
    
    bar <- function(n, x, y)
        replicate(n, rfun(x = x, y = y))

    z = bar(iter, x = x, y = y)
    p1 <- mean(abs(D) <= abs(z))
    p2 <- mean(D <= z)
    p2 <- min(p2, 1-p2)
    pv = switch(alternative, two.sided = p1, 
                one.sided = p2, stop("Wrong test type!"))

    RVAL <- list(statistic = c(D = D),
                 p.value = pv, p1=p1, p2=p2,
                 method = "Permutation test for NGS data", 
                 data.name = nam)
    class(RVAL) <- "htest"
    return(RVAL)
}

###  Fitting Copula 09/29/2017 ##################################
## Creates the grid counts from a bivariate data set X
## over an equally-spaced set of grid points
## contained in "gpoints" using the linear
## binning strategy. Note that the FORTRAN subroutine
## "lbtwod" is called.
.bin2D <- function(X, gpoints1, gpoints2)
{
    n <- nrow(X)
    X <- c(X[, 1L], X[, 2L])
    M1 <- length(gpoints1)-1
    M2 <- length(gpoints2)-1
    g1 <- gpoints1
    g2 <- gpoints2
    out <- .Fortran(.F_bintwod, as.double(X), as.integer(n),
                    as.double(g1[-1]), as.double(g2[-1]), 
                    as.integer(M1), as.integer(M2),
                    M = double(M1*M2))
    matrix(out$M, M1, M2)
}

.binshrink <- function(xg,xc,k=5){
    while(length(xc) > k){
        l <- length(xc)
        m <- which(xc == min(xc))[1L]
        if(m == 1){
            xc[2] <- xc[2] + xc[m];
            xg <- xg[-(m+1)]
        }else if(m == l){
            xc[m-1] <- xc[m-1] + xc[m]
            xg <- xg[-m]
        }else{
            if(xc[m-1] > xc[m+1]){
                xc[m+1] <- xc[m+1] + xc[m]
                xg <- xg[-(m+1)]
            }else{
                xc[m-1] <- xc[m-1] + xc[m]
                xg <- xg[-m]
            }
        }
        xc <- xc[-m]
    }
    list(grid=xg, counts=xc)
}

.mixnorm2d <- function(Fx,Fy, Psi){
    stopifnot(Psi > 0)
    Sxy <- 1 + (Fx+Fy)*(Psi-1)
    if(Psi==1){
        Hxy <- Fx * Fy
    }else{
        Hxy <- 0.5 * (Sxy-sqrt(Sxy^2-4*Psi*(Psi-1)*Fx*Fy))/(Psi-1);
    }    
    Hxy
}

fit.nmix.copula <- function(x,y,mle.large=FALSE){
    ## must turn mixed off otherwise will have trouble to compute the
    ## copula -- f(0,0) needs to be redefined.
    N <- length(x)
    stopifnot(length(y)==N)
    xfit <- fit.NGS(x)
    yfit <- fit.NGS(y)
    
    nbinx <- length(xfit$xhist$xhist$counts)
    nbiny <- length(yfit$xhist$xhist$counts)
    xgrid <- xfit$xhist$xhist$breaks
    ygrid <- yfit$xhist$xhist$breaks
    xcount <- xfit$xhist$xhist$counts
    ycount <- yfit$xhist$xhist$counts
    ## transform raw data such that the zero's will be grouped
    ## correctly.
    lx <- log(x); ly <- log(y)
    xmin <- min(c(x[x>0],y[y>0]))*.25
    if(xfit$p0>0){
        lx[x==0] <- log(xmin)
    }
    if(yfit$p0>0){
        ly[y==0] <- log(xmin)
    }

    ngroup <- 5
    mat2d.large <- .bin2D(cbind(lx,ly), xgrid, ygrid)
    xgrid2 <- .binshrink(xgrid,xcount,k=ngroup)$grid
    ygrid2 <- .binshrink(ygrid,ycount,k=ngroup)$grid
    mat2d.small <- .bin2D(cbind(lx,ly), xgrid2, ygrid2)

    ## rough estimates of Psi
    nr <- nrow(mat2d.small)
    nc <- ncol(mat2d.small)
    psis <- NULL
    for(i in 1:(nr-1)){
        a <- sum(mat2d.small[1:i,1:i])
        b <- sum(mat2d.small[1:i,(i+1):nc])
        c <- sum(mat2d.small[(i+1):nr,1:i])
        d <- sum(mat2d.small[(i+1):nr,(i+1):nc])
        tmp <- a * c
        if(tmp==0){
            out <- 1000
        }else{
            out <- a*d/tmp
        }
        psis <- c(psis, out)
    }

    ##  Compute and save fx, fy, Fx and Fy 
    Fx <- rep(0, ngroup+1)
    Fy <- rep(0, ngroup+1)
    for(i in 1:(ngroup+1)){
        Fx[i] <- pmixnorm(xgrid2[i],xfit$fit$p,xfit$fit$mu,xfit$fit$s)
        Fy[i] <- pmixnorm(ygrid2[i],yfit$fit$p,yfit$fit$mu,yfit$fit$s)
    }

    ## grid search for Psi 
    psis1 <- seq(0.1*max(0.01,min(psis)), 10*max(psis), length=100)
    psis2 <- seq(11*max(psis), 1000*max(psis), length=200)
    psis <- c(psis1, psis2)
    Chi2 <- NULL
    LLK <- NULL
    for(psi in psis){
        ecounts <- matrix(0,nrow=ngroup,ncol=ngroup)
        for(i in 1:ngroup){
            for(j in 1:ngroup){
                a <- Fx[i]
                b <- Fx[i+1]
                c <- Fy[j]
                d <- Fy[j+1]
                Hbd <- .mixnorm2d(b,d,psi)
                Had <- .mixnorm2d(a,d,psi)
                Hbc <- .mixnorm2d(b,c,psi)
                Hac <- .mixnorm2d(a,c,psi)
                ecounts[i,j] <- Hbd - Had - Hbc + Hac
            }
        }
        llk <- sum(mat2d.small*log(ecounts), na.rm=TRUE)
        ecounts <- N*ecounts
        
        chistat <- sum((mat2d.small-ecounts)^2/ecounts)
        Chi2 <- c(Chi2, chistat)
        LLK <- c(LLK, llk)
    }
    df0 <- ngroup^2 - 3*(xfit$fit$ng+yfit$fit$ng) - 2
    if(df0<1) df0 <- 1

    sele1 <- which(Chi2 == min(Chi2))[1]
    sele2 <- which(LLK == max(LLK))[1]
    ## Compute MLE using the large matrix.  This need to be tuned off.
    ## Otherwise, we need to compute using Fortran or C code for large
    ## data analysis.
    if(mle.large){
        ##  Compute and save fx, fy, Fx and Fy
        ng1 <- nrow(mat2d.large)
        ng2 <- ncol(mat2d.large)
        Fx <- rep(0, ng1+1)
        Fy <- rep(0, ng2+1)
        for(i in 1:(ng1+1)){
            Fx[i] <- pmixnorm(xgrid[i],xfit$fit$p,xfit$fit$mu,xfit$fit$s)
        }
        for(i in 1:(ng2+1)){
            Fy[i] <- pmixnorm(ygrid[i],yfit$fit$p,yfit$fit$mu,yfit$fit$s)
        }
        
        LLK2 <- NULL
        for(psi in psis){
            ecounts <- matrix(0,nrow=ng1,ncol=ng2)
            for(i in 1:ng1){
                for(j in 1:ng2){
                    a <- Fx[i]
                    b <- Fx[i+1]
                    c <- Fy[j]
                    d <- Fy[j+1]
                    Hbd <- .mixnorm2d(b,d,psi)
                    Had <- .mixnorm2d(a,d,psi)
                    Hbc <- .mixnorm2d(b,c,psi)
                    Hac <- .mixnorm2d(a,c,psi)
                    ecounts[i,j] <- Hbd - Had - Hbc + Hac
                }
            }
            llk <- sum(mat2d.large*log(ecounts), na.rm=TRUE)
            LLK2 <- c(LLK2, llk)
        }
        sele3 <- which(LLK2 == max(LLK2))[1]
        psi.mle.large <- psis[sele3]
    }else{
        psi.mle.large <- NULL
        LLK2 <- NULL
    }
    
    structure(list(xfit=xfit, yfit=yfit,
                   Psi.chisq = psis[sele1],
                   Psi.mle.small = psis[sele2],
                   Psi.mle.large = psi.mle.large,
                   p.value = 1-pchisq(Chi2[sele1], df0),
                   Psis = psis,
                   ChiSq = Chi2,
                   df=df0,
                   LLK.small=LLK,
                   LLK.large=LLK2,
                   mat2large=mat2d.large,
                   mat2small=mat2d.small),
              class='NGS.nmix.copula')
}

print.NGS.nmix.copula <- function(x,...){
    cat("\nFinite normal mixture model for 'x'\n")
    print(x$xfit)
    cat("\nFinite normal mixture model for 'y'\n")
    print(x$yfit)
    cat("\nPsi estimates:\n")
    cat("\tChiSq est=",x$Psi.chisq, "\n\tMLE=",x$Psi.mle.small,"\n")
}

plot.NGS.nmix.copula <- function(x,...){
    persp(x$mat2large)
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

    
.lnorm <- function(parms, x) {
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

fit.lognormal <- function(x,x.limits){
    xnam = deparse(substitute(x))
    x <- as.numeric(x)
    nx <- length(x)
    nx0 <- sum(x==0)
    if(any(x<0)) stop("negative 'x' found")
    ## initialize the output variables
    x0 <- seq(min(x),max(x), length=401)
    ## unadjusted mean and sd of sample
    mu0 <- mean(x); s0 <- sd(x)

    if(missing(x.limits)){
        x.lmts <- x.limits <- rep(0, nx)
    }else{
        if(any(x.limits < 0)) stop("non-positive 'x.limits' found")
        stopifnot(nx == length(x.limits))
        sele <- x == 0
        p0 <- mean(sele)
        x.lmts <- x.limits
        x.lmts[sele] <- x.lmts[sele] * (1-p0)
    }
        
    ##adjusted mean and sd of sample
    yhat <- x + 0.5*x.lmts
    y.mean <- mean(yhat)
    y.median <- median(yhat)
    mu1 <- y.mean; s1 <- sd(yhat)
    
    ## perform an exact test using binomial distribution to test
    ## whether the positive measures are random (by chance) as results
    ## of measurement errors and others.
    pv <- pbinom(nx0,size=nx, prob=0.5, lower.tail = FALSE)

    ## for the zero measures (after adjusted for the measurement
    ## errors due to sensitivity of the platform/technology), we add
    ## random positive noise and obtain MLE for lognormal
    ## distribution.  If a gene is noweakly expressed, we approximate
    ## using such MLEs.  If a gene is positively measured, we use such
    ## MLEs as initial values to fit lognormal distribution based on
    ## pre-binned data.
    sele <- yhat <= 0
    n0 <- sum(sele)
    if(n0 > 0){
        yhat[sele] <- runif(n0, 0, 0.1) # 0.1 is chosen arbitrarily
    }
    ly <- log(yhat)
    mu.mle <- mean(ly)
    sd.mle <- mean((ly-mu.mle)^2)
    ## initialize the estimates
    mu.est <- mu.mle
    sd.est <- sd.mle
    
    ## we setup a cutoff of 0.001.  If pv<0.001, we conclude the gene
    ## expression is weak (undetectable, not significantly different
    ## from zero, ...) and we simply estimate the lognoraml parameters
    ## using the MMEs.

    if(pv >= 0.005/nx){
        tmp <- optim(c(mu.mle, sd.mle), .lnorm, NULL,
                     method="L-BFGS-B",
                     lower=c(0.1*mu.mle,0.1*sd.mle),
                     upper=c(10*mu.mle,10*sd.mle),
                     x = x)
        mu.est <- tmp$par[1];
        sd.est <- tmp$par[2];
    }

    ## goodness-of-fit test using chisq-test
    xbreaks <- qlnorm(c(0.2,0.4,0.6,0.8),
                      meanlog=mu.est, sdlog=sd.est)
    xbreaks <- c(0,xbreaks,Inf)
    xhist <- hist(x, breaks=xbreaks,plot=FALSE)
    x1 <- xhist$counts
    y1 <- rep(round(nx*0.2),5)
    p.gof <- chisq.test(cbind(x1,y1))$p.value

    y0 <- dlnorm(x0,meanlog=mu.est,sdlog=sd.est,log=FALSE)
    MU <- exp(mu.est+0.5*sd.est^2)
    SD <- sqrt((exp(sd.est^2)-1)*exp(2*mu.est+sd.est^2))

    structure(list(
        data = x,
        x.name = xnam,
        x=x0, y=y0, size=c(nx,nx0),
        meanlog = mu.est, sdlog = sd.est,
        p.value = pv, p.gof = p.gof,
        xbar=mu0, s=s0,
        Mean = MU, SD = SD,
        xbar.adj=mu1, s.adj=s1),
        class='NGS.lognormal')
}

plot.NGS.lognormal <- function(x,...){
    hist(x$data, nclass=30, prob=TRUE)
    lines(x, ...)
}

print.NGS.lognormal <- function(x,...){
    cat("\nFitting lognormal distribution to data '",x$x.name,"'",sep='')
    cat("\n\nSummaries:\n")
    tmp1 <- c(x$Mean,x$SD)
    tmp2 <- c(x$xbar,x$s)
    tmp3 <- c(x$xbar.adj,x$s.adj)
    tmp <- data.frame(MLE=tmp1, raw=tmp2, adj=tmp3)
    row.names(tmp) <- c("Mean","Std.Dev")
    print(tmp, digits=3)
    cat("\n n.obs =", x$size[1], ", n.zero =", x$size[2])        
    cat("\n Parameter estimates: meanlog =",
        round(x$meanlog,3), ", sdlog =",
        round(x$sdlog,3), "\n GOF-test, p-value =",
        round(x$p.gof,4))
    cat("\n mean>0 (detectable expression level), p-value =",
        round(x$p.value,4),"\n")
}


### DEG test

.perm.test.NGS <- function(x,y,size){
    n <- length(x)
    sele <- sample(1:n, size)
    x2 <- x[sele]; y2 <- x[-sele]
    nx <- size; ny <- n-nx
    lmts.x <- y[sele]; lmts.y <- y[-sele]
    fit.x <- fit.lognormal(x2, lmts.x)
    fit.y <- fit.lognormal(y2, lmts.y)
    ## t.test
    dmu <- fit.x$meanlog-fit.y$meanlog
    se <- sqrt((fit.x$sdlog)^2/nx+(fit.y$sdlog)^2/ny)
    d.f. <- nx+ny-2
    tstat <- dmu/se
}

.find.array.limit <- function(x) min(x[x>0])

deg.NGS <- function(x,y, gene,iter=1001){
    stopifnot(is.numeric(iter))
    if(iter<100) iter <- 100
    gnames <- row.names(x)
    tmp <- row.names(y)
    if(!any(gnames == tmp))
        stop("gene names are different in the two data files")
    xp.names <- names(x)
    yp.names <- names(y)
    
    xnam = deparse(substitute(x))
    ynam = deparse(substitute(y))
    
    x1 <- as.matrix(x)
    y1 <- as.matrix(y)
    gsize <- nrow(x1)  # number of genes
    nx <- ncol(x1)
    ny <- ncol(y1)
    ## choose the gene to be tested.
    if(missing(gene)){
        sele <- 1
    }else{
        if(is.numeric(gene)){
            sele <- round(gene)
            stopifnot(sele>0 && sele<gsize)
        }else{
            tmp <- tolower(as.character(gene))
            sele <- match(tmp, tolower(gnames))
            if(is.na(sele))
                stop("gene not found")
        }
    }
    gene.name <- gnames[sele]
    ## the smalest non-zero measurement for each profile
    lmts.x <- as.numeric(apply(x1,2,.find.array.limit))
    lmts.y <- as.numeric(apply(y1,2,.find.array.limit))

    ## now we impute the log-expression by fitting a log-normal
    ## distribution based on binned data (so we don't need to
    ## log-transform the zero measures). For each gene, we do the
    ## following: fit log-normal distributions for the control data
    ## and treatment data, then save, (1) size.x, (2) meanlog.x, (3)
    ## sdlog.x, (4) size.y, (5) meanlog.y, (6) sdlog.y, (7)
    ## normal.test, (8) permutation test.
    x2 <- x1[sele,]
    y2 <- y1[sele,]
    p0.x <- mean(x2==0) #percent of zeros
    p0.y <- mean(y2==0) #percent of zeros
    if(length(unique(x2))>30){
        xfit <- fit.NGS(x2)
    }else{
        xfit <- NULL
    }
    if(length(unique(y2))>30){
        yfit <- fit.NGS(y2)
    }else{
        yfit <- NULL
    }

    Control <- x2
    Treat <- y2
    fit.x <- fit.lognormal(Control, lmts.x)
    fit.y <- fit.lognormal(Treat, lmts.y)
    ## t.test
    dmu <- fit.x$meanlog-fit.y$meanlog
    se <- sqrt((fit.x$sdlog)^2/nx+(fit.y$sdlog)^2/ny)
    d.f. <- nx+ny-2
    tstat <- dmu/se
    pv1 <- pt(-abs(tstat),df=d.f.)*2.
    ## perform a permutation test
    x <- c(x2,y2)
    y <- c(lmts.x, lmts.y)
    ##    tmp <- replicate(iter, .perm.test.NGS(x=x,y=y,size=nx))
    ##    pv2 <- mean(abs(tmp)>abs(tstat))
    ##    pv3 <- mean(tstat < tmp)
    ##    pv3 <- min(pv3, 1-pv3)
    tmp <- perm.test.NGS(x=x,y=y,iter=iter)
    pv2 <- tmp$p1; pv3 <- tmp$p2;
    structure(list(gene=gene.name,
                   x.name = xnam, y.name = ynam,
                   fit.x=fit.x, fit.y=fit.y,
                   xfit=xfit,yfit=yfit,
                   p.t.test=pv1, t.stat=tstat, t.df=d.f.,
                   p.perm.test=pv2,
                   p.perm.test.1t=pv3),
              class='NGS.test')
}

print.NGS.test <- function(x,...){
    cat("\nGene name:", x$gene,"\n")
    print(x$fit.x)
    if(!is.null(x$xfit)){
        cat(" Homogeneity test, p-value =",
            round(x$xfit$p.hetero,4), "\n")
    }
    print(x$fit.y)
    if(!is.null(x$xfit)){
        cat(" Homogeneity test, p-value =",
            round(x$yfit$p.hetero,4),"\n")
    }
    
    cat("\nDEG test results:\n")
    if(x$fit.x$p.value >= 0.0001){
        if(x$fit.y$p.value >= 0.0001){
            cat("expressed in both", x$x.name,
                "and", x$y.name,"\n")
            tmp <- c('t.test'= x$p.t.test,
                     'perm (2-sided)' = x$p.perm.test,
                     'perm (1-sided)'= x$p.perm.test.1t)
            print(tmp, digits=4)
        }else{
            cat("expressed in", x$x.name,
                "but not expressed in", x$y.name,
                "(undetectable)\n")
        }
    }else{
        if(x$fit.y$p.value >= 0.0001){
            cat("expressed in", x$y.name,
                "but not expressed in", x$x.name,
                "(undetectable)\n")
        }else{
            cat("not expressed in both", x$x.name,
                "and", x$y.name,"\n")
        }
    }
}

plot.NGS.test <- function(x,main=NULL,...){
    if(is.null(main)) main <- x$gene
    if(!is.null(x$xfit)){
        plot(x$xfit,main=main,...)
        if(!is.null(x$yfit)){
            lines(x$yfit, lty=2,...)
        }
    }else{
        if(!is.null(x$yfit)){
            plot(x$yfit,main=main,...)
        }else{
            warning("check raw data -- all zero measures?")
        }
    }
}
