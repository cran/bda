## 2018/01/11: Fit NGS data using raw data (without binning, but we
## will handle the zeros and extreme high expressions differently).

## to fit a mixed normal mixture with three components (j=1,2,3) and
## (j=0) to a single gene expression profile.  Return the fitted
## components.  Also, we return the weights w(i,j) for gene i.

EM.NGS <- function(x, mixed=TRUE){
    tol <- .Machine$double.eps
    if(any(is.na(x)))
        stop("missing value(s) found")
    if(any(x<0))
        stop("negative value(s) not allowed")
    if(mixed){
        p0 <- mean(x==0) # component 0
        y <- log(x[x > 0]) # log-transformation
    }else{
        xmin <- min(x[x>0])
        y <- log(x + 0.25*xmin)
        p0 <- 0
    }
    ## exclude a 5% of the extremely large data in fitting the model
    n <- length(y)
    m <- round(n * 0.95)
    y <- sort(y)[1:m]
    ## initialize the estimates.  In case some values have very high
    ## frequencies, we bin the data to group the data into three
    ## classes.
    tmp <- NULL
    l <- hist(y, plot=FALSE,nclass=50)$counts
    ncum <- cumsum(l)
    n1 <- ncum[which(ncum>m/3)[1]]
    n2 <- ncum[which(ncum>2*m/3)[1]]
    p1 <- n1/m; p2 <- (n2-n1)/m; p3 <- 1-p1-p2
    mu1 <- mean(y[1:n1]);s1 <- sd(y[1:n1])
    mu2 <- mean(y[(n1+1):n2]);s2 <- sd(y[(n1+1):n2])
    mu3 <- mean(y[(n2+1):m]);s3 <- sd(y[(n2+1):m])
    
    pars <- c(mu1,mu2,mu3,s1,s2,s3,p1,p2)
    res <- .Fortran(.F_em3,
                    iter=as.integer(length(y)), 
                    as.double(y),
                    pars=as.double(pars),
                    as.double(tol))
    tmp <- 1- sum(res$pars[7:8])
    mu1 <- res$pars[1]
    mu2 <- res$pars[2]
    mu3 <- res$pars[3]
    s1 <- res$pars[4]
    s2 <- res$pars[5]
    s3 <- res$pars[6]
    p1 <- res$pars[7]
    p2 <- res$pars[8]
    p3 <- 1-p1-p2

    x0 <- seq(min(y), max(y),length=401)
    f1 <- dnorm(x0,mu1,s1)
    f2 <- dnorm(x0,mu2,s2)
    f3 <- dnorm(x0,mu3,s3)
    f <- p1*f1+p2*f2+p3*f3
    
    out <- data.frame(
        mean = res$pars[1:3],
        sd   = res$pars[4:6],
        p    = c(res$pars[7:8],tmp)
    )
    structure(list(x=x0, y=f, p0=p0,
                   iter=res$iter,
                   fit=out),
              class='NGS.MM')
}

print.NGS.MM <- function(x,...){
    if(x$iter<50000)
        cat("\nIteration:", x$iter)
    else
        cat("\nConverge failed")
    cat("\nComponent of zeros:", x$p0,"\n")
    print(x$fit)
}

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
    m <- length(unique(x))
    n.first <- min(m, n.first)
    
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
    ybreaks <- yhist$breaks
    tmp.breaks <- c(xq[1]-abs(xq[1]*.1),xq[2:k],ybreaks[-1])    
    xhist <- hist(x, breaks=tmp.breaks, plot=FALSE)
    xhist
}

## We assume all data are NGS data and without any transformation such
## as log- or glog-transformation.  We will transform the data within
## the fit.NGS function.

## 2018/01/21: if the same size is too small or when there are two few
## distinct values, we fit a lognormal distribution.

fit.NGS <- function(x, iter.max=30,reflect=TRUE){
    LogNormal <- TRUE

    if(class(x) == 'histogram'){
        xhist0 <- x;
        p0 <- 0;
        sLog <- -1
        p.signal <- NA
        x <- structure(list(xhist=xhist0, 
                            sLog = sLog,
                            p0 = p0),
                       class='NGS.hist')
        
    }else if(class(x) == 'NGS.hist'){
        xhist0 <- x$xhist
        p0 <- x$p0;
        sLog <- x$sLog
        n <- sum(xhist0$counts)
        n0 <- round(n*p0)
        p.signal <- pbinom(n-n0,n,0.5, lower.tail=FALSE)
    }else{
        n <- length(x)
        n0 <- sum(x==0,na.rm=TRUE)
        p.signal <- pbinom(n-n0,n,0.5, lower.tail=FALSE)
        tmp <- unique(x)
        if(length(tmp) > 30 || p.signal < 0.05){
            LogNormal <- FALSE
            x <- bin.NGS(x, sLog=0)
            xhist0 <- x$xhist
            p0 <- x$p0;
            sLog <- x$sLog
        }else{
            LogNormal <- TRUE
            p0 <- mean(x==0,na.rm=TRUE)
        }
    }

    if(LogNormal){
        out <- fit.mlnorm(x, method='lognormal',k=1)
        xhist <- hist(x, plot=FALSE)

        fit0 <- structure(list(xhist=xhist, 
                               ng=1,p=1, mu=out$meanlog,
                               s=out$sdlog,
                               llk = NULL, nl = 0, nu = 0,
                               iter = NA, ifault = NA,
                               "AIC"=NA, "BIC"=NA,"AICc"=NA,
                               lognormal = FALSE, c.glog=0,
                               x.range = range(x), N = NA,K=NA,
                               call = match.call(),
                               trunc = FALSE),
                          class='nmix')
        
        res <- structure(list(xhist=xhist, 
                              fitted.model=fit0,
                              fnmm = NULL,
                              p0 = p0, p.hetero = 1.0,
                              p.signal = p.signal,
                              sLog = 0),
                         class='NGS.nmix')
    }else{
        ## reflection is optional.  We force reflection=FALSE if p0 is
        ## small enough.
        
        if(p0 < 0.25) reflect <- FALSE
        
        ## To stablize the mixture model estimation, we
        ## repeat the procedue and carefully choose initial values.
        
        if(reflect)
            xhist0$breaks[1] <- 2*xhist0$breaks[1]-xhist0$breaks[2]
        
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
        
        res <- structure(list(xhist=x, 
                              fitted.model=fit0,
                              fnmm = fnmm, p0 = p0,
                              p.hetero = p.hetero,
                              p.signal = p.signal,
                              sLog = sLog),
                         class='NGS.nmix')
    }
    res
}

.fitMM.NGS <- function(xhist, iter.max,k=5){
    res <- NULL
    llk <- NULL;N <- NULL;K <- NULL
    if(missing(iter.max)) iter.max <- 10
    if(iter.max<10) iter.max <- 10
    for(i in 1:iter.max){
        tmp <- fit.mixnorm(xhist, k=k)
        x.ord <- order(tmp$mu)
        tmp$mu <- tmp$mu[x.ord]
        tmp$s <- tmp$s[x.ord]
        tmp$p <- tmp$p[x.ord]
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
    if(is.null(x$fnmm)){
        a <- x$fit$x.range[1]
        b <- x$fit$x.range[2]+1
        x0 <- seq(a, b, length=401)
        y0 <- dlnorm(x0, x$fit$mu, x$fit$s)
        out <- list(x=x0, y=y0)
    }else{
        out <- pdf.NGS.nmix(x)
    }
    
    if(hist){
        plot(x$xhist,main=main,freq=FALSE,...)
        lines(out)
    }else{
        plot(out, main=main,type='l',...)
    }
}

lines.NGS.nmix <- function(x,...){
    if(is.null(x$fnmm)){
        a <- x$fit$x.range[1]
        b <- x$fit$x.range[2]+1
        x0 <- seq(a, b, length=401)
        y0 <- dlnorm(x0, x$fit$mu, x$fit$s)
        out <- list(x=x0, y=y0)
    }else{
        out <- pdf.NGS.nmix(x)
    }
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
            fx <- fit.mlnorm(xtmp)
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
    D <- max(abs(as.numeric(Fy)-as.numeric(Fx)))
}



### DEG test

.perm.test.NGS <- function(x,y,size){
    n <- length(x)
    sele <- sample(1:n, size)
    x2 <- x[sele]; y2 <- x[-sele]
    nx <- size; ny <- n-nx
    lmts.x <- y[sele]; lmts.y <- y[-sele]
    fit.x <- fit.mlnorm(x2)
    fit.y <- fit.mlnorm(y2)
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
    fit.x <- fit.mlnorm(Control)
    fit.y <- fit.mlnorm(Treat)
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
    tmp <- .permtest0(x=x,y=y,iter=iter)
    pv2 <- tmp$pv; pv3 <- NA;
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

