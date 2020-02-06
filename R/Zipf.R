ZipfPlot <- function(x, x0, plot.new=TRUE, fitted=TRUE,weights,...){
    if(class(x)=="bdata"){
        out <- .zipfbin(x$breaks,x$counts,x0,
                        plot.new=plot.new,
                        fitted=fitted,
                        weights=weights, ...)
    }else if(class(x)=="histogram"){
        out <- .zipfbin(x$breaks,x$counts,x0,
                        plot.new=plot.new,
                        fitted=fitted,
                        weights=weights, ...)
    }else if(is.numeric(x)){
        out <- .Zipf(x)
    }else
        stop("data type not supported")
    out
}

ZPlot <- function(x, plot.new=TRUE,...){
    x <- as.numeric(x)
    n <- length(x)
    if(any(x<0)) warning("negative value(s) in 'x'")
    tmp <- table(x)
    Counts <- as.numeric(names(tmp))
    Freq <- as.numeric(tmp)
    x <- log(rev(Counts))
    y <- log(cumsum(rev(Freq))/n)
    if(plot.new){
        plot(x, y, ...)
    }else{
        lines(x,y,...)
    }
}

.zipfbin <- function(xbrks, fn, x0, plot.new=TRUE,
                     weights, ...){
    nclass <- length(fn)
    if(!missing(x0)){
        i <- which(xbrks >= x0)[1]
        if(nclass - i + 1 < 3){
            stop("too few data points")
        }else{
            fn <- fn[2:nclass]
            xbrks <- xbrks[2:(nclass+1)]
        }
    }
    
    N <- sum(fn)
    r <- N - cumsum(fn) + 1
    xr <- xbrks[-1]
    k <- length(fn)
    if(!is.finite(xr[k])){
        xr <- xr[-k]
        r <- r[-k]
    }
    Rank <- r
    X <- xr
    if(missing(weights)){
        wts <- rep(1, length(X))
    }else{
        if(length(weights) != length(X)){
            warning("lower boundaries weighted by counts")
            wts <- fn[-1]
            ##print(X)
            ##print(wts)           
        }else{
            wts <- weights
        }
        if(any(weights < 0))
            stop("negative weight(s) not allowed")
    }
    
    if(any(X<0))
        stop("negative 'x' values not allowed")
    if(plot.new)
        plot(log(X), log(Rank), ...)
    else
        lines(log(X), log(Rank), ...)
    
    lmout <- lm(log(Rank)~log(X),weights=wts)
    out <- summary(lmout) 
    alpha <- -out$coef[2,1]
    b <- out$coef[1,1]
    sea <- out$coef[2,2]
    seb <- out$coef[1,2]
    t0 <- qt(0.975,out$df[2])
    abline(a=b, b=-alpha, ...)
    xm <- exp((b-log(N))/alpha)
    a1 <- alpha - t0*sea
    a2 <- alpha + t0*sea
    b1 <- b-t0*seb
    b2 <- b+t0*seb
    xm1 <- exp((b1-log(N))/alpha)
    xm2 <- exp((b2-log(N))/alpha)
    res <- data.frame(est=c(alpha,xm),
                      ll.95=c(a1,xm1),
                      ul.95=c(a2,xm2)
                      )
    rownames(res) <- c("alpha","xm")

    ##print(cbind(Rank,X))
    k <- length(Rank)
    alpha <- NULL
    ll <- NULL
    ul <- NULL
    radj <- NULL
    M <- NULL
    
    for(i in 1:(k-2)){
        sele <- X >= X[i]
        M <- c(M, sum(sele))
        y2 <- Rank[sele]
        x2 <- X[sele]
        wts2 <- wts[sele]
        lmout <- lm(log(y2)~log(x2),weights=wts2)
        out <- summary(lmout)
        a <- out$coef[2,1]
        alpha <- c(alpha, a)
        sea <- out$coef[2,2]
        t0 <- qt(0.975,out$df[2])
        a1 <- a - t0*sea
        a2 <- a + t0*sea
        ll <- c(ll, a1)
        ul <- c(ul, a2)
        radj <- c(radj, out$adj)
    }
    
    out <- data.frame(
        nclass=M,
        P.index=round(-alpha,3),
        ll.95=round(-ul,3),
        ul.95=round(-ll,3),
        R.adj=round(radj,3)
    )
    
    list(pars=res, best=out,y=log(Rank),x=log(X))
}

.zipf <- function(x, x0,plot.new=TRUE, fitted=TRUE,weights,...){
    if(missing(x0))
        x0 <- 0

    if(missing(weights)){
        wts <- FALSE
    }else{
        wts <- TRUE
    }

    if(x0 < 0)
        stop("'x0' must be non-negative")

    if(x0==0){
        sele <- x > x0
    }else{
        sele <- x >= x0
    }
    y <- x[sele]
    
    ##if(any(y<0)) stop("negative 'x' values not allowed")
    tmp <- table(y)
    N <- length(y)
    Counts <- as.numeric(names(tmp))
    Freq <- as.numeric(tmp)
    X <- rev(Counts)
    Rank <- cumsum(rev(Freq))
    if(plot.new){
        plot(log(X), log(Rank), ...)
        lines(log(X), log(Rank), ...)
    }else{
        points(log(X), log(Rank), ...)
        lines(log(X), log(Rank), ...)
    }
    ##lmout <- lm(log(Rank)~log(X))
    if(wts){
        ri <- rank(y, ties.method='min')
        ri <- length(y) - ri + 1
        lmout <- lm(log(ri)~log(y))
        y2 <- log(ri)
        x1 <- log(y)
    }else{
        lmout <- lm(log(Rank)~log(X))
        y2 <- log(Rank)
        x1 <- log(X)
    }
    x2 <- x1*x1
    x3 <- x2*x1
    lmout2 <- lm(y2~x1+x2+x3)
    x0 <- seq(min(x1),max(x1)+1,length=100)
    y0 <- lmout2$coef[[1]]+lmout2$coef[[2]]*x0 +
        lmout2$coef[[3]]*x0*x0 + lmout2$coef[[4]]*x0*x0*x0
    out <- summary(lmout) 
    ##print(out)
    alpha <- -out$coef[2,1]
    b <- out$coef[1,1]
    sea <- out$coef[2,2]
    seb <- out$coef[1,2]
    t0 <- qt(0.975,out$df[2])
    if(fitted){
        abline(a=b, b=-alpha, ...)
    }
    xm <- exp((b-log(N))/alpha)
    a1 <- alpha - t0*sea
    a2 <- alpha + t0*sea
    b1 <- b-t0*seb
    b2 <- b+t0*seb
    xm1 <- exp((b1-log(N))/alpha)
    xm2 <- exp((b2-log(N))/alpha)
    res <- data.frame(est=c(alpha,xm),
                      ll.95=c(a1,xm1),
                      ul.95=c(a2,xm2)
                      )
    rownames(res) <- c("alpha","xm")
    list(x=x0,y=y0,coef=res,model1=lmout,model2=lmout2)
}

.Zipf <- function(x){
    x <- as.numeric(x)
    tmp <- table(x)
    Counts <- as.numeric(names(tmp))
    Freq <- as.numeric(tmp)
    mydat <- data.frame(Counts=Counts, Freq=Freq)
    log.X <- log(rev(Counts))
    log.R <- log(cumsum(rev(Freq)))
    list(y=log.R, x=log.X)
}

.FindDmax <- function(pars, x, y, cutoff){
    gamma <- pars[1]
    scaler <- pars[2]
    stopifnot(scaler > 0)
    cutoff <- exp(cutoff)

    x0 <- x[x > cutoff]
    y0 <- sort((y/scaler)^gamma, decreasing=TRUE)
    y0 <- y0[1:length(x0)]
    xy.ord <- order(c(x0,y0))
    xy0 <- c(rep(1,length(x0)), rep(-1, length(y0)))
    tmp <- xy0[xy.ord]
    max(abs(cumsum(tmp)))
    
#    x[x<=0] <- 0.25 # used as baseline
#    y[y<=0] <- 0.25
#    x0 <- x[x > cutoff]
#    out <- binning(log(x0), nclass=50)
#    xbrks <- out$breaks
#    xgrid <- c(-Inf,xbrks[-c(1,length(xbrks))], Inf)
#    xt <- binning(log(x), breaks=xgrid)
#    y <- (y/scaler)^gamma
#    yt <- binning(log(y), breaks=xgrid)
#    Fx <- cumsum(xt$counts)/sum(xt$counts)
#    Fy <- cumsum(yt$counts)/sum(yt$counts)
#    max(abs(Fx-Fy))
}

.ZNorm <- function(x, y, cutoff=6, optim=FALSE){
    xnam = deparse(substitute(x))
    ynam = deparse(substitute(y))

    out <- .Zipf(x);lx <- out$x;rx <- out$y
    out <- .Zipf(y);ly <- out$x;ry <- out$y

    sele <- rx < cutoff;
    lmx <- lm(rx[sele]~lx[sele])
    res <- as.numeric(lmx$coef)

    sele <- ry < cutoff;
    sum(sele);mean(sele)
    lmy <- lm(ry[sele]~ly[sele])
    res2 <- as.numeric(lmy$coef)
    res <- rbind(res, res2)

    r2x <- summary(lmx)$r.squared
    r2y <- summary(lmy)$r.squared
    if(r2x < 0.85 || r2y < 0.85)
        warning("'cutoff' might be too small")
       
    ## power transformation
    gamma <- res[2,2]/res[1,2]
    y3 <- y^gamma

    out <- .Zipf(y3);ly <- out$x;ry <- out$y
    sele <- ry < cutoff;
    sum(sele);mean(sele)
    lmy <- lm(ry[sele]~ly[sele])
    res2 <- as.numeric(lmy$coef)
    res <- rbind(res, res2)

    Y <- c(lx,ly);X <- c(rx,ry);
    G <- c(rep("x",length(lx)), rep("y",length(ly)))
    sele <- is.finite(Y)
    Y <- Y[sele]; X <- X[sele];G <- G[sele]
    lmout <- lm(Y~X+G)
    r <- exp(lmout$coef[[3]])
    y2 <- y/r
    out <- .Zipf(y2);ly <- out$x;ry <- out$y

    sele <- ry < cutoff;
    sum(sele);mean(sele)
    lmy <- lm(ry[sele]~ly[sele])
    res2 <- as.numeric(lmy$coef)
    res <- rbind(res, res2)

    res <- as.data.frame(res)
    rownames(res) <- c(xnam, ynam,
                       paste(ynam,".",sep='pwr.transform'),
                       paste(ynam,".scaling",sep='')
                       )
    names(res) <- c("y-intercept","slope")
    gamma <- as.numeric(gamma)
    r <- as.numeric(r)

    scaler.optim <- r
    power.optim <- gamma
    mat.optim <- NULL

    if(optim){
        J <- 30
        Gammas <- seq(0.9*gamma,1.1*gamma, length=J);
        K <- 30
        Alphas <- seq(0.9*r,1.1*r, length=K);
        PARS <- cbind(sort(rep(Gammas,K)), rep(Alphas, J))
        D <- apply(PARS, 1, FUN=.FindDmax, x=x,y=y,cutoff=cutoff)
        Dmax <- matrix(D, nrow=J, ncol=K)
        tmp <- PARS[which(D==min(D))[1],]
        scaler.optim <- tmp[2]
        power.optim <- tmp[1]
        mat.optim <- list(x=Gammas, y=Alphas, z=Dmax)
    }
    y <- (y/scaler.optim)^power.optim
    
    list(x=x, y=y,
         scaler=r, power=gamma,
         scaler.optim=scaler.optim,
         power.optim=power.optim,
         mat.optim=mat.optim,
         coef=res)
}

Zipf.Normalize <- function(x, y, cutoff=6, optim=FALSE,method){
    xnam = deparse(substitute(x))
    ynam = deparse(substitute(y))
    if(missing(method)){
        out <- .ZNorm(x=x,y=y,cutoff=cutoff,optim=optim)
    }else{
        out <- .ZNorm2(x=x,y=y,cutoff=cutoff,optim=optim)
    }
    invisible(out)
}

.ZNorm2 <- function(x, y, cutoff=6, optim=FALSE){
    xnam = deparse(substitute(x))
    ynam = deparse(substitute(y))

    out <- .Zipf(x);lx <- out$x;rx <- out$y
    out <- .Zipf(y);ly <- out$x;ry <- out$y

    sele <- rx < cutoff;
    lmx <- lm(rx[sele]~lx[sele])
    res <- as.numeric(lmx$coef)

    sele <- ry < cutoff;
    sum(sele);mean(sele)
    lmy <- lm(ry[sele]~ly[sele])
    res2 <- as.numeric(lmy$coef)
    res <- rbind(res, res2)

    Y <- c(lx,ly);X <- c(rx,ry);
    G <- c(rep("x",length(lx)), rep("y",length(ly)))
    sele <- is.finite(Y)
    Y <- Y[sele]; X <- X[sele];G <- G[sele]
    lmout <- lm(Y~X+G)
    r <- exp(lmout$coef[[3]])
    y2 <- y/r
    out <- .Zipf(y2);ly <- out$x;ry <- out$y

    sele <- ry < cutoff;
    sum(sele);mean(sele)
    lmy <- lm(ry[sele]~ly[sele])
    res2 <- as.numeric(lmy$coef)
    res <- rbind(res, res2)

    res <- as.data.frame(res)
    rownames(res) <- c(xnam, ynam,
                       paste(ynam,".scaling",sep='')
                       )
    names(res) <- c("y-intercept","slope")
    gamma <- 1
    r <- as.numeric(r)

    scaler.optim <- r
    power.optim <- gamma

    if(optim){
        Alphas <- seq(0.9*r,1.1*r, length=100);
        Gammas <- rep(1,100)
        PARS <- cbind(Gammas, Alphas)
        D <- apply(PARS, 1, FUN=.FindDmax, x=x,y=y,cutoff=cutoff)
        tmp <- Alphas[which(D==min(D))[1]]
        scaler.optim <- tmp
        mat.optim <- list(x=Alphas, y=D)
    }
    y <- y/scaler.optim
    mat.optim <- NULL
    
    list(x=x, y=y,
         scaler=r, power=gamma,
         scaler.optim=scaler.optim,
         power.optim=power.optim,
         mat.optim=mat.optim,
         coef=res)
}
