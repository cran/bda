## to compute the distribution of the difference of two finite Gaussian mixtures
## parameters: mu1, mu2,sig1,sig2,p1,p2,x.range

fnm <- function(p1,p2,mu1,mu2,sig1,sig2,from,to){
    if(missing(p1)){
        stop("mixing coefficient(s) missing for treatment arm")
    }else{
        if(any(is.na(p1)))
            stop("missing value(s) in 'p1'")
        if(any(p1<=0|p1>1))
            stop("invalid mixing coefficient(s) in 'p1'")
        if(sum(p1) != 1)
            stop("sum of 'p1' must be 100%")
        k1 <- length(p1)
    }
    
    if(missing(p2)){
        stop("mixing coefficient(s) missing for control arm")
    }else{
        if(any(is.na(p2)))
            stop("missing value(s) in 'p2'")
        if(any(p2<=0|p2>1))
            stop("invalid mixing coefficient(s) in 'p2'")
        if(sum(p2) != 1)
            stop("sum of 'p2' must be 100%")
        k2 <- length(p2)
    }

    if(missing(mu1)){
        stop("mean values missing for treatment arm")
    }else{
        if(any(is.na(mu1)))
            stop("missing value(s) in 'mu1'")
        if(any(!is.finite(mu1)))
            stop("infinite value(s) in 'mu1'")
        if(length(mu1) != k1)
            stop("length of 'mu1' differs from 'p1'")
    }
            
    if(missing(mu2)){
        stop("mean values missing for control arm")
    }else{
        if(any(is.na(mu2)))
            stop("missing value(s) in 'mu2'")
        if(any(!is.finite(mu2)))
            stop("infinite value(s) in 'mu2'")
        if(length(mu2) != k2)
            stop("length of 'mu2' differs from 'p2'")
    }

    if(missing(sig1)){
        stop("SD values missing for treatment arm")
    }else{
        if(any(is.na(sig1)))
            stop("missing value(s) in 'sig1'")
        if(any(!is.finite(sig1)))
            stop("infinite value(s) in 'sig1'")
        if(any(sig1<=0))
            stop("invalid value(s) in 'sig1'")
        if(length(sig1) != k1)
            stop("length of 'sig1' differs from 'p1'")
    }

    if(missing(sig2)){
        stop("SD values missing for control arm")
    }else{
        if(any(is.na(sig2)))
            stop("missing value(s) in 'sig2'")
        if(any(!is.finite(sig2)))
            stop("infinite value(s) in 'sig2'")
        if(any(sig2<=0))
            stop("invalid value(s) in 'sig2'")
        if(length(sig2) != k2)
            stop("length of 'sig2' differs from 'p2'")
    }

    mu.max <- max(c(mu1,mu2))
    mu.min <- min(c(mu1,mu2))
    s.max <- max(c(sig1,sig2))
    x.max <- mu.max + 2*s.max
    x.min <- mu.min - 2*s.max
    
    if(missing(from)) from <- x.min
    if(missing(to)) to <- x.max
    stopifnot(to > from)
    xgrid <- seq(from, to, length=401)
    ##fhat <- .Fortran(.F_lpsmooth,
    ##fx = as.double(x), as.integer(n),
    ##as.double(x),y0=as.double(y), as.integer(n),
    ##bw = as.double(bw1), as.integer(0),
    ##as.double(c(from, to)), as.integer(0),
    ##ellx = double(n), kappa=double(1))

    if(k2==1){
        fx <- dmnorm(xgrid,p=p1,mean=mu1-mu2,sd=sqrt(sig1^2+sig2^2))
        Fx <- pmnorm(xgrid,p=p1,mean=mu1-mu2,sd=sqrt(sig1^2+sig2^2))
    }else{
        fx <- rep(1/diff(range(xgrid)),length(xgrid))
        Fx <- cumsum(fx)
    }
    list(x=xgrid,y=fx,F=Fx)
}


tkde <- function(x,y,group,type="absolute",x.max=10,conf.level = 0.95){
    stopifnot(x.max>0)
    stopifnot(conf.level>0&&conf.level<1)
    if(missing(x))
        stop("'x' missing")
    else{
        if(any(x<=0,na.rm=TRUE))
            stop("'x' must be positive")
        if(any(x > x.max,na.rm=TRUE))
            stop("'x' cannot exceed 'x.max'")
    }

    if(missing(y))
        stop("'y' missing")
    else{
        if(any(y<=0,na.rm=TRUE))
            stop("'y' must be positive")
        if(any(y > x.max,na.rm=TRUE))
            stop("'y' cannot exceed 'x.max'")
    }

    if(missing(group))
        stop("'group' missing")
    else{
        group <- as.factor(group)
        if(nlevels(group) != 2)
            stop("'group' must have two levels")
    }

    if(length(x) != length(y))
        stop("lengths differ for 'x' and 'y'")
    if(length(x) != length(group))
        stop("lengths differ for 'x' and 'group'")

    sele <- is.na(x) | is.na(y) 
    if(any(sele)){
        x <- x[!sele]
        y <- y[!sele]
        group <- group[!sele]
    }

    type <- match.arg(tolower(type),
                      c("absolute","relative","percent","mean",
                        "average"))
    options(warn=-4)
    if(type=="absolute"||type=="mean"||type=="average"){
        d <- x-y #absolute changes
        s <- sd(d,na.rm=TRUE)
        a <- min(d); b <- max(d)
    }else{
        d <- 1-y/x #absolute changes
        s <- sd(d,na.rm=TRUE)
        a <- 1-10/min(x); b <- 1
    }
    a0 <- a-s; b0 <- b+s;
    sele <- group==levels(group)[1]
    n1 <- sum(sele)
    n2 <- sum(!sele)
    fx <- density(d[sele],from=a0,to=b0)
    f1 <- .funRR(fx,a,b)
    fx <- density(d[!sele],from=a0,to=b0)
    f2 <- .funRR(fx,a,b)
    F1 <- cumsum(f1$y); f1$F <-1- F1/max(F1)
    F2 <- cumsum(f2$y); f2$F <-1- F2/max(F2)
    tmp <- cbind(n2,n1,f2$F,f1$F)
    out <- apply(tmp,1,.bincb,conf.level=conf.level)
    tmp <- cbind(n1,f1$F)
    rr1 <- apply(tmp,1,.bincb1,conf.level=conf.level)
    tmp <- cbind(n2,f2$F)
    rr2 <- apply(tmp,1,.bincb1,conf.level=conf.level)
    options(warn=0)
    risk <- list(x=f1$x,y=f2$y/f1$y)
    responder <- list(y=f2$F-f1$F,x=f1$x,ll=out[1,],ul=out[2,],
                      pv=out[3,],rr1=rr1,rr2=rr2)
    list(g1=f1,g2=f2,risk=risk,responder=responder)    
}

.bincb1 <- function(x,conf.level=0.95){
    n <- x[1];
    p <- x[2]
    alpha <- 0.5*(1-conf.level)
    ll <- qbinom(alpha,n,p)/n
    ul <- qbinom(1-alpha,n,p)/n
    c(ll,p,ul)
}

.bincb <- function(x,conf.level=0.95){
    n1 <- x[1];
    n2 <- x[2]
    p1 <- x[3]
    p2 <- x[4]
    tmp <- prop.test(x=c(n1*p1,n2*p2),
                     n=c(n1,n2),
                     conf.level=conf.level)
    c(as.numeric(tmp$conf.int),tmp$p.value)
}

.funRR <- function(fx,a,b){
    x <- fx$x
    y <- fx$y
    sele1 <- x < a
    sele2 <- x > b
    sele <- (!sele1) & (!sele2)
    x0 <- x[sele]
    y0 <- y[sele]
    k <- sum(sele1)
    y0[1:k] <- rev(y[sele1]) + y0[1:k]
    k <- sum(sele2); n <- length(y0)
    y0[(n-k+1):n] <- rev(y[sele2]) + y0[(n-k+1):n]
    list(x=x0,y=y0)
}
