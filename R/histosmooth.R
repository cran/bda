### binning:  2014/04/17  To manually bin the data

## either 'x' or 'freq' is given, not both/none
## if 'breaks' is not given, 'bw' should be given
## if 'from' is not given, take 'from = min(x)-runif(1)*bw'

binning <- function(x, freq, breaks, bw)
{
    top.coded <- 0; #0=none, -1=left/lower,1=right/upper, 2=both
    if(missing(x)){
        if(missing(freq))
            stop("Both 'x' and 'freq' are missing")
        name <- deparse(substitute(freq))
        nclass <- length(freq)
        stopifnot(is.numeric(freq))
        if(any(is.na(freq)))
            stop("Missing value(s) in 'freq'")
        if(any(freq < 0))
            stop("Negative value(s) in 'freq'")
        tmp <- freq - round(freq)
        if(any(tmp !=0))
            stop("Invalid 'freq' value(s)")
        if(any(!is.finite(freq)))
            stop("Infinite frequeny value(s)")
           
        if(missing(breaks)){
            stop("'breaks' not specified, or stating point missing")
        }else{
            stopifnot(is.numeric(breaks))
            if(length(breaks) == nclass + 1){
                if(any(diff(breaks) <= 0))
                    stop("Invalid break point(s)")
                bw <- diff(breaks)
                if(any(bw != diff(range(breaks))/nclass))
                    equidist <- FALSE
                else{
                    equidist <- TRUE
                }
                if(!is.finite(breaks[1])){
                    bw1 <- 5 * freq[1]*bw[2]/freq[2]
                    breaks[1] <- breaks[2] - bw1
                    top.coded <- -1;
                    equidist <- FALSE
                    bw <- diff(breaks)
                }
                if(!is.finite(breaks[nclass+1])){
                    bw1 <- 5 * freq[nclass]*bw[nclass-1]/freq[nclass-1]
                    breaks[nclass+1] <- breaks[nclass] + bw1
                    top.coded <- ifelse(top.coded==0, 1, 2);
                    equidist <- FALSE
                    bw <- diff(breaks)
                }
                mids <- 0.5*(breaks[-1]+breaks[-(nclass+1)])
            }else{ # breaks not all specified, just use the first one as the starting points
                x0 <- breaks[1]
                stopifnot(is.finite(x0))
                stopifnot(!missing(bw))
                stopifnot(is.numeric(bw))
                stopifnot(is.finite(bw))
                stopifnot(length(bw)==1)
                stopifnot(bw>0)
                breaks <- x0 + c(0:nclass) * bw
                mids <- x0 + (c(1:nclass) - 0.5)*bw
                equidist <- TRUE
            }
        }
        mu <- sum(freq * mids)/sum(freq)
        sig2 <- (sum(mids^2*freq)-sum(freq)*mu^2)/sum(freq)
        sig <- sqrt(sig2)
    }else{
        name <- deparse(substitute(x))
        n <- length(x)
        if(missing(freq))
            freq <- rep(1,n)
        stopifnot(is.numeric(freq))
        if(any(is.na(freq)))
            stop("Missing value(s) in 'freq'")
        if(any(freq < 0))
            stop("Negative value(s) in 'freq'")
        tmp <- freq - round(freq)
        if(any(tmp !=0))
            stop("Invalid 'freq' value(s)")
        if(any(!is.finite(freq)))
            stop("Infinite frequency value(s)")
        if(n != length(freq))
            stop("'x' and 'freq' have difference lengths")
        x.na <- is.na(x)
        if(any(x.na)){
            x <- x[!x.na]; freq <- freq[!x.na]
        }
        if(any(!is.finite(x)))
            stop("Infinite 'x' value(s)")
        
        x <- rep(x,freq)
        mu <- mean(x)
        sig <- sd(x)
        
        if(missing(breaks)){
            stopifnot(!missing(bw))
            stopifnot(is.numeric(bw))
            stopifnot(length(bw)==1)
            stopifnot(bw>0)
            stopifnot(is.finite(bw))
            x0 <- min(x) - runif(1)*bw
            nclass <- ceiling((max(x)-x0)/bw)
            breaks <- x0 + c(0:nclass) * bw
            mids <- x0 + (c(1:nclass) - 0.5)*bw
            equidist <- TRUE
        }else if(length(breaks)==1){
            stopifnot(!missing(bw))
            stopifnot(is.numeric(bw))
            stopifnot(is.finite(bw))
            stopifnot(length(bw)==1)
            stopifnot(bw>0)
            x0 <- breaks
            nclass <- ceiling((max(x)-x0)/bw)
            breaks <- x0 + c(0:nclass) * bw
            mids <- x0 + (c(1:nclass) - 0.5)*bw
            equidist <- TRUE
        }else{
            nclass <- length(breaks) - 1
            if(any(breaks[1]>x)||any(breaks[nclass+1]<x))
                stop("some 'x' not counted; maybe 'breaks' do not span range of 'x'")
            if(any(diff(breaks) <= 0))
                stop("Invalid break point(s)")
            
            bw <- diff(breaks)
            if(any(bw != diff(range(breaks))/nclass))
                equidist <- FALSE
            else{
                equidist <- TRUE
            }
        }
        res <- hist(x, breaks=breaks,plot=FALSE)
        freq <- res$counts
           
        if(!is.finite(breaks[1])){
            bw1 <- 5 * freq[1]*bw[2]/freq[2]
            breaks[1] <- breaks[2] - bw1
            top.coded <- -1;
            equidist <- FALSE
            bw <- diff(breaks)
        }
        
        if(!is.finite(breaks[nclass+1])){
            bw1 <- 5 * freq[nclass]*bw[nclass-1]/freq[nclass-1]
            breaks[nclass+1] <- breaks[nclass] + bw1
            top.coded <- ifelse(top.coded==0, 1, 2);
            equidist <- FALSE
            bw <- diff(breaks)
        }
        mids <- 0.5*(breaks[-1]+breaks[-(nclass+1)])
        
    }

    h <-  (1/(4*pi))^(1/10)*(243/(35*sum(freq)))^(1/5)*sig*15^(1/5)
    

    structure(
        list(breaks = breaks,
             counts = freq,
             density = freq/(sum(freq)*bw),
             mids = mids,
             xname = name,
             equidist=FALSE,
             equalbw = equidist,
             top.coded = top.coded,
             bw = bw, h = max(h,mean(bw)),
             mean = mu, stdev = sig,
             call = match.call()),
        class="histogram")  
}



## by default, split the histogram into k*nclass bins
.split <- function(x,k)
{
    if(missing(k)) k <- 2
    k <- round(k)
    if(k > 3)
        warning("Number of sub-classes 'f' might be too large")
    if(k <= 1){
        out <- x
    }else{
        a <- min(x$breaks)
        b <- max(x$breaks)
        F <- x$counts
        nbin <- length(F)
        M <- k*nbin + 1
        bw <- diff(x$breaks)/k
        bws <- rep(bw, rep(k, nbin))
        breaks <- a + c(0,cumsum(bws))  # good for unequal bins
        
        out <- .bootkde(x)
        fx <- out$y;
        x0 <- out$x
        Fx <- cumsum(fx); Fx <- Fx/max(Fx)
        Fb <- approx(x0, Fx, breaks)$y
        m <- length(Fb)
        if(is.na(Fb[m])) Fb[m] <- 1
        if(is.na(Fb[1])) Fb[1] <- 0
        mx <- diff(Fb)
        counts <- NULL
        for(i in 1:nbin){
            total <- sum(mx[((i-1)*k+1):(i*k)])
            if(F[i]==0 || total ==0){
                counts <- c(counts, rep(0,k))
            }else{
                tmp2 <- 0; tmp <- 0
                for(j in 1:(k-1)){
                    p <- mx[k*(i-1)+j]/total
                    tmp <- round(p*F[i])
                    counts <- c(counts, tmp)
                    tmp2 <- tmp2 + tmp
                }
                counts <- c(counts, max(0,F[i]-tmp2))
            }
        }
        
        nclass <- length(breaks)-1
        mids <- 0.5 * (breaks[-1]+breaks[-(nclass+1)])
        bw <- diff(breaks)
        
        out <- structure(
            list(breaks = breaks,
                 counts = counts,
                 density = counts/(sum(counts)*bw),
                 mids = mids,
                 bw = bw, h=x$h,
                 xname = x$xname,
                 equalbw = x$equalbw,
                 equidist=FALSE,
                 top.coded = x$top.coded,
                 mean = x$mean, stdev = x$stdev,
                 call = x$call),
            class="histogram")
    }
    out
}

##  histospline applies only to 'histogram'

.histospline <- function(x, from, to, gridsize=512L)
    {
        stopifnot(class(x)=='histogram')
        nbin <- length(x$mids)
        stopifnot(gridsize > 10)
        a <- ifelse(missing(from), x$breaks[1] - x$h, from)
        b <- ifelse(missing(to), x$breaks[nbin+1] + x$h, to)


        out <- spline(x$mids, x$density, n=gridsize,
                      xmin = a, xmax = b)
        f0 <- out$y; f0[f0<0] <- 0;
        x0 <- out$x;

        xdiff <- diff(x0)[1]
        mu <- sum(f0*x0)*xdiff
        F0 <- cumsum(f0)
        F0 <- (F0-min(F0))/diff(range(F0))
        med.x <- approx(F0, x0, 0.5)$y
        return(structure(list(y = f0, x = x0,
                              ucb=NULL, lcb=NULL,
                              bw = NULL,
                              mean = mu,
                              median = med.x,
                              data.name = x$xname,
                              call = match.call()),
                         class='smooth'))
    }

.smkde <- function(x, from, to, gridsize=512L)
    {
        if(!x$equalbw)
            stop("Smooth KDE is not applicable for unequal width bins")

        nbin <- length(x$mids)
        stopifnot(gridsize > 10)
        a <- ifelse(missing(from), x$breaks[1] - x$h, from)
        b <- ifelse(missing(to), x$breaks[nbin+1] + x$h, to)
        gpoints <- seq(a, b, length = gridsize)
        ## to approximate the initial estimate of f(x) by simulation
        F <- x$counts; X <- x$mids; n <- sum(F)
        A <- x$breaks[1:nbin]; B <- x$breaks[-1]
        x0 <- runif(n,rep(A,F),rep(B,F))
        f0 <- density(x0, n=gridsize)$y
        bw.2 <- x$bw[1]*0.5
        ## iterate to find a smoothed KDE
        iter <- 1000;
        out <- .Fortran(.F_smoothkde,
                        fx=as.double(f0), as.double(gpoints),
                        as.integer(gridsize), as.double(X), as.double(F),
                        as.integer(nbin), as.double(bw.2),
                        bw = as.double(x$h),iter=as.integer(iter))
        y = out$fx
        if(out$iter>=iter) y = rep(0, gridsize)
        sele1 = is.na(y) | !is.finite(y)
        if(any(sele1)) y[sele1] = 0.0
        tot.mass <- sum(y)*(b-a)/(gridsize-1)

        f0 <- y/tot.mass;
        x0 <- gpoints;

        xdiff <- diff(x0)[1]
        mu <- sum(f0*x0)*xdiff
        F0 <- cumsum(f0)
        F0 <- (F0-min(F0))/diff(range(F0))
        med.x <- approx(F0, x0, 0.5)$y

        return(structure(list(y = f0,
                              x = gpoints,
                              ucb=NULL, lcb=NULL,
                              bw = out$bw,
                              mean = mu,
                              median = med.x,
                              data.name = x$xname,
                              call = match.call()),
                         class='smooth'))
    }

print.smooth <- function (x, digits = NULL, ...)
{
    cat("Call:  ", deparse(x$call), sep = "")

    if(!is.null(x$bw))
        cat(";\tbw =", round(x$bw,3), sep = "")
    
    if(!is.null(x$dist))
        cat("\n  Dist := '", x$dist, sep="")
    if(!is.null(x$pars)){
        cat("'\n\tPara := [", sep = "")
        k <- length(x$pars)
        if(k==1){
            cat(x$pars, "]\n",sep='')
        }else{
            for(i in 1:(k-1)){
                cat(x$pars[i], ", ",sep='')
            }
            cat(x$pars[k], "]\n", sep='')
        }
    }
    if(!is.null(x$mean))
        cat("\nMean = ", round(x$mean,3), sep="")
    if(!is.null(x$median))
        cat(",\tMedian = ", round(x$median,3),"\n", sep='')

    cat("\n")
  invisible(x)
}


plot.smooth <- function (x, col=1, lwd=1, lty=1,
                         shade,border="gray",scb=FALSE,...)
{
    if(length(col)==1){
        col1 <- col; col2 <- col
    }else{
        col1 <- col[1]; col2 <- col[2]
    }
    if(length(lwd)==1){
        lwd1 <- lwd; lwd2 <- lwd
    }else{
        lwd1 <- lwd[1]; lwd2 <- lwd[2]
    }
    if(length(lty)==1){
        lty1 <- lty; lty2 <- lty
    }else{
        lty1 <- lty[1]; lty2 <- lty[2]
    }
    
    plot(x$x, x$y, col=col1, lty=lty1,lwd=lwd1,...)

    if(!is.null(x$ucb)&&!is.null(x$lcb)&&scb){
        if(missing(shade)){
            lines(x$ucb~x$x,col=col2,lty=lty2,lwd=lwd2,...)
            lines(x$lcb~x$x,col=col2,lty=lty2,lwd=lwd2,...)
        }else{
            y0 <- c(x$ucb, rev(x$lcb))
            x0 <- c(x$x, rev(x$x))
            polygon(x0, y0, col=shade, border=border,...)
            lines(x$x, x$y, col=col1,lty=lty1,lwd=lwd1,...)
        }
    }
    
    invisible(x)
}

lines.smooth <- function (x, col=1, lwd=1, lty=1,
                         shade,border="gray",scb=FALSE,...)
{
    if(length(col)==1){
        col1 <- col; col2 <- col
    }else{
        col1 <- col[1]; col2 <- col[2]
    }
    if(length(lwd)==1){
        lwd1 <- lwd; lwd2 <- lwd
    }else{
        lwd1 <- lwd[1]; lwd2 <- lwd[2]
    }
    if(length(lty)==1){
        lty1 <- lty; lty2 <- lty
    }else{
        lty1 <- lty[1]; lty2 <- lty[2]
    }
    
    if(!is.null(x$ucb)&&!is.null(x$lcb)&&scb){
        if(missing(shade)){
            lines(x$ucb~x$x,col=col2,lty=lty2,lwd=lwd2,...)
            lines(x$lcb~x$x,col=col2,lty=lty2,lwd=lwd2,...)
        }else{
            y0 <- c(x$ucb, rev(x$lcb))
            x0 <- c(x$x, rev(x$x))
            polygon(x0, y0, col=shade, border=border,...)
        }
    }
    lines(x$x, x$y, col=col1,lty=lty1,lwd=lwd1,...)
    
    invisible(x)
}


.bootkde <- function(x, from, to, gridsize=512L,
                    method='quantile', conf.level)
{
    stopifnot(class(x)=='histogram')

    method <- match.arg(tolower(method),c("z.score","quantile","cdf"))

    X <- x$mids;
    N <- length(X)
    F <- x$counts;
    B <- diff(x$breaks)/2;
    A <- -B;
    a <- x$breaks[1];
    b <- rev(x$breaks)[1]
    ucb=NULL; lcb=NULL;
    ## Set default bandwidth
    h0 = bw.nrd(rep(X,F)+runif(sum(F),-B,B))
    h = .Fortran(.F_hbmise,as.double(X), as.double(F),as.double(2*B),
        as.integer(N), hopt=as.double(h0))$hopt

    a <- ifelse(missing(from), x$breaks[1] - x$h, from)
    b <- ifelse(missing(to), x$breaks[N+1] + x$h, to)
    gpoints <- seq(a, b, length = gridsize)
    
    y <- .Fortran(.F_ofcpdf,
                  as.double(X), as.double(F),as.double(-B),as.double(B),
                  as.integer(N), y=as.double(gpoints), as.integer(gridsize),
                  para=as.double(h))$y
    y <- cumsum(y)

    if(missing(conf.level)){
        conf.level <- NULL
        y <- .bootemp(sum(F),x=gpoints,y=X,f=F,lb=A,ub=B,Fx=y,from=a,to=b)
    }else{
        stopifnot(is.numeric(conf.level))
        stopifnot(conf.level>0&&conf.level<1)
        alpha <- 1.0 - conf.level
        iter <- 1000
        tmp <- as.matrix(rep(sum(F),iter), ncol=1);
        ## simulate the rounding errors by a pilot estimate of F(x) based on a BME
        out <- apply(tmp,1, .bootemp,x=gpoints,y=X,f=F,lb=A,ub=B,
                     Fx=y,from=from,to=to)
        
        if(method=="z.score"){
            sigs <- apply(out,1,sd);
            z0 <- qnorm(1-alpha/2.);
            ym <- apply(out,1,median);
            lcb <- ym-sigs*z0;
            lcb[lcb<0] <- 0;
            ucb <- ym+sigs*z0;
            ucb[ucb>1] <- 1
            y <- out[,1];
            type <- "Density"
        }else if(method=="quantile"){
            ym <- apply(out,1,median);
            lcb <- as.numeric(apply(out,1,quantile, alpha/2));
            lcb[lcb<0] <- 0;
            ucb <- as.numeric(apply(out,1,quantile, 1-alpha/2));
            ucb[ucb>1] <- 1
            y <- out[,1];
            type <- "Density"
        }else{
            out <- apply(out,2,cumsum)*(b-a)/gridsize
            ym <- apply(out,1,median);
            lcb <- as.numeric(apply(out,1,quantile, alpha/2));
            lcb[lcb<0] <- 0;
            ucb <- as.numeric(apply(out,1,quantile, 1-alpha/2));
            ucb[ucb>1] <- 1
            y <- out[,1];
            type <- "Probability"
        }
    }
    
    f0 <- y
    x0 <- gpoints;
    xdiff <- diff(x0)[1]
    mu <- sum(f0*x0)*xdiff
    F0 <- cumsum(f0)
    F0 <- (F0-min(F0))/diff(range(F0))
    med.x <- approx(F0, x0, 0.5)$y

    return(structure(list(y=y,x=gpoints,
                          mean = mu, median = med.x,
                          ucb=ucb,lcb=lcb, conf.level=conf.level,
                          call = match.call(), data.name = x$xname
                          ), class = "smooth"))
}

.bootemp <- function(n,x,y,f,lb,ub,Fx,from,to){
    M=length(x);u=runif(n);N=length(f);

    smpl = .Fortran(.F_remp, as.integer(N), as.double(y), as.double(f),
        as.double(lb), as.double(ub), as.integer(M), as.double(Fx),
        as.double(x), smpl=as.double(u), as.integer(n))$smpl

    density(smpl,n=M,from=from,to=to)$y
}





##################################################################
### Function:  bde

##  bdde: binned data density estimation
## create an internal R object "bindata" with
## f -> frequencies, a -> lower limits, b -> upper limits

## 2014/04/26: use regression method to estimate the parameters, also
## search the neiborhood of the initial estimates to find the MLEs.
## The method applies only to the top-coded data only.

## 2014/04/29:  add the other histosmooth algorithms to bde


bde <- function(x, breaks,freq, bw, type="weibull",
                from,to,gridsize=512L, conf.level)
    UseMethod("bde")

bde.default <- function(x,breaks,freq,bw,type="weibull",
                        from,to,gridsize=512L, conf.level)
{
    f.call <- match.call()
    xhist <- binning(x=x,freq=freq, breaks=breaks,bw=bw)
    out <- bde(xhist,type=type, from=from,to=to,
               gridsize=gridsize, conf.level=conf.level)
    out$call <- f.call
    out
}

bde.histogram <- function(x,breaks,freq,bw,type="weibull",
                          from,to,gridsize=512L, conf.level)
{
    f.call <- match.call()
    ## support Dagum and Weibull only
    type <- match.arg(tolower(type),
                      c('ewd',
                        'dagum',
                        'weibull',
                        'histospline','spline',
                        'smkde','smoothkde',
                        'fnmm','nmix','normmix','nm',
                        'lpr','npr','root-unroot',
                        'bootkde'))
    out <- switch(type,
                  'ewd' = .bdeEWD(x,from=from,to=to,gridsize=gridsize),
                  'dagum' = .bdeDagum(x,from=from,to=to,gridsize=gridsize),
                  'weibull' = .bdeWeibull(x,from=from,to=to,gridsize=gridsize),
                  'spline' = .histospline(x,from=from,to=to,gridsize=gridsize),
                  'histospline' = .histospline(x,from=from,to=to,
                      gridsize=gridsize),
                  'smkde' = .smkde(x,from=from,to=to,gridsize=gridsize),
                  'smoothkde' = .smkde(x,from=from,to=to,gridsize=gridsize),
                  'lpr' = .histonpr(x,from=from,to=to,gridsize=gridsize,
                      conf.level=conf.level),
                  'npr' = .histonpr(x,from=from,to=to,gridsize=gridsize,
                      conf.level=conf.level),
                  'root-unroot' = .histonpr(x,from=from,to=to,
                      gridsize=gridsize, conf.level=conf.level),
                  'bootkde' = .bootkde(x,from=from,to=to,gridsize=gridsize,
                      conf.level=conf.level),
                  "fnmm" = bnmm(x,from=from,to=to,gridsize=gridsize)
                  )
    out$call <- f.call
    out
}

.bdeWeibull <- function(x.hist, from, to, gridsize=512L){
    breaks <- x.hist$breaks
    freq <- x.hist$counts
    nclass <- length(freq)
    if(breaks[1] < 0)
        stop("Negative value(s) in 'x'")
    Fn <- cumsum(freq)/sum(freq)
    x <- breaks[-1]
    nu <- 0;  #frequency for the upper bin if top.coded
    if(Fn[nclass] >= 1) Fn[nclass] <- 0.99999
    if(x.hist$top.coded ==1 ||x.hist$top.coded == 2){
        Fn <- Fn[-nclass]
        x <- x[-nclass]
        nu <- freq[nclass]
        freq <- freq[-nclass]
    }
    if(Fn[1]<=0){
        Fn <- Fn[-1]
        x <- x[-1]
        freq <- freq[-1]
    }
    
    ## define fine grid point
    stopifnot(gridsize > 10)
    ## if top/bottom coded, compute an upper/lower bound.
    ub <- ifelse(is.finite(breaks[nclass+1]), breaks[nclass+1],
                 2* breaks[nclass] - breaks[nclass-1]);
    if(missing(from)) from <- 0
    if(missing(to)) to <- ub
    gpoints <- seq(from, to, length=gridsize)
    bw <- (to-from)/gridsize
    ##  parameter estimation
    pars <- rep(0,3) # 2 for Weibulll  3 parameters for dagum and ewd
    idist <- 0
    res <- .Fortran(.F_bdregmle,
                    as.double(Fn), as.double(x), as.double(freq),
                    as.integer(nu),
                    as.integer(length(x)), as.integer(idist),
                    pars = as.double(pars))
        
    ##  compute fx(gpoints, pars)
    x0 <- gpoints;
    f0 <- .bdweibull(x0, res$pars[1], res$pars[2])
    f0[is.na(f0)] <- 0
        
    xdiff <- diff(x0)[1]
    mu <- sum(f0*x0)*xdiff
    F0 <- cumsum(f0)
    F0 <- (F0-min(F0))/diff(range(F0))
    med.x <- approx(F0, x0, 0.5)$y
        
    out <- structure(list(y = f0, x = gpoints,
                          pars = res$pars[-3],
                          mean = mu,
                          median = med.x,
                          call = match.call(),
                          data.name = x.hist$xname,
                          type = "Weibull"
                          ),
                     class = "smooth")
    invisible(out)
}

.bdeEWD <- function(x.hist, from, to, gridsize=512L){
    breaks <- x.hist$breaks
    freq <- x.hist$counts
    nclass <- length(freq)
    if(breaks[1] < 0)
        stop("Negative value(s) in 'x'")
    Fn <- cumsum(freq)/sum(freq)
    x <- breaks[-1]
    nu <- 0;  #frequency for the upper bin if top.coded
    if(Fn[nclass] >= 1) Fn[nclass] <- 0.99999
    if(x.hist$top.coded ==1 ||x.hist$top.coded == 2){
        Fn <- Fn[-nclass]
        x <- x[-nclass]
        nu <- freq[nclass]
        freq <- freq[-nclass]
    }
    if(Fn[1]<=0){
        Fn <- Fn[-1]
        x <- x[-1]
        freq <- freq[-1]
    }
    
    ## define fine grid point
    stopifnot(gridsize > 10)
    ## if top/bottom coded, compute an upper/lower bound.
    ub <- ifelse(is.finite(breaks[nclass+1]), breaks[nclass+1],
                 2* breaks[nclass] - breaks[nclass-1]);
    if(missing(from)) from <- 0
    if(missing(to)) to <- ub
    gpoints <- seq(from, to, length=gridsize)
    bw <- (to-from)/gridsize
    ##  parameter estimation
    pars <- rep(0,3) # 2 for Weibulll  3 parameters for dagum and ewd
    idist <- 1
    res <- .Fortran(.F_bdregmle,
                    as.double(Fn), as.double(x), as.double(freq),
                    as.integer(nu),
                    as.integer(length(x)), as.integer(idist),
                    pars = as.double(pars))
        
    ##  compute fx(gpoints, pars)
    x0 <- gpoints;
    f0 <- .bdewd(x0,res$pars[1], res$pars[2],res$pars[3])
    f0[is.na(f0)] <- 0
        
    xdiff <- diff(x0)[1]
    mu <- sum(f0*x0)*xdiff
    F0 <- cumsum(f0)
    F0 <- (F0-min(F0))/diff(range(F0))
    med.x <- approx(F0, x0, 0.5)$y
        
    out <- structure(list(y = f0, x = gpoints,
                          pars = res$pars,
                          mean = mu,
                          median = med.x,
                          call = match.call(),
                          data.name = x.hist$xname,
                          type = "EWD"
                          ),
                     class = "smooth")
    invisible(out)
}

.bdeDagum <- function(x.hist, from, to, gridsize=512L){
    breaks <- x.hist$breaks
    freq <- x.hist$counts
    nclass <- length(freq)
    if(breaks[1] < 0)
        stop("Negative value(s) in 'x'")
    Fn <- cumsum(freq)/sum(freq)
    x <- breaks[-1]
    nu <- 0;  #frequency for the upper bin if top.coded
    if(Fn[nclass] >= 1) Fn[nclass] <- 0.99999
    if(x.hist$top.coded ==1 ||x.hist$top.coded == 2){
        Fn <- Fn[-nclass]
        x <- x[-nclass]
        nu <- freq[nclass]
        freq <- freq[-nclass]
    }
    if(Fn[1]<=0){
        Fn <- Fn[-1]
        x <- x[-1]
        freq <- freq[-1]
    }
    
    ## define fine grid point
    stopifnot(gridsize > 10)
    ## if top/bottom coded, compute an upper/lower bound.
    ub <- ifelse(is.finite(breaks[nclass+1]), breaks[nclass+1],
                 2* breaks[nclass] - breaks[nclass-1]);
    if(missing(from)) from <- 0
    if(missing(to)) to <- ub
    gpoints <- seq(from, to, length=gridsize)
    bw <- (to-from)/gridsize
    ##  parameter estimation
    pars <- rep(0,3) # 2 for Weibulll  3 parameters for dagum and ewd
    idist <- 2
    res <- .Fortran(.F_bdregmle,
                    as.double(Fn), as.double(x), as.double(freq),
                    as.integer(nu),
                    as.integer(length(x)), as.integer(idist),
                    pars = as.double(pars))
        
    ##  compute fx(gpoints, pars)
    x0 <- gpoints;
    f0 <- .bddagum(x0,res$pars[1], res$pars[2],res$pars[3])
    f0[is.na(f0)] <- 0
        
    xdiff <- diff(x0)[1]
    mu <- sum(f0*x0)*xdiff
    F0 <- cumsum(f0)
    F0 <- (F0-min(F0))/diff(range(F0))
    med.x <- approx(F0, x0, 0.5)$y
        
    out <- structure(list(y = f0, x = gpoints,
                          pars = res$pars,
                          mean = mu,
                          median = med.x,
                          call = match.call(),
                          data.name = x.hist$xname,
                          type = "Dagum"
                          ),
                     class = "smooth")
    invisible(out)
}


.bddagum <- function(x, a,b,p){
    a*p/x*(x/b)^(a*p)/((x/b)^a+1)^(p+1)
}

.bdewd <- function(x,kappa,lambda,alpha){
    alpha*kappa/lambda*(x/lambda)^(kappa-1)*(1-exp(-(x/lambda)^kappa))^(alpha-1)*exp(-(x/lambda)^kappa)
}

.bdweibull <- function(x,kappa,lambda){
    kappa/lambda*(x/lambda)^(kappa-1)*exp(-(x/lambda)^kappa)
}
######################################################
######################################################
.betabde <- function(f,a,b,iter,top.coded)
{
  alphas <- rep(1,iter)
  tmp <- apply(as.matrix(alphas,ncol=1),1,.betaplot,
               f=f,a=a,b=b,top.coded=top.coded)
  tmp[-1,which(tmp[1,]==max(tmp[1,]))[[1]]]
}

.betaplot <- function(x,f,a,b,top.coded)
  {
    n <- sum(f); nbin <- length(f)
    x <- runif(n,rep(a,f),rep(b,f))
    mu <- mean(x); s2 <- var(x)
    alpha <- mu*(mu*(1-mu)/s2-1)
    pbeta <- alpha/mu-alpha
    Fa <- pbeta(a,alpha,pbeta)
    if(top.coded) Fa <- c(Fa, 1)
    else Fa <- c(Fa, pbeta(b[nbin],alpha,pbeta))
    llk <- sum(log(diff(Fa)))
    c(llk,alpha, pbeta)
  }

######################################################
######################################################

.dagumbde <- function(f,a,b,iter,top.coded){
  alphas <- c(runif(iter), runif(iter,1,5))
  tmp <- apply(as.matrix(alphas,ncol=1),1,.dagumplot,
               f=f,a=a,b=b)
  tmp[1:3,which(tmp[4,]==max(tmp[4,]))[[1]]]
}

.dagumplot <- function(alpha,f,a,b)
  {
    tol <- 0.000001
    n <- sum(f)
    x <- runif(n,rep(a,f),rep(b,f))
    Fhat <- (rank(x)-0.3)/(n+0.4)
    ly <- log(-1-log(Fhat)/alpha+tol)
    lx <- log(x+tol)
    out <- lm(ly~lx)
    kappa <- -out$coef[[2]]
    lambda <- exp(out$coef[[1]]/kappa)
    llk <- .dagumllk(f,a,b,alpha,kappa,lambda)
    pars <- c(kappa, lambda, alpha,llk)
  }

.dagumllk <- function(f,a,b,alpha,kappa,lambda){
  if(is.na(lambda)||is.na(kappa)){
    llk <- -9999999999999.99
  }else if(lambda<=0||kappa<=0){
    llk <- -9999999999999.99
  }else{
    tol <- 0.0000001
    Fa <- .pdagum(a, kappa,lambda,alpha)
    Fb <- .pdagum(b, kappa,lambda,alpha)
    llk <- sum(f*log(Fb-Fa+tol))
  }
}

.ewdbde <- function(f,a,b,iter,top.coded,fixed){
  if(fixed)
    alphas <- rep(1,iter)
  else
    alphas <- c(runif(iter), runif(iter,1,5))
  tmp <- apply(as.matrix(alphas,ncol=1),1,.ewdplot,
               f=f,a=a,b=b)
  tmp[1:3,which(tmp[4,]==max(tmp[4,]))[[1]]]
}

.ewdplot <- function(alpha,f,a,b)
  {
    tol <- 0.000001
    n <- sum(f)
    x <- runif(n,rep(a,f),rep(b,f))
    Fhat <- (rank(x)-0.3)/(n+0.4)
    ly <- log(-log(1-Fhat^(1/alpha))+tol)
    lx <- log(x+tol)
    out <- lm(ly~lx)
    kappa <- out$coef[[2]]
    lambda <- exp(-out$coef[[1]]/kappa)
    llk <- .ewdllk(f,a,b,alpha,kappa,lambda)
    pars <- c(kappa, lambda, alpha,llk)
  }

.ewdllk <- function(f,a,b,alpha,kappa,lambda){
  if(is.na(lambda)||is.na(kappa)){
    llk <- -9999999999999.99
  }else if(kappa<=0||lambda<=0){
    llk <- -999999999999
  }else{
    tol <- 0.00000000001
    Fa <- .pewd(a, kappa,lambda,alpha)
    Fb <- .pewd(a, kappa,lambda,alpha)
    llk <- sum(f*log(Fb-Fa+tol))
  }
  llk
}


##################################################
##  Weibull(x,lambda,kappa)
##  f(x) = k/l*(x/l)^(k-1)*exp(-(x/l)^k)
##  in R Weibill(x,kappa,lambda)

## for binned data from a Weibull distribution
.bweibullplot <- function(f,a,b,top.coded){
  nbin <- length(f)
  n <- sum(f)
  y <- runif(n,rep(a,f),rep(b,f))
  pars <- .weibullplot(y)
  Fa <- pweibull(a,pars[[2]],pars[[1]])
  if(top.coded) Fa <- c(Fa, 1)
  else Fa <- c(Fa, pweibull(b[nbin],pars[[2]],pars[[1]]))
  llk <- sum(log(diff(Fa)))
  list(pars=pars,llk=llk)
}


## for raw data from a Weibull distribution
.weibullplot <- function(x){
  tol <- 0.0000001
  if(any(x<0))
    stop("Negative value not allowed in Weibull distribution")
  n <- length(x)
  Fhat <- (rank(x)-0.3)/(n+0.4)
  ly <- log(-log(1-Fhat))
  lx <- log(x+tol)
  out <- lm(ly~lx)
  kappa <- out$coef[[2]]
  lambda <- exp(-out$coef[[1]]/kappa)
  c(lambda, kappa)
}

.weibullbde <- function(f,a,b){
  nbin <- length(f)
  if(!is.finite(b[nbin])) b[nbin] <- a[nbin]*2-a[nbin-1]
  n <- sum(f)
  y <- runif(n,rep(a,f),rep(b,f))
  pars <- .weibullplot(y)
  ##.bdmle(f,a,b,dist=0,pars=pars)
}

.bdmle <- function(f,a,b,dist,pars){
  nbin <- length(f); npar <- length(pars)
  .Fortran(.F_BDMLE, as.double(f), as.double(a),
           as.double(b),as.integer(nbin),
           pars = as.double(pars), as.integer(npar),
           as.integer(dist))$pars
}

#############################################################


.pewd <- function(x,kappa,lambda,alpha)
  {
    n <- length(x)
    res <- rep(0,n)
    sele1 <- x > 0
    sele2 <- is.finite(x); 
    sele3 <- is.na(x);
    if(any(sele3)) res[sele3] <- NA
    if(any(!sele3&!sele2)) res[!sele3&!sele2] <- 1.0
    x <- x[sele1&sele2]
    stopifnot(kappa>0&&lambda>0&&alpha>0)
    tmp <- (1-exp(-(x/lambda)^kappa))^alpha
    res[sele1&sele2] <- tmp
    res
  }

.dewd <- function(x,kappa,lambda,alpha)
  {
    n <- length(x)
    res <- rep(0,n)
    sele1 <- x > 0
    sele2 <- is.finite(x); 
    sele3 <- is.na(x);
    if(any(sele3)) res[sele3] <- NA
    x <- x[sele1&sele2]
    stopifnot(kappa>0&&lambda>0&&alpha>0)
    tmp <- alpha * (1-exp(-(x/lambda)^kappa))^(alpha-1) *
      exp(-(x/lambda)^kappa) * kappa * (x/lambda)^(kappa-1)/lambda
    res[sele1&sele2] <- tmp
    res
  }

.pdagum <- function(x,a,b,p)
  {
    n <- length(x)
    res <- rep(0,n)
    sele1 <- x > 0
    sele2 <- is.finite(x); 
    sele3 <- is.na(x);
    if(any(sele3)) res[sele3] <- NA
    if(any(!sele3&!sele2)) res[!sele3&!sele2] <- 1.0
    x <- x[sele1&sele2]
    stopifnot(a>0&&b>0&&p>0)
    tmp <- (1+(x/b)^(-a))^(-p)
    res[sele1&sele2] <- tmp
    res
  }

.ddagum <- function(x,a,b,p)
  {
    n <- length(x)
    res <- rep(0,n)
    sele1 <- x > 0
    sele2 <- is.finite(x); 
    sele3 <- is.na(x);
    if(any(sele3)) res[sele3] <- NA
    x <- x[sele1&sele2]
    stopifnot(a>0&&b>0&&p>0)
    tmp <- (a*p)/x*(x/b)^(a*p)*((x/b)^a+1)^(-p-1)
    res[sele1&sele2] <- tmp
    res
  }


.histonpr <- function(x,from,to,gridsize,conf.level=0.95)
    {
        if(!x$equalbw)
            stop("Smooth KDE is not applicable for unequal width bins")
        if(missing(conf.level)) conf.level <- 0.95
        
        xhist <- x
        mcounts <- mean(x$counts)
        if(mcounts > 50)
            x <- .split(x,k=2)
        else if(mcounts > 100)
            x <- .split(x, k=3)
        
        F <- x$counts;
        X <- x$mids;
        nbin <- length(X);
        n <- sum(F)
        a <- x$breaks[1]; b <- x$breaks[nbin+1];

        ##        if(missing(bw)){
        lscv <- 1
        A <- x$breaks[1:nbin]; B <- x$breaks[-1]
        x0 <- runif(n,rep(A,F),rep(B,F))
        bw <- bw.nrd(x0)
        adaptive = 1  # use adaptive bandwidth selector
        ##        }else{
        ##            lscv <- 0
        ##            adaptive = 0  # don't use adaptive bandwidth selector
        ##        }

        if(missing(from)) from <- a - 3*bw
        if(missing(to)) to <- b + 3*bw
        stopifnot(to > from)

        if(missing(gridsize)) gridsize <- 512L
        stopifnot(gridsize > 10)
        ##  root-transformation
        Y = sqrt(nbin/n*(F+0.25))
        ##print(x$counts)
        ## call npr(...) to estimate the nonparametric regression function
        out <- npr(y=Y, x=X, bw=bw, sd.y= 0.5*sqrt(nbin/n),
                   from=from, to=to, gridsize=gridsize,
                   conf.level=conf.level)
      
        y <- out$y; x <- out$x
        ll <- out$lcb; ul <- out$ucb;
        y <- y*y; tot.mass <- sum(y)*(to-from)/gridsize;
        ll[ll<0] <- 0 # density can not be negative
        ll <- ll*ll/tot.mass
        ul <- ul*ul/tot.mass
        
        f0 <- y/tot.mass
        x0 <- x;
        xdiff <- diff(x0)[1]
        mu <- sum(f0*x0)*xdiff
        F0 <- cumsum(f0)
        F0 <- (F0-min(F0))/diff(range(F0))
        med.x <- approx(F0, x0, 0.5)$y

        structure(list(y = y/tot.mass, x = x, cv=out$cv,
                       conf.level = out$conf.level,
                       pars = c(out$bw, out$kappa),
                       mean = mu, median = med.x,
                       ucb=ul, lcb=ll,
                       call = match.call()),
                  class = 'smooth')
    }


