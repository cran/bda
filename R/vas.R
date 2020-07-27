## VAS & MEM #########################################################
## non-parametric regression for weighted data -- an alternative way
## to handle measurement errors. Remarks: (1) a single valued 'bw' is
## used to avoid heavy computational burdens. (2) LSCV is used to
## select the optimal bandwidth.

.vasbin2 <- function(x,y){
    x <- as.matrix(x)
    y <- as.numeric(y)
    
    nx <- apply(!is.na(x),1,sum)
    wx <- matrix(1/nx, nrow=nrow(x),ncol=ncol(x))
    y0 <- matrix(y,nrow=length(y), ncol=ncol(x))
    
    wx[is.na(x)] <- 0
    x <- as.numeric(x)
    w <- as.numeric(wx)
    y1 <- as.numeric(y0)
    
    sele <- w > 0
    x <- x[sele]
    w <- w[sele]
    y <- y1[sele]
    list(x=x,w=w,y=y)
}

VAS.npr <- function(y,x,bw,s.x,from,to,gridsize,conf.level=0.95)
{
    ## y: dependent variable
    ## x: predictor
    ## wlpsmooth if x is a matrix
    ## lpsmooth if x is a vector

    if(missing(s.x)){
        s.x <- 0
    }else{
        stopifnot(is.numeric(s.x))
        stopifnot(length(s.x)==1)
        stopifnot(s.x>=0)
    }

    y <- as.matrix(y)
    x <- as.matrix(x)
    stopifnot(nrow(x)==nrow(y))
    stopifnot(ncol(y)==1)
    if(ncol(x)==1){
        y <- as.numeric(y)
        x <- as.numeric(x)
        sele <- is.na(x) | is.na(y)
        if(any(sele)){
            if(sum(!sele)<5)
                stop("too many missing values")
            x <- x[!sele]
            y <- y[!sele]            
        }
        out <- lpsmooth(y=y,x=x,bw=bw,
                        from=from,to=to,gridsize=gridsize,
                        conf.level=conf.level)
    }else{
        res <- .vasbin2(x,y)
        out <- wlpsmooth(y=res$y,x=res$x,w=res$w,s.x=s.x,bw=bw,
                         from=from,to=to,gridsize=gridsize,
                         conf.level=conf.level)
    }
    out
}

## VAS.pdf estimate the density function using nonparametric
## regression (local polynomial with df=1) via binning. 'x' will be
## binned automatically. An existing R function lpsommoth will be
## called with the following settings: (1) lscv=TRUE (to refine the
## bandwidth estimate); (2) adaptive=FALSE (enhance the speed).
VAS.pdf <- function(x,y,bw,s.x, gridsize = 401L,type="absolute",
                    range.x, truncate = FALSE){
    type <- match.arg(tolower(type),
                      c("absolute","percent","relative"))
    if(missing(s.x)){
        s.x <- 0
    }else{
        stopifnot(is.numeric(s.x))
        stopifnot(length(s.x)==1)
        stopifnot(s.x>0)
    }
    
    wx <- .vasbin(x=x,y=y,type=type)
    x <- wx$x
    w <- wx$w
    n <- sum(w)
    w <- w/n
    mu <- sum(w*x)
    vx <- sum(w*(x-mu)^2)*n/(n-1)
    ## density estimation with RR
    if(wx$ME){
        if(s.x > 0){
            out <- wdekde(x=wx$x,w=wx$w,bandwidth=bw,s.x=s.x)
        }else{
            out <- wkde(x=wx$x,w=wx$w,bandwidth=bw)
        }
    }else{
        out <- density(wx$x)
    }
    
    x0 <- out$x; x1 <- x0
    f0 <- out$y; f1 <- f0
    if(!missing(range.x)){
        stopifnot(is.numeric(range.x))
        stopifnot(length(range.x)==2)
        a <- range.x[[1]]
        b <- range.x[[2]]
        ## apply RR to left end
        sele <- x0 < a
        if(any(sele)){
            k <- sum(sele)
            x1[sele] <- NA
            f1[(k+1):(2*k)] <- f1[(k+1):(2*k)]+rev(f1[sele])
            f1[sele] <- 0
        }
        ## apply RR to right end
        sele <- x0 > b
        if(any(sele)){
            k <- sum(sele)
            n <- length(x1)
            x1[sele] <- NA
            f1[(n-2*k+1):(n-k)] <- f1[(n-2*k+1):(n-k)]+rev(f1[sele])
            f1[sele] <- 0
        }
        sele <- f1 > 0
        x1 <- x1[sele]
        f1 <- f1[sele]
    }
    
    list(x=x1,y=f1,x0=x0,y0=f0,
         data=wx$x, wt=wx$w,
         mean=mu,s=sqrt(vx),bw=out$bw)
}

.pc <- function(x,y) (x-y)/x
.vaschange <- function(x,nx,ny,type){
    x0 <- x[1:nx]
    x1 <- x[(nx+1):(nx+ny)]
    if(type=="absolute"){
        out <- as.numeric(outer(x0,x1,"-"))
    }else{
        out <- as.numeric(outer(x0,x1,FUN=.pc))
    }
    out
}

.vasbin <- function(x,y,type){
    x <- as.matrix(x)
    if(ncol(x)==1) ME <- FALSE
    else ME <- TRUE
    if(missing(y)){
        nx <- apply(!is.na(x),1,sum)
        wx <- matrix(1/nx, nrow=nrow(x),ncol=ncol(x))
        wx[is.na(x)] <- 0
        x <- as.numeric(x)
        w <- as.numeric(wx)
        sele <- w > 0
        x <- x[sele]
        w <- w[sele]
    }else{
        y <- as.matrix(y)
        if(ncol(x)==1&&ncol(y)==1) ME <- FALSE
        else ME <- TRUE
        
        stopifnot(nrow(x)==nrow(y))
        nx <- ncol(x); ny <- ncol(y)
        xy <- cbind(x,y)
        x <- apply(xy, 1, .vaschange,nx=nx,ny=ny,type=type)
        x <- as.matrix(x)
        nx <- apply(!is.na(x),1,sum)
        wx <- matrix(1/nx, nrow=nrow(x),ncol=ncol(x))
        wx[is.na(x)] <- 0
        x <- as.numeric(x)
        w <- as.numeric(wx)
        sele <- w > 0
        x <- x[sele]
        w <- w[sele]
    }
    list(x=x,w=w,ME=ME)
}

.vasdenest <- function(x,nclass=20,bw){
    ## 'x' can be a matrix if repeated measures are available
    if(any(is.na(x))){
        x <- x[!is.na(x)]
        warning("missing value(s) removed")
    }
    stopifnot(is.numeric(nclass))
    stopifnot(length(nclass)==1)
    nclass <- round(nclass)
    stopifnot(nclass>5)
    k <- length(x)/3
    if(nclass > k) nclass <- k
    
    a <- min(x); b <- max(x)
    ## determine the number of classes
    breaks <- seq(a,b,length=nclass)
    xhist <- hist(x, breaks=breaks,plot=FALSE)
    y0 <- sqrt(0.25+xhist$counts)
    d <- diff(xhist$breaks)
    x0 <- xhist$breaks[-1] - d/2
    x1 <- (x0 - a)/(b-a)
    out <- lpsmooth(y=y0, x=x1, sd.x=0,lscv=FALSE,bw=bw)
    fhat <- (out$y)^2
    x0 <- out$x
    d <- diff(x0)
    fh <- (fhat[-1]+fhat[-length(fhat)])/2
    fTotal <- sum(fh*d)
    f0 <- fhat/fTotal/(b-a)
    x0 <- x0  * (b-a) + a
    list(y=f0, x=x0, pars=out$pars)
}

#.rootunroot <- function(x,breaks,conf.level=0.95){
#    name <- deparse(substitute(x))
#    if(any(is.na(x)))
#        stop("missing value(s) not allowed in 'x'")
#    if(any(x<0))
#        stop("negative value(s) not allowed in 'x'")
#    n <- length(x)
#    if(length(breaks)-n != 1)
#        stop("invalid length of 'breaks'")
#    if(any(diff(breaks)<=0))
#        stop("invalid values in 'breaks'")    
#}


VAS.ecdf <- function(x,w,alpha=0.05){
    if(missing(w)){
        if(any(is.na(x)))
            x <- x[!is.na(x)]
        stopifnot(alpha>0 && alpha<1)
        n <- length(x)
        xord <- sort(x)
        delta <- 0.05*diff(range(xord))
        xord <- c(xord[1]-delta,xord)
        Fn <- (0:n)/n
        en <- sqrt(0.5 * log(2/alpha))/n
        Ux <- Fn + en
        Ux[Ux>1] <- 1
        Lx <- Fn - en
        Lx[Lx<0] <- 0
    }else{
        n <- length(x)
        stopifnot(length(w)==n)
        stopifnot(all(w>=0))
        sele <- is.na(x) | is.na(w)
        if(any(sele)){
            x <- x[!sele]
            w <- w[!sele]
            n <- length(x)
        }
        n0 <- sum(w)
        xord <- order(x)
        x.ord <- x[xord]
        w.ord <- w[xord]
        delta <- 0.05*diff(range(x.ord))
        Fn <- cumsum(w.ord)/n0
        xord <- c(x.ord[1]-delta, x.ord)
        Fn <- c(0,Fn)
        en <- sqrt(0.5 * log(2/alpha))/n0
        Ux <- Fn + en
        Ux[Ux>1] <- 1
        Lx <- Fn - en
        Lx[Lx<0] <- 0        
    }
    
    structure(list(x=xord,
                   y=Fn,
                   lb=Lx,
                   ub=Ux,
                   alpha=alpha,
                   data=x,
                   ecdf=TRUE),
              class="VAS")
}

print.VAS <- function(x,...){
    print(summary(x$data))
}

.drawECDF <- function(x,y,ecdf=TRUE,...){
    if(ecdf){
        k <- length(x)
        xn <- 2*x[k] - x[k-1]
        x <- sort(c(x[1],x[-1],x[-1],xn))
        y <- sort(c(y,y))
    }
    lines(x,y,...)
}

plot.VAS <- function(x,cb=FALSE,type='n',...){
    type <- 'n'
    plot(x$y~x$x,type=type,...)
    lines(x,cb=cb,...)
}

lines.VAS <- function(x,cb=FALSE,...){
    .drawECDF(x$x, x$y,ecdf=x$ecdf,...)
    if(cb){
        .drawECDF(x$x, x$lb,col='gray',lty=2)
        .drawECDF(x$x, x$ub,col='gray',lty=2)
    }
}
