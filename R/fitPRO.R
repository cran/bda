#####################################################################
## Created on Oct 11, 2023 by Bin Wang
##Laste updated on Oct 14, 2023

## if 'x' is a matrix or data.frome with 2 columns, assume it is
## interval-valuded PRO; otherwise single value PRO measures (VAS or
## NRS), and estimate via binning.

fit.PRO <- function(x,dist,x.range,nclass){
    if(missing(dist)) dist <- "norm"
    else dist <- as.character(dist)

    if(inherits(x,"data.frame")) x <- as.matrix(x)
    
    if(inherits(x,"matrix")){
        if(ncol(x) == 1){
            xb <- binning(x[,1])
        }else if(ncol(x) == 2){
            xb <- .mybinPRO(ll=x[,1],ul=x[,2],
                            x.range=x.range,nclass=nclass)
        }else
            stop("'x' must have two columns")
    }else if(inherits(x,'histogram')){
        xb <- binning(x)
    }else if(inherits(x,'bdata')){
        xb <- x
        x0 <- xb$freq
        xbrk <- xb$breaks
        sele <- x0 == 0
        if(any(sele)){
            xb$freq <- x0[-which(sele)]
            xb$breaks <- xbrk[-(which(sele)+1)]
        }
        
        x.fit <- .fit.FSD1(xb,dist=dist)
        y.fit <- NULL
        x0 <- x.fit$x
        sele <- is.finite(x0)
        x0 <- x0[sele]
        y0 <- x.fit$y[sele]
        z0 <- NULL
        out <- structure(
            list(x.fit=x.fit,
                 y.fit = y.fit,
                 Psi = NULL,
                 Psi.rng = NULL,
                 x=x0,y=y0,z=z0),
            class="FSD")
    }else{
        if(length(x)<10){
            warning("too few data points")
            out <- NULL
        }else{
            xt <- table(x)
            if(length(xt)<3){
                warning("too few distinct data values")
                out <- NULL
            }else{
                xb <- hist(x,plot=FALSE)
                xb <- binning(xb)
                x0 <- xb$freq
                xbrk <- xb$breaks
                sele <- x0 == 0
                if(any(sele)){
                    xb$freq <- x0[-which(sele)]
                    xb$breaks <- xbrk[-(which(sele)+1)]
                }
                
                x.fit <- .fit.FSD1(xb,dist=dist)
                y.fit <- NULL
                x0 <- x.fit$x
                sele <- is.finite(x0)
                x0 <- x0[sele]
                y0 <- x.fit$y[sele]
                z0 <- NULL
                out <- structure(
                    list(x.fit=x.fit,
                         y.fit = y.fit,
                         Psi = NULL,
                         Psi.rng = NULL,
                         x=x0,y=y0,z=z0),
                    class="FSD")
            }
        }
    }
    out
}

.mybinPRO <- function(ll,ul,x.range,nclass){
    if(missing(nclass)){
        nclass <- 10
    }else if(!is.numeric(nclass)){
        nclass <- 10
    }else{
        nclass <- round(nclass)
        if(nclass>15 || nclass < 3)
            stop("too few classes")
    }
    
    if(missing(x.range)){
        a <- min(c(ll,ul))
        b <- max(c(ll,ul))
        ll[ll < a] <- a; ll[ll > b] <- b
        ul[ll < a] <- a; ul[ll > b] <- b
    }else{
        k <- length(x.range)
        if(k != 2) stop("'x.range' must have length 2")
        a <- x.range[1]
        b <- x.range[2]
        if(a >= b) stop("invalid 'x.range'")
    }
    if(any(ll > ul))
        stop("lower limit can not be larger than upper limit")
       
    res <- .Fortran(.F_probin,
                    as.double(ll),
                    as.double(ul),
                    as.integer(length(ll)),
                    as.double(a),
                    as.double(b),
                    as.integer(nclass),
                    freq = as.double(rep(0,nclass)))
    
    xbrk <- seq(a,b,length=nclass+1)
    xtmp <- rnorm(100)
    xh <- hist(xtmp, plot=FALSE)
    structure(
        list(ll = xbrk[-(nclass+1)],
             ul = xbrk[-1],
             breaks=xbrk,
             freq = res$freq,
             xhist = xh,
             xZipf=NULL),
        class="bdata")
}
