
binning <- function(x, counts, breaks, lower.limit, upper.limit)
    UseMethod("binning")

binning.histogram <-
    function(x, counts, breaks, lower.limit, upper.limit)
{
    freq <- x$counts
    nclass <- length(freq)
    ll <- x$breaks[-(1+nclass)]
    ul <- x$breaks[-1]
    if(any(ll<0)){
        xZipf <- NULL
    }else{
        lx <- log(ul)
        Fn <- 1-cumsum(freq)/sum(freq)
        ly <- log(Fn[-nclass]); lx <- lx[-nclass]
        xZipf <- list(x=lx,y=ly)
    }
    
    structure(
        list(ll = ll,
             ul = ul,
             breaks=x$breaks,
             freq = freq,
             xhist = x,
             xZipf=xZipf),
        class="bdata")
}

binning.default <- function(x, counts, breaks, lower.limit, upper.limit)
{
    if(missing(x)){
        if(missing(counts))
            stop("'x' and 'counts' cannot be both missing")        
        n <- length(counts)
        if(any(counts < 0))
            stop("negative counts not allowed")
        icounts <- round(counts)
        if(any(icounts != counts))
            stop("'counts' must be integers")

        if(missing(breaks)){
            if(missing(lower.limit))
                stop("'breaks' and 'lower.limit' cannot be both missing")
            if(length(lower.limit) != n)
                stop("lengths of 'counts' and 'lower.limit' differ")
            tmp <- diff(lower.limit)
            if(any(tmp <= 0))
                stop("invalid 'lower.limit' value(s)")
            ll <- lower.limit
            if(missing(upper.limit)){
                ul <- c(ll[-1],Inf)
            }else{
                if(length(upper.limit) != n)
                    stop("lengths of 'counts' and 'upper.limit' differ")
                tmp <- diff(upper.limit)
                if(any(tmp <= 0))
                    stop("invalid 'upper.limit' value(s)")
                ul <- upper.limit
            }
            ll2 <- ll[-1]
            ul2 <- ul[-n]
            breaks <- c(ll[1],(ll2+ul2)/2,ul[n])
        }else{#lower and upper limits won't be used
            if(length(breaks) - n != 1)
                stop("wrong length of 'breaks'")
            tmp <- diff(breaks)
            if(any(tmp <= 0))
                stop("invalid break value(s)")
            ll <- breaks[-(n+1)]
            ul <- breaks[-1]
        }

        if(any(ll<0)){
            xZipf <- NULL
        }else{
            lx <- log(ul)
            Fn <- 1-cumsum(counts)/sum(counts)
            ly <- log(Fn[-n]); lx <- lx[-n]
            xZipf <- list(x=lx,y=ly)
        }

        x <- rnorm(100)
        xhist <- hist(x, plot=FALSE)
        brk2 <- (ll[-1]+ul[-n])/2
        if(!is.finite(ll[1])){
            if(!is.finite(ul[n])){
                xbrks <- brk2
                xcnts <- counts[-c(1,n)]
            }else{
                xbrks <- c(brk2,ul[n])
                xcnts <- counts[-1]
            }
        }else{
            if(!is.finite(ul[n])){
                xbrks <- c(ll[1],brk2)
                xcnts <- counts[-n]
            }else{
                xbrks <- c(ll[1],brk2,ul[n])
                xcnts <- counts
            }
        }
                    
        xhist$breaks <- xbrks
        xhist$counts <- xcnts
        bw <- diff(xbrks)
        dens <- xcnts/sum(counts)/bw
        xhist$density <- dens
        xhist$mids <- xbrks[-1] - bw/2
        if(all(bw==bw[1])){
            xhist$equidist <- TRUE
        }else{
            xhist$equidist <- FALSE
        }
        
        out <- structure(
            list(ll = ll,
                 ul = ul,
                 breaks = breaks,
                 freq = counts,
                 xhist = xhist,
                 xZipf=xZipf),
            class="bdata")
        
    }else{
        if(inherits(x,"matrix")){
            if(nrow(x) == 2){
                x <- t(x)
            }
            if(ncol(x) == 2){
                n.na <- apply(is.na(x),1,sum)
                if(any(n.na > 0)){
                    x <- x[n.na == 0,]
                    if(nrow(x) < 10)
                        stop("too few data points")
                }
                out <- .bin2d(x,breaks=breaks)
            }else{
                out <- NULL
            }
        }else{
            sele <- is.na(x)
            if(any(sele)){
                if(mean(sele)==1)
                    stop("all values are missing")
                x <- x[!sele]
            }
            sele <- !is.finite(x)
            if(any(sele)){
                if(mean(sele)==1)
                    stop("all values are infinite")
                x <- x[!sele]
            }
            if(missing(breaks)){
                res <- hist(x, plot=FALSE)
            }else{
                res <- hist(x, breaks=breaks, plot=FALSE)
            }
            out <- binning(res)
        }
    }
    out
}

.bin2d <- function(x,breaks){
    ## 'x' has been checked previously
    ## 'breaks' needs to be checked
    if(missing(breaks)){
        xh1 <- binning(x[,1])
        xh2 <- binning(x[,2])
        brks <- list(x=xh1$breaks,y=xh2$breaks)
    }else{
        if(is.null(breaks$x)){
            xh1 <- binning(x[,1])
            if(is.null(breaks$y)){
                xh2 <- binning(x[,2])
                brks <- list(x=xh1$breaks,y=xh2$breaks)
            }else{
                if(any(diff(breaks$y)<=0))
                    stop("invalid value(s) in 'breaks' for y-variable")
                if(min(breaks$y) > min(x[,2]))
                    stop("invalid value(s) in 'breaks' for y-variable")
                if(max(breaks$y) < max(x[,2]))
                    stop("invalid value(s) in 'breaks' for y-variable")
                brks <- list(x=xh1$breaks,y=breaks$y)
            }
        }else{
            if(any(diff(breaks$x)<=0))
                    stop("invalid value(s) in 'breaks' for x-variable")
            if(min(breaks$x) > min(x[,1]))
                stop("invalid value(s) in 'breaks' for x-variable")
            if(max(breaks$x) < max(x[,1]))
                stop("invalid value(s) in 'breaks' for x-variable")
            if(is.null(breaks$y)){
                xh2 <- binning(x[,2])
                brks <- list(x=breaks$x,y=xh2$breaks)
            }else{
                if(any(diff(breaks$y)<=0))
                    stop("invalid value(s) in 'breaks' for y-variable")
                if(min(breaks$y) > min(x[,2]))
                    stop("invalid value(s) in 'breaks' for y-variable")
                if(max(breaks$y) < max(x[,2]))
                    stop("invalid value(s) in 'breaks' for y-variable")
                brks <- list(x=breaks$x,y=breaks$y)
            }
        }
    }
    ## find a 2-d matrix xy using 'x'[nx2] and 'breaks'
    xy <- .counts2d(x,breaks=brks)
    list(xy=xy,breaks=brks)
}

.counts2d <- function(x, breaks){
    nx <- nrow(x)
    x1 <- x[,1]; x2 <- x[,2]
    x1brk <- breaks$x; nbrk1 <- length(x1brk)
    x2brk <- breaks$y; nbrk2 <- length(x2brk)
    counts <- rep(0, (nbrk1-1)*(nbrk2-1))
    res <- .Fortran(.F_bin2d,
                    as.double(x1),
                    as.double(x2),
                    as.integer(nx),
                    as.double(x1brk),
                    as.integer(nbrk1),
                    as.double(x2brk),
                    as.integer(nbrk2),
                    xy=as.double(counts))
    xy <- matrix(res$xy, ncol=(nbrk1-1),nrow=(nbrk2-1))
}
print.bdata <- function(x,...){
    out <- data.frame(Frequency=x$freq)
    rnames <- paste(x$ll," < x <= ",x$ul,sep='')
    rownames(out) <- rnames
    print(out)
}


plot.bdata <- function(x,type="hist",...){
    type <- match.arg(tolower(type),c("histogram","zipf"))
    if(type=='histogram'){
        X <- x$xhist
        plot(X,...)
    }else{
        Log.x <- x$xZipf$x
        Log.Rank <- x$xZipf$y
        plot(Log.x, Log.Rank,...)
    }
    invisible(NULL)
}

