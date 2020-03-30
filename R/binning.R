
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
        if(any(counts <= 0))
            stop("zeros and negative counts not allowed")
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
    out
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

