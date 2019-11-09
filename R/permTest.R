## updated on 01/21/2018: write c code to speed up the test.

## updated on 01/22/2018: don't use grid. no binning.

## updated on 01/24/2018: write a new function perm.test to handle
## vectors, matrices, and histograms.
.mean <- function(x,y) mean(x,na.rm=TRUE) - mean(y,na.rm=TRUE)

perm.test <- function(x,y,tstat,alternative = "two.sided",iter=999)
    UseMethod("perm.test")

perm.test.default <- function(x,y,tstat,alternative = "two.sided",iter=999){
    ## the only function that user-defined test statistics are
    ## allowed.
    xnam <- deparse(substitute(x))
    ynam <- deparse(substitute(y))
    nam <- paste(xnam,"(",length(x),") versus ", ynam,"(",length(y),").")
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
                            c("less","greater","two.sided"))
    if(missing(tstat)){
        out <- .permtest0(x,y,iter=iter)
        D <- out$D
        pv <- out$pv
    }else{
        fun <- tstat
        D = fun(x, y)
        rfun <- function(x, y) {
            nx = length(x)
            ny = length(y)
            n = nx + ny
            n0 = sample(n)
            xy = c(x, y)
            fun(xy[n0[1:nx]], xy[n0[(nx + 1):n]])
        }
        z <- replicate(iter, rfun(x,y))
        p1 <- mean(abs(D) <= abs(z))
        p2 <- mean(D > z)
        p3 <- mean(D < z)
        pv = switch(alternative, two.sided = p1, 
                    less = p2,
                    greater = p3)
    }
    
    names(D) <- "D"
    RVAL <- list(statistic = D,
                 p.value = pv,
                 method = "Permutation test", 
                 data.name = nam)
    class(RVAL) <- "htest"
    return(RVAL)
}

perm.test.histogram <- function(x,y,tstat,alternative = "two.sided",iter=999){
    ## to test distributional differences between two samples --
    ## univariate analysis only.
    ## one-sided test only -- 'alternative' will be omitted.
    xbreaks <- x$breaks;
    xcounts <- x$counts;
    nx <- length(xcounts)
    xlo <- xbreaks[1:nx]
    xhi <- xbreaks[-1]
    ybreaks <- y$breaks;
    ycounts <- y$counts;
    ny <- length(ycounts)
    ylo <- ybreaks[1:ny]
    yhi <- ybreaks[-1]

    xnam <- x$xname
    ynam <- y$xname
    nam <- paste(xnam,"(",sum(xcounts),") versus ",
                 ynam,"(",sum(ycounts),").")
    
    if(!is.finite(xlo[1])){
        if(!is.finite(ylo[1])){
            xlo[1] <- ylo[1] <- min(xlo[2],ylo[2])-1
        }else{
            xlo[1] <- min(xlo[2],ylo[1])-1
        }
    }else{
        if(!is.finite(ylo[1])){
            ylo[1] <- min(xlo[1],ylo[2])-1
        }
    }
    if(!is.finite(xhi[nx])){
        if(!is.finite(yhi[ny])){
            xhi[nx] <- yhi[ny] <- max(xlo[nx],ylo[ny])+1
        }else{
            xhi[nx] <- max(xlo[nx],yhi[ny])+1
        }
    }else{
        if(!is.finite(yhi[ny])){
            yhi[ny] <- max(xhi[nx],ylo[ny])+1
        }
    }
    x <- as.matrix(cbind(xcounts, xlo, xhi))
    y <- as.matrix(cbind(ycounts, ylo, yhi))
    
    out <- .permtest3(x,y,iter=iter)
    
    D <- 000
    names(D) <- "D"
    
    RVAL <- list(statistic = D,
                 p.value = out$pv, 
                 method = "Permutation test for distributional difference", 
                 data.name = nam)
    class(RVAL) <- "htest"
    return(RVAL)
}

perm.test.bdata <- function(x,y,tstat,alternative = "two.sided",iter=999){
    ## to test distributional differences between two samples --
    ## univariate analysis only.
    ## one-sided test only -- 'alternative' will be omitted.
    xbreaks <- x$breaks;
    xcounts <- x$counts;
    nx <- length(xcounts)
    xlo <- xbreaks[1:nx]
    xhi <- xbreaks[-1]
    ybreaks <- y$breaks;
    ycounts <- y$counts;
    ny <- length(ycounts)
    ylo <- ybreaks[1:ny]
    yhi <- ybreaks[-1]

    xnam <- x$name
    ynam <- y$name
    nam <- paste(xnam,"(",sum(xcounts),") versus ",
                 ynam,"(",sum(ycounts),").")

    if(!is.finite(xlo[1])){
        if(!is.finite(ylo[1])){
            xlo[1] <- ylo[1] <- min(xlo[2],ylo[2])-1
        }else{
            xlo[1] <- min(xlo[2],ylo[1])-1
        }
    }else{
        if(!is.finite(ylo[1])){
            ylo[1] <- min(xlo[1],ylo[2])-1
        }
    }
    if(!is.finite(xhi[nx])){
        if(!is.finite(yhi[ny])){
            xhi[nx] <- yhi[ny] <- max(xlo[nx],ylo[ny])+1
        }else{
            xhi[nx] <- max(xlo[nx],yhi[ny])+1
        }
    }else{
        if(!is.finite(yhi[ny])){
            yhi[ny] <- max(xhi[nx],ylo[ny])+1
        }
    }
    x <- as.matrix(cbind(xcounts, xlo, xhi))
    y <- as.matrix(cbind(ycounts, ylo, yhi))
    
    out <- .permtest3(x,y,iter=iter)
    
    D <- 111
    names(D) <- "D"
    
    RVAL <- list(statistic = D,
                 p.value = out$pv, 
                 method = "Permutation test for distributional difference", 
                 data.name = nam)
    class(RVAL) <- "htest"
    return(RVAL)
}

perm.test.matrix <- function(x,y,tstat,alternative = "two.sided",iter=999){
    ## to test distributional differences between two samples --
    ## univariate analysis only. Both 'x' and 'y' should have three
    ## columns: col1=frequencies, col2=lower limits, and col3=upper
    ## limits. Infinite values are allowed.  Missing values should be
    ## excluded.
    ## one-sided test only -- 'alternative' will be omitted.

    ## data validation:
    stopifnot(is.matrix(x)) # classes
    stopifnot(is.matrix(y))
    stopifnot(ncol(x)==3) # columns
    stopifnot(ncol(y)==3)
    if(any(round(x[,1])<=0))
        stop("Invalid count(s) in 'x'")
    if(any(round(y[,1])<=0))
        stop("Invalid count(s) in 'y'")
    if(any(x[,2]>x[,3]))
        stop("Invalid interval(s) in 'x'")
    if(any(y[,2]>y[,3]))
        stop("Invalid interval(s) in 'y'")
        
    xnam <- deparse(substitute(x))
    ynam <- deparse(substitute(y))
    nam <- paste(xnam,"(",sum(x[,1]),") versus ",
                 ynam,"(",sum(y[,1]),").")

    out <- .permtest3(x,y,iter=iter)
    
    D <- 111
    names(D) <- "D"
    
    RVAL <- list(statistic = D,
                 p.value = out$pv, 
                 method = "Permutation test for distributional difference", 
                 data.name = nam)
    class(RVAL) <- "htest"
    return(RVAL)
}

##perm.test.NGS <- function(x,y,tstat,alternative = "two.sided",iter=999){
    ## not for general data.frame.  Just for NGS which is also a
    ## data.frame. To handle ks-type test with very large iternations
    ## -- coded in C. 
    ## one-sided test only -- 'alternative' will be omitted.
##    warning("to be developed")
##}

.PermTestNGS <- function(x,y, iter = 1001, algorithm=0){
    out <- NULL
    if(algorithm==1){
        out <- .permtest1(x=x,y=y,iter=iter)
    }else if(algorithm==2){
        out <- .permtest2(x=x,y=y,iter=iter)
    }else{
        out <- .permtest0(x=x,y=y,iter=iter)
    }
    out
}


.permtest0 <- function(x,y, iter = 1001){
    xnam <- deparse(substitute(x))
    ynam <- deparse(substitute(y))
    nam <- paste(xnam, "(", length(x), ") vs ",
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

    nx <- length(x); ny <- length(y)
    options(warn = -1)
    D <- ks.test(x,y)$statistic # observed t-score
    options(warn = 0)
    xy <- table(c(x,y))
    M <- length(xy)
    F <- as.numeric(xy)
    out <- .Fortran(.F_permtest2, pv = as.double(D),
                    as.integer(M), as.integer(nx),
                    as.integer(ny), as.integer(F),
                    as.integer(iter))
    list(pv=out$pv, D=D)
}

.permtest1 <- function(x,y, iter = 1001){
    xnam <- deparse(substitute(x))
    ynam <- deparse(substitute(y))
    nam <- paste(xnam, "(", length(x), ") vs ",
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
    options(warn = -1)
    nx <- length(x); ny <- length(y)
    D <- ks.test(x,y)$statistic # observed t-score
    xy <- c(x,y)
    xmin <- min(xy)
    xmax <- max(xy)

    out <- .Fortran(.F_permtest,
                    as.double(x), as.integer(nx),
                    as.double(y), as.integer(ny),
                    as.double(xmin), as.double(xmax),
                    D=as.double(D), pv=as.double(0),
                    as.integer(iter))
  
    RVAL <- list(statistic = D,
                 p.value = out$pv, 
                 method = "Permutation test for NGS data", 
                 data.name = nam)
    class(RVAL) <- "htest"
    options(warn = 0)
    return(RVAL)
}

.permtest2 <- function(x,y, iter = 1001){
    xnam <- deparse(substitute(x))
    ynam <- deparse(substitute(y))
    nam <- paste(xnam, "(", length(x), ") vs ",
                 ynam, "(", length(y), ")")
    if (!is.numeric(x) && !is.complex(x) && !is.logical(x)) {
        warning("'x' is not numeric or logical: returning NA")
        return(NA_real_)
    }
    x <- x[!is.na(x)]
    nx <- length(x)
    if (!is.numeric(y) && !is.complex(y) && !is.logical(y)) {
        warning("'y' is not numeric or logical: returning NA")
        return(NA_real_)
    }
    y <- y[!is.na(y)]
    
    options(warn = -1)
    R <- iter
    z <- c(x,y) # pooled sample
    K <- 1:length(z)
    reps <- numeric(R)  # storage for replicates
    D <- ks.test(x,y)$statistic # observed t-score
    for(i in 1:R){
        k <- sample(K, size=nx, replace=FALSE)
        x1 <- z[k]
        y1 <- z[-k]
        reps[i] <- ks.test(x1, y1)$statistic
    }
    ## this is a one-tailed test.  D = sup|Fn(X) - Fn(Y)|
    p <- mean(c(D, reps) >= D)

    RVAL <- list(statistic = D,
                 p.value = p, 
                 method = "Permutation test for NGS data", 
                 data.name = nam)
    class(RVAL) <- "htest"
    options(warn = 0)
    return(RVAL)
}

.permtest3 <- function(x,y, iter = 999){
    cat("\nData X :\n")
    print(x)
    cat("\nData Y :\n")
    print(y)
#    out <- .Fortran(.F_permtest,
#                    as.double(x), as.integer(nx),
#                    as.double(y), as.integer(ny),
#                    as.double(xmin), as.double(xmax),
#                    D=as.double(D), pv=as.double(0),
#                    as.integer(iter))
    cat("\n\nThis program is incomplete. The outputs are faked!")
    list(pv=runif(1),D=runif(1))
}

perm.mean <- function(x,y, iter = 9999){
    ## call a permutation test to compare the mean values. Return
    ## three p-values: p.lt, p.gt and p.ne
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

    nx <- length(x); ny <- length(y)
    D <- mean(x) - mean(y)
    xy <- c(x,y)
    x.sum <- sum(xy)
    out <- .Fortran(.F_permtest3, as.double(xy),
                    as.integer(nx), as.integer(ny),
                    pv = as.double(c(D,x.sum,0)),
                    as.integer(iter))
    list(pv=out$pv)
}
