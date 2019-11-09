###  Lognormal distribution
#####################################################################
## Created on Feb 8, 2019 by Bin Wang
## Laste updated on Feb 8, 2019

fit.lognormal <- function(x,method='mle'){
    method <- match.arg(tolower(method),
                        c("mle","percentile"))

    if(class(x)=="histogram"){
        if(method=='mle'){
            x <- .rmzero(x)
            res <- .fitFSDbin(x$breaks, x$counts,1,1)
            out <- c(res$pars$mean, res$pars$sd)
        }else{
            out <- .lnorm.lm(x)
        }
        res <- structure(
            list(name=x$name,
                 data = x,
                 pars = out),
            class='lognormal')
    }else if(class(x)=="bdata"){
        if(method=='mle'){
            x <- .rmzero(x)
            res <- .fitFSDbin(x$breaks, x$counts,1,1)
            out <- c(res$pars$mean, res$pars$sd)
        }else{
            out <- .lnorm.lm(x)
        }
        res <- structure(
            list(name=x$name,
                 data = x,
                 pars = out),
            class='lognormal')
    }else if(class(x)=="numeric"){
        dname <- deparse(substitute(x))
        out <- fit.mlnorm(x=x, method="lognormal")
        res <- structure(
            list(name=dname,
                 data = x,
                 pars = c(out$meanlog, out$sdlog)),
            class='lognormal')
    }else
        stop("data type not supported")
    
    res
}

plot.lognormal <- function(x,xlim=NULL,ngrid=1000,...){
    if(is.null(xlim)){
        if(class(x$data)=='bdata'){
            xbrks <- x$data$breaks
            k <- length(xbrks)
            if(is.finite(xbrks[1]))
                xmin <- xbrks[1]
            else
                xmin <- xbrks[2]
            if(is.finite(xbrks[k]))
                xmax <- xbrks[k]
            else
                xmax <- xbrks[k-1]
            xlim <- c(xmin, xmax)
        }else{
            xlim <- range(x$data)
        }
    }
    x0 <- seq(xlim[1], xlim[2], length=ngrid)
    f0 <- dlnorm(x0,x$pars[1],x$pars[2])
    
    if(class(x$data)=='bdata'){
        plot(x$data,xlim=xlim,...)
        lines(x0, f0, ...)
    }else{
        plot(x0, f0, xlim=xlim, ...)
    }
}

lines.lognormal <- function(x,xlim=NULL,ngrid=1000,...){
    if(is.null(xlim)){
        if(class(x$data)=='bdata'){
            xbrks <- x$data$breaks
            k <- length(xbrks)
            if(is.finite(xbrks[1]))
                xmin <- xbrks[1]
            else
                xmin <- xbrks[2]
            if(is.finite(xbrks[k]))
                xmax <- xbrks[k]
            else
                xmax <- xbrks[k-1]
            xlim <- c(xmin, xmax)
        }else{
            xlim <- range(x$data)
        }
    }
    x0 <- seq(xlim[1], xlim[2], length=ngrid)
    f0 <- dlnorm(x0,x$pars[1],x$pars[2])
    lines(x0, f0, ...)
}

print.lognormal <- function(x,...){
    cat("Lognormal distribution:")
    cat("\n  location :=", round(x$pars[1],4)) 
    cat("\n  scale :=", round(x$pars[2],4)) 
    cat("\n")
}

