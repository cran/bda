## VAS & MEM

VAS.ecdf <- function(x,alpha=0.05){
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
