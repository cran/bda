
##################
pmlnorm <- function(q,p,mean,sd){
    sapply(q,FUN=.pmlnorm, p=p,mu=mean,s=sd)
}

.pmlnorm <- function(x,p,mu,s){
    k <- length(p)
    stopifnot(length(x) == 1)
    stopifnot(length(mu) == k)
    stopifnot(length(s) == k)
    stopifnot(all(p>0))
    p <- p/sum(p)
    .myfun <- function(x){
        if(is.finite(x[2])){
            out <- plnorm(x[1],meanlog = x[2], sdlog = x[3])
        }else{
            out <- 1.0
        }
    }
    x0 <- rep(x, k)
    tmp <- cbind(x0,mu,s)
    sum(p*apply(tmp, 1, .myfun))
}

dmlnorm <- function(x,p,mean,sd){
    sapply(x,FUN=.dmlnorm, p=p,mu=mean,s=sd)
}

.dmlnorm <- function(x,p,mu,s){
    k <- length(p)
    stopifnot(length(x) == 1)
    stopifnot(length(mu) == k)
    stopifnot(length(s) == k)
    stopifnot(all(p>0))
    .myfun <- function(x){
        out <- 0
        if(is.finite(x[2])){
            out <- dlnorm(x[1],meanlog = x[2], sdlog = x[3])
        }
        out
    }
    x0 <- rep(x, k)
    tmp <- cbind(x0,mu,s)
    sum(p*apply(tmp, 1, .myfun))
}

rmlnorm <- function(n,p,mean,sd){
    mu <- mean
    s <- sd
    k <- length(p)
    stopifnot(length(n) == 1)
    n <- round(n)
    stopifnot(n>0)
    stopifnot(length(mu) == k)
    stopifnot(length(s) == k)
    stopifnot(all(p>0))
    .myfun <- function(x){
        if(is.finite(x[2])){
            out <- rlnorm(x[1],meanlog = x[2], sdlog = x[3])
        }else{
            out <- rep(0, x[1])
        }
        out
    }
    x0 <- rep(n, k)
    tmp <- cbind(n,mu,s)
    out <-  apply(tmp, 1, .myfun)
    gid <- sample(1:k,size=n,replace=TRUE,prob=p)
    gid <- cbind(1:n,gid)
    out[gid]
}


qmlnorm <- function(prob,p,mean,sd){
    sapply(prob,FUN=.qmlnorm, p=p,mu=mean,s=sd)
}

.qmlnorm <- function(x,p,mu,s){
    k <- length(p)
    stopifnot(length(x) == 1)
    stopifnot(length(mu) == k)
    stopifnot(length(s) == k)
    stopifnot(all(p>0))
    p <- p/sum(p)
    if(x <= 0){
        res <- 0
    }else if(x==1){
        res <- Inf
    }else{
        ## find xa and xb such that F(xa) > x and F(xb) <= x
        xa <- 0.1
        xb <- 10
        Fxa <- pmlnorm(xa,p,mu,s)
        Fxb <- pmlnorm(xb,p,mu,s)
        if(Fxa > x){
            xb <- xa
            for(i in 1:10000){
                xa <- 0.5 * xb
                Fxa <- pmlnorm(xa,p,mu,s)
                if(Fxa > x){
                    xb <- xa
                    xa <- 0.5 * xb
                }else{
                    ##cat("\n iter1=", i, "xa=",xa, "xb=",xb)
                    break
                }
            }
        }else if(Fxb < x){
            xa <- xb
            for(i in 1:10000){
                xb <- 2 * xa
                Fxb <- pmlnorm(xb,p,mu,s)
                if(Fxb < x){
                    xa <- xb
                    xb <- 2.0 * xb
                }else{
                    ##cat("\n iter1=", i, "xa=",xa, "xb=",xb)
                    break
                }
            }
        }
        for(i in 1:10000){
            xc <- (xa+xb)*0.5
            Fxc <- pmlnorm(xc,p,mu,s)
            if(Fxc < x){
                xa <- xc
            }else{
                xb <- xc
            }
            if(is.na(xb - xa)){
                res <- xa
                break
            }else if(xb - xa < 0.000001){
                res <- xc
                break
            }
        }
    }
    return(res)
}
