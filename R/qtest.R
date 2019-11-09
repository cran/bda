## 2018/06/23
### p.test to test the goodness of fit
q.test <- function(x, p, q)
    UseMethod("q.test")

q.test.NGS.lognormal <- function(x, p, q)
{
    xnam <- x$x.name
    n <- x$n
    x.raw <- x$data
    p0 <- x$n0/n

    xt <- table(x.raw)
    Fx <- cumsum(as.numeric(xt))
    x0 <- as.numeric(names(xt))

    if(missing(p) && missing(q)){
        if(p0 > 0.1){
            p <- p0
            q <- x0[2] 
        }else{
            sele <- which(Fx/n >= 0.3)[1]
            q <- x0[sele+1]
            p <- Fx[sele]/n
        }
        cat("\n p=", round(p,3), ", q=", round(q,3),"\n")
    }

    if(p <= 0 || p >= 1)
        stop("Invalid 'p' value")

        
    fx0 <- dmlnorm(q,  x$p, x$meanlog, x$sdlog)
    qp0 <- qmlnorm(p,  x$p, x$meanlog, x$sdlog)

    z <- sqrt(n/p/(1-p))*fx0*(qp0-q)
    p.homo <- pnorm(-abs(z))*2.0
    
    names(z) <- "Z"
    RVAL <- list(statistic = z,
                 p.value = p.homo,
                 method = "Q-test for Log-Normal distribution", 
                 data.name = xnam)
    class(RVAL) <- "htest"
    return(RVAL)
}

q.test.default <- function(x, p, q)
{
    xnam <- deparse(substitute(x))
    x <- x[!is.na(x)]
    n <- length(x)
    
    if(p <= 0 || p >= 1)
        stop("Invalid 'p' value")
    
    mu <- mean(x, na.rm=TRUE)
    s <- sd(x, na.rm=TRUE)

    phat <- p
    qhat <- q
    
    fx0 <- dnorm(qhat, mu, s)
    qp0 <- qnorm(phat,  mu, s)

    z <- sqrt(n/phat/(1-phat))*fx0*(qp0-qhat)
    p.homo <- pnorm(-abs(z))*2.0
    
    names(z) <- "Z"
    RVAL <- list(statistic = z,
                 p.value = p.homo,
                 method = "Q-test for Normal distribution", 
                 data.name = xnam)
    class(RVAL) <- "htest"
    return(RVAL)
}

