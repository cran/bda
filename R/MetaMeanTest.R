#################################################
## MCAT scores

## under H0, all years have the same mean MCAT scores.  Because we
## don't have the std. dev. values for each year, we cannot
## standardize the mean scores.  We can assume homogenity and estimate
## the SD accordingly.

## Assuming homogenity. We want to estimate the common variance. For
## each group, the size, mean are n.i and xbar.i.  When n.i is large
## enough, xbar.i follows a normal distribution with mean mu and
## variance sigma^2/n.i. We can estimate the mean mu.bar=sum(n.i *
## xbar.i)/sum(n.i), and then centerize x.bar as x0.i =
## sqrt(n.i)*(xbar.i - mu.bar) which follows a normal distribtution
## with mean 0 and variance sigma^2.  Hence we can estimate
## sigma^2=sd(x0.i).

.mysd <- function(x,n){
    k <- length(x)
    stopifnot(k>1)
    mu <- sum(x*n)/sum(n)
    x0 <- sqrt(n)*(x-mu)
    ##sd(x0)
    sqrt(sum(x0*x0)/(k-1))
}

.mymean <- function(x,n){
    mu <- sum(x*n)/sum(n)
}

metamean.test <- function(x, n, group, s, mu, sigma,years,
                          alternative = c("two.sided", "less", "greater"),
                          conf.level=0.95, sizable=FALSE,legend.pos=2
                       ){
    ## sample sizes cannot be decimal values, cannot be zero or less.
    n <- round(n)
    if(any(n<1)) stop("negative size not allowed")
    ## number of years
    l <- length(n)
    ## all vectors should have the same lenth
    stopifnot(length(x)==l)
    stopifnot(length(group)==l)
    ## missing values are not allowed in 'x' and 'n', but not 's'
    if(any(is.na(x))) stop("missing value(s) in 'x'")
    if(any(is.na(n))) stop("missing value(s) in 'n'")
    ## only two levels are allowed for 'group'.
    group <- as.factor(as.character(group))
    lvlgroup <- levels(group)
    stopifnot(length(lvlgroup)==2)
    ## total size and weights
    nsum <- tapply(n,group,sum)
    ##    size1 <- size[1]; size2 <- size[2]
    grp1 <- group == lvlgroup[1]
    grp2 <- group == lvlgroup[2]
    wts <- rep(1,l) #initialize the weights to 1
    wts[grp1] <- n[grp1]/nsum[1]
    wts[grp2] <- -n[grp2]/nsum[2]
    tmp <- wts * x
    xbar1 <- sum(tmp[grp1])
    xbar2 <- abs(sum(tmp[grp2]))
    ## mu and sigma are for the underlying populations.  If they are
    ## available, we use them to correct the sample means, and use
    ## sigma to better estimate the standard deviations.
    if(missing(mu)){
        mu <- rep(0,l)
    }else{
        ## it doesn't make sense to pass one single population mean
        ## for analysis.
        stopifnot(length(mu)==l)
    }
    ## if both sigma and s are missing, we estimate the stdev using
    ## the means only.  Otherwise, compute the weighted variance using
    ## s or sigma (preferred if both are available).  Missing values
    ## are allowed in both 's' and 'sigma'.  We impute the missing
    ## values using the rest available information.
    shat <- .mysd(x,n) # estimate stdev based on the sample means ONLY.
    if(missing(sigma)){
        sigma <- rep(NA, l)
        if(missing(s)){
            s <- rep(NA, l)
        }else{
            if(length(s)==1){
                if(is.na(s)){
                    s <- rep(NA, l)
                }else{
                    s <- rep(s, l)
                }
            }else{
                stopifnot(length(s) == l)
            }
        }
    }else{
        if(length(sigma)==1){
            if(is.na(sigma)){
                sigma <- rep(NA, l)
            }else{
                sigma <- rep(sigma, l)
            }
        }else{
            stopifnot(length(sigma) == l)
        }
    }
    ## now we have two vectors of 'sigma' and 's'
    ## if 'sigma' not missing, use them
    sigma.na <- is.na(sigma)
    s.na <- is.na(s)
    sele <- !sigma.na
    if(any(sele)){
        s[sele] <- sigma[sele]
    }
    ## now if any 's' is missing, incompute if non-NA values
    ## available.  Otherwise use 'shat'
    s.na <- is.na(s)
    if(any(s.na)){
        if(sum(!s.na)==0) s <- rep(shat, l)
        else{ #impute
            n.tmp <- n[!s.na]
            s.tmp <- s[!s.na]
            sp2 <- sum((n.tmp-1)*s.tmp^2)/sum(n.tmp-1)
            s[s.na] <- sqrt(sp2)
        }
    }
    if(any(s <= 0)) stop("Non-positive SD value(s)")

    alternative <- match.arg(tolower(alternative),
                             c("two.sided", "less", "greater"))
    ##delta <- sum((x-mu)*wts)
    ##sdelta <- sqrt(sum(wts^2*s^2/n))

    delta0 <- sum((x-mu)*wts)
    delta <- sum((x-mu)/sigma*wts)
    sdelta <- sqrt(sum(1/nsum))

    TS <- delta/sdelta
    pv <- 2*pnorm(-abs(TS))
    pv = switch(alternative,
                two.sided = 2*pnorm(-abs(TS)),
                less = 1.0 - pnorm(TS),
                greater = pnorm(TS)
                )

    #list(statistic=TS, std.dev=s,means=c(xbar1,xbar2),
    #     p.value=pv)

    method <- paste("Z-test for equal means (", alternative, ")") 
    dname <- paste("means between", lvlgroup[1], "and", lvlgroup[2])

    names(TS) <- "Z"
    RVAL <- list(statistic = TS,
                 parameter = c(mean1=xbar1,
                               mean2=xbar2,
                               diff.adj=delta0),
                 method = method,
                 p.value= pv,
                 data.name=dname)
    
    class(RVAL) <- "htest"
    ## return(RVAL)
    if(missing(years)) years <- c(1:l)
    else if(length(years)==1){
        years <- as.numeric(years) + c(1:l) - 1
    }else{
        stopifnot(length(years)==l)
        years <- as.numeric(years)
    }
    
    structure(
        list(test=RVAL,
             years=years,
             xbar=x, s=s,
             mu=mu, sigma=sigma,
             size=n, group=group,
             conf.level=conf.level,
             sizable = sizable,
             legend.pos = legend.pos
             ),
        class="MetaMeans")
}

plot.MetaMeans <- function(x,ylim=NULL,type=NULL,...){
    y <- x$mu
    s <- x$s
    n <- x$size
    group <- x$group
    y2 <- x$xbar
    conf.level <- x$conf.level
    sizable <- x$sizable
    legend.pos <- x$legend.pos
    
    x <- x$years
    
    k <- length(n)
    stopifnot(conf.level>0 && conf.level <1)
    tmp <- min(conf.level, 1 - conf.level)
    uls <- y + abs(qnorm(0.5*tmp)) * s/sqrt(n)
    lls <- y - abs(qnorm(0.5*tmp)) * s/sqrt(n)
    ymax <- 1.05 * max(uls)
    ymin <- 0.95 * min(lls)
    if(is.null(ylim)){
        ylim <- c(ymin, ymax)
    }else{
        ymin <- min(ymin,ylim[1])
        ymax <- max(ymax,ylim[2])
        ylim <- c(ymin, ymax)
    }
    if(is.null(type)) type <- "b"
    cols <- c(rep("red",k))
    group <- as.factor(group)
    sele <- group == levels(group)[1]
    cols[sele] <- "blue"
    if(is.null(col)) col <- cols

    plot(y ~ x, ylim=ylim,type=type,...)
    x.diff <- diff(x)*0.25;
    xd <- c(x.diff[1], x.diff, x.diff[k-1])
    for(i in 1:k){
        segments(x[i]-xd[i],lls[i],x[i]+xd[i+1],lls[i], col=cols[i])
        segments(x[i]-xd[i],uls[i],x[i]+xd[i+1],uls[i], col=cols[i])
        segments(x[i],lls[i],x[i],uls[i], col=cols[i],lty=2)
    }

    pch2 <- rep(20,k)
    sele <- y2 > y
    pch2[sele] <- 24
    sele <- y2 < y
    pch2[sele] <- 25
    if(sizable){
        cexs <- abs(y2-y)/s+0.1
        points(x, y2, pch=pch2,cex=cexs/min(cexs))
    }else{
        points(x, y2, pch=pch2)
    }

    if(legend.pos==1 || legend.pos==2){
        y0 <- ymax
        yj <- 1
    }else{
        y0 <- ymin
        yj <- 0
    }
    if(legend.pos==1 || legend.pos==4){
        x0 <- max(x)
        xj <- 1
    }else{
        x0 <- min(x)
        xj <- 0
    }
    if(legend.pos > 0){
        legend(x0, y0, xjust=xj, yjust=yj,
               legend=c("Years without TBL","Years with TBL",
                        "National Average",
                        "Above National Average",
                        "Equals to National Average",
                        "Below National Average"),
               lty=c(1,1,0,0,0,0),
               col=c("red","blue","black","black","black","black"),
               pch=c(NA,NA,21,24,20,25))
    }

}

lines.MetaMeans <- function(x,...){
    y <- x$mu
    s <- x$s
    n <- x$size
    group <- x$group
    y2 <- x$xbar
    conf.level <- x$conf.level
    sizable <- x$sizable
    legend.pos <- x$legend.pos
    
    x <- x$years
    
    k <- length(n)
    stopifnot(conf.level>0 && conf.level <1)
    tmp <- min(conf.level, 1 - conf.level)
    uls <- y + abs(qnorm(0.5*tmp)) * s/sqrt(n)
    lls <- y - abs(qnorm(0.5*tmp)) * s/sqrt(n)
    cols <- c(rep("red",k))
    group <- as.factor(group)
    sele <- group == levels(group)[1]
    cols[sele] <- "blue"
    if(is.null(col)) col <- cols

    lines(y ~ x, ...)
    x.diff <- diff(x)*0.25;
    xd <- c(x.diff[1], x.diff, x.diff[k-1])
    for(i in 1:k){
        segments(x[i]-xd[i],lls[i],x[i]+xd[i+1],lls[i], col=cols[i])
        segments(x[i]-xd[i],uls[i],x[i]+xd[i+1],uls[i], col=cols[i])
        segments(x[i],lls[i],x[i],uls[i], col=cols[i],lty=2)
    }

    pch2 <- rep(20,k)
    sele <- y2 > y
    pch2[sele] <- 24
    sele <- y2 < y
    pch2[sele] <- 25
    if(sizable){
        cexs <- abs(y2-y)/s+0.1
        points(x, y2, pch=pch2,cex=cexs/min(cexs))
    }else{
        points(x, y2, pch=pch2)
    }
}

print.MetaMeans <- function(x,...){
    x$test
}
