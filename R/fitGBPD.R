###  Fit Mixture of Generalized Beta and Pareto (Type I) Distribution
#####################################################################
## Created on Mar 18, 2020 by Bin Wang
##Laste updated on Mar 18, 2020

## try a mixture of two Gaussian distribution first: there must be at
## least four classes (4 upper boundaries) so we can compute the
## parameters using different quantiles (add 1 more class with
## frequency 1). If the data is top-coded, five classes are needed.

## Created on Mar 18, 2020 by Bin Wang
##Laste updated on Mar 18, 2020


fit.GBP <- function(x,breaks){
    pars <- NULL

    if(class(x)=="histogram"){
        freq <- x$counts
        breaks <- x$breaks
    }else if(class(x)=="bdata"){
        freq <- x$freq
        breaks <- x$breaks
    }else if(is.numeric(x)){
        if(missing(breaks)){
            y <- hist(x, plot=FALSE)
            freq <- y$counts
            breaks <- y$breaks
        }else{
            if(length(x)+1 == length(breaks)){
                x <- round(x)
                stopifnot(all(x>0))
                d <- diff(breaks)
                stopifnot(all(d>0))
                stopifnot(breaks[1]>=0)
                freq <- x
            }else
                stop("'x' and 'breaks' have different lengths")
        }
    }
    ## find initial estimates using MOP
    out <- .estGBP(freq, breaks)
    ## estimate the mixing coefficient
    p1 <- plnorm(breaks, meanlog=out[1,1],sdlog=out[1,2])
    p2 <- plnorm(breaks, meanlog=out[2,1],sdlog=out[2,2])
    p1 <- diff(p1); p2 <- diff(p2)
    P1 <- sum(p1/(p1+p2)*freq)/sum(freq)
    pars <- data.frame(prop=c(P1,1-P1),
                       Mean=c(out[1,1],out[2,1]),
                       SD=c(out[1,2],out[2,2]))
    rownames(pars) <- c("Comp1","Comp2")
    pars
}

.estGBP <- function(freq, breaks){
    nclass <- length(freq)
    if(is.finite(rev(breaks)[1])){
        stopifnot(nclass>3)
        tmp <- c(freq,1); n <- sum(tmp)
        Fn <- cumsum(tmp)/n
        Fn <- Fn[-(nclass+1)]
        qtl1 <- breaks[2:3]
        lvl1 <- Fn[1:2]
        par1 <- .mopGBP(qtl1,lvl1)
        qtl2 <- breaks[nclass:(nclass+1)]
        lvl2 <- Fn[(nclass-1):nclass]
        par2 <- .mopGBP(qtl2,lvl2)
        
    }else{
        stopifnot(nclass>4)
        Fn <- cumsum(freq)/sum(freq)
        qtl1 <- breaks[2:3]
        lvl1 <- Fn[1:2]
        par1 <- .mopGBP(qtl1,lvl1)

        qtl2 <- breaks[(nclass-1):nclass]
        lvl2 <- Fn[(nclass-2):(nclass-1)]
        par2 <- .mopGBP(qtl2,lvl2)
    }
    rbind(par1,par2)
}

.mopGBP <- function(qtl,lvl){
    qtl <- log(qtl)
    q1 <- qtl[1]; q2 <- qtl[2]
    Q1 <- qnorm(lvl[1]);Q2 <- qnorm(lvl[2])
    if(Q2 != 0){
        mu <- (q1-Q1/Q2*q2)/(1-Q1/Q2)
        sigma <- (q1-q2)/(Q1-Q2)
    }else{
        mu <- q2
        sigma <- (q1-q2)/Q1
    }
    c(mu,sigma)
}

pGBP <- function(x,pars){
    p1 <- plnorm(x,pars$Mean[1],pars$SD[1])
    p2 <- plnorm(x,pars$Mean[2],pars$SD[2])
    p <- p1*pars$prop[1]+p2*pars$prop[2]
 }

dGBP <- function(x,pars){
    d1 <- dlnorm(x,pars$Mean[1],pars$SD[1])
    d2 <- dlnorm(x,pars$Mean[2],pars$SD[2])
    d <- d1*pars$prop[1]+d2*pars$prop[2]
}
