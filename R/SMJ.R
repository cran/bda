
#############################################################
## Find the best fit

## Inputs: 1) counts, 2) breaks; 3) weights

fitGLD2bdata <- function(x, weights){
    if(class(x)=="histogram"||class(x)=="bdata"){
        x <- .rmzero(x)
        counts <- x$counts;
        breaks <- x$breaks
        stopifnot(length(breaks)-length(counts) == 1)
        stopifnot(length(counts) >= 4)
        ## weights can be activated later
        if(missing(weights)){
            weights <- rep(1, length(counts)-1)
        }else{
            stopifnot(all(weights>=0 & weights<=1))
        }
        stopifnot(length(breaks)-length(weights) == 2)
        ## find different combinations
        xbrks <- breaks
        w <- weights
        x <- counts
        k <- length(xbrks)-2; 
        pmat <- NULL
        for(i1 in 1:(k-4)){
            for(i2 in (i1+1):(k-3)){
                for(i3 in (i2+1):(k-2)){
                    for(i4 in (i3+1):(k-1)){
                        for(i5 in (i4+1):k){
                            y2 <- c(i1,i2,i3,i4,i5)
                            pmat <- rbind(pmat, y2)
                        }
                    }
                }
            }
        }
        
        Fn <- cumsum(x)/sum(x)
        k <- length(Fn)
        F1 <- log(1.0 - Fn[-k])     
        xhist <- binning(breaks=xbrks, counts=x)
        x0 <- xbrks[-c(1,length(xbrks))];
        lx0 <- log(x0)
        
        parGLD <- NULL
        MSE <- NULL
        Dmax <- 1e9
        parbest <- NULL
        ibest <- NULL
        gldfit <- NULL
        zipf <- NULL
        for(i in 1:nrow(pmat)){
            pts <- pmat[i,]
            fit0 <- fit.GLD.FMKL(xhist,
                                 qtl=xbrks[pts+1], qtl.levels=Fn[pts],
                                 percentile='exact', mle=FALSE)
            parGLD <- rbind(parGLD, fit0$pars)
            Fgld <- pgld(x0, fit0$pars)
            F2 <- log(1.0-Fgld)
            wts <- w/sum(w)
            mse <- sum((F1-F2)^2*wts)
            MSE <- c(MSE, mse)
            if(mse < Dmax){
                ibest <- i
                Dmax <- mse
                parbest <- fit0$pars
                gldfit <- fit0
                zipf <- list(x=lx0,y=F2)
            }
        }
        ## fit Pareto (type I)
        lmout <- lm(F1~lx0)
        pidx <- -coef(lmout)[[2]]
        xm <- exp(coef(lmout)[[1]]/pidx)
        ##cat("\nBest fit with percentiles:", pmat[ibest,])
        ##cat("\nFitted FMKL GLD parameters:\n\t", parbest,"\n")
        out <- list(pMat=pmat,
                    parMat=parGLD,
                    dMat=MSE,
                    Zipf = list(x=lx0,y=F1),
                    parsPareto = list(xm=xm,p.index=pidx),
                    ZipfGLD = zipf,
                    fitGLD=gldfit,
                    parsGLD = parbest,
                    percentiles=pmat[ibest,],
                    D = Dmax
                    )
    }else
        stop("data type not supported")
    out
}

.lnorm.lm <- function(x){
    xc <- x$counts
    bk <- x$breaks
    stopifnot(all(bk>=0))
    k <- length(bk)
    Fn <- cumsum(xc)/sum(xc);
    bk <- bk[-1]
    sele <- which(diff(Fn)==0)
    if(any(sele)){
        Fn <- Fn[-c(sele+1)]
        bk <- bk[-c(sele+1)]
    }
    sele <- which(Fn==1)
    if(any(sele)){
        Fn <- Fn[-sele]
        bk <- bk[-sele]
    }
    z <- qnorm(Fn)
    lx <- log(bk)
    lmout <- lm(lx~z)
    coef(lmout)
}
#(lnout <- myest(xhist)$par)

.rmzero <- function(x){
    sele <- x$count==0
    xc <- x$count
    bk <- x$breaks
    xbrks <- bk
    if(any(sele)){
        sele2 <- which(sele)
        xc <- xc[-sele2]
        bk <- bk[-(sele2+1)]
        if(any(sele2==length(xbrks)-1)){
            bk[length(bk)] <- Inf
        }
    }
    binning(counts=xc, breaks=bk)
}
