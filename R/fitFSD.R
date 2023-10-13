###  Firm Size and/or Age Distribution
#####################################################################
## Created on Feb 17, 2020 by Bin Wang
##Laste updated on Feb 17, 2020

## 'x' can be a vector or a matrix. If 'x' is a matrix, both
## 'x.breaks' and 'y.breaks' must be provided. If 'x' is a vector,
## only 'x.breaks' is required and 'y.breaks' won't be used.

##  Part of R package BDA
##  Copyright (C) 2009-2010 Bin Wang
##
##  Unlimited use and distribution (see LICENCE).

##  Created on 2020/03/09

## GPD parameters: 1) scale psi>0; 2) shape kappa.
## If kappa>0 -> range=(0,psi/kappa)
## If kappa<0 -> range=(0,Inf)


fit.FSD <- function(x,breaks,dist){
    Psi <- NULL
    Psi.rng <- NULL
    if(is.matrix(x)){
        if(missing(breaks))
            stop("'breaks' cannot be missing for bivariate distribution")
        if(missing(breaks))
            stop("'breaks' cannot be missing")
        x.breaks <- breaks$age
        if(any(diff(x.breaks)<=0))
            stop("invalid value(s) in 'breaks' for x-variable")
        y.breaks <- breaks$size
        if(any(diff(y.breaks)<=0))
            stop("invalid value(s) in 'breaks' for y-variable")
        
        xy <- x
        if(any(is.na(xy)))
            stop("missing value(s) in 'x'")
        if(any(!is.finite(xy)))
            stop("infinite value(s) in 'x'")
        
        x <- apply(xy,1,sum)
        y <- apply(xy,2,sum)
        
        if(missing(dist)){
            x.dist <- "GLD2"
            y.dist <- "GLD2"
        }else{
            dist <- as.character(dist)
            if(length(dist)==1){
                x.dist <- dist
                y.dist <- dist
            }else{
                x.dist <- dist[1]
                y.dist <- dist[2]
            }
        }
                
        xb <- binning(counts=x,breaks=x.breaks)
        x.fit <- .fit.FSD1(xb,dist=x.dist)

        yb <- binning(counts=y,breaks=y.breaks)
        y.fit <- .fit.FSD1(yb,dist=y.dist)

        Psi.rng <- .est.Psi(xy)
        
        options(warn=-1)
        ##res <- optimize(.X2Copula,
        res <- optimize(.llkCopula,
                        interval=c(Psi.rng$min,Psi.rng$max),
                        Fx=x.fit$y,Fy=y.fit$y,cnts=xy)
        options(warn=0)
        Psi <- res$minimum
        x0 <- x.fit$x
        sele.x <- is.finite(x0)
        y0 <- y.fit$x
        sele.y <- is.finite(y0)
        x0 <- x0[sele.x]; y0 <- y0[sele.y]
        z0 <- .pdfCopula(Psi,x.fit$y[sele.x],y.fit$y[sele.y],
                         x.fit$y2[sele.x],y.fit$y2[sele.y])
    }else{
        
        if(missing(dist)) dist <- "GLD2"
        else dist <- as.character(dist)
        
        if(inherits(x,'histogram')||inherits(x,'bdata')){
            x.fit <- .fit.FSD1(x,dist=dist)
        }else{
            if(missing(breaks)){
                if(length(x)<10)
                    stop("too few data points")
                xt <- table(x)
                if(length(xt)<3)
                    stop("too few distinct data values")
                xb <- hist(x,plot=FALSE)
            }else{
                xb <- binning(counts=x,breaks=breaks)
            }
            x.fit <- .fit.FSD1(xb,dist=dist)
        }
        y.fit <- NULL
        x0 <- x.fit$x
        sele <- is.finite(x0)
        x0 <- x0[sele]
        y0 <- x.fit$y[sele]
        z0 <- NULL
    }
    out <- structure(
        list(x.fit=x.fit,
             y.fit = y.fit,
             Psi = Psi,
             Psi.rng = Psi.rng$summ,
             x=x0,y=y0,z=z0),
        class="FSD")
}

fit.Copula <- function(x,y,xy){
    if(!inherits(x,"FSD"))
        stop("'x' must be a fitted FSD")
    if(!inherits(y,"FSD"))
        stop("'y' must be a fitted FSD")
    if(any(is.na(xy)))
        stop("missing value(s) in 'x'")
    if(any(!is.finite(xy)))
        stop("infinite value(s) in 'x'")

    Psi.rng <- .est.Psi(xy)
    x.fit <- x$x.fit
    y.fit <- y$x.fit
    options(warn=-1)

    ##print(c(length(x.fit$y),length(y.fit$y)))
          
    res <- optimize(.llkCopula,
                    interval=c(Psi.rng$min,Psi.rng$max),
                    Fx=x.fit$y,Fy=y.fit$y,cnts=xy)
    options(warn=0)
    Psi <- res$minimum
    x0 <- x.fit$x
    sele.x <- is.finite(x0)
    y0 <- y.fit$x
    sele.y <- is.finite(y0)
    x0 <- x0[sele.x]; y0 <- y0[sele.y]
    z0 <- .pdfCopula(Psi,x.fit$y[sele.x],y.fit$y[sele.y],
                     x.fit$y2[sele.x],y.fit$y2[sele.y])
    out <- structure(
        list(x.fit = x.fit,
             y.fit = y.fit,
             Psi = Psi,
             Psi.rng = Psi.rng$summ,
             x=x0,y=y0,z=z0),
        class="FSD")
}

.est.Psi <- function(x){
    x[is.na(x)] <- 0
    .funPsi <- function(ind,xmat,I,J){
        i <- ind[1]; j <- ind[2]
        nI <- sum(xmat[(i+1):I,(j+1):J])
        nIII <- sum(xmat[1:i,1:j])
        nIV <- sum(xmat[1:i,(j+1):J])
        nII <- sum(xmat[(i+1):I,1:j])
        nI <- ifelse(nI==0,1,nI)
        nII <- ifelse(nII==0,1,nII)
        nIII <- ifelse(nIII==0,1,nIII)
        nIV <- ifelse(nIV==0,1,nIV)
        (nIII*nII)/(nI*nIV)
    }
    tmp <- dim(x)
    I <- tmp[1];J <- tmp[2]
    ij <- expand.grid(1:(I-1),1:(J-1))
    out <- apply(ij,1,.funPsi,xmat=x,I=I,J=J)
    res <- quantile(out,c(0,.025,.5,.0975,1))
    list(min=min(out), max=max(out),summ=res)
}

lines.FSD <- function(x,grid.size,...){
    out <- NULL
    if(is.null(x$y.fit)){
        y <- x$x.fit
        if(missing(grid.size)){
            x0 <- as.numeric(y$lx)
            y0 <- as.numeric(y$ly)
            sele <- y0 > 0
            x0 <- x0[sele]
            y0 <- y0[sele]
        }else{
            n <- as.numeric(grid.size)[1]
            if(n<10) stop("grid size too small")
            out1 <- .mypdf(y,n)
            x0 <- as.numeric(out1$x)
            y0 <- as.numeric(out1$y)
            sele <- y0 > 0
            x0 <- x0[sele]
            y0 <- y0[sele]
        }
        lines(x0, y0,...)
        out <- list(x=x0, y=y0)
    }
    invisible(out)
}
    
plot.FSD <- function(x,grid.size,Psi,plot.new=TRUE,...){
    if(is.null(x$y.fit)){
        if(missing(grid.size)){
            x0 <- as.numeric(x$x.fit$x)
            y0 <- as.numeric(x$x.fit$y)
            sele <- y0>0
            x0 <- x0[sele]
            y0 <- y0[sele]
        }else{
            n <- as.numeric(grid.size)[1]
            if(n<10) stop("grid size too small")
            out1 <- .mypdf(x$x.fit,n)
            x0 <- as.numeric(out1$x)
            y0 <- as.numeric(out1$y)
            sele <- y0>0
            x0 <- x0[sele]
            y0 <- y0[sele]
        }
        if(plot.new){
            plot(x0, y0,...)
        }else{
            lines(x0, y0,...)
        }
        
        out <- list(x=x0, y=y0)
    }else{
        if(missing(grid.size)){
            if(plot.new){
                contour(x=x$x,y=x$y,z=x$z,...)
            }else{
                contour(x=x$x,y=x$y,z=x$z,add=TRUE,...)
            }
            out <- list(x=x$x, y=x$y,z=x$z)
        }else{
            n <- as.numeric(grid.size)[1]
            if(n<10) stop("grid size too small")
            outx <- .mycdf(x$x.fit,n)
            outy <- .mycdf(x$y.fit,n)
            if(missing(Psi)){
                Psi <- x$Psi
            }else{
                stopifnot(Psi>0||!is.finite(Psi))
            }
            z0 <- .pdfCopula(Psi,outx$f,outy$f,outx$F,outy$F)
            
            x0 <- outx$x
            n <- length(x0)
            if(!is.finite(x0[n]))
                x0[n] <- 2*x0[n-1]
                
            y0 <- outy$x
            n <- length(y0)
            if(!is.finite(y0[n]))
                y0[n] <- 2*y0[n-1]
            
            if(plot.new){
                contour(x=x0,y=y0,z=z0,...)
            }else{
                contour(x=x0,y=y0,z=z0,add=TRUE,...)
            }
            out <- list(x=x0, y=y0,z=z0)
        }
    }
    invisible(out)
}

.mypdf <- function(x,n){
    xbrks <- x$x
    k <- length(xbrks)
    m <- ceiling(n/(k-1))
    ## define grid #############
    a <- xbrks[1]
    if(a<=0) a <- 0.1
    b <- xbrks[k]
    if(!is.finite(b))
        b <- xbrks[k-1]*2 - xbrks[k-2] 
    x0 <- exp(seq(log(a),log(b),length=n))
    
    if(x$dist=="weibull")
        y0 <- dweibull(x0,x$pars[1],x$pars[2])
    else if(x$dist=="gpd")
        y0 <- .dGPD(x0,x$pars[1],x$pars[2],x$pars[3])
    else if(x$dist=="ewd")
        y0 <- .dEWD(x0,x$pars[1],x$pars[2],x$pars[3])
    else if(x$dist=="exp")
        y0 <- dexp(x0,x$pars[1])
    else if(x$dist=="ep"||x$dist=="mep"||x$dist=="mixed")
        y0 <- .dMEP(x0,x$pars[1],x$pars[2],x$pars[3])
    else if(x$dist=="pd"||x$dist=="pareto")
        y0 <- .dPD(x0,x$pars[1],x$pars[2])
    else if(x$dist=="kde"||x$dist=="bkde")
        y0 <- .dKDE(x0,x$pars[1],x$breaks,x$counts)
    else if(x$dist=="spline")
        y0 <- .dspline(x0,x$breaks,x$counts)
    else if(x$dist=="lnorm"||x$dist=="lognormal")
        y0 <- dlnorm(x0,x$pars[1],x$pars[2])
    else
        y0 <- .dgld(x0,x$pars)
    
    y0[is.na(y0)] <- 0
    sele <- is.finite(y0)&(y0>0)
    x0 <- x0[sele]; y0 <- y0[sele]
    list(x=x0,y=y0)
}

.mycdf <- function(x,n){
    xbrks <- x$x
    k <- length(xbrks)
    m <- ceiling(n/(k-1))
    x0 <- NULL
    M <- xbrks[k]
    if(!is.finite(M))
        xbrks[k] <- xbrks[k-1]*5 - xbrks[k-2] 
    for(i in 1:(k-1)){
        delta <- (xbrks[i+1]-xbrks[i])/(m-1)
        tmp <- seq(xbrks[i],xbrks[i+1]-delta,length=m)
        x0 <- c(x0,tmp)
    }
    if(max(x0)<M) x0 <- c(x0,M)
    if(x$dist=="weibull"){
        y0 <- dweibull(x0,x$pars[1],x$pars[2])
        y1 <- pweibull(x0,x$pars[1],x$pars[2])
    }else if(x$dist=="gpd"){
        y0 <- .dGPD(x0,x$pars[1],x$pars[2],x$pars[3])
        y1 <- .pGPD(x0,x$pars[1],x$pars[2],x$pars[3])
    }else if(x$dist=="exp"){
        y0 <- dexp(x0,x$pars[1])
        y1 <- pexp(x0,x$pars[1])
    }else if(x$dist=="ewd"){
        y0 <- .dEWD(x0,x$pars[1],x$pars[2],x$pars[3])
        y1 <- .pEWD(x0,x$pars[1],x$pars[2],x$pars[3])
    }else if(x$dist=="pd"||x$dist=="pareto"){
        y0 <- .dPD(x0,x$pars[1],x$pars[2])
        y1 <- .pPD(x0,x$pars[1],x$pars[2])
    }else if(x$dist=="kde"||x$dist=="bkde"){
        y0 <- .dKDE(x0,x$pars[1],x$breaks,x$counts)
        y1 <- .pKDE(x0,x$pars[1],x$breaks,x$counts)
    }else if(x$dist=="lnorm"||x$dist=="lognormal"){
        y0 <- dlnorm(x0,x$pars[1],x$pars[2])
        y1 <- plnorm(x0,x$pars[1],x$pars[2])
    }else if(x$dist=="spline"){
        y0 <- .dspline(x0,x$breaks,x$counts)
        y1 <- .pspline(x0,x$breaks,x$counts)
    }else{
        y0 <- .dgld(x0,x$pars)
        y1 <- .pgld(x0,x$pars)
    }
    
    y0[is.na(y0)] <- 0
    y1[is.na(y1)] <- 0
    
    list(x=x0,f=y0,F=y1)
}

print.FSD <- function(x,...){
    cat("Fitting Firm Size Distribution(s)...")
    res <- x$x.fit
    cat("\n x: ", res$dist)
    cat("\n fitted parameters: ",
        signif(as.numeric(res$pars), digits=3))
    if(!is.null(res$aux)){
        cat("\nGoodness-of-fit:\n")
        print(signif(res$aux,digits=3))
    }
    
    if(!is.null(x$y.fit)){
        res <- x$y.fit
        cat("\n\n y: ", res$dist)
        cat("\n fitted parameters: ",
            signif(as.numeric(res$pars)),digits=3)
        if(!is.null(res$aux)){
            cat("\nGoodness-of-fit:\n")
            print(signif(res$aux,digits=3))
        }
        
        cat("\n Psi :=", x$Psi)
        cat("\n Summary of Psi :\n")
        print(x$Psi.rng)
    }
    cat("\n\n")
}

.fit.FSD1 <- function(x, breaks, dist="weibull"){
    res <- NULL
    dist <- match.arg(tolower(dist),
                      c("weibull","pareto","gpd","pd",
                        "ewd","gld","gld1","gld2",
                        "mep","ep", "mixed","exp",
                        "lnorm","lognormal","kde",
                        "bkde","spline","norm"))
    if(inherits(x,"histogram")){
        brks <- x$breaks
        cnts <- x$counts
        xh <- x
    }else if(inherits(x,"bdata")){
        brks <- x$breaks
        cnts <- x$freq
        xh <- x$xhist
    }else{
        xh <- hist(x, plot=FALSE)
        brks <- xh$breaks
        cnts <- xh$counts
    }
    
    if(dist=="weibull")
        res <- .fit.Weibull(freq=cnts,breaks=brks)
    else if(dist=="gpd")
        res <- .fit.GPD(freq=cnts,breaks=brks)
    else if(dist=="ewd")
        res <- .fit.EWD(freq=cnts,breaks=brks)
    else if(dist=="gld1")
        res <- .fit.GLD(freq=cnts,breaks=brks,type=1)
    else if(dist=="exp")
        res <- .fit.EXP(freq=cnts,breaks=brks)
    else if(dist=="gld2"||dist=="gld")
        res <- .fit.GLD(freq=cnts,breaks=brks,type=2)
    else if(dist=="pareto"||dist=="pd")
        res <- .fit.PD(freq=cnts,breaks=brks)
    else if(dist=="kde"||dist=="bkde")
        res <- .fit.KDE(freq=cnts,breaks=brks)
    else if(dist=="ep"||dist=="mep"||dist=="mixed")
        res <- .fit.MEP(freq=cnts,breaks=brks)
    else if(dist=="lnorm"||dist=="lognormal")
        res <- .fit.lnorm(freq=cnts,breaks=brks)
    else if(dist=="norm")
        res <- .fit.norm(freq=cnts,breaks=brks)
    else if(dist=="spline")
        res <- .fit.spline(freq=cnts,breaks=brks)
    else
        stop("distribution type not supported")
    
    list(xhist = xh,
         dist = dist,
         size = sum(cnts),
         breaks = brks,
         counts = cnts,
         pars = as.numeric(res$par),
         Dn.Zipf = res$Dn.Zipf,
         Dn = res$Dn,
         aux = res$aux,
         AIC = res$AIC,
         BIC = res$BIC,
         AICc = res$AICc,
         x=res$x,y=res$y,y2=res$y2,
         lx=res$lx,ly=res$ly)
}

.fitFSD <- function(x,type){
    if(type=="exp"){
        out <- .fitFirmAge(x)
    }else if(type=="pareto"){
        out <- .fitFirmSize(x)
    }else{
        out <- .fitFirmSize(x)
    }
    out
}

.fitFirmAge <- function(x){
    .funG <- function(l,x,freq){
        k <- length(x)
        f <- freq
        x <- x[-1]
        F <- pexp(x, l)
        dF <- diff(c(0,F,1))
        ldF <- log(dF)
        ldF[!is.finite(ldF)] <- -1000
        -sum(f*ldF)
    }
    lmd <- optimize(.funG,c(0.001,0.3),
                    x = x$ll,
                    freq = x$freq)$minimum
    list(lambda=lmd)
}

.fitFirmSize <- function(x){
    if(x$ll[1] > 0)
        xm <- x$ll[1]
    else
        xm <- x$ll[2]
    out <- fit.Pareto(x,xm=xm)
}

####################################################
.fit.norm <- function(freq,breaks){
    k <- length(freq)
    n <- sum(freq)
    Fn <- cumsum(freq)/(n+1)
    ## find initial estimates using LSE
    if(is.finite(breaks[k+1])){
        xi <- breaks[-1]
        zi <- qnorm(Fn)
    }else{
        xi <- breaks[-c(1,k+1)]
        zi <- qnorm(Fn[-k])
    }
    lm0 <- lm(xi~zi)
    shat <- lm0$coef[[2]]
    muhat <- lm0$coef[[1]]
    par1 <- c(muhat,shat)
    pars <- par1

    ##print(pars)
    
    x0 <- breaks[-c(1,k+1)]
    Fx <- pnorm(x0, par1[1], par1[2])
    llk <- sum(log(diff(c(0,Fx,1)))*freq)
    AIC1 <- .AIC(llk,2,n)
    Dn.Zipf1 <- .funCOD(freq,Fx)
    L1 <- .funL2(freq,Fx)
    ## find MLE
    options(warn=-1)
    out <- optim(par1, .funNORM, 
                 freq=freq, brks=breaks)
    if(out$conv==0){
        par2 <- out$par
        pars <- par2
        x0 <- breaks[-c(1,k+1)]
        Fx <- pnorm(x0, par2[1], par2[2])
        llk <- sum(log(diff(c(0,Fx,1)))*freq)
        AIC2 <- .AIC(llk,2,n)
        Dn.Zipf2 <- .funCOD(freq,Fx)
        L2 <- .funL2(freq,Fx)
    }else{
        par2 <- c(NA, NA)
        AIC2 <- list(AIC=NA, BIC=NA, AICc=NA)
        Dn.Zipf2 <- NA
        L2 <- NA
    }
    options(warn=0)

    x <- breaks
    k <- length(x)
    
    y <- dnorm(x, pars[1],pars[2])
    y2 <- pnorm(x, pars[1],pars[2])
    a <- x[1];b <- rev(x)[1]
    if(a<=0) a <- 0.1
    if(!is.finite(b)) b <- 2*rev(x)[2]-rev(x)[3]
    lx <- exp(seq(log(a),log(b),length=401))
    ly <- dnorm(lx, pars[1],pars[2])
                 
    Fx <- y2[-c(1,k)]
    llk <- sum(log(diff(c(0,Fx,1)))*freq)
    AIC0 <- .AIC(llk,1,n)
    Dn.Zipf0 <- .funCOD(freq,Fx)
    L0 <- .funL2(freq,Fx)
    
    
    tmp <- data.frame(mean.log = c(par1[1],par2[1]),
                      sd.log = c(par1[2],par2[2]),
                      Dn.Zipf = c(Dn.Zipf1,Dn.Zipf2),
                      Dn = c(L1,L2),
                      AIC = c(AIC1$AIC,AIC2$AIC),
                      BIC = c(AIC1$BIC,AIC2$BIC),
                      AICc = c(AIC1$AICc,AIC2$AICc))
    
    rownames(tmp) <- c("LSE","MLE")
  
    list(pars=pars, aux=tmp,
         Dn.Zipf=Dn.Zipf0,Dn=L0,
         AIC=AIC0$AIC,
         BIC=AIC0$BIC,
         AICc=AIC0$AICc,
         x=x,y=y,y2=y2,lx=lx,ly=ly)
}

.funNORM <- function(pars,freq,brks){
    n <- length(brks)
    stopifnot(length(freq)+1==n)
    brks <- brks[-c(1,n)]
    Fx = pnorm(brks, pars[1], pars[2])
    if(any(is.na(Fx))){
        llk <- 1.0e+9
    }else{
        Px <- diff(c(0,Fx,1))
        sele = Px <= 1.0e-6
        delta = 0 #penalty
        if(any(sele)){
            Px = Px[!sele]
            freq = freq[!sele]
            delta = 1.0e+6 * sum(freq[sele])
        }
        llk <- delta-sum(log(Px)*freq)
    }
    if(!is.finite(llk)) llk <- .Machine$double.xmax
    llk
}

####################################################
.fit.lnorm <- function(freq,breaks){
    k <- length(freq)
    n <- sum(freq)
    Fn <- cumsum(freq)/(n+1)
    ## find initial estimates using LSE
    if(is.finite(breaks[k+1])){
        xi <- log(breaks[-1])
        zi <- qnorm(Fn)
    }else{
        xi <- log(breaks[-c(1,k+1)])
        zi <- qnorm(Fn[-k])
    }
    lm0 <- lm(xi~zi)
    shat <- lm0$coef[[2]]
    muhat <- lm0$coef[[1]]
    par1 <- c(muhat,shat)
    pars <- par1

    ##print(pars)
    
    x0 <- breaks[-c(1,k+1)]
    Fx <- plnorm(x0, par1[1], par1[2])
    llk <- sum(log(diff(c(0,Fx,1)))*freq)
    AIC1 <- .AIC(llk,2,n)
    Dn.Zipf1 <- .funCOD(freq,Fx)
    L1 <- .funL2(freq,Fx)
    ## find MLE
    options(warn=-1)
    out <- optim(par1, .funLNORM, 
                 freq=freq, brks=breaks)
    if(out$conv==0){
        par2 <- out$par
        pars <- par2
        x0 <- breaks[-c(1,k+1)]
        Fx <- plnorm(x0, par2[1], par2[2])
        llk <- sum(log(diff(c(0,Fx,1)))*freq)
        AIC2 <- .AIC(llk,2,n)
        Dn.Zipf2 <- .funCOD(freq,Fx)
        L2 <- .funL2(freq,Fx)
    }else{
        par2 <- c(NA, NA)
        AIC2 <- list(AIC=NA, BIC=NA, AICc=NA)
        Dn.Zipf2 <- NA
        L2 <- NA
    }
    options(warn=0)

    x <- breaks
    k <- length(x)
    
    y <- dlnorm(x, pars[1],pars[2])
    y2 <- plnorm(x, pars[1],pars[2])
    a <- x[1];b <- rev(x)[1]
    if(a<=0) a <- 0.1
    if(!is.finite(b)) b <- 2*rev(x)[2]-rev(x)[3]
    lx <- exp(seq(log(a),log(b),length=401))
    ly <- dlnorm(lx, pars[1],pars[2])
                 
    Fx <- y2[-c(1,k)]
    llk <- sum(log(diff(c(0,Fx,1)))*freq)
    AIC0 <- .AIC(llk,1,n)
    Dn.Zipf0 <- .funCOD(freq,Fx)
    L0 <- .funL2(freq,Fx)
    
    
    tmp <- data.frame(mean.log = c(par1[1],par2[1]),
                      sd.log = c(par1[2],par2[2]),
                      Dn.Zipf = c(Dn.Zipf1,Dn.Zipf2),
                      Dn = c(L1,L2),
                      AIC = c(AIC1$AIC,AIC2$AIC),
                      BIC = c(AIC1$BIC,AIC2$BIC),
                      AICc = c(AIC1$AICc,AIC2$AICc))
    
    rownames(tmp) <- c("LSE","MLE")
  
    list(pars=pars, aux=tmp,
         Dn.Zipf=Dn.Zipf0,Dn=L0,
         AIC=AIC0$AIC,
         BIC=AIC0$BIC,
         AICc=AIC0$AICc,
         x=x,y=y,y2=y2,lx=lx,ly=ly)
}

.funLNORM <- function(pars,freq,brks){
    n <- length(brks)
    stopifnot(length(freq)+1==n)
    brks <- brks[-c(1,n)]
    Fx = plnorm(brks, pars[1], pars[2])
    if(any(is.na(Fx))){
        llk <- 1.0e+9
    }else{
        Px <- diff(c(0,Fx,1))
        sele = Px <= 1.0e-6
        delta = 0 #penalty
        if(any(sele)){
            Px = Px[!sele]
            freq = freq[!sele]
            delta = 1.0e+6 * sum(freq[sele])
        }
        llk <- delta-sum(log(Px)*freq)
    }
    if(!is.finite(llk)) llk <- .Machine$double.xmax
    llk
}

####################################################
.fit.EXP <- function(freq,breaks){
    n <- sum(freq)
    Fn <- cumsum(freq)/(n+1)
    k <- length(Fn)
    pars <- mean(-log(1-Fn[-k])/breaks[-c(1,k+1)])
    
    x <- breaks
    k <- length(x)

    y <- dexp(x, pars)
    y2 <- pexp(x, pars)
    a <- x[1];b <- rev(x)[1]
    if(a<=0) a <- 0.1
    if(!is.finite(b)) b <- 2*rev(x)[2]-rev(x)[3]
    lx <- exp(seq(log(a),log(b),length=401))
    ly <- dexp(lx, pars)

    Fx <- y2[-c(1,k)]
    llk <- sum(log(diff(c(0,Fx,1)))*freq)
    AIC0 <- .AIC(llk,1,n)
    Dn.Zipf0 <- .funCOD(freq,Fx)
    L0 <- .funL2(freq,Fx)
    
    
    tmp <- data.frame(lambda = pars[1],
                      Dn.Zipf=Dn.Zipf0,
                      Dn = L0,
                      AIC = AIC0$AIC,
                      BIC = AIC0$BIC,
                      AICc = AIC0$AICc)
    
    rownames(tmp) <- c("LS")
  
    list(pars=pars, aux=tmp,
         Dn.Zipf=Dn.Zipf0,Dn=L0,
         AIC=AIC0$AIC,
         BIC=AIC0$BIC,
         AICc=AIC0$AICc,
         x=x,y=y,y2=y2,lx=lx,ly=ly)
}

####################################################
.dKDE <- function(x0, h.adj,x,f){
    k <- length(x)
    lb <- x[-k]
    ub <- x[-1]
    xc <- (lb+ub)/2
    h <- (ub - lb)*h.adj
    n <- length(x0)
    if(any(h<=0))
        stop("invalid break point(s)")
    res <- .Fortran(.F_dKDE,
                    y=as.double(x0),
                    as.integer(n),
                    as.double(xc),
                    as.double(h),
                    as.double(f),
                    as.integer(k-1))
    res$y
}

.pKDE <- function(x0, h.adj,x,f){
    k <- length(x)
    lb <- x[-k]
    ub <- x[-1]
    xc <- (lb+ub)/2
    h <- (ub - lb)*h.adj
    n <- length(x0)
    if(any(h<=0))
        stop("invalid break point(s)")
    res <- .Fortran(.F_pKDE,
                    y=as.double(x0),
                    as.integer(n),
                    as.double(xc),
                    as.double(h),
                    as.double(f),
                    as.integer(k-1))
    res$y
}

.fit.KDE <- function(freq,breaks){
    x <- breaks
    k <- length(x)
    if(!is.finite(x[1]))
        stop("the lower boundary of the first class must be finite")
    if(!is.finite(x[k]))
        stop("the upper boundary of the last class must be finite")
    pars <- 1.0
    y <- .dKDE(x, pars,breaks,freq)
    y2 <- .pKDE(x, pars,breaks,freq)
    
    a <- x[1];b <- rev(x)[1]
    lx <- seq(a,b,length=401)
    ly <- .dKDE(lx, pars,breaks,freq)

    Fx <- y2[-c(1,k)]
    llk <- sum(log(diff(c(0,Fx,1)))*freq)
    n <- sum(freq)
    AIC0 <- .AIC(llk,1,n)
    Dn.Zipf0 <- .funCOD(freq,Fx)
    L0 <- .funL2(freq,Fx)
    
    
    tmp <- data.frame(h.adj = pars,
                      Dn.Zipf=Dn.Zipf0,
                      Dn = L0,
                      AIC = AIC0$AIC,
                      BIC = AIC0$BIC,
                      AICc = AIC0$AICc)
    
    rownames(tmp) <- c("KDE")
  
    list(pars=pars, aux=tmp,
         Dn.Zipf=Dn.Zipf0,
         Dn=L0,
         AIC=AIC0$AIC,
         BIC=AIC0$BIC,
         AICc=AIC0$AICc,
         x=x,y=y,y2=y2,lx=lx,ly=ly)
}

####################################################
.fit.GLD <- function(freq,breaks,type=2){
    xh <- binning(counts=freq, breaks=breaks)
    n <- sum(freq)
    Fn <- cumsum(freq)/n
    k <- length(Fn)
    i1 <- c(1:4,(k-1))
    i2 <- c(1,(k-4):(k-1))
    i3 <- c(1,2,floor(k/2), (k-2),(k-1))

    isele <- i1
    qtls <- breaks[isele+1]
    qlvls <- Fn[isele]
    lbound <- breaks[1]
    out <- fit.GLD(xh, qtl=qtls, qtl.levels=qlvls,lbound=lbound)
    gld1.mop <- .funGLD(out$pars,breaks,freq)
    gld1.mle <- .funGLD(out$aux,breaks,freq)
    isele <- i2
    qtls <- breaks[isele+1]
    qlvls <- Fn[isele]
    lbound <- breaks[1]
    out <- fit.GLD(xh, qtl=qtls, qtl.levels=qlvls,lbound=lbound)
    gld2.mop <- .funGLD(out$pars,breaks,freq)
    gld2.mle <- .funGLD(out$aux,breaks,freq)
    isele <- i3
    qtls <- breaks[isele+1]
    qlvls <- Fn[isele]
    lbound <- breaks[1]
    out <- fit.GLD(xh, qtl=qtls, qtl.levels=qlvls,lbound=lbound)
    gld3.mop <- .funGLD(out$pars,breaks,freq)
    gld3.mle <- .funGLD(out$aux,breaks,freq)

    i <- 1
    lmd1 <- c(gld1.mop$pars[i],gld1.mle$pars[i],
              gld2.mop$pars[i],gld2.mle$pars[i],
              gld3.mop$pars[i],gld3.mle$pars[i])
    i <- 2
    lmd2 <- c(gld1.mop$pars[i],gld1.mle$pars[i],
              gld2.mop$pars[i],gld2.mle$pars[i],
              gld3.mop$pars[i],gld3.mle$pars[i])
    i <- 3
    lmd3 <- c(gld1.mop$pars[i],gld1.mle$pars[i],
              gld2.mop$pars[i],gld2.mle$pars[i],
              gld3.mop$pars[i],gld3.mle$pars[i])
    i <- 4
    lmd4 <- c(gld1.mop$pars[i],gld1.mle$pars[i],
              gld2.mop$pars[i],gld2.mle$pars[i],
              gld3.mop$pars[i],gld3.mle$pars[i])
    Dn.Zipf <- c(gld1.mop$Dn.Zipf,gld1.mle$Dn.Zipf,
             gld2.mop$Dn.Zipf,gld2.mle$Dn.Zipf,
             gld3.mop$Dn.Zipf,gld3.mle$Dn.Zipf)
    Dn <- c(gld1.mop$Dn,gld1.mle$Dn,
            gld2.mop$Dn,gld2.mle$Dn,
            gld3.mop$Dn,gld3.mle$Dn)
    AIC <- c(gld1.mop$AIC,gld1.mle$AIC,
             gld2.mop$AIC,gld2.mle$AIC,
             gld3.mop$AIC,gld3.mle$AIC)
    BIC <- c(gld1.mop$BIC,gld1.mle$BIC,
             gld2.mop$BIC,gld2.mle$BIC,
             gld3.mop$BIC,gld3.mle$BIC)
    AICc <- c(gld1.mop$AICc,gld1.mle$AICc,
              gld2.mop$AICc,gld2.mle$AICc,
              gld3.mop$AICc,gld3.mle$AICc)
    ## choose the best estimate using BIC
    isele <- which(BIC==min(BIC))[1]
    if(isele==1)
        out <- gld1.mop
    else if(isele==2)
        out <- gld1.mle
    else if(isele==3)
        out <- gld2.mop
    else if(isele==4)
        out <- gld2.mle
    else if(isele==5)
        out <- gld3.mop
    else
        out <- gld3.mle
  
    tmp <- data.frame(lmd1 = lmd1,
                      lmd2 = lmd2,
                      lmd3 = lmd3,
                      lmd4 = lmd4,
                      Dn.Zipf = Dn.Zipf,
                      Dn = Dn,
                      AIC = AIC,
                      BIC = BIC,
                      AICc = AICc)
    rownames(tmp) <- c("MOP.lo","MLE.lo",
                       "MOP.hi","MLE.hi",
                       "MOP.med","MLE.med")
  
    res <- list(pars=out$pars,aux=tmp,
                Dn.Zipf=Dn.Zipf[isele],Dn=Dn[isele],
                AIC=AIC[isele],
                BIC=BIC[isele],
                AICc=AICc[isele],
                x=out$x, y=out$y,y2=out$y2,
                lx=out$lx, ly=out$ly)
}

.funGLD <- function(pars, breaks, freq){
    n <- sum(freq)
    x <- breaks
    k <- length(x)
    y <- .dgld(x, pars)
    y2 <- .pgld(x, pars)
    a <- x[1];b <- rev(x)[1]
    if(a<=0) a <- 0.1
    if(!is.finite(b)) b <- 2*rev(x)[2]-rev(x)[3]
    lx <- exp(seq(log(a),log(b),length=401))
    ly <- .dgld(lx, pars)

    Fx <- y2[-c(1,k)]
    llk <- sum(log(diff(c(0,Fx,1)))*freq)
    AIC0 <- .AIC(llk,4,n)
    Dn.Zipf0 <- .funCOD(freq,Fx)
    L0 <- .funL2(freq,Fx)
    list(pars=pars,x=x,y=y,y2=y2,lx=lx,ly=ly,
         Dn.Zipf=Dn.Zipf0, Dn=L0,AIC=AIC0$AIC,
         BIC=AIC0$BIC,AICc=AIC0$AICc)
}

####################################################
.pGPD <- function(x,mu,sigma,xi){
  stopifnot(is.numeric(x))
  if(any(is.na(x)))
    stop("missing value(s) in 'x'")
  sapply(x,FUN=.subpGPD,mu=mu,sigma=sigma,xi=xi)
}

.subpGPD <- function(x,mu,sigma,xi){
  res = 0
  stopifnot(sigma>0)
  if(is.finite(x)){
    z <- (x-mu)/sigma
    if(z>0){
      if(xi==0){
        res = 1-exp(-z)
      }else{
        res = 1-(1+xi*z)^(-1/xi)        
      }
    }
  }else
  res = 1
  if(!is.finite(res)) res <- 0
  if(is.na(res)) res <- 0
  
  res
}

.dGPD <- function(x,mu,sigma,xi){
  stopifnot(is.numeric(x))
  if(any(is.na(x)))
    stop("missing value(s) in 'x'")
  sapply(x,FUN=.subdGPD,mu=mu,sigma=sigma,xi=xi)
}

.subdGPD <- function(x,mu,sigma,xi){
  res = 0
  stopifnot(sigma>0)
  if(is.finite(x)){
    z <- (x-mu)/sigma
    if(z>0){
      if(xi==0){
        res = exp(-z)
      }else{
        res = (1+xi*z)^(-1/xi-1)/sigma        
      }
    }
  }else
  res = 0
  if(!is.finite(res)) res <- 0
  if(is.na(res)) res <- 0
  res
}

.funGPD <- function(pars,freq,brks){
  n <- length(brks)
  stopifnot(length(freq)+1==n)
  brks <- brks[-c(1,n)]
  Fx = .pGPD(brks, pars[1], pars[2],pars[3])
  if(any(is.na(Fx))){
    llk <- 1.0e+9
  }else{
    Px <- diff(c(0,Fx,1))
    sele = Px <= 1.0e-6
    delta = 0 #penalty
    if(any(sele)){
      Px = Px[!sele]
      freq = freq[!sele]
      delta = 1.0e+6 * log(sum(freq[sele]))
    }
    llk <- delta-sum(log(Px)*freq)
  }
  if(!is.finite(llk)) llk <- .Machine$double.xmax
  llk
}

.funGPD2 <- function(pars,freq,brks,p1){
  n <- length(brks)
  stopifnot(length(freq)+1==n)
  brks <- brks[-c(1,n)]
  ##print(c(p1,pars))
  Fx = .pGPD(brks, p1, pars[1],pars[2])
  u <- runif(1,.8,.99)
  tol <- 9e15 #.Machine$double.xmax
  if(any(is.na(Fx))){
      llk <- tol*u
  }else{
      Px <- diff(c(0,Fx,1))
      sele = Px <= 1.0e-6
      delta = 0 #penalty
      if(any(sele)){
          if(mean(sele)<1){
              Px = Px[!sele]
              freq = freq[!sele]
              delta = 1.0e+6 * log(sum(freq[sele]))
              llk <- delta-sum(log(Px)*freq)
          }else{
              llk <- tol*u
          }
      }else{
          llk <- delta-sum(log(Px)*freq)
      }
  }
  if(!is.finite(llk)) llk <- tol*u
  llk
}


.GPDBIC <- function(xi,brks,counts,mu){
    out <- 9e100
    N <- sum(counts)
    x <- brks[-1]
    k <- length(counts)
    Fn <- cumsum(counts)/(N+1)
    if(xi<0){
        sigs <- xi*(x-mu)/((1-Fn)^(-xi)-1)
    }else{
        sigs <- xi*(x-mu)*((1-Fn)^(xi)-1)
    }
    print(sigs)
    sele <- is.finite(sigs)
    if(any(sele)){
        sigs <- sigs[sele]
        if(any(sigs>0)){
            sig <- mean(sigs[sigs>0])
            par0 <- c(mu,sig,xi)
            Fx <- .pGPD(x[-k], par0[1], par0[2], par0[3])
            llk <- sum(log(diff(c(0,Fx,1)))*counts)
            AIC0 <- .AIC(llk,2,N)
            out <- AIC0$BIC
        }
    }
    out
}

.fit.GPD <- function(freq,breaks){
    n <- length(freq)
    if(length(breaks) != n +1)
        stop("length of the breaks points not match")
    brks <- breaks[-c(1,n+1)]
    N <- sum(freq)

    ##xis <- runif(100,-2.5,-1)
    ##res <- sapply(xis,FUN=.GPDBIC,brks=breaks,counts=freq,mu=breaks[1])
    ##print(res)
    tmp <- NULL
    brks2 <- breaks[-1]
    Fn <- cumsum(freq)/(N+1)
    k <- length(brks2)
    if(!is.finite(brks2[k])){
        brks2 <- brks2[-k]
        Fn <- Fn[-k]
        k <- length(brks2)
    }
    xk <- brks2[k]
    Fk <- Fn[k]
    yi <- NULL
    xi <- NULL
    for(i in 1:(k-1)){
        tmp <- log((1-Fn[i])/(1-Fk))
        yi <- c(yi,tmp)
        tmp <- log(brks2[i]/xk)
        xi <- c(xi,-tmp)
    }
    ##print(cbind(tmp,ks))
    ##plot(tmp~ks)
    lmout <- lm(xi~yi)
    xi <- lmout$coef[[2]]
    ##print(lmout$coef)
    if(xi<0){
        sig <- xi*xk/((1-Fk)^(-xi)-1)
    }else{
        sig <- xi*xk/(1/((1-Fk)^(xi))-1)
    }
    
    par0 <- c(NA,NA,NA)
    AIC0 <- list(AIC=NA,BIC=NA,AICc=NA)
    Dn.Zipf0 <- NA
    L0 <- NA

    options(warn=-1)
    phats <- c(sig,xi)



    print(phats)


    
    out <- optim(phats,.funGPD2, 
                 method= "L-BFGS-B",
                 freq=freq, 
                 brks=breaks,
                 ##lower=c(0.00001, abs(xi)*.5),
                 ##upper=c(2*sig,abs(xi)*2),
                 p1=breaks[1])

    
    if(out$conv==0){
        par0 <- c(breaks[1],out$par)
        Fx <- .pGPD(brks, par0[1], par0[2], par0[3])
        llk <- sum(log(diff(c(0,Fx,1)))*freq)
        AIC0 <- .AIC(llk,2,N)
        Dn.Zipf0 <- .funCOD(freq,Fx)
        L0 <- .funL2(freq,Fx)

    }
    options(warn=0)

    ##m1 <- min(0.01,b0/10)
    ##m2 <- b0
    ## finding MLE
    ##phats <- c((m1+m2)*.5,1.2,1.3)
    ##out <- optim(phats,.funGPD,
    ##method= "L-BFGS-B",
    ##freq=freq, 
    ##brks=breaks,
    ##lower=c(m1,0.00001, -5),
    ##upper=c(m2,25,10))
    ##if(out$conv==0){
    ##par0 <- out$par
    ##Fx <- .pGPD(brks, par0[1], par0[2], par0[3])
    ##llk <- sum(log(diff(c(0,Fx,1)))*freq)
    ##AIC0 <- .AIC(llk,2,N)
    ##Dn.Zipf0 <- .funCOD(freq,Fx)
    ##L0 <- .funL2(freq,Fx)

    x <- breaks
    k <- length(x)
    y <- .dGPD(x, par0[1], par0[2], par0[3])
    y2 <- .pGPD(x, par0[1], par0[2], par0[3])
    a <- x[1];b <- rev(x)[1]
    if(a<=0) a <- 0.1
    if(!is.finite(b)) b <- 2*rev(x)[2]-rev(x)[3]
    lx <- exp(seq(log(a),log(b),length=401))
    ly <- .dGPD(lx, par0[1], par0[2], par0[3])

    tmp <- data.frame(
        mu = c(par0[1]),
        sigma = c(par0[2]),
        xi = c(par0[3]),
        Dn.Zipf=c(Dn.Zipf0),
        Dn = c(L0),
        AIC = c(AIC0$AIC),
        BIC = c(AIC0$BIC),
        AICc = c(AIC0$AICc))
  
    rownames(tmp) <- c("MLE")

    
    list(pars=par0,Dn.Zipf=Dn.Zipf0,Dn=L0,
         AIC=AIC0$AIC,
         BIC=AIC0$BIC,
         AICc=AIC0$AICc,
         x=x,y=y,y2=y2,aux=tmp,
         lx=lx,ly=ly)
}


#######################################################
.pEWD <- function(x,kappa,lambda,alpha){
  stopifnot(is.numeric(x))
  if(any(is.na(x)))
    stop("missing value(s) in 'x'")
  sapply(x,FUN=.subpEWD,kappa,lambda,alpha)
}

.subpEWD <- function(x,kappa,lambda,alpha){
  res = 0
  stopifnot(kappa>0)
  stopifnot(lambda>0)
  stopifnot(alpha>0)
  if(is.finite(x)){
    res = (1-exp(-(x/lambda)^kappa))^alpha        
  }else
    res = 1
  
  res
}

.dEWD <- function(x,kappa,lambda,alpha){
  stopifnot(is.numeric(x))
  if(any(is.na(x)))
    stop("missing value(s) in 'x'")
  sapply(x,FUN=.subdEWD,kappa,lambda,alpha)
}

.subdEWD <- function(x,kappa,lambda,alpha){
  res = 0
  stopifnot(kappa>0)
  stopifnot(lambda>0)
  stopifnot(alpha>0)
  if(is.finite(x)){
      a <- exp(-(x/lambda)^kappa)
      b <- alpha*kappa/lambda*(x/lambda)^(kappa-1)     
      res <- b*(1-a)^(alpha-1)*a
  }else
      res = 0
  
  res
}

.funEWD <- function(pars,freq,brks){
  n <- length(brks)
  stopifnot(length(freq)+1==n)
  brks <- brks[-c(1,n)]
  Fx = .pEWD(brks, pars[1], pars[2],pars[3])
  if(any(is.na(Fx))){
    llk <- 1.0e+9
  }else{
    Px <- diff(c(0,Fx,1))
    sele = Px <= 1.0e-6
    delta = 0 #penalty
    if(any(sele)){
      Px = Px[!sele]
      freq = freq[!sele]
      delta = 1.0e+6
    }
    llk <- delta-sum(log(Px)*freq)
  }
  if(!is.finite(llk)) llk <- .Machine$double.xmax
  llk
}

.funEWD2 <- function(pars,freq,brks,alp){
    n <- length(brks)
    stopifnot(length(freq)+1==n)
    brks <- brks[-c(1,n)]
    Fx = .pEWD(brks, pars[1], pars[2],alp)
    if(any(is.na(Fx))){
        llk <- 1.0e+9
    }else{
        Px <- diff(c(0,Fx,1))
        sele = Px <= 1.0e-6
        delta = 0 #penalty
        if(any(sele)){
            Px = Px[!sele]
            freq = freq[!sele]
            delta = 1.0e+6
        }
        llk <- delta-sum(log(Px)*freq)
    }
    if(!is.finite(llk)) llk <- .Machine$double.xmax
    llk
}

.fit.EWD <- function(freq,breaks){
  options(warn=-1)
  cnts <- freq
  brks <- breaks
  n <- length(cnts)
  if(length(brks) != n +1)
    stop("length of the breaks points not match")
  
  brks <- brks[-c(1,n+1)]
  N <- sum(cnts)
  Fhat <- cumsum(cnts)/(N+1)
  Fhat <- Fhat[-n]
  lx <- log(brks)
  ## fix alpha and search alpha over a fine grid. For given aplha,
  ## find initial estimates of the other two parameters using Weibull
  ## plot methods, and search for MLE numerically.

  alp1 <- seq(0.5,0.999, length=100)
  alp2 <- seq(1.001,2, length=100)
  alps <- c(alp1, alp2)

  alp0 <- 1
  ly <- log(-log(1-Fhat^(1/alp0)))
  sele <- is.finite(lx) & is.finite(ly)
  lx <- lx[sele]; ly <- ly[sele]
  out <- lm(ly~lx)
  kappa <- out$coef[[2]]
  lambda <- exp(-out$coef[[1]]/kappa)
  pars  = c(kappa, lambda, 1)
  Fx <- .pEWD(brks, pars[1], pars[2], pars[3])
  llk <- sum(log(diff(c(0,Fx,1)))*freq)
  BIC0 <- .AIC(llk,3,N)$BIC
  par0 <- pars
  
  for(alp0 in alps){
      ly <- log(-log(1-Fhat^(1/alp0)))
      sele <- is.finite(lx) & is.finite(ly)
      lx <- lx[sele]; ly <- ly[sele]
      out <- lm(ly~lx)
      kappa <- out$coef[[2]]
      lambda <- exp(-out$coef[[1]]/kappa)
      parhat  = c(kappa, lambda)
      phats <- c(parhat)
      out <- optim(phats,.funEWD2, 
                   freq=freq, 
                   brks=breaks,
                   lower=c(0.5*kappa,0.5*lambda),
                   upper=c(2*kappa,2*lambda),
                   alp=alp0)
      if(out$conv==0){
          pars <- c(out$par, alp0)
          Fx <- .pEWD(brks, pars[1], pars[2], pars[3])
          llk <- sum(log(diff(c(0,Fx,1)))*freq)
          BIC1 <- .AIC(llk,3,N)$BIC
          if(BIC1 < BIC0){
              par0 <- pars
              BIC0 <- BIC1
          }
      }
  }
  options(warn=0)
  x <- breaks
  k <- length(x)
  y <- .dEWD(x, par0[1], par0[2], par0[3])
  y2 <- .pEWD(x, par0[1], par0[2], par0[3])
  a <- x[1];b <- rev(x)[1]
  if(a<=0) a <- 0.1
  if(!is.finite(b)) b <- 2*rev(x)[2]-rev(x)[3]
  lx <- exp(seq(log(a),log(b),length=401))
  ly <- .dEWD(lx, par0[1], par0[2], par0[3])

  Fx <- .pEWD(brks, par0[1], par0[2], par0[3])
  llk <- sum(log(diff(c(0,Fx,1)))*freq)
  AIC0 <- .AIC(llk,3,N)
  Dn.Zipf0 <- .funCOD(freq,Fx)
  L0 <- .funL2(freq,Fx)


  tmp <- data.frame(kappa = par0[1],
                    lambda=par0[2],
                    alpha = par0[3],
                    Dn.Zipf=Dn.Zipf0,
                    Dn = L0,
                    AIC = AIC0$AIC,
                    BIC = AIC0$BIC,
                    AICc = AIC0$AICc)
  
  rownames(tmp) <- c("MLE")
  
  list(pars=par0, aux=tmp,Dn.Zipf=Dn.Zipf0,Dn=L0,
       AIC=AIC0$AIC,
       BIC=AIC0$BIC,
       AICc=AIC0$AICc,
       x=x,y=y,y2=y2,lx=lx,ly=ly)
}

####################################################
.fit.Weibull <- function(freq,breaks)
{
  tol <- 1e-10
  delta <- 0.0
  cnts <- freq
  brks <- breaks
  
  ## to estimate the initial values using Weibull plot method
  options(warn=-1)
  n <- length(cnts)
  if(length(brks) != n +1)
    stop("length of the breaks points not match")
  brks <- brks[-c(1,n+1)]
  N <- sum(cnts)
  Fhat <- cumsum(cnts)/N
  Fhat <- Fhat[-n]
  ly <- log(-log(1-Fhat))
  lx <- log(brks)
  sele <- is.finite(lx) & is.finite(ly)
  lx <- lx[sele]; ly <- ly[sele]
  out <- lm(ly~lx)
  kappa <- out$coef[[2]]
  lambda <- exp(-out$coef[[1]]/kappa)
  parhat  = c(kappa, lambda)
  
  .weibull.bllk <- function(phats,breaks,counts){
    tmp <- pweibull(breaks, phats[1], phats[2])
    tmp <- diff(c(0,tmp,1))
    sele <- tmp <= 1.0e-6
    if(any(sele)){
      tmp <- tmp[!sele]
      counts <- counts[!sele]
      delta <- 9.0e+6
    }
    res <- delta-sum(counts * log(tmp))
    res
  }

  out <- optim(parhat, .weibull.bllk, breaks=brks, counts=cnts,
               method="L-BFGS-B", lower=c(tol,tol))
  options(warn=0)
  x <- breaks
  k <- length(x)
  y <- dweibull(x, out$par[1], out$par[2])
  y2 <- pweibull(x, out$par[1], out$par[2])

  a <- x[1];b <- rev(x)[1]
  if(a<=0) a <- 0.1
  if(!is.finite(b)) b <- 2*rev(x)[2]-rev(x)[3]
  lx <- exp(seq(log(a),log(b),length=401))
  ly <- dweibull(lx, out$par[1], out$par[2])

  par0 <- as.numeric(out$par)
  Fx <- pweibull(brks, par0[1], par0[2])
  llk <- sum(log(diff(c(0,Fx,1)))*freq)
  AIC0 <- .AIC(llk,2,N)
  Dn.Zipf0 <- .funCOD(freq,Fx)
  L0 <- .funL2(freq,Fx)

  tmp <- data.frame(lambda = par0[1],
                    kappa=par0[2],
                    Dn.Zipf=Dn.Zipf0,
                    Dn = L0,
                    AIC = AIC0$AIC,
                    BIC = AIC0$BIC,
                    AICc = AIC0$AICc)
  rownames(tmp) <- c("MLE")

  list(pars=par0, aux=tmp,
       Dn.Zipf=Dn.Zipf0, Dn=L0,
       AIC=AIC0$AIC,
       BIC=AIC0$BIC,
       AICc=AIC0$AICc,
       x=x,y=y,y2=y2,lx=lx,ly=ly)
}
#####################################################################
## the location parameter xm needs to be estimated differently as the
## first class will dominate the lieklihood. If the lower boundary of
## the first class, b0, is zero, we estimate the MLE using the data
## without class 1. Otherwise, if b0>0, we find the MLE using the
## whole dataset. The LS-estimate will be used to find the initial
## estimate. It will also provide rough estimate of xm so we can
## determine the range to seach for the MLE of xm.
.pPD <- function(x,xm,alpha){
  stopifnot(is.numeric(x))
  if(any(is.na(x)))
    stop("missing value(s) in 'x'")
  sapply(x,FUN=.subpPD,xm,alpha)
}

.subpPD <- function(x,xm,alpha){
  res = 0
  stopifnot(xm>0)
  stopifnot(alpha>0)
  if(is.finite(x)){
    if(x>xm)
      res = 1-(xm/x)^alpha        
  }else
    res = 1
  
  res
}

.dPD <- function(x,xm,alpha){
  stopifnot(is.numeric(x))
  if(any(is.na(x)))
    stop("missing value(s) in 'x'")
  sapply(x,FUN=.subdPD,xm,alpha)
}

.subdPD <- function(x,xm,alpha){
  res = 0
  stopifnot(xm>0)
  stopifnot(alpha>0)
  if(is.finite(x)){
    if(x>xm)
      res = alpha*(xm/x)^alpha/x        
  }else
    res = 0
  
  res
}

.funPD <- function(pars,freq,brks){
  n <- length(brks)
  stopifnot(length(freq)+1==n)
  brks <- brks[-c(1,n)]
  Fx = .pPD(brks, pars[1], pars[2])
  if(any(is.na(Fx))){
    llk <- 1.0e+9
  }else{
    Px <- diff(c(0,Fx,1))
    sele = Px <= 1.0e-6
    delta = 0 #penalty
    if(any(sele)){
      Px = Px[!sele]
      freq = freq[!sele]
      delta = 1.0e+6 * sum(freq[sele])
    }
    llk <- delta-sum(log(Px)*freq)
  }
  if(!is.finite(llk)) llk <- .Machine$double.xmax
  llk
}

.funPD1 <- function(alp,freq,brks,xm){
    n <- length(brks)
    stopifnot(length(freq)+1==n)
    brks <- brks[-c(1,n)]
    Fx = .pPD(brks, xm, alp)
    if(any(is.na(Fx))){
        llk <- 1.0e+9
    }else{
        Px <- diff(c(0,Fx,1))
        sele = Px <= 1.0e-6
        delta = 0 #penalty
        if(any(sele)){
            Px = Px[!sele]
            freq = freq[!sele]
            delta = 1.0e+6 * sum(freq[sele])
        }
        llk <- delta-sum(log(Px)*freq)
    }
    if(!is.finite(llk)) llk <- .Machine$double.xmax
    llk
}

.funPD2 <- function(xm,freq,brks,alp){
    n <- length(brks)
    stopifnot(length(freq)+1==n)
    brks <- brks[-c(1,n)]
    Fx = .pPD(brks, xm, alp)
    if(any(is.na(Fx))){
        llk <- 1.0e+9
    }else{
        Px <- diff(c(0,Fx,1))
        sele = Px <= 1.0e-6
        delta = 0 #penalty
        if(any(sele)){
            Px = Px[!sele]
            freq = freq[!sele]
            delta = 1.0e+6 * sum(freq[sele])
        }
        llk <- delta-sum(log(Px)*freq)
    }
    if(!is.finite(llk)) llk <- .Machine$double.xmax
    llk
}

## we fit the two parameters (xm and alpha) based on grouped data. The
## range to search for xm has to be from (0, b0=breaks[1])
.fit.PD <- function(freq,breaks){
  cnts <- freq
  brks <- breaks
  n <- length(cnts)
  if(length(brks) != n +1)
    stop("length of the breaks points not match")
  brks <- brks[-c(1,n+1)]
  N <- sum(cnts)
  Fhat <- cumsum(cnts)/N
  Fhat <- Fhat[-n]
  lF <- log(1-Fhat)
  lx <- log(brks)
  sele <- is.finite(lx) & is.finite(lF)
  lx <- lx[sele]; lF <- lF[sele]
  lmout <- lm(lF~lx)
  options(warn=-1)
  
  p2 <- as.numeric(-lmout$coef[2])
  p1 <- as.numeric(exp(lmout$coef[1]/p2))
  if(p1<0) p1 <- max(breaks[1], 0.5*breaks[2])
  par5 <- c(p1,p2) # Least-squares estimates
  
  Fx <- .pPD(brks, par5[1], par5[2])
  llk <- sum(log(diff(c(0,Fx,1)))*freq)
  AIC5 <- .AIC(llk,2,N)
  Dn.Zipf5 <- .funCOD(freq,Fx)
  L5 <- .funL2(freq,Fx)
  
  b0 <- breaks[1]; b1 <- breaks[2]
  if(b0 <= 0 || p1 > b1){
      m1 <- min(b1/10,0.01)
      m2 <- b1
  }else{
      m1 <- min(b0/10,0.01)
      m2 <- b0
  }
  ybrks <- breaks
  ycnts <- freq
  phats <- c(p1,p2)
  out3 <- optim(phats,.funPD, 
                freq=ycnts, brks=ybrks,
                lower=c(m1,max(p2-0.5,0.5)),
                upper=c(m2,min(p2+.5,1.4)))
  if(out3$conv==0){
      par3 <- out3$par
      x0 <- ybrks[-c(1,length(ybrks))]
      Fx <- .pPD(x0, par3[1], par3[2])
      llk <- sum(log(diff(c(0,Fx,1)))*ycnts)
      AIC3 <- .AIC(llk,2,N)
      Dn.Zipf3 <- .funCOD(ycnts,Fx)
      L3 <- .funL2(ycnts,Fx)
  }else{
      par3 <- c(NA, NA)
      AIC3 <- list(AIC=NA, BIC=NA, AICc=NA)
      Dn.Zipf3 <- NA
      L3 <- NA
  }

  
  ## estimate xm with alpha=LS estimate
  out4 <- optimize(.funPD2,
                   interval=c(m1,m2),
                   freq=ycnts,brks=ybrks,alp=p2)
  xmhat <- out4$min
  if(xmhat < 0)
      xmhat <- max(breaks[1], p1)
  
  par4 <- c(xmhat,p2)      
  x0 <- ybrks[-c(1,length(ybrks))]
  Fx <- .pPD(x0, par4[1], par4[2])
  llk <- sum(log(diff(c(0,Fx,1)))*ycnts)
  AIC4 <- .AIC(llk,1,N)
  Dn.Zipf4 <- .funCOD(ycnts,Fx)
  L4 <- .funL2(ycnts,Fx)

      
  if(par5[1]<breaks[2]){
      out2 <- optimize(.funPD1,
                       interval=c(max(p2-0.5,0.5),min(p2+.5,1.4)),
                       freq=ycnts,brks=ybrks,xm=par5[1])
      par2 <- c(par5[1], out2$min)
      Fx <- .pPD(brks, par2[1], par2[2])
      llk <- sum(log(diff(c(0,Fx,1)))*freq)
      AIC2 <- .AIC(llk,1,N)
      Dn.Zipf2 <- .funCOD(freq,Fx)
      L2 <- .funL2(freq,Fx)
  }else{
      par2 <- c(NA,NA)
      AIC2 <- list(AIC=NA, BIC=NA, AICc=NA)
      Dn.Zipf2 <- NA
      L2 <- NA
  }

  ## fit with xm=l[1]
  #m1 <- breaks[1]
  #if(m1<=0) m1 <- min(1, breaks[2]/4)
  #out1 <- optimize(.funPD1,
  #                interval=c(max(p2-.5,.1),min(p2+1,10)),
  #                freq=freq,brks=breaks,xm=m1)
  #par1 <- c(m1,out1$min)

  #Fx <- .pPD(brks, par1[1], par1[2])
  #llk <- sum(log(diff(c(0,Fx,1)))*freq)
  #AIC1 <- .AIC(llk,2,N)
  #Dn.Zipf1 <- .funCOD(freq,Fx)

  ## fit with xm=u[1]
  #m2 <- breaks[2]
  #if(m2<=0) m2 <- m1+0.01
  #out2 <- optimize(.funPD1,
  #                interval=c(max(p2-.5,.1),min(p2+1,10)),
  #                freq=freq,brks=breaks,xm=m2)
  #par2 <- c(m2,out2$min)

  #Fx <- .pPD(brks, par2[1], par2[2])
  #llk <- sum(log(diff(c(0,Fx,1)))*freq)
  #AIC2 <- .AIC(llk,2,N)
  #Dn.Zipf2 <- .funCOD(freq,Fx)
  
  ##if(breaks[1]>0){# fix alpha and search for xm
  ##    p1 <- min(p1,breaks[1])
  ##    lF <- log(1-Fhat)
  ##    lx <- log(p1/brks)
  ##    sele <- is.finite(lx) & is.finite(lF)
  ##    lx <- lx[sele]; lF <- lF[sele]
  ##    p2 <- lm(lF~lx-1)$coef[1]
  ##}
  
  ## search for both parameters
  
  options(warn=0) #################################

  ## choose the best fit in order
  par0 <- NA
  Dn.Zipfopt <- NA
  Dnopt <- NA
  aicopt <- NA
  bicopt <- NA
  aiccopt <- NA
  
  if(!any(is.na(par4))){
      par0 <- par4
      if(par0[1]<breaks[2]&par0[1]>0){
          Dn.Zipfopt <- Dn.Zipf4
          Dnopt <- L4
          aicopt <- AIC4$AIC
          bicopt <- AIC4$BIC
          aiccopt <- AIC4$AICc
      }else{
          par0 <- NA
      }
  }
  
  if(!any(is.na(par2))&any(is.na(par0))){
      par0 <- par2
      Dn.Zipfopt <- Dn.Zipf2
      Dnopt <- L2
      aicopt <- AIC2$AIC
      bicopt <- AIC2$BIC
      aiccopt <- AIC2$AICc
  }

  if(!any(is.na(par5))&any(is.na(par0))){
      par0 <- par5
      Dn.Zipfopt <- Dn.Zipf5
      Dnopt <- L5
      aicopt <- AIC5$AIC
      bicopt <- AIC5$BIC
      aiccopt <- AIC5$AICc
  }

  if(!any(is.na(par3))&any(is.na(par0))){
      par0 <- par3
      Dn.Zipfopt <- Dn.Zipf3
      Dnopt <- L3
      aicopt <- AIC3$AIC
      bicopt <- AIC3$BIC
      aiccopt <- AIC3$AICc
  }

  p1 <- par0[1]; p2 <- par0[2]
  x <- breaks
  k <- length(x)
  y <- .dPD(x, p1, p2)
  y2 <- .pPD(x, p1, p2)
  
  tmp <- data.frame(
      Xm = c(par5[1],par4[1],par3[1],par2[1]),
      alpha = c(par5[2],par4[2],par3[2],par2[2]),
      Dn.Zipf=c(Dn.Zipf5,Dn.Zipf4,Dn.Zipf3,Dn.Zipf2),
      Dn = c(L5,L4,L3,L2),
      AIC = c(AIC5$AIC,AIC4$AIC,AIC3$AIC,AIC2$AIC),
      BIC = c(AIC5$BIC,AIC4$BIC,AIC3$BIC,AIC2$BIC),
      AICc = c(AIC5$AICc,AIC4$AICc,AIC3$AICc,AIC2$AICc))
  
  rownames(tmp) <- c(
      "LSE",
      "MLE(xm)",
      "MLE",
      "MLE(alpha)")
      
  a <- x[1];b <- rev(x)[1]
  if(a<=0) a <- 0.1
  if(!is.finite(b)) b <- 2*rev(x)[2]-rev(x)[3]
  lx <- exp(seq(log(a),log(b),length=401))
  ly <- .dPD(lx, p1, p2)

  list(pars=par0,Dn.Zipf=Dn.Zipfopt, Dn=Dnopt,
       AIC=aicopt,BIC=bicopt,AICc=aiccopt,
       x=x,y=y,y2=y2,aux=tmp,
       lx=lx,ly=ly)
}

## sub-routines for Copula estimation ###########
.Fxy <- function(psi,Fx, Fy){
  stopifnot(psi>0)
  if(psi==1){
    Hxy <- Fx * Fy
  }else{
    Sxy <- 1 + (Fx+Fy)*(psi-1)
    Hxy <- 0.5*(Sxy-sqrt(Sxy^2-4*psi*(psi-1)*Fx*Fy))/(psi-1)
  }
  Hxy
}

.fxy <- function(psi,fx,fy,Fx,Fy){
  stopifnot(psi>0)
  Sxy <- 1 + (Fx+Fy)*(psi-1)
  fn <- psi*fx*fy*(1+(psi-1)*(Fx+Fy-2*Fx*Fy))
  fd <- (Sxy^2-4*psi*(psi-1)*Fx*Fy)^1.5
  fn/fd
}

.pdfCopula <- function(psi,fx,fy,Fx,Fy){
  nx <- length(Fx)
  ny <- length(Fy)
  Fxy <- matrix(0, nrow=nx, ncol=ny)
  for(i in 1:nx){
    for(j in 1:ny){
      Fxy[i,j] <- .fxy(psi,fx[i],fy[j],Fx[i],Fy[j])
    }
  }
  Fxy[is.na(Fxy)] <- 0
  Fxy[Fxy<0] <- 0
  Fxy
}

.cdfCopula <- function(psi,Fx,Fy){
  nx <- length(Fx)
  ny <- length(Fy)
  Fxy <- matrix(0, nrow=nx, ncol=ny)
  for(i in 1:nx){
    for(j in 1:ny){
      Fxy[i,j] <- .Fxy(psi,Fx[i],Fy[j])
    }
  }
  Fxy
}

.llkCopula <- function(psi, Fx, Fy, cnts){

  ## Computing P(cell_ij)
  Fxy <- .cdfCopula(psi,Fx, Fy)

  ##cat("\nDim of 'Fxy'\n")
  ##print(dim(Fxy))

  Mxy <- matrix(0, nrow=nrow(Fxy)-1,
                ncol=ncol(Fxy)-1)

  ##cat("\nDim of 'Mxy'\n")
  ##print(dim(Mxy))

  
  for(i in 1:nrow(Mxy)){
    for(j in 1:ncol(Mxy)){
      Mxy[i,j] <- Fxy[i+1,j+1]+Fxy[i,j]-Fxy[i+1,j]-Fxy[i,j+1]
    }
  }
  Mxy[is.na(Mxy)] <- 0
  tmp <- min(Mxy[Mxy>0])
  Mxy[Mxy<=0] <- tmp/100

  ##cat("\nDim of 'Mxy'\n")
  ##print(dim(Mxy))

  ##cat("\nDim of 'cntx (xy)'\n")
  ##print(dim(cnts))

  
  res2 <- sum(log(Mxy)*cnts)
  #list(LLK=res2,LL0=res1)
  -res2
}

.X2Copula <- function(psi, Fx, Fy, cnts){

  ## Computing P(cell_ij)
  Fxy <- .cdfCopula(psi,Fx, Fy)
  Mxy <- matrix(0, nrow=nrow(Fxy)-1,
                ncol=ncol(Fxy)-1)
  for(i in 1:nrow(Mxy)){
    for(j in 1:ncol(Mxy)){
      Mxy[i,j] <- Fxy[i+1,j+1]+Fxy[i,j]-Fxy[i+1,j]-Fxy[i,j+1]
    }
  }
  #res1 <- sum(log(Mxy[Mxy>0]))
  # cat("\nCell with zero probabilities :=", sum(Mxy<=0))
  #tmp <- min(Mxy[Mxy>0])
  #cat("\nMinimal positive probability :=", tmp,"\n")
  Mxy[is.na(Mxy)] <- 0
  Mxy[Mxy<0] <- 0
  Ex <- Mxy*sum(cnts); Ex[Ex==0] <- 1
  Ox <- cnts; Ox[Ox==0] <- 1
  ##sum((Ex-Ox)^2/Ox)
  out1 <- sum((Ex-Ox)^2/Ex)
  out2 <- sum((Ex-Ox)^2/Ox)
  ##print(c(psi,out1,out2))
  out1
}

.AIC <- function(llk,npar,n){
    AIC <- 2*npar-2*llk
    AICc <- AIC + 2*npar*(npar+1)/(n-npar-1)
    BIC <- log(n)*npar-2*llk
    list(AIC=AIC,BIC=BIC,AICc=AICc)
}

.funCOD <- function(cnts,Fhat){
    Sn <- 1-cumsum(cnts)/(1+sum(cnts))
    k <- length(Sn)
    Sn <- Sn[-k]
    lFn <- log(Sn)
    lFx <- log(1-Fhat)
#    sele <- is.finite(lFn)&is.finite(lFx)
#    penalty <- sum(!sele)
#    if(penalty > 0){
#        lFn <- lFn[sele]
#        lFx <- lFx[sele]
#    }
#    SSE <- sum((lFn-lFx)^2)
#    SSR <- sum((lFx-mean(lFn))^2)
#    CoD <- SSR/(SSR+SSE)
#    CoD^(1+penalty)
    out <- 1

    ##print(rbind(lFn,lFx))



    sele <- is.finite(lFn)&is.finite(lFx)
    if(sum(sele)>0){
        lFn <- lFn[sele]; lFx <- lFx[sele]
        sele <- !is.na(lFn)&!is.na(lFx)
        if(sum(sele)>0){
            lFn <- lFn[sele]; lFx <- lFx[sele]
            out <- max(abs(lFn-lFx))
        }
    }
    out
}

.funL2 <- function(cnts,Fhat){
    Fn <- cumsum(cnts)/(sum(cnts)+1)
    an <- Fn[-length(Fn)]
    bn <- Fhat
    max(abs(an-bn)) # Dn statistic
}

ZipfPlot.FSD <- function(x, x0, plot=FALSE,plot.new=TRUE, weights,...){
    x <- x$x.fit
    if(length(x$breaks) < 3)
        stop("too few classes...")
    bi <- x$breaks[-1]
    n <- sum(x$counts)
    Fn <- cumsum(x$counts)/(n+1)
    ri <- n*(1-Fn)
        
    x0 <- log(bi)
    y0 <- log(ri)
    sele <- is.finite(x0)&is.finite(y0)
    x0 <- x0[sele]; y0 <- y0[sele]
        
    x1 <- log(x$x)
    y1 <- log(n*(1-x$y2)+1)
    sele <- is.finite(x1)&is.finite(y1)
    x1 <- x1[sele]; y1 <- y1[sele]
        
    if(plot.new){
        plot(x0, y0, ...)
        lines(x1, y1,...)
    }else{
        lines(x1,y1,...)
    }
    invisible(NULL)
}

.fit.MEP <- function(freq,breaks){
    n <- length(freq)
    if(length(breaks) != n +1)
        stop("length of the breaks points not match")
    brks <- breaks[-c(1,n+1)]
    N <- sum(freq)
    Fhat <- cumsum(freq)/N
    F1 <- Fhat[1]
    Fhat <- Fhat[-n]
    lF <- log(1-Fhat)
    lx <- log(brks)
    sele <- is.finite(lx) & is.finite(lF)
    lx <- lx[sele]; lF <- lF[sele]
    lmout <- lm(lF~lx)
    b1 <- breaks[2]
  
    alpha <- -lmout$coef[2]
    xm <- b1#exp(lmout$coef[1]/alpha)
    lambda.hat <- -log(1-F1)/b1
    par0 <- c(lambda.hat, xm, alpha)

    p1 <- par0[1];p2 <- par0[2];p3 <- par0[3]
    x <- breaks; k <- length(x)
    y <- .dMEP(x,p1,p2,p3)
    y2 <- .pMEP(x,p1,p2,p3)
    ##print(cbind(x,y2))
    Fx <- y2[-c(1,k)]
    llk <- sum(log(diff(c(0,Fx,1)))*freq)
    AIC5 <- .AIC(llk,3,N)
    Dn.Zipf5 <- .funCOD(freq,Fx)
    L5 <- .funL2(freq,Fx)

  
    tmp <- data.frame(
        lambda = p1,
        xm = p2,
        alpha = p3,
        Dn.Zipf =Dn.Zipf5,
        Dn = L5,
        AIC = AIC5$AIC,
        BIC = AIC5$BIC,
        AICc = AIC5$AICc)
  
    rownames(tmp) <- "LS"
      
    a <- x[1];b <- rev(x)[1]
    if(a<=0) a <- 0.1
    if(!is.finite(b)) b <- 2*rev(x)[2]-rev(x)[3]
    lx <- exp(seq(log(a),log(b),length=401))
    ly <- .dMEP(lx, p1, p2,p3)

    list(pars=par0,Dn.Zipf=Dn.Zipf5, Dn=L5,
         AIC=AIC5$AIC,BIC = AIC5$BIC,
         AICc = AIC5$AICc,
         x=x,y=y,y2=y2,aux=tmp,
         lx=lx,ly=ly)
}
    
.pMEP <- function(x,lambda,xm,alpha){
    stopifnot(is.numeric(x))
    stopifnot(lambda>0)
    stopifnot(xm>0)
    stopifnot(alpha>0)
    if(any(is.na(x)))
        stop("missing value(s) in 'x'")
    sapply(x,FUN=.subpMEP,lambda,xm,alpha)
}

.subpMEP <- function(x,lambda,xm,alpha){
    res = 0
    if(is.finite(x)){
        if(x>xm){
            F1 <- 1-exp(-lambda*xm)
            F2 <- 1-(xm/x)^alpha
            res = F1 + (1-F1)*F2
        }else{
            res <- 1-exp(-lambda*x)
        }
    }else
        res = 1

    res[res>1] <- 1; res[res<0] <- 0
    res
}

.dMEP <- function(x,lambda,xm,alpha){
    stopifnot(is.numeric(x))
    stopifnot(lambda>0)
    stopifnot(xm>0)
    stopifnot(alpha>0)
    if(any(is.na(x)))
        stop("missing value(s) in 'x'")
    sapply(x,FUN=.subdMEP,lambda,xm,alpha)
}


.subdMEP <- function(x,lambda,xm,alpha){
    res = 0
    if(is.finite(x)){
        if(x>xm){
            res = alpha/x*(xm/x)^alpha
        }else{
            res <- lambda*exp(-lambda*x)
        }
    }else
        res = 0
    
    res
}

ImportFSD <- function(x, type, year){
    options(warn=-1)
    
    brks.size <- c(0,5,10,20,50,100,250,500,1000,2500,5000,10000,Inf)
    if(missing(year)){#import from a 2-way table
        year <- 2014
        year <- 2014; cutoff <- 2014-1977+1; cutoff
        brks.age <- c(0,1,2,3,4,5,6,11,16,21,26,cutoff,Inf)
        ## get the XY table
        xy <- NULL
        k <- ncol(x)
        for(i in 1:k){
            tmp <- as.character(x[,i])
            tmp <- as.numeric(tmp)
            tmp[is.na(tmp)] <- 0
            xy <- cbind(xy,tmp)
        }
        xy <- matrix(c(as.numeric(xy)), ncol=k)
        xcnts <- apply(xy,1,sum)
        ycnts <- apply(xy,2,sum)
        xh <- binning(counts=xcnts, breaks=brks.age)
        yh <- binning(counts=ycnts, breaks=brks.size)
    }else{
        year <- round(year)
        if(year <1977 || year > 2014)
            stop("data for this year are not available")
        cutoff <- year-1977+1; 
        brks.age <- c(0,1,2,3,4,5,6,11,16,21,26,cutoff,Inf)
            
        xy <- NULL
        if(missing(type)){
            type <- "size"
        }else{
            type <- match.arg(tolower(type), c("size","age"))
        }
        i <- round(year - 1976)
        if(i<1 || i > nrow(x))
            stop("data for this year are not available")
        
        cnts <- as.numeric(as.character(x[i,]))
        cnts[is.na(cnts)] <- 0
        cnts[cnts<0] <- 0
        sele <- which(cnts > 0)
        cnts <- cnts[sele]

        if(type == 'size'){
            brks.size <- brks.size[c(1,sele+1)]
            yh <- binning(counts=cnts, breaks=brks.size)
            xh <- NULL
        }else{
            brks.age <- brks.age[c(1,sele+1)]
            xh <- binning(counts=cnts, breaks=brks.age)
            yh <- NULL
        }
    }
    options(warn=2)

    breaks <- list(age=brks.age, size=brks.size)
    list(xy=xy, breaks=breaks,age=xh,size=yh)
}

####################################################
.fit.spline <- function(freq,breaks){
    pars <- 1
    n <- sum(freq)
    x <- breaks
    k <- length(x)
    y <- .dspline(x,breaks,freq)
    y2 <- .pspline(x,breaks,freq)
    a <- x[1];b <- rev(x)[1]
    if(a<=0) a <- 0.1
    if(!is.finite(b)) b <- 2*rev(x)[2]-rev(x)[3]
    lx <- exp(seq(log(a),log(b),length=401))
    ly <- .dspline(lx,breaks,freq)

    Fx <- y2[-c(1,k)]
    llk <- sum(log(diff(c(0,Fx,1)))*freq)
    AIC0 <- .AIC(llk,1,n)
    Dn.Zipf0 <- .funCOD(freq,Fx)
    L0 <- .funL2(freq,Fx)
    
    
    tmp <- data.frame(lambda = pars[1],
                      Dn.Zipf=Dn.Zipf0,
                      Dn = L0,
                      AIC = AIC0$AIC,
                      BIC = AIC0$BIC,
                      AICc = AIC0$AICc)
    
    rownames(tmp) <- c("Spline")
  
    list(pars=pars, aux=tmp,
         Dn.Zipf=Dn.Zipf0,Dn=L0,
         AIC=AIC0$AIC,
         BIC=AIC0$BIC,
         AICc=AIC0$AICc,
         x=x,y=y,y2=y2,lx=lx,ly=ly)
}


.dspline <- function(x, breaks,freq)
    {
        F <- freq;
        n <- sum(F)
        x0 <- breaks
        k <- length(x0)
        if(!is.finite(x0[k]))
            x0[k] <- 10*x0[k-1]
        dx <- diff(x0)
        X <- x0[-1]-0.5*dx;
        Y <- freq/sum(freq*dx)
        A <- x0[-k]; B <- x0[-1]
        x0 <- runif(n,rep(A,F),rep(B,F))
        h <-  (1/(4*pi))^(1/10)*(243/(35*n))^(1/5)*sqrt(var(x0))*15^(1/5)
        selex <- is.finite(x)
        X <- c(x0[1],X)
        Y <- c(0,Y)
        out <- spline(X, Y, xout=x[selex])
        f0 <- x 
        f0[selex] <- out$y
        f0[!selex] <- 0
        f0[f0<0] <- 0
        f0
    }

.pspline <- function(x, breaks,freq)
    {
        F <- freq;
        n <- sum(F)
        x0 <- breaks
        k <- length(x0)
        if(!is.finite(x0[k]))
            x0[k] <- 10*x0[k-1]
        X <- x0;
        Y <- c(0, cumsum(freq)/(n+1))
        A <- x0[-k]; B <- x0[-1]
        x0 <- runif(n,rep(A,F),rep(B,F))
        h <-  (1/(4*pi))^(1/10)*(243/(35*n))^(1/5)*sqrt(var(x0))*15^(1/5)
        selex <- is.finite(x)
        
        out <- spline(X, Y, xout=x[selex])
        f0 <- x 
        f0[selex] <- out$y
        f0[!selex] <- 0
        f0[f0<0] <- 0
        f0[f0>1] <- 1
        f0
    }
