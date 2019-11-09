
###  Fitting Copula 09/29/2017 ##################################
## Creates the grid counts from a bivariate data set X
## over an equally-spaced set of grid points
## contained in "gpoints" using the linear
## binning strategy. Note that the FORTRAN subroutine
## "lbtwod" is called.
.bin2D <- function(X, gpoints1, gpoints2)
{
    n <- nrow(X)
    X <- c(X[, 1L], X[, 2L])
    M1 <- length(gpoints1)-1
    M2 <- length(gpoints2)-1
    g1 <- gpoints1
    g2 <- gpoints2
    out <- .Fortran(.F_bintwod, as.double(X), as.integer(n),
                    as.double(g1[-1]), as.double(g2[-1]), 
                    as.integer(M1), as.integer(M2),
                    M = double(M1*M2))
    matrix(out$M, M1, M2)
}

.binshrink <- function(xg,xc,k=5){
    while(length(xc) > k){
        l <- length(xc)
        m <- which(xc == min(xc))[1L]
        if(m == 1){
            xc[2] <- xc[2] + xc[m];
            xg <- xg[-(m+1)]
        }else if(m == l){
            xc[m-1] <- xc[m-1] + xc[m]
            xg <- xg[-m]
        }else{
            if(xc[m-1] > xc[m+1]){
                xc[m+1] <- xc[m+1] + xc[m]
                xg <- xg[-(m+1)]
            }else{
                xc[m-1] <- xc[m-1] + xc[m]
                xg <- xg[-m]
            }
        }
        xc <- xc[-m]
    }
    list(grid=xg, counts=xc)
}

.mixnorm2d <- function(Fx,Fy, Psi){
    stopifnot(Psi > 0)
    Sxy <- 1 + (Fx+Fy)*(Psi-1)
    if(Psi==1){
        Hxy <- Fx * Fy
    }else{
        Hxy <- 0.5 * (Sxy-sqrt(Sxy^2-4*Psi*(Psi-1)*Fx*Fy))/(Psi-1);
    }    
    Hxy
}

fit.NGS.nmix.copula <- function(x,y,mle.large=FALSE,gridsize=101){
    ## must turn mixed off otherwise will have trouble to compute the
    ## copula -- f(0,0) needs to be redefined.
    N <- length(x)
    stopifnot(length(y)==N)
    xfit <- fit.NGS(x)
    yfit <- fit.NGS(y)
    
    nbinx <- length(xfit$xhist$xhist$counts)
    nbiny <- length(yfit$xhist$xhist$counts)
    xgrid <- xfit$xhist$xhist$breaks
    ygrid <- yfit$xhist$xhist$breaks
    xcount <- xfit$xhist$xhist$counts
    ycount <- yfit$xhist$xhist$counts
    ## transform raw data such that the zero's will be grouped
    ## correctly.
    lx <- log(x); ly <- log(y)
    xmin <- min(c(x[x>0],y[y>0]))*.3
    if(any(x == 0)){
        lx[x==0] <- log(xmin)
    }
    
    if(any(y == 0)){
        ly[y==0] <- log(xmin)
    }

    ngroup <- 5
    mat2d.large <- .bin2D(cbind(lx,ly), xgrid, ygrid)
    xgrid2 <- .binshrink(xgrid,xcount,k=ngroup)$grid
    ygrid2 <- .binshrink(ygrid,ycount,k=ngroup)$grid
    mat2d.small <- .bin2D(cbind(lx,ly), xgrid2, ygrid2)

    ## rough estimates of Psi
    nr <- nrow(mat2d.small)
    nc <- ncol(mat2d.small)
    psis <- NULL
    for(i in 1:(nr-1)){
        a <- sum(mat2d.small[1:i,1:i])
        b <- sum(mat2d.small[1:i,(i+1):nc])
        c <- sum(mat2d.small[(i+1):nr,1:i])
        d <- sum(mat2d.small[(i+1):nr,(i+1):nc])
        if(a==0) a <- 0.1
        if(b==0) b <- 0.1
        if(c==0) c <- 0.1
        if(d==0) d <- 0.1
        
        out <- exp(log(a)+log(d)-log(b)-log(c))

        psis <- c(psis, out)
    }

    ##  Compute and save fx, fy, Fx and Fy 
    Fx <- rep(0, ngroup+1)
    Fy <- rep(0, ngroup+1)
    for(i in 1:(ngroup+1)){
        Fx[i] <- pmixnorm(xgrid2[i],xfit$fit$p,xfit$fit$mu,xfit$fit$s)
        Fy[i] <- pmixnorm(ygrid2[i],yfit$fit$p,yfit$fit$mu,yfit$fit$s)
    }

    ## grid search for Psi 
    psis1 <- seq(0.1*max(0.01,min(psis)), 10*max(psis), length=100)
    psis2 <- seq(11*max(psis), 1000*max(psis), length=200)
    psis <- c(psis1, psis2)
    Chi2 <- NULL
    LLK <- NULL
    for(psi in psis){
        ecounts <- matrix(0,nrow=ngroup,ncol=ngroup)
        for(i in 1:ngroup){
            for(j in 1:ngroup){
                a <- Fx[i]
                b <- Fx[i+1]
                c <- Fy[j]
                d <- Fy[j+1]
                Hbd <- .mixnorm2d(b,d,psi)
                Had <- .mixnorm2d(a,d,psi)
                Hbc <- .mixnorm2d(b,c,psi)
                Hac <- .mixnorm2d(a,c,psi)
                ecounts[i,j] <- Hbd - Had - Hbc + Hac
            }
        }
        llk <- sum(mat2d.small*log(ecounts), na.rm=TRUE)
        ecounts <- N*ecounts
        
        chistat <- sum((mat2d.small-ecounts)^2/ecounts)
        Chi2 <- c(Chi2, chistat)
        LLK <- c(LLK, llk)
    }
    df0 <- ngroup^2 - 3*(xfit$fit$ng+yfit$fit$ng) - 2
    if(df0<1) df0 <- 1

    sele1 <- which(Chi2 == min(Chi2))[1]
    sele2 <- which(LLK == max(LLK))[1]
    ## Compute MLE using the large matrix.  This need to be tuned off.
    ## Otherwise, we need to compute using Fortran or C code for large
    ## data analysis.
    if(mle.large){
        ##  Compute and save fx, fy, Fx and Fy
        ng1 <- nrow(mat2d.large)
        ng2 <- ncol(mat2d.large)
        Fx <- rep(0, ng1+1)
        Fy <- rep(0, ng2+1)
        for(i in 1:(ng1+1)){
            Fx[i] <- pmixnorm(xgrid[i],xfit$fit$p,xfit$fit$mu,xfit$fit$s)
        }
        for(i in 1:(ng2+1)){
            Fy[i] <- pmixnorm(ygrid[i],yfit$fit$p,yfit$fit$mu,yfit$fit$s)
        }
        
        LLK2 <- NULL
        for(psi in psis){
            ecounts <- matrix(0,nrow=ng1,ncol=ng2)
            for(i in 1:ng1){
                for(j in 1:ng2){
                    a <- Fx[i]
                    b <- Fx[i+1]
                    c <- Fy[j]
                    d <- Fy[j+1]
                    Hbd <- .mixnorm2d(b,d,psi)
                    Had <- .mixnorm2d(a,d,psi)
                    Hbc <- .mixnorm2d(b,c,psi)
                    Hac <- .mixnorm2d(a,c,psi)
                    ecounts[i,j] <- Hbd - Had - Hbc + Hac
                }
            }
            llk <- sum(mat2d.large*log(ecounts), na.rm=TRUE)
            LLK2 <- c(LLK2, llk)
        }
        sele3 <- which(LLK2 == max(LLK2))[1]
        psi.mle.large <- psis[sele3]
    }else{
        psi.mle.large <- NULL
        LLK2 <- NULL
    }

    tmp <- structure(list(xfit=xfit, yfit=yfit,
                          Psi.chisq = psis[sele1],
                          Psi.mle.small = psis[sele2],
                          Psi.mle.large = psi.mle.large,
                          p.value = 1-pchisq(Chi2[sele1], df0),
                          Psis = psis,
                          ChiSq = Chi2,
                          df=df0,
                          type="NGS.nmix",
                          LLK.small=LLK,
                          LLK.large=LLK2,
                          mat2large=mat2d.large,
                          mat2small=mat2d.small),
                     class='copula')

    out <- .pdf.nmix.copula(tmp, ngrid=gridsize)

    structure(list(x=out$x, y=out$y,z=out$z,
                   type="NGS.nmix",
                   xfit=xfit, yfit=yfit,
                   Psi.chisq = psis[sele1],
                   Psi.mle.small = psis[sele2],
                   Psi.mle.large = psi.mle.large,
                   p.value = 1-pchisq(Chi2[sele1], df0),
                   Psis = psis,
                   ChiSq = Chi2,
                   df=df0,
                   LLK.small=LLK,
                   LLK.large=LLK2,
                   mat2large=mat2d.large,
                   mat2small=mat2d.small),
              class='copula')
}

print.copula <- function(x,...){
    cat("Copula estimate of f(x,y). Type=",x$type)
    cat("\n\tfitting 'x'...\n")
    print(x$xfit)
    cat("\n\tfitting 'y'...\n")
    print(x$yfit)
    cat("\nPsi estimates:\n")
    cat("\tChiSq est=",x$Psi.chisq, "\n\tMLE=",x$Psi.mle.small,"\n")
}

print.copula2 <- function(x,...){
    cat("Copula estimate of f(x,y). Type=",x$type)
    cat("\n\tfitting 'x'...\n")
    print(x$xfit)
    cat("\n\tfitting 'y'...\n")
    print(x$yfit)
    cat("\nPsi estimates:",x$Psi, "\n")
}

persp.copula <- function(x,xlab,ylab,zlab,...){
    if(missing(xlab)) xlab <- "x"
    if(missing(ylab)) ylab <- "y"
    if(missing(zlab)) zlab <- "f(x,y)"
    persp(x=x$x, y=x$y,z=x$z,xlab=xlab,ylab=ylab,zlab=zlab,...)
}

.pdf.nmix.copula <- function(x,ngrid=101){
    ngrid <- round(ngrid)
    stopifnot(ngrid>3)
    tmp <- x$xfit$xhist$xhist$breaks
    x0 <- seq(min(tmp), max(tmp),length=ngrid)
    tmp <- x$yfit$xhist$xhist$breaks
    y0 <- seq(min(tmp), max(tmp),length=ngrid)
    Fx <- rep(0, ngrid); fx <- Fx;
    Fy <- rep(0, ngrid); fy <- Fy; 
    fxy <- matrix(0, nrow=ngrid,ncol=ngrid)
    for(i in 1:ngrid){
        Fx[i] <- pmixnorm(x0[i],x$xfit$fit$p,x$xfit$fit$mu,x$xfit$fit$s)
        Fy[i] <- pmixnorm(y0[i],x$yfit$fit$p,x$yfit$fit$mu,x$yfit$fit$s)
        fx[i] <- dmixnorm(x0[i],x$xfit$fit$p,x$xfit$fit$mu,x$xfit$fit$s)
        fy[i] <- dmixnorm(y0[i],x$yfit$fit$p,x$yfit$fit$mu,x$yfit$fit$s)
    }
    Psi <- x$Psi.mle.small
    for(i in 1:ngrid){
        for(j in 1:ngrid){
            Sxy <- 1+(Fx[i]+Fy[j])*(Psi - 1)
            ftop <- Psi*fx[i]*fy[j]*(1+(Psi-1)*(Fx[i]+Fy[j]-2*Fx[i]*Fy[j]));
            fbottom <- (Sxy^2-4*Psi*(Psi-1)*Fx[i]*Fy[j])^1.5
            fxy[i,j] <- ftop/fbottom
        }
    }
    list(x=x0,y=y0,z=fxy)
}

.pdf.gld.copula <- function(x,ngrid=101){
    ngrid <- round(ngrid)
    stopifnot(ngrid>3)
    tmp <- x$xfit$xhist$breaks
    k <- length(tmp)
    if(is.finite(tmp[1])){
        xmin <- tmp[1]
    }else{
        xmin <- tmp[2]
    }
    if(is.finite(tmp[k])){
        xmax <- tmp[k]
    }else{
        xmax <- tmp[k-1]
    }
    
    x0 <- seq(xmin, xmax,length=ngrid)

    tmp <- x$yfit$xhist$breaks
    k <- length(tmp)
    if(is.finite(tmp[1])){
        xmin <- tmp[1]
    }else{
        xmin <- tmp[2]
    }
    if(is.finite(tmp[k])){
        xmax <- tmp[k]
    }else{
        xmax <- tmp[k-1]
    }
    
    y0 <- seq(xmin, xmax,length=ngrid)
    
    Fx <- rep(0, ngrid); fx <- Fx;
    Fy <- rep(0, ngrid); fy <- Fy; 
    fxy <- matrix(0, nrow=ngrid,ncol=ngrid)
    for(i in 1:ngrid){
        Fx[i] <- pgld(x0[i],x$xfit$pars)
        Fy[i] <- pgld(y0[i],x$yfit$pars)
        fx[i] <- dgld(x0[i],x$xfit$pars)
        fy[i] <- dgld(y0[i],x$yfit$pars)
    }
    Psi <- x$Psi.mle.small
    for(i in 1:ngrid){
        for(j in 1:ngrid){
            Sxy <- 1+(Fx[i]+Fy[j])*(Psi - 1)
            ftop <- Psi*fx[i]*fy[j]*(1+(Psi-1)*(Fx[i]+Fy[j]-2*Fx[i]*Fy[j]));
            fbottom <- (Sxy^2-4*Psi*(Psi-1)*Fx[i]*Fy[j])^1.5
            fxy[i,j] <- ftop/fbottom
        }
    }
    list(x=x0,y=y0,z=fxy)
}

.pdf.mlnorm.copula <- function(x,ngrid=101){
    ngrid <- round(ngrid)
    stopifnot(ngrid>3)
    tmp <- x$xfit$data;
    xmin <- log(min(tmp[tmp>0])*.25)
    xmax <- log(max(tmp)+1)
    tmp <- x$yfit$data;
    ymin <- log(min(tmp[tmp>0])*.25)
    ymax <- log(max(tmp)+1)

    x0 <- seq(xmin, xmax,length=ngrid)
    y0 <- seq(ymin, ymax,length=ngrid)
    Fx <- rep(0, ngrid); fx <- Fx;
    Fy <- rep(0, ngrid); fy <- Fy; 
    fxy <- matrix(0, nrow=ngrid,ncol=ngrid)
    for(i in 1:ngrid){
        Fx[i] <- pmixnorm(x0[i],x$xfit$p,x$xfit$meanlog,x$xfit$sdlog)
        Fy[i] <- pmixnorm(y0[i],x$yfit$p,x$yfit$meanlog,x$yfit$sdlog)
        fx[i] <- dmixnorm(x0[i],x$xfit$p,x$xfit$meanlog,x$xfit$sdlog)
        fy[i] <- dmixnorm(y0[i],x$yfit$p,x$yfit$meanlog,x$yfit$sdlog)
    }
    Psi <- x$Psi.mle.small
    for(i in 1:ngrid){
        for(j in 1:ngrid){
            Sxy <- 1+(Fx[i]+Fy[j])*(Psi - 1)
            ftop <- Psi*fx[i]*fy[j]*(1+(Psi-1)*(Fx[i]+Fy[j]-2*Fx[i]*Fy[j]));
            fbottom <- (Sxy^2-4*Psi*(Psi-1)*Fx[i]*Fy[j])^1.5
            fxy[i,j] <- ftop/fbottom
        }
    }
    list(x=x0,y=y0,z=fxy)
}

fit.nmix.copula <- function(x,y,mle.large=FALSE,gridsize=101){
    fit.NGS.nmix.copula(exp(x), exp(y), mle.large = mle.large,
                        gridsize = gridsize)
}


fit.GLD.copula <- function(x,y,mle.large=FALSE,gridsize=101){
    ## must turn mixed off otherwise will have trouble to compute the
    ## copula -- f(0,0) needs to be redefined.
    N <- length(x)
    stopifnot(length(y)==N)
    xfit <- fit.GLD(x)
    yfit <- fit.GLD(y)

    nbinx <- length(xfit$xhist$counts)
    xgrid <- xfit$xhist$breaks
    xcount <- xfit$xhist$counts
    
    nbiny <- length(yfit$xhist$counts)
    ygrid <- yfit$xhist$breaks
    ycount <- yfit$xhist$counts

    ngroup <- 5
    mat2d.large <- .bin2D(cbind(x,y), xgrid, ygrid)
    xgrid2 <- .binshrink(xgrid,xcount,k=ngroup)$grid
    ygrid2 <- .binshrink(ygrid,ycount,k=ngroup)$grid
    mat2d.small <- .bin2D(cbind(x,y), xgrid2, ygrid2)

    ## rough estimates of Psi
    nr <- nrow(mat2d.small)
    nc <- ncol(mat2d.small)
    psis <- NULL
    for(i in 1:(nr-1)){
        a <- sum(mat2d.small[1:i,1:i])
        b <- sum(mat2d.small[1:i,(i+1):nc])
        c <- sum(mat2d.small[(i+1):nr,1:i])
        d <- sum(mat2d.small[(i+1):nr,(i+1):nc])
        if(a==0) a <- 0.1
        if(b==0) b <- 0.1
        if(c==0) c <- 0.1
        if(d==0) d <- 0.1
        
        out <- exp(log(a)+log(d)-log(b)-log(c))

        psis <- c(psis, out)
    }

    ##  Compute and save fx, fy, Fx and Fy 
    Fx <- rep(0, ngroup+1)
    Fy <- rep(0, ngroup+1)
    for(i in 1:(ngroup+1)){
        Fx[i] <- pgld(xgrid2[i],xfit$pars)
        Fy[i] <- pgld(ygrid2[i],yfit$pars)
    }

    ## grid search for Psi 
    psis1 <- seq(0.1*max(0.01,min(psis)), 10*max(psis), length=100)
    psis2 <- seq(11*max(psis), 1000*max(psis), length=200)
    psis <- c(psis1, psis2)
    Chi2 <- NULL
    LLK <- NULL
    for(psi in psis){
        ecounts <- matrix(0,nrow=ngroup,ncol=ngroup)
        for(i in 1:ngroup){
            for(j in 1:ngroup){
                a <- Fx[i]
                b <- Fx[i+1]
                c <- Fy[j]
                d <- Fy[j+1]
                Hbd <- .mixnorm2d(b,d,psi)
                Had <- .mixnorm2d(a,d,psi)
                Hbc <- .mixnorm2d(b,c,psi)
                Hac <- .mixnorm2d(a,c,psi)
                ecounts[i,j] <- Hbd - Had - Hbc + Hac
            }
        }
        llk <- sum(mat2d.small*log(ecounts), na.rm=TRUE)
        ecounts <- N*ecounts
        
        chistat <- sum((mat2d.small-ecounts)^2/ecounts)
        Chi2 <- c(Chi2, chistat)
        LLK <- c(LLK, llk)
    }
    df0 <- ngroup^2 - 12
    if(df0<1) df0 <- 1

    sele1 <- which(Chi2 == min(Chi2))[1]
    sele2 <- which(LLK == max(LLK))[1]
    ## Compute MLE using the large matrix.  This need to be tuned off.
    ## Otherwise, we need to compute using Fortran or C code for large
    ## data analysis.
    if(mle.large){
        ##  Compute and save fx, fy, Fx and Fy
        ng1 <- nrow(mat2d.large)
        ng2 <- ncol(mat2d.large)
        Fx <- rep(0, ng1+1)
        Fy <- rep(0, ng2+1)
        for(i in 1:(ng1+1)){
            Fx[i] <- pgld(xgrid[i],xfit$pars)
        }
        for(i in 1:(ng2+1)){
            Fy[i] <- pgld(ygrid[i],yfit$pars)
        }
        
        LLK2 <- NULL
        for(psi in psis){
            ecounts <- matrix(0,nrow=ng1,ncol=ng2)
            for(i in 1:ng1){
                for(j in 1:ng2){
                    a <- Fx[i]
                    b <- Fx[i+1]
                    c <- Fy[j]
                    d <- Fy[j+1]
                    Hbd <- .mixnorm2d(b,d,psi)
                    Had <- .mixnorm2d(a,d,psi)
                    Hbc <- .mixnorm2d(b,c,psi)
                    Hac <- .mixnorm2d(a,c,psi)
                    ecounts[i,j] <- Hbd - Had - Hbc + Hac
                }
            }
            llk <- sum(mat2d.large*log(ecounts), na.rm=TRUE)
            LLK2 <- c(LLK2, llk)
        }
        sele3 <- which(LLK2 == max(LLK2))[1]
        psi.mle.large <- psis[sele3]
    }else{
        psi.mle.large <- NULL
        LLK2 <- NULL
    }

    tmp <- structure(list(xfit=xfit, yfit=yfit,
                          Psi.chisq = psis[sele1],
                          Psi.mle.small = psis[sele2],
                          Psi.mle.large = psi.mle.large,
                          p.value = 1-pchisq(Chi2[sele1], df0),
                          Psis = psis,
                          ChiSq = Chi2,
                          df=df0,
                          type="NGS.nmix",
                          LLK.small=LLK,
                          LLK.large=LLK2,
                          mat2large=mat2d.large,
                          mat2small=mat2d.small),
                     class='copula')

    out <- .pdf.gld.copula(tmp, ngrid=gridsize)

    structure(list(x=out$x, y=out$y,z=out$z,
                   type="NGS.nmix",
                   xfit=xfit, yfit=yfit,
                   Psi.chisq = psis[sele1],
                   Psi.mle.small = psis[sele2],
                   Psi.mle.large = psi.mle.large,
                   p.value = 1-pchisq(Chi2[sele1], df0),
                   Psis = psis,
                   ChiSq = Chi2,
                   df=df0,
                   LLK.small=LLK,
                   LLK.large=LLK2,
                   mat2large=mat2d.large,
                   mat2small=mat2d.small),
              class='copula')
}

.bin.mlnorm.copula <- function(x){
    stopifnot(all(x>=0))
    if(any(x==0)){
        x1 <- x[x>0]
        stopifnot(length(x1) > 30)
        out <- hist(log(x1), plot=FALSE)
        xmin <- min(x1)
        brks <- out$breaks
        x2 <- x
        x2[x2 == 0] <- xmin*.3
        brks <- c(log(xmin*.25), log(xmin*.999),brks[-1])
        out <- hist(log(x2), breaks=brks, plot=FALSE)
    }else{
        stopifnot(length(x) > 30)
        out <- hist(log(x), plot=FALSE)
    }
        out
    }

.fitmlnormk <- function(x,k){
    np <- sum(x > 0)
    
    if(np > 20){
        fx <- fit.mlnorm(x, k=k, method='fnmm',optim=FALSE)
    }else{
        xt <- table(x)
        if(length(xt) <= 2){
            fx <- fit.mlnorm(x, method='mop')
        }else{
            fx <- fit.mlnorm(x, k=1, method='lognormal')
        }
    }
    fx
}

fit.mlnorm.copula <- function(x,y,mle.large=FALSE,gridsize=101,k=2){
    ## assuming we have large data.  otherwise, stop! We treat the
    ## zeros separately because they are not allowed in mlnorm.
    N <- length(x)
    stopifnot(length(y)==N)
    xfit <- .fitmlnormk(x,k=k)
    yfit <- .fitmlnormk(y,k=k)

    xhist <- .bin.mlnorm.copula(x)
    yhist <- .bin.mlnorm.copula(y)
    
    nbinx <- length(xhist$counts)
    nbiny <- length(yhist$counts)
    xgrid <- xhist$breaks
    ygrid <- yhist$breaks
    xcount <- xhist$counts
    ycount <- yhist$counts
    ## transform raw data such that the zero's will be grouped
    ## correctly.
    lx <- log(x); ly <- log(y)
    xmin <- min(c(x[x>0],y[y>0]))*.25
    if(any(x == 0)){
        lx[x==0] <- log(xmin)
    }
    
    if(any(y == 0)){
        ly[y==0] <- log(xmin)
    }

    ngroup <- 5
    mat2d.large <- .bin2D(cbind(lx,ly), xgrid, ygrid)
    xgrid2 <- .binshrink(xgrid,xcount,k=ngroup)$grid
    ygrid2 <- .binshrink(ygrid,ycount,k=ngroup)$grid
    mat2d.small <- .bin2D(cbind(lx,ly), xgrid2, ygrid2)

    ## rough estimates of Psi
    nr <- nrow(mat2d.small)
    nc <- ncol(mat2d.small)
    psis <- NULL
    for(i in 1:(nr-1)){
        a <- sum(mat2d.small[1:i,1:i])
        b <- sum(mat2d.small[1:i,(i+1):nc])
        c <- sum(mat2d.small[(i+1):nr,1:i])
        d <- sum(mat2d.small[(i+1):nr,(i+1):nc])
        if(a==0) a <- 0.1
        if(b==0) b <- 0.1
        if(c==0) c <- 0.1
        if(d==0) d <- 0.1
        
        out <- exp(log(a)+log(d)-log(b)-log(c))
        psis <- c(psis, out)
    }

    ##  Compute and save fx, fy, Fx and Fy 
    Fx <- rep(0, ngroup+1)
    Fy <- rep(0, ngroup+1)
    for(i in 1:(ngroup+1)){
        Fx[i] <- pmixnorm(xgrid2[i],xfit$p,xfit$meanlog,xfit$sdlog)
        Fy[i] <- pmixnorm(ygrid2[i],yfit$p,yfit$meanlog,yfit$sdlog)
    }

    ## grid search for Psi 
    psis1 <- seq(0.1*max(0.01,min(psis)), 10*max(psis), length=100)
    psis2 <- seq(11*max(psis), 1000*max(psis), length=200)
    psis <- c(psis1, psis2)
    Chi2 <- NULL
    LLK <- NULL
    for(psi in psis){
        ecounts <- matrix(0,nrow=ngroup,ncol=ngroup)
        for(i in 1:ngroup){
            for(j in 1:ngroup){
                a <- Fx[i]
                b <- Fx[i+1]
                c <- Fy[j]
                d <- Fy[j+1]
                Hbd <- .mixnorm2d(b,d,psi)
                Had <- .mixnorm2d(a,d,psi)
                Hbc <- .mixnorm2d(b,c,psi)
                Hac <- .mixnorm2d(a,c,psi)
                ecounts[i,j] <- Hbd - Had - Hbc + Hac
            }
        }
        llk <- sum(mat2d.small*log(ecounts), na.rm=TRUE)
        ecounts <- N*ecounts
        
        chistat <- sum((mat2d.small-ecounts)^2/ecounts)
        Chi2 <- c(Chi2, chistat)
        LLK <- c(LLK, llk)
    }
    df0 <- ngroup^2 - xfit$npara - yfit$npara
    if(df0<1) df0 <- 1

    sele1 <- which(Chi2 == min(Chi2))[1]
    sele2 <- which(LLK == max(LLK))[1]
    ## Compute MLE using the large matrix.  This need to be tuned off.
    ## Otherwise, we need to compute using Fortran or C code for large
    ## data analysis.
    if(mle.large){
        ##  Compute and save fx, fy, Fx and Fy
        ng1 <- nrow(mat2d.large)
        ng2 <- ncol(mat2d.large)
        Fx <- rep(0, ng1+1)
        Fy <- rep(0, ng2+1)
        for(i in 1:(ng1+1)){
            Fx[i] <- pmixnorm(xgrid2[i],xfit$p,xfit$meanlog,xfit$sdlog)
        }
        for(i in 1:(ng2+1)){
            Fy[i] <- pmixnorm(ygrid2[i],yfit$p,yfit$meanlog,yfit$sdlog)
        }
        
        LLK2 <- NULL
        for(psi in psis){
            ecounts <- matrix(0,nrow=ng1,ncol=ng2)
            for(i in 1:ng1){
                for(j in 1:ng2){
                    a <- Fx[i]
                    b <- Fx[i+1]
                    c <- Fy[j]
                    d <- Fy[j+1]
                    Hbd <- .mixnorm2d(b,d,psi)
                    Had <- .mixnorm2d(a,d,psi)
                    Hbc <- .mixnorm2d(b,c,psi)
                    Hac <- .mixnorm2d(a,c,psi)
                    ecounts[i,j] <- Hbd - Had - Hbc + Hac
                }
            }
            llk <- sum(mat2d.large*log(ecounts), na.rm=TRUE)
            LLK2 <- c(LLK2, llk)
        }
        sele3 <- which(LLK2 == max(LLK2))[1]
        psi.mle.large <- psis[sele3]
    }else{
        psi.mle.large <- NULL
        LLK2 <- NULL
    }

    tmp <- structure(list(xfit=xfit, yfit=yfit,
                          Psi.chisq = psis[sele1],
                          Psi.mle.small = psis[sele2],
                          Psi.mle.large = psi.mle.large,
                          p.value = 1-pchisq(Chi2[sele1], df0),
                          Psis = psis,
                          ChiSq = Chi2,
                          df=df0,
                          type="NGS.nmix",
                          LLK.small=LLK,
                          LLK.large=LLK2,
                          mat2large=mat2d.large,
                          mat2small=mat2d.small),
                     class='copula')

    out <- .pdf.mlnorm.copula(tmp, ngrid=gridsize)

    structure(list(x=out$x, y=out$y,z=out$z,
                   type="mlnorm",
                   xfit=xfit, yfit=yfit,
                   Psi.chisq = psis[sele1],
                   Psi.mle.small = psis[sele2],
                   Psi.mle.large = psi.mle.large,
                   p.value = 1-pchisq(Chi2[sele1], df0),
                   Psis = psis,
                   ChiSq = Chi2,
                   df=df0,
                   LLK.small=LLK,
                   LLK.large=LLK2,
                   mat2large=mat2d.large,
                   mat2small=mat2d.small),
              class='copula')
}


fit.GLD2bin <- function(x,y,data,xlim,ylim,ngrid){
    if(class(x)=="copula2"){
        res <- x
        out <- .pdf.gld.copula2(res, xlim=xlim,ylim=ylim,ngrid=ngrid)
        res$x <- out$x
        res$y <- out$y
        res$z <- out$z
    }else{
        if(class(x)!="FMKL"&&class(x)!="FSD")
            stop("'x' must be a fitted model: 'FKML-GLD' or 'FSD'")
        if(class(y)!="FMKL"&&class(y)!="FSD")
            stop("'y' must be a fitted model: 'FKML-GLD' or 'FSD'")
        xfit <- x
        yfit <- y
        xn <- sum(xfit$xhist$counts)
        yn <- sum(yfit$xhist$counts)
        N <- sum(data)
        if(N != xn)
            stop("'x' is not a fitted GLD on 'data'")
        if(N != yn)
            stop("'y' is not a fitted GLD on 'data'")
        
        nbinx <- length(xfit$xhist$counts)
        xgrid <- xfit$xhist$breaks
        xcount <- xfit$xhist$counts
        
        nbiny <- length(yfit$xhist$counts)
        ygrid <- yfit$xhist$breaks
        ycount <- yfit$xhist$counts
        
        mat2d.large <- data # only (aggregated) data
        
        ## rough estimates of Psi
        nr <- nrow(mat2d.large)
        nc <- ncol(mat2d.large)
        psis <- NULL
        for(i in 1:(min(nr,nc)-1)){
            a <- sum(mat2d.large[1:i,1:i])
            b <- sum(mat2d.large[1:i,(i+1):nc])
            c <- sum(mat2d.large[(i+1):nr,1:i])
            d <- sum(mat2d.large[(i+1):nr,(i+1):nc])
            if(a==0) a <- 0.1
            if(b==0) b <- 0.1
            if(c==0) c <- 0.1
            if(d==0) d <- 0.1
            
            out <- exp(log(a)+log(d)-log(b)-log(c))
            
            psis <- c(psis, out)
        }
        
        ## grid search for Psi 
        psis1 <- seq(0.1*max(0.01,min(psis)), 10*max(psis), length=100)
        psis2 <- seq(11*max(psis), 1000*max(psis), length=200)
        psis <- c(psis1, psis2)
        
        ##  Compute and save fx, fy, Fx and Fy
        ng1 <- nrow(mat2d.large)
        ng2 <- ncol(mat2d.large)
        Fx <- rep(0, ng1+1)
        Fy <- rep(0, ng2+1)
        for(i in 1:(ng1+1)){
            if(xfit$method=="GLD"||class(xfit)=="FMKL"){
                Fx[i] <- pgld(xgrid[i],xfit$pars)
            }else{
                Fx[i] <- pmlnorm(xgrid[i], xfit$pars$p,
                                  xfit$pars$mean,
                                  xfit$pars$sd)
            }
        }
        
        for(i in 1:(ng2+1)){
            if(yfit$method=="GLD"||class(yfit)=="FMKL"){
                Fy[i] <- pgld(ygrid[i],yfit$pars)
            }else{
                Fy[i] <- pmlnorm(ygrid[i], yfit$pars$p,
                                  yfit$pars$mean,
                                  yfit$pars$sd)
            }
        }
        
        LLK2 <- NULL
        for(psi in psis){
            ecounts <- matrix(0,nrow=ng1,ncol=ng2)
            for(i in 1:ng1){
                for(j in 1:ng2){
                    a <- Fx[i]
                    b <- Fx[i+1]
                    c <- Fy[j]
                    d <- Fy[j+1]
                    Hbd <- .mixnorm2d(b,d,psi)
                    Had <- .mixnorm2d(a,d,psi)
                    Hbc <- .mixnorm2d(b,c,psi)
                    Hac <- .mixnorm2d(a,c,psi)
                    tmp <- Hbd - Had - Hbc + Hac
                    if(tmp<0) tmp <- 0
                    ecounts[i,j] <- tmp
                }
            }
            llk <- sum(mat2d.large*log(ecounts), na.rm=TRUE)
            LLK2 <- c(LLK2, llk)
        }
        sele3 <- which(LLK2 == max(LLK2))[1]
        psi.mle.large <- psis[sele3]
        
        res <- structure(list(xfit=xfit, yfit=yfit,
                              Psi = psi.mle.large),
                         class='copula2')
        
        out <- .pdf.gld.copula2(res, xlim=xlim,ylim=ylim,ngrid=ngrid)
        res$x <- out$x
        res$y <- out$y
        res$z <- out$z
    }
    res
    
}


.pdf.gld.copula2 <- function(x,xlim, ylim,ngrid=101){
    ngrid <- round(ngrid)
    stopifnot(ngrid>3)
    if(missing(xlim)){
        tmp <- x$xfit$xhist$breaks
        k <- length(tmp)
        if(is.finite(tmp[1])){
            xmin <- tmp[1]
        }else{
            xmin <- tmp[2]
        }
        if(is.finite(tmp[k])){
            xmax <- tmp[k]
        }else{
            xmax <- tmp[k-1]
        }
    }else{
        xmin <- xlim[1]; xmax <- xlim[2]
    }
    x0 <- seq(xmin, xmax,length=ngrid)
        
    if(missing(ylim)){
        tmp <- x$yfit$xhist$breaks
        k <- length(tmp)
        if(is.finite(tmp[1])){
            xmin <- tmp[1]
        }else{
            xmin <- tmp[2]
        }
        if(is.finite(tmp[k])){
            xmax <- tmp[k]
        }else{
            xmax <- tmp[k-1]
        }
    }else{
        xmin <- ylim[1]; xmax <- ylim[2]
    }
    y0 <- seq(xmin, xmax,length=ngrid)
    
    Fx <- rep(0, ngrid); fx <- Fx;
    Fy <- rep(0, ngrid); fy <- Fy; 
    fxy <- matrix(0, nrow=ngrid,ncol=ngrid)
    for(i in 1:ngrid){
        if(class(x$xfit)=="FMKL"||x$xfit$method=="GLD"){
            Fx[i] <- pgld(x0[i],x$xfit$pars)
            fx[i] <- dgld(x0[i],x$xfit$pars)
        }else{
            Fx[i] <- pmlnorm(x0[i],x$xfit$pars$p,
                             x$xfit$pars$mean,x$xfit$pars$sd)
            fx[i] <- dmlnorm(x0[i],x$xfit$pars$p,
                             x$xfit$pars$mean,x$xfit$pars$sd)
        }
        if(class(x$yfit)=="FMKL"||x$yfit$method=="GLD"){
            Fy[i] <- pgld(y0[i],x$yfit$pars)
            fy[i] <- dgld(y0[i],x$yfit$pars)
        }else{
            Fy[i] <- pmlnorm(y0[i], x$yfit$pars$p,
                             x$yfit$pars$mean,x$yfit$pars$sd)
            fy[i] <- dmlnorm(y0[i], x$yfit$pars$p,
                             x$yfit$pars$mean,x$yfit$pars$sd)
        }
    }
    Psi <- x$Psi
    for(i in 1:ngrid){
        for(j in 1:ngrid){
            Sxy <- 1+(Fx[i]+Fy[j])*(Psi - 1)
            ftop <- Psi*fx[i]*fy[j]*(1+(Psi-1)*(Fx[i]+Fy[j]-2*Fx[i]*Fy[j]));
            fbottom <- (Sxy^2-4*Psi*(Psi-1)*Fx[i]*Fy[j])^1.5
            fxy[i,j] <- ftop/fbottom
        }
    }
    list(x=x0,y=y0,z=fxy)
}
