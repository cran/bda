
pro.test <- function(x,y,group,cutoff,x.range,type){
    if(missing(x)){
        out <- .prottest1(x=y,group=group,x.range=x.range)
    }else{
        if(missing(type)){
            out <- .prottest2(x=x,y=y,group=group,cutoff=cutoff,x.range=x.range)
        }else{
            res <- .prottest2(x=x,y=y,group=group,cutoff=cutoff,x.range=x.range)
            group <- as.factor(group)
            glvs <- levels(group)
            n <- tapply(group,group,length); #print(n)
            mu0 <- res$T0[1]
            mut <- res$T1.trt[1]
            muc <- res$T1.ctrl[1]
            s0 <- res$T0[2]
            st <- res$T1.trt[2]
            sc <- res$T1.trt[2]
            Z1 <- (mut-muc)/sqrt(st^2/n[1]+sc^2/n[2]);
            pv1 <- pnorm(-abs(Z1))
            res$Z1 <- c(Z1,pv1)

            if(is.numeric(type)){
                if(length(type) == 2){
                    if(any(type < -1) || any(type > 1))
                        stop("invalid correlation coefficient")
                    vt <- st^2 + s0^2 - 2*type[2]*st*s0
                    vc <- sc^2 + s0^2 - 2*type[2]*sc*s0
                    Z2 <- (mut-muc)/sqrt(vc/n[1]+vt/n[2]);
                    pv2 <- pnorm(-abs(Z2))
                    res$Z2 <- c(Z2,pv2)
                }else{
                    sele <- group == glvs[1]
                    vc <- var(x[sele]-y[sele],na.rm=TRUE)
                    vt <- var(x[!sele]-y[!sele],na.rm=TRUE)
                    Z2  <- (mut - muc)/sqrt(vt/n[2]+vc/n[1])
                    pv2 <- pnorm(-abs(Z2))
                    res$Z2 <- c(Z2,pv2)
                }
            }else{
                if(nlevels(group) != 2)
                    stop("this test support RCT with 2 arms only")
                mu <- c(res[1,1],res[1,2],res[1,3])
                sele <- group == glvs[1]
                s0 <- .bdvar(x[sele],y[sele],mu=mu[1]-mu[2],type=type)
                sele <- group == glvs[2]
                s1 <- .bdvar(x[sele],y[sele],mu=mu[1]-mu[3],type=type)
                Z2  <- (s0[1] - s1[1])/sqrt(s0[2]^2+s1[2]^2)
                pv2 <- pnorm(-abs(Z2))
                res$Z2 <- c(Z2,pv2)
            }
            out <- res
        }
    }
    out
}

.bdvar <- function(x,y,mu,type="VAS"){
    n <- length(x)
    stopifnot(is.numeric(x))
    stopifnot(is.numeric(x))
    if(length(y) != n)
        stop("'x' and 'y' length not match")
    sele <- is.na(x) | is.na(y)
    if(any(sele)){
        if(sum(!sele)<5)
            stop("too few data points")
        x <- x[!sele]
        y <- y[!sele]
        n <- length(x)
    }
    x.limit <- mean(x) + 3 * sd(x)
    if(x.limit <= 10) x.limit <- 10 + sd(x)
    y.limit <- mean(y) + 3 * sd(y)
    if(y.limit <= 10) y.limit <- 10 + sd(y)
    if(missing(mu))
        mu <- mean(x) - mean(y)
    v0 <- var(x-y)
    
    type <- match.arg(tolower(type),c("vas","wbf","nrs"))
    out <- NULL
    if(type=='vas'){
        ## treat '0' as 0, and 10 as >=10
        x10 <- x == 10
        y10 <- y == 10
        ## d2 if x0 or x1 = 10; d1 for x0 or x1 <10
        d1 <- NULL
        d2 <- NULL
        sele <- !x10 & !y10
        if(any(sele))
            d1 <- x[sele] - y[sele]

        sele <- !x10 & y10
        if(any(sele)){
            ll <- -999; #x[sele] - y.limit
            ul <- x[sele] - 9.5
            d2 <- data.frame(ll=ll,ul=ul)
        }
        
        sele <- x10 & !y10
        if(any(sele)){
            ll <- 9.5 - x[sele]
            ul <- 999; #x.limit - x[sele]
            tmp <- data.frame(ll=ll,ul=ul)
            d2 <- rbind(d2,tmp)
        }

        sele <- x10 & y10
        if(any(sele)){
            ll <- 9.5 - y.limit
            ul <- x.limit - 9.5
            tmp <- data.frame(ll=rep(ll,sum(sele)),
                              ul=rep(ul,sum(sele)))
            d2 <- rbind(d2,tmp)
        }
        if(is.null(d2)){
            out <- var(d1)
        }else if(is.null(d1)){
            out <- .dbvmle1(d2,mu0=mu,var0=v0)
        }else{
            out <- .dbvmle2(d1,d2,mu0=mu,var0=v0)
        }
    }else if(type=='nrs'){
        ## treat '0' as 0, and 10 as >=10
        x10 <- x == 10
        y10 <- y == 10
        ## d2 if x0 or x1 = 10; d1 for x0 or x1 <10
        d2 <- NULL
        sele <- !x10 & !y10
        if(any(sele)){
            tmp <- x[sele] - y[sele]
            d2 <- data.frame(ll=tmp-1,ul=tmp+1)
        }

        sele <- !x10 & y10
        if(any(sele)){
            ll <- -999; #x[sele] -0.5 - y.limit
            ul <- x[sele] - 9.0
            tmp <- data.frame(ll=ll,ul=ul)
            d2 <- rbind(d2,tmp)
        }
        
        sele <- x10 & !y10
        if(any(sele)){
            ll <- 9.5 - x[sele] - 0.5
            ul <- 999; #x.limit - x[sele] + 0.5
            tmp <- data.frame(ll=ll,ul=ul)
            d2 <- rbind(d2,tmp)
        }

        sele <- x10 & y10
        if(any(sele)){
            ll <- 9.5 - y.limit
            ul <- x.limit - 9.5
            tmp <- data.frame(ll=rep(ll,sum(sele)),
                              ul=rep(ul,sum(sele)))
            d2 <- rbind(d2,tmp)
        }
        out <- .dbvmle1(d2,mu0=mu,var0=v0)
    }else if(type=='wbf'){
        ## treat '0' as 0, and 10 as >=10
        x10 <- x == 10
        y10 <- y == 10
        ## d2 if x0 or x1 = 10; d1 for x0 or x1 <10
        d2 <- NULL
        sele <- !x10 & !y10
        if(any(sele)){
            tmp <- x[sele] - y[sele]
            d2 <- data.frame(ll=tmp-2,ul=tmp+2)
        }

        sele <- !x10 & y10
        if(any(sele)){
            ll <- -999; #x[sele] - 1 - y.limit
            ul <- x[sele] - 8 #+ 1 - 9.0
            tmp <- data.frame(ll=ll,ul=ul)
            d2 <- rbind(d2,tmp)
        }
        
        sele <- x10 & !y10
        if(any(sele)){
            ll <- 9 - x[sele] - 1
            ul <- 999 #x.limit - x[sele] + 1
            tmp <- data.frame(ll=ll,ul=ul)
            d2 <- rbind(d2,tmp)
        }

        sele <- x10 & y10
        if(any(sele)){
            ll <- 9 - y.limit
            ul <- x.limit - 9
            tmp <- data.frame(ll=rep(ll,sum(sele)),
                              ul=rep(ul,sum(sele)))
            d2 <- rbind(d2,tmp)
        }
        out <- .dbvmle1(d2,mu0=mu,var0=v0)
    }else{
        mu <- mean(x-y,na.rm=TRUE)
        s <- sd(x-y,na.rm=TRUE)
        out <- c(mu,s/sqrt(n))
    }
    out
}

.funmle1 <- function(paras, x){
    mu <- paras[1]; s <- paras[2]
    dF <- pnorm(x$ul, mu,s) - pnorm(x$ll, mu,s)
    sele <- dF > 0
    sum(log(dF[sele]))
}

.funmle2 <- function(paras, x,y){
    mu <- paras[1]; s <- paras[2]
    dF <- pnorm(y$ul, mu,s) - pnorm(y$ll, mu,s)
    sele <- dF > 0
    mle0 <- sum(log(dF[sele]))
    df <- dnorm(x,mu,s)
    sele <- df > 0
    mle0 <- mle0 + sum(log(df[sele]))
}


.dbvmle1 <- function(d2,mu0,var0){
    s0 <- sqrt(var0)
    n <- length(d2$ll)
    
    res <- .Fortran(.F_fitdpro1,
                    as.double(d2$ll),
                    as.double(d2$ul),
                    as.integer(n),
                    mu=as.double(mu0),
                    s=as.double(s0))
    c(res$mu,res$s/sqrt(n))
}

.dbvmle2 <- function(d1,d2,mu0,var0){
    s0 <- sqrt(var0)
    n2 <- length(d2$ll)
    n1 <- length(d1)
    res <- .Fortran(.F_fitdpro2,
                    as.double(d2$ll),
                    as.double(d2$ul),
                    as.integer(n2),
                    as.double(d1),
                    as.integer(n1),
                    mu=as.double(mu0),
                    s=as.double(s0))
    ##mu <- seq(.8*mu0, 1.2*mu0, length=50)
    ##s <- seq(s0*.9, 4*s0, length=50)
    ##musd <- expand.grid(mu=mu,sig=s)
    ##mles <- apply(musd,1,.funmle2,x=d1,y=d2)
    ##isele <- which(mles == max(mles))[1]
    ##mu <- musd[isele,1]
    ##n <- length(d2$ll) + length(d1)
    ##s <- musd[isele,2]/sqrt(n)
    ##c(mu,s)
    n <- n1+n2
    c(res$mu,res$s/sqrt(n))
}

.pairedSD2 <- function(x,y,b){
    if(length(x)!=length(y))
        stop("'x' and 'y' lengths differ")
    if(any(x > b) || any(y > b))
        stop("value(s) fall beyond data range")
    bs <- matrix(b,nrow=21,ncol=1)
    res <- apply(bs,1,.funpsd2,x=x,y=y)
    apply(res,2,median)
    #median(res)
}

.pairedSD <- function(x,y,b){
    n <- length(y)
    if(length(x) != n)
        stop("'x' and 'y' lengths differ")
    if(any(x > b) || any(y > b))
        stop("value(s) fall beyond data range")
    mux <- 1.5 * (x == b)
    x <- x + mux
    muy <- 1.5 * (y == b)
    y <- y + muy
    sx <- 0.25 + 0.5* abs(x-round(x))
    sy <- 0.25 + 0.5* abs(y-round(y))
    iter <- 201
    sig <- rep(0,iter)
    rho <- rep(0,iter)
    out <- .Fortran(.F_bootsd, as.integer(n),
                    as.double(x), as.double(y),
                    as.double(sx), as.double(sy),
                    as.integer(iter),
                    sig=as.double(sig),
                    rho=as.double(rho))
    c(median(out$sig),median(out$rho))
}

.funpsd2 <- function(b,x,y){
    n <- length(x)
    mux <- 1. * (x == b)
    muy <- 1. * (y == b)
    sx <- 0.125 + 0.25* abs(x-round(x))
    sy <- 0.125 + 0.25* abs(y-round(y))
    dx <- rnorm(n,mux,sx)
    dy <- rnorm(n,muy,sy)
    tmp <- (y+dy) - (x+dx)
    c(mean(tmp),sd(tmp))
}

.funpsd <- function(b,x,y){
    selex <- round(x)==x
    if(mean(selex) < 1.0){
        seleb <- x == b
        if(any(seleb)){
            dx <- rnorm(sum(seleb),0,0.5)
            x[seleb] <- x[seleb] + abs(dx)
            selea <- selex & !seleb
            if(any(selea)){
                dx <- rnorm(sum(selea),0,0.25)
                x[selea] <- x[selea] + abs(dx)
            }
        }else{
            dx <- rnorm(sum(selex),0,0.25)
            x[selex] <- x[selex] + dx
        }
    }else{
        n <- length(x)
        xt <- table(x)
        xi <- as.numeric(names(xt))
        xd <- mean(diff(xi))
        xsd <- xd * 0.25
        dx <- rnorm(n,0,xsd)
        x <- abs(x + dx)
    }
    xe <- x

    x <- y
    selex <- round(x)==x
    if(mean(selex) < 1.0){
        seleb <- x == b
        if(any(seleb)){
            dx <- rnorm(sum(seleb),0,0.5)
            x[seleb] <- x[seleb] + abs(dx)
            selea <- selex & !seleb
            if(any(selea)){
                dx <- rnorm(sum(selea),0,0.25)
                x[selea] <- x[selea] + abs(dx)
            }
        }else{
            dx <- rnorm(sum(selex),0,0.25)
            x[selex] <- x[selex] + dx
        }
    }else{
        n <- length(x)
        xt <- table(x)
        xi <- as.numeric(names(xt))
        xd <- mean(diff(xi))
        xsd <- xd * 0.25
        dx <- rnorm(n,0,xsd)
        x <- abs(x + dx)
    }
    ye <- x
    
    sd(xe-ye)
}

.PROnum1 <- function(x,y,group,cutoff,x.range){
    nclass <- length(cutoff) + 1
    class.names <- paste("State-",1:nclass,sep='')
    n <- length(group)
    stopifnot(length(x) == n)
    stopifnot(length(y) == n)
    if(!missing(x.range)){
        if(!is.numeric(x.range))
            stop("'x.range' must be numeric")
        if(length(x.range) != 2)
            stop("'x.range' must have length 2")
        a <- x.range[1]; b <- x.range[2]
        if(a >= b) stop("invalid range")
        if(!is.finite(a)||!is.finite(b))
            stop("Infinite range for PRO")
    }else{
        a <- min(c(x,y)); b <- max(c(x,y))
    }
    breaks <- c(a,cutoff,b)
    if(any(diff(breaks) <= 0)) stop("invalid breaks value(s)")
    if(any(x<a) || any(x>b)) stop("'x' out of range")
    if(any(y<a) || any(y>b)) stop("'y' out of range")
    tmp <- paste("x <=",cutoff)
    tmp[-1] <- paste(cutoff[-1],"<",tmp[-1])
    tmp <- c(tmp,paste(">", rev(cutoff)[1]))
  
    State.tbl <- data.frame(state=class.names, range=tmp)

    out <- .mytestnum(x=x,y=y,group=group,cutoff=cutoff,
                      class.names=class.names,
                      x.range=x.range)    
    res <- .mytestiv3(x=x,y=y,
                      breaks=breaks,group=group,
                      class.names=class.names)
    out$Results <- res
    out$Test <- .proivtest(res)
    out$states <- State.tbl
    out
}

.PROnum2 <- function(x,y,group,x.range){
    n <- length(group)
    if(length(y) != n)
            stop("'y' and 'group' have different lengths")
    if(length(x) != n)
        stop("'group' and 'x' have different lengths")
    sele1 <- is.na(x)
    sele2 <- is.na(y)
    sele3 <- is.na(group)
    sele <- sele1 | sele2 | sele3
    if(any(sele)){
        x <- x[!sele]
        y <- y[!sele]
        group <- group[!sele]
    }
    fgrp <- as.factor(as.character(group))
    if(nlevels(fgrp) != 2)
        stop("'group' must have two levels (arms)")
    sele <- fgrp == levels(fgrp)[1]
    x0 <- x[sele]
    x1 <- x[!sele]
    y0 <- y[sele]
    y1 <- y[!sele]
    fx0 <- fit.PRO(x=x0,x.range=x.range,dist='norm')
    fx1 <- fit.PRO(x=x1,x.range=x.range,dist='norm')
    fy0 <- fit.PRO(x=y0,x.range=x.range,dist='norm')
    fy1 <- fit.PRO(x=y1,x.range=x.range,dist='norm')
    if(missing(x.range)){
        b <- 10
    }else{
        b <- x.range[2]
    }
    res <- .pairedSD(x=x0,y=y0,b=b)
    sigc <- res[1]; rhoc <- res[2]
    res <- .pairedSD(x=x1,y=y1,b=b)
    sigt <- res[1]; rhot <- res[2]

    if(is.null(fx0)){
        mu0c <- NA
        s0c <- NA
    }else{
        mu0c <- fx0$x.fit$pars[1]
        s0c <- fx0$x.fit$pars[2]
    }
    if(is.null(fx1)){
        mu0t <- NA
        s0t <- NA
    }else{
        mu0t <- fx1$x.fit$pars[1]
        s0t <- fx1$x.fit$pars[2]
    }
    
    mu0 <- (mu0c+mu0t)/2
    s0 <- sqrt(0.5*(s0c^2+s0t^2))

    if(is.null(fy0)){
        mu1c <- NA
        s1c <- NA
    }else{
        mu1c <- fy0$x.fit$pars[1]
        s1c <- fy0$x.fit$pars[2]
    }
    if(is.null(fy1)){
        mu1t <- NA
        s1t <- NA
    }else{
        mu1t <- fy1$x.fit$pars[1]
        s1t <- fy1$x.fit$pars[2]
    }
    
    Mean <- c(mu0,mu1c,mu1t)
    SD <- c(s0,s1c,s1t)
    out <- as.data.frame(rbind(Mean,SD))
    names(out) <- c("T0","T1.ctrl","T1.trt")
    out
}


.prottest2 <- function(x,y,group,cutoff,x.range){
    options(warn=-2)
    group <- as.factor(as.character(group))
    if(nlevels(group) != 2) stop("only 2-arm trial data are supported")
    n <- length(group) # group cannot be missing
    if(inherits(x,"data.frame")) x <- as.matrix(x)
    if(inherits(y,"data.frame")) y <- as.matrix(y)
    
    if(any(is.na(x)) || any(is.na(y)) || any(is.na(group)))
        stop("missing value not allowed in 'x', 'y' and 'group'")
    
    if(inherits(x,"character") || inherits(y,"character")){
        stopifnot(length(x) == n)
        stopifnot(length(y) == n)
        x <- as.factor(as.character(x))
        y <- as.factor(as.character(y))
        out <- .mytestchar(x=x,y=y,group=group)
    }else if(inherits(x,"factor") || inherits(y,"factor")){
        stopifnot(length(x) == n)
        stopifnot(length(y) == n)
        out <- .mytestchar(x=x,y=y,group=group)
    }else if(inherits(x,"matrix") && inherits(y,"matrix")){
        nclass <- length(cutoff) + 1
        class.names <- paste("State-",1:nclass,sep='')
        stopifnot(nrow(x) == n)
        stopifnot(nrow(y) == n)
        if(ncol(x) != 2 && ncol(y) != 2)
            stop("'x' and 'y' must have two columns.")
        if(!missing(x.range)){
            if(!is.numeric(x.range))
                stop("'x.range' must be numeric")
            if(length(x.range) != 2)
                stop("'x.range' must have length 2")
            a <- x.range[1]; b <- x.range[2]
            if(a >= b) stop("invalid range")
            if(!is.finite(a)||!is.finite(b))
                stop("Infinite range for PRO")
        }else{
            a <- min(c(x,y)); b <- max(c(x,y))
        }
        breaks <- c(a,cutoff,b)
        if(any(diff(breaks) <= 0)) stop("invalid breaks value(s)")
        if(any(x<a) || any(x>b)) stop("'x' out of range")
        if(any(y<a) || any(y>b)) stop("'y' out of range")

        tmp <- paste("x <=",cutoff)
        tmp[-1] <- paste(cutoff[-1],"<",tmp[-1])
        tmp <- c(tmp,paste(">", rev(cutoff)[1]))
  
        if(missing(class.names)){
            class.names <- paste("state-",1:length(tmp),sep='')
        }else{
            stopifnot(length(class.names) == length(tmp))    
        }
  
        State.tbl <- data.frame(state=class.names, range=tmp)

        out <- .mytestiv2(x=x,y=y,breaks=breaks,group=group,
                          class.names=class.names)
        out$Test <- .proivtest(out)
        out$states <- State.tbl

    }else if(inherits(x,"numeric") && inherits(y,"numeric")){
        if(missing(cutoff)){
            if(!missing(x.range)){
                if(length(x.range)==2){
                    out <- .PROnum2(x=x,y=y,group=group,x.range=x.range)
                }else{
                    out <- .PROnum3(x=x,y=y,group=group,b=x.range[1])
                }
            }else{
                out <- .PROnum2(x=x,y=y,group=group,x.range=x.range)
            }
        }else{
            out <- .PROnum1(x=x,y=y,group=group,cutoff=cutoff,
                            x.range=x.range)
        }
    }else{
        stop("'x' and 'y' data types differ")
    }
    options(warn=1)
    out
}

.proivtest <- function(x){
    pv <- NULL
    m0 <- as.matrix(x$xtab0)
    m1 <- as.matrix(x$xtab1)
    n0 <- sum(m0)
    n1 <- sum(m1)
    ## improving test --------------------------
    k0 <- sum(m0[lower.tri(m0)])
    k1 <- sum(m1[lower.tri(m1)])
    tmp <- .ptestnormal(n0,n1,k0,k1)
    pv <- rbind(pv,tmp)
    ## no-change test --------------------------
    k0 <- sum(diag(m0))
    k1 <- sum(diag(m1))
    tmp <- .ptestnormal(n0,n1,k0,k1)
    pv <- rbind(pv,tmp)
    ## worsening test --------------------------
    k0 <- sum(m0[upper.tri(m0)])
    k1 <- sum(m1[upper.tri(m1)])
    tmp <- .ptestnormal(n0,n1,k0,k1)
    pv <- rbind(pv,tmp)
    ## CMH test --------------------------
    
    r0 <- apply(m0==0,1,mean)
    r1 <- apply(m1==0,1,mean)
    sele <- (r0 == 1) & (r1 == 1)
    if(any(sele)){
        m0 <- m0[-which(sele),]
        m1 <- m1[-which(sele),]
    }
    r0 <- apply(m0==0,2,mean)
    r1 <- apply(m1==0,2,mean)
    sele <- (r0 == 1) & (r1 == 1)
    if(any(sele)){
        m0 <- m0[,-which(sele)]
        m1 <- m1[,-which(sele)]
    }
    mts <- array(0, dim = c(nrow(m0), ncol(m0), 2))
    mts[,,1] <- m0
    mts[,,2] <- m1
    tmp <- mantelhaen.test(mts)
    pv <- rbind(pv,c(NA,NA,round(tmp$p.value,4)))
    
    res <- as.data.frame(pv)
    names(res) <- c("p0(%)","p1(%)","p.value")
    rownames(res) <- c("Improving","No-change","Worsening","CMH")
    res
}

.ptestnormal <- function(n0,n1,k0,k1){
    p0 <- k0/n0
    p1 <- k1/n1
    v0 <- p0*(1-p0)/n0
    v1 <- p1*(1-p1)/n1
    sp <- sqrt(v0+v1)
    if(sp > 0){
        Z <- (p1-p0)/sp
        out <- pnorm(-abs(Z))
    }else{
        out <- 1.0
    }
    c(round(100*p0,2),round(100*p1,2),round(out,4))
}

.mytestiv2 <- function(x,y,breaks,group,class.names){
    sele <- group == levels(group)[1]
    tbl1 <- .binivpro2(x=x[sele,],y=y[sele,],breaks=breaks,
                       class.names=class.names)
    tbl2 <- .binivpro2(x=x[!sele,],y=y[!sele,],breaks=breaks,
                       class.names=class.names)
    list(xtab0=tbl1, xtab1=tbl2)
}

.mytestiv3 <- function(x,y,breaks,group,class.names){
    sele <- group == levels(group)[1]
    tbl1 <- .binivpro3(x=x[sele],y=y[sele],breaks=breaks,
                       class.names=class.names)
    tbl2 <- .binivpro3(x=x[!sele],y=y[!sele],breaks=breaks,
                       class.names=class.names)
    list(xtab0=tbl1, xtab1=tbl2)
}

.mytestiv1 <- function(x,y,breaks,group,class.names){
    tbl1 <- .binivpro2(x=x,y=y,breaks=breaks,
                       class.names=class.names)
    list(xtab1=tbl1)
}

.binivpro2 <- function(x,y,breaks,class.names){
    n <- nrow(x)
    k <- length(breaks) - 1
    xdist <- matrix(0,nrow=k,ncol=k)
    for(i in 1:n){
        xi <- .binivpro1(x[i,],breaks)
        yi <- .binivpro1(y[i,],breaks)
        xyi <- xi %*% t(yi)
        xdist <- xdist + xyi
    }
    xdist <- as.data.frame(xdist)
    names(xdist) <- class.names
    rownames(xdist) <- class.names
    xdist
}

.binivpro3 <- function(x,y,breaks,class.names){
    n <- length(x)
    k <- length(breaks) - 1
    xdist <- matrix(0,nrow=k,ncol=k)
    for(i in 1:n){
        xi <- .binivpro1(x[i],breaks)
        yi <- .binivpro1(y[i],breaks)
        xyi <- xi %*% t(yi)
        xdist <- xdist + xyi
    }
    xdist <- as.data.frame(xdist)
    names(xdist) <- class.names
    rownames(xdist) <- class.names
    xdist
}

.binivpro1 <- function(x, breaks){
    k <- length(breaks) - 1
    out <- rep(0,k)
    if(length(x) == 1){
        j <- which(breaks >= x)
        out[j[1]-1]  <- 1
    }else if(length(x) == 2){
        l <- x[1]; u <- x[2]
        w <- u - l
        if(l>u){
            stop("invalid interval limits")
        }else if(l == u){
            j <- which(breaks >= l)
            out[j[1]-1]  <- 1
        }else{
            for(i in 1:k){
                if(l < breaks[i+1]){
                    if(u <= breaks[i+1]){
                        out[i] <- out[i] + (u-l)/w
                        break;
                    }else{
                        out[i] <- out[i] + (breaks[i+1]-l)/w
                        l <- breaks[i+1]
                    }
                }
            }
        }
    }
    out
}

.mytestchar <- function(x,y,group){
  if(nlevels(x)>5 || nlevels(y)>5)
    stop("too many states/levels")
  Outcome <- table(x,y,group)
  out.tbl <- ftable(Outcome)
  out.CMH <- mantelhaen.test(Outcome)
  out <- list(Results=out.tbl, Test=out.CMH)
}

.mytestnum <- function(x,y,group,cutoff,class.names,x.range){
  
  if(any(x<0 | y<0))
    stop("negative values found in 'x' or 'y'")

  ## Absolute changes ######################################
  dx <- x - y
  sele <- is.na(dx) | is.na(group)
  dx <- dx[!sele]
  g <- group[!sele]
  ## t-test
  tmp1 <- tapply(dx, g,mean,na.rm=TRUE)
  tmp2 <- t.test(dx~g)$p.value
  tout <- c(tmp1, tmp2,NA)
  ## Mann-Whitney U-test
  tmp1 <- tapply(dx, g,median,na.rm=TRUE)
  tmp2 <- wilcox.test(dx~g)  
  tout <- rbind(tout, c(tmp1,tmp2$p.value,NA))
  ## KS-test
  tmp1 <- .subedf(dx, g)
  tout <- rbind(tout, tmp1)
  ## Relative changes ######################################
  dx <- 1 - y/x
  sele <- is.na(dx) | is.na(group)
  dx <- dx[!sele]
  g <- group[!sele]
  ## t-test
  tmp1 <- tapply(dx, g,mean,na.rm=TRUE)
  tmp2 <- t.test(dx~g)$p.value
  tout <- rbind(tout, c(tmp1, tmp2,NA))
  ## Mann-Whitney U-test
  tmp1 <- tapply(dx, g,median,na.rm=TRUE)
  tmp2 <- wilcox.test(dx~g)  
  tout <- rbind(tout, c(tmp1,tmp2$p.value,NA))
  ## KS-test
  tmp1 <- .subedf(dx, g)
  tout <- rbind(tout, tmp1)
  tout <- as.data.frame(tout)
  names(tout)[3] <- "p-value"
  names(tout)[4] <- "Ref"
  row.names(tout) <- c("t-test (abs)","U-test (abs)","KS-test (abs)",
                       "t-test (rel)", "U-test (rel)","KS-test (rel)")
  #print(tout)
  
  
  out <- .subMyTest(x=x,y=y,group=group,
                    cutoff=cutoff,class.names=class.names)
  
#  if(!missing(conf.level)){
#    stopifnot(is.numeric(conf.level))  
#    stopifnot(length(conf.level) == 1)
#    stopifnot(conf.level<1 && conf.level > 0)
#    alpha <- min(conf.level,1-conf.level)
#    tmp <- .subMyCI(x=x,y=y,group=group,cutoff=cutoff,
#                    class.names=class.names,
#                    x.range=x.range,delta=delta,
#                    alpha=alpha)    
#    out.test <- data.frame(estimate=c(out$Test$p.value,out$RA),
#                           t(tmp))
#    names(out.test) <- c("estimate", "median", "lcl","ucl")
#    rownames(out.test)[1] <- "CMH.test" 
#    rownames(out.test)[4] <- "2p.test"
#    rownames(out.test)[2] <- paste("resp.rate",rownames(out.test)[2])
#    rownames(out.test)[3] <- paste("resp.rate",rownames(out.test)[3])
#    out <- list(Multi.State.Analysis=out$Results,
#                conf.level=1-alpha,Test=out.test)
#  }
  
  out$others <- tout[,-4]
  out
}

.subedf <- function(x, group){
  x0 <- unique(sort(x))
  n <- length(x0)
  g <- as.factor(as.character(group))
  sele <- g==levels(g)[1]
  x1 <- x[sele]; n1 <- length(x1)
  x2 <- x[!sele]; n2 <- length(x2)
  pv <- ks.test(x1,x2)$p.value
  Fx1 <- Fx2 <- rep(0,n)
  for(i in 1:n){
    Fx1[i] <- sum(x1<x0[i])
    Fx2[i] <- sum(x2<x0[i])
  }
  Fx1 <- Fx1/n1
  Fx2 <- Fx2/n2
  dFx <- abs(Fx1-Fx2)
  isele <- which(dFx==max(dFx))
  if(length(isele)>1){
    warning("maximum differences reached at multiple levels")
    #print(x0[isele])
    isele <- isele[1]    
  }
  c(Fx1[isele],Fx2[isele],pv,x0[isele])
}

  
.mytest <- function(cutoff,x,y,group,class.names){
  out <- .mytestnum(x=x,y=y,group=group,cutoff=cutoff,class.names=class.names)
  c(CMH.test=out$Test$p.value,out$RA)
}


.subMyTest <- function(x,y,group,cutoff,class.names){  
  if(missing(cutoff)){
    cutoff <- median(x)
    if(missing(class.names)){
      class.names <- c("1-low","2-high")    
    }else{
      if(length(class.names) ==2)
        class.names <- paste(1:2,"-", class.names, sep='')
      else
        class.names <- c("1-low","2-high")    
    }
  }else{
    if(any(diff(cutoff)<=0))
      stop("'cutoff' needs to be sorted in increasing order")
    xy <- c(x,y)
    if(min(xy) >= min(cutoff) || max(xy) <= max(cutoff))
      stop("invalid cutoff values")
    k <- length(cutoff)
    if(missing(class.names)){
      class.names <- paste("State-",c(1:(k+1)), sep='')
    }else{
      if(length(class.names) == k+1){
        class.names <- paste(1:(k+1),"-", class.names, sep='')
      }
      else{
        class.names <- paste("State-",c(1:(k+1)), sep='')
        warning("incorrect class names")        
      }
    }
  }
  
  .labelclass <- function(x,x0,z){
    k <- length(x0)
    tmp <- rep(z[1],length(x))
    for(i in 1:(k-1)){
      sele <- x>x0[i] & x<=x0[i+1]
      tmp[sele] <- z[i+1]
    }
    sele <- x>x0[k]
    tmp[sele] <- z[k+1]
    tmp
  }
  
  Before <- .labelclass(x,cutoff,class.names)
  After <- .labelclass(y,cutoff,class.names)
  
  ## Show summaries by arm/group
  #require(gmodels)
  #out.tbl2 <- NULL
  #for(i in levels(group)){
  #  sele <- group == i
    ##print(CrossTable(Before[sele], After[sele]))
    #out.tbl2 <- CrossTable(Before[sele], After[sele])
  #}
  
  ## print 3-way table for comparisons
  Outcome <- table(Before,After,group)
  ## print(ftable(Outcome))
  
  out.tbl <- ftable(Outcome)
  ## Cochran-Mantel-Haenszel chi-squared test of the null hypothesis 
  ## that two nominal variables are conditionally independent in each 
  ## stratum, assuming that there is no three-way interaction. 
  ## x is a 3 dimensional contingency table, where the last dimension 
  ## refers to the strata.
  
  out.CMH <- mantelhaen.test(Outcome)
  ##print(out.CMH)

  out <- list(Results=out.tbl, Test=out.CMH)
  
  ## Loglinear Models: https://www.statmethods.net/stats/frequencies.html
  #mydat <- data.frame(group=group,Before=Before,After=After)
  #tmp <- xtabs(~Before+After+group, data=mydat)
  #print(tmp)
  #out.loglin <- loglin(~Before+After+group, tmp) #mutual independend
  #print(out.loglin)
  #group is independent of A*B  
  #out.loglin2 <- loglin(~group+Before+After+Before*After, tmp) 
  #print(out.loglin2)
  
  ## if there are only two arms
  #ndim <- dim(Outcome)
  #if(ndim[3]==2){
  #  gnames <- colnames(Outcome[1,,])
  #  xcount <- rep(0,ndim[3])
  #  n1 <- sum(Outcome[,,1])
  #  n2 <- sum(Outcome[,,2])
  #  for(k in 1:ndim[3]){
  #    inames <- rownames(Outcome[,,k])
  #    jnames <- colnames(Outcome[,,k])
  #    xt <- Outcome[,,k]
  #    for(i in 1:ndim[1]){
  #      ilvl <- as.numeric(substr(inames[i],1,1))
  #      for(j in 1:ndim[2]){
  #        jlvl <- as.numeric(substr(jnames[j],1,1))
  #        if(ilvl > jlvl) xcount[k] <- xcount[k] + xt[i,j]
  #      }
  #    }
  #  }
    
    #cat("\nRespnder rate:")
    #cat("\n  ", gnames[1], ":", xcount[1], "/", n1,"=",
    #    round(xcount[1]/n1*100,2),"%")
    #cat("\n  ", gnames[2], ":", xcount[2], "/", n2,"=",
    #    round(xcount[2]/n2*100,2),"%")
    #ptest <- prop.test(x=xcount,n=c(n1,n2))
    #cat("\n   p-value = ", ptest$p.value)
    #cat("\n\t", ptest$method,"\n")
    
    #out.responder <- c(Rate1=xcount[1]/n1,Rate2=xcount[2]/n2,p.value=ptest$p.value)
    #names(out.responder) <- c(gnames[1], gnames[2],"p-value")
    #out <- list(Results=out.tbl, Test=out.CMH)#, Responder.Analysis=out.responder)
    #}
  
  invisible(out)
}

.PROnum3 <- function(x,y,group,b){
    n <- length(group)
    if(length(y) != n)
            stop("'y' and 'group' have different lengths")
    if(length(x) != n)
        stop("'group' and 'x' have different lengths")
    sele1 <- is.na(x)
    sele2 <- is.na(y)
    sele3 <- is.na(group)
    sele <- sele1 | sele2 | sele3
    if(any(sele)){
        x <- x[!sele]
        y <- y[!sele]
        group <- group[!sele]
    }
    fgrp <- as.factor(as.character(group))
    if(nlevels(fgrp) != 2)
        stop("'group' must have two levels (arms)")
    sele <- fgrp == levels(fgrp)[1]
    x0 <- x[sele]
    x1 <- x[!sele]
    y0 <- y[sele]
    y1 <- y[!sele]

    mu0 <- mean(x, na.rm=TRUE)
    s0 <- sd(x,na.rm=TRUE)
    mu10 <- mean(y0, na.rm=TRUE)
    mu11 <- mean(y1, na.rm=TRUE)
    s10 <- sd(y0, na.rm=TRUE)
    s11 <- sd(y1, na.rm=TRUE)
    n0 <- sum(sele); n1 <- sum(!sele)
    
    fx <- .fitPROhist(x=x0,y=x1,b=b)
    fy <- .fitPROhist(x=y0,y=y1,b=b)
    mud <- fy$mu - fx$mu
    sigd <- sqrt(fy$sig*fy$sig/n0+fx$sig*fx$sig/n1)

    zd <- mud/sigd
    pvd <- pnorm(zd)
    
    Mean <- c(mu0,mu10,mu10,mud,zd)
    SD <- c(s0,s10,s11,sigd,pvd)
    out <- as.data.frame(rbind(Mean,SD))
    names(out) <- c("T0","T1.ctrl","T1.trt","D","Z")
    out
}

.fitPROhist <- function(x,y,b){
    xb <- x >= b
    yb <- y >= b
    dxy <- y - x
    case1 <- xb & !yb
    case2 <- !xb & yb
    case0 <- !case1 & !case2
    xy0 <- dxy[case0]
    bxy <- binning(xy0)
    xbrk <- bxy$breaks
    xcnt <- bxy$freq
    nclass <- length(xcnt)
    sele <- case1
    if(any(sele)){
        z <- dxy[sele]
        k <- sum(sele)
        for(i in 1:k){
            j <- which(xbrk > z[i])[1]
            if(!is.na(j)){
                if(j==1){
                    xcnt[1] <- xcnt[1] + 1
                }else{
                    xcnt[1:(j-1)] <- xcnt[1:(j-1)] + 1/(j-1)
                }
            }
        }
    }
    sele <- case2
    if(any(sele)){
        z <- dxy[sele]
        k <- sum(sele)
        for(i in 1:k){
            j <- which(xbrk > z[i])[1]
            if(!is.na(j)){
                if(j>nclass){
                    xcnt[nclass] <- xcnt[nclass] + 1
                }else{
                    xcnt[j:nclass] <- xcnt[j:nclass] + 1/(nclass-j+1)
                }
            }else{
                xcnt[nclass] <- xcnt[nclass] + 1
            }
        }
    }
    print(xcnt)
    bxy$freq <- xcnt
    sele <- xcnt == 0
    if(any(sele)){
        bxy$freq <- xcnt[-which(sele)]
        bxy$breaks <- xbrk[-(which(sele)+1)]
    }
    print(bxy)
    x.fit <- .fit.FSD1(bxy,dist="norm")
    list(mu=x.fit$pars[1],sig=x.fit$pars[2])
}

.prottest1 <- function(x,group,x.range){
    options(warn=-2)
    group <- as.factor(as.character(group))
    if(nlevels(group) != 2) stop("only 2-arm trial data are supported")
    n <- length(group) # group cannot be missing
    
    if(any(is.na(x)) || any(is.na(group)))
        stop("missing value not allowed in 'x' and 'group'")

    if(missing(x.range)){
        stop("'x.range' must be provided")
    }else if(length(x.range) != 2){
        stop("'x.range' must be have length of 2")
    }

    if(inherits(x,"numeric")){
        out <- .PROnum4(x=x,group=group,x.range=x.range)
    }else{
        stop("'x' must be a vector")
    }
    options(warn=1)
    out
}

.PROnum4 <- function(x,group,x.range){
    n <- length(group)
    if(length(x) != n)
        stop("'group' and 'x' have different lengths")
    sele1 <- is.na(x)
    sele3 <- is.na(group)
    sele <- sele1 | sele3
    if(any(sele)){
        x <- x[!sele]
        group <- group[!sele]
    }
    fgrp <- as.factor(as.character(group))
    if(nlevels(fgrp) != 2)
        stop("'group' must have two levels (arms)")
    sele <- fgrp == levels(fgrp)[1]
    x0 <- x[sele]
    x1 <- x[!sele]
    n0 <- length(x0)
    n1 <- length(x1)
    fx0 <- fit.PRO(x=x0,x.range=x.range,dist='norm')
    fx1 <- fit.PRO(x=x1,x.range=x.range,dist='norm')

    if(is.null(fx0)){
        mu0c <- NA
        s0c <- NA
    }else{
        mu0c <- fx0$x.fit$pars[1]
        s0c <- fx0$x.fit$pars[2]
    }
    if(is.null(fx1)){
        mu0t <- NA
        s0t <- NA
    }else{
        mu0t <- fx1$x.fit$pars[1]
        s0t <- fx1$x.fit$pars[2]
    }
    
    mud <- mu0t-mu0c
    sdd <- sqrt(s0c^2/n0+s0t^2/n1)

    Z2 <- mud/sdd
    pv1 <- pnorm(Z2)
    
    Mean <- c(mu0c,mu0t,mud,Z2)
    SD <- c(s0c,s0t,sdd,pv1)
    out <- as.data.frame(rbind(Mean,SD))
    names(out) <- c("ctrl","trt","D","Z")
    out
}
