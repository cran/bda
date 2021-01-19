.bootPRO.mean <- function(x,iter=999,conf.level=.95,mcid=NULL){
    .rmean <- function(x)
    {
        y <- as.numeric(x[-1])
        n <- length(y)
        n1 <- round(n/2)
        x1 <- y[1:n1]
        x2 <- y[(n1+1):n]
        res <- NA
        if(any(!is.na(x1)) && any(!is.na(x2))){
            x1 <- x1[!is.na(x1)]
            x2 <- x2[!is.na(x2)]
            xbar1 <- sample(x1,size=1,replace=TRUE)
            xbar2 <- sample(x2,size=1,replace=TRUE)
            res <- xbar1 - xbar2
        }
        res
    }

    .diff.means <- function(data,indices)
    {
        d <- data[indices,]
        tmp <- apply(d,1,.rmean)
        rd <- as.numeric(tmp)
        dbar <- tapply(rd, d[,1], mean, na.rm=TRUE)
        delta <- diff(dbar)
        as.numeric(c(delta,dbar))
    }

    .rmean2 <- function(x)
    {
        mcid <- x[1]
        y <- as.numeric(x[-c(1,2)])
        n <- length(y)
        n1 <- round(n/2)
        x1 <- y[1:n1]
        x2 <- y[(n1+1):n]
        res <- NA
        if(any(!is.na(x1)) && any(!is.na(x2))){
            x1 <- x1[!is.na(x1)]
            x2 <- x2[!is.na(x2)]
            xbar1 <- sample(x1,size=1,replace=TRUE)
            xbar2 <- sample(x2,size=1,replace=TRUE)
            tmp <- xbar1 - xbar2
            if(tmp >= mcid)
                res <- 1
            else
                res <- 0
        }
        res
    }
    
    .diff.means2 <- function(data,indices)
    {
        d <- data[indices,]
        tmp <- apply(d,1,.rmean2)
        rd <- as.numeric(tmp)
        dbar <- tapply(rd, d[,2], mean, na.rm=TRUE)
        delta <- diff(dbar)
        as.numeric(c(delta,dbar))
    }
    
    if(is.null(mcid)){
        out1 <- boot::boot(x, .diff.means, R = round(iter),
                     stype = "i", strata = x[,1])
    }else{
        x2 <- cbind(mcid,x)
        out1 <- boot::boot(x2, .diff.means2, R = round(iter),
                     stype = "i", strata = x2[,2])
    }
    
    out2 <- boot::boot.ci(out1, conf.level=conf.level)
    n <- ncol(x)
    xb1 <- apply(x[,c(2:((1+n)/2))],1,mean,na.rm=TRUE)
    xb2 <- apply(x[,c(((1+n)/2)+1):n],1,mean,na.rm=TRUE)
    d <- xb1 - xb2
    ##Normal <- c(ll=out2$normal[2],ul=out2$normal[3])
    ##Basic <- c(ll=out2$basic[4],ul=out2$basic[5])
    ##Studentized <- c(ll=out2$student[4],ul=out2$student[5])
    ##Percentile <- c(ll=out2$percent[4],ul=out2$percent[5])
    ##BCa <- c(ll=out2$bca[4],ul=out2$bca[5])
    ##Test <- data.frame(Normal,Basic,Studentized, Percentile,BCa);
    if(is.null(mcid)){
        Mean <- as.numeric(tapply(d,x[,1],mean,na.rm=TRUE))
        Mean <- c(diff(Mean),Mean)
        Bias <- apply(out1$t,2,mean)
        SE <- apply(out1$t,2,sd)
        se0 <- tapply(d,x[,1],sd,na.rm=TRUE)/sqrt(tapply(d,x[,1],length))
        se0 <- c(sqrt(sum(se0^2)),se0)
        Stat <- data.frame(mean.orig=Mean,se.orig=se0,
                           mean.boot=Bias,se.boot=SE) 
        rownames(Stat) <- c("Diff",levels(x[,1]))
        cat("\nAbsolute changes (t0-t1):\n")
        print(round(Stat,3))
        if(nlevels(x[,1])==2){
            pv1 <- t.test(d~x[,1])$p.value
            pv2 <- wilcox.test(d~x[,1])$p.value
            
            Diff <- data.frame(t.test=round(pv1,4),
                               Wilcox.test=round(pv2,4),
                               BCa.ll=round(out2$bca[4],2),
                               BCa.ul=round(out2$bca[5],2))
            rownames(Diff) <- paste(levels(x[,1])[1],"-",levels(x[,1])[2],sep='')
            cat("\nComparisons:\n")
            print(Diff)
        }
    }else{## responder analysis
        d <- 1.0 * (d > mcid)
        Mean <- as.numeric(tapply(d,x[,1],mean,na.rm=TRUE))
        Mean <- c(diff(Mean),Mean)
        Bias <- apply(out1$t,2,mean)
        SE <- apply(out1$t,2,sd)
        Stat <- data.frame(rate.orig=Mean,
                           rate.boot=Bias,se.boot=SE) 
        rownames(Stat) <- c("Diff",levels(x[,1]))
        cat("\nResponder analysis (absolute changes t0-t1):\n")
        print(round(Stat,3))
        if(nlevels(x[,1])==2){
            sele <- !is.na(d)
            n <- tapply(d[sele],x[sele,1],length)
            fx <- tapply(d[sele],x[sele,1],sum,na.rm=TRUE)
            ##print(cbind(fx,n))
            pv1 <- prop.test(x=fx,n=n)$p.value
            TeaTasting <-
                matrix(c(fx,n-fx), nrow = 2)
            pv2 <- fisher.test(TeaTasting)$p.value

            Diff <- data.frame(prop.test=round(pv1,4),
                               fisher.test=round(pv2,4),
                               BCa.ll=round(out2$bca[4],2),
                               BCa.ul=round(out2$bca[5],2))
            rownames(Diff) <- paste(levels(x[,1])[1],"-",
                                    levels(x[,1])[2],sep='')
            cat("\nComparisons:\n")
            print(Diff)
        }

    }
                       
    ##cat("\n", out2$bca[1]*100,"% Confidence Intervals:\n",sep="")
    ##print(round(Test,3))
    ##out <- list(Statistics=Stat, CI=Test)
    invisible(NULL)
}

.bootPRO.prop <- function(x,iter=999,conf.level=.95,mcid=NULL){
    .rmean <- function(x)
        {
            y <- as.numeric(x[-1])
            n <- length(y)
            n1 <- round(n/2)
            x1 <- y[1:n1]
            x2 <- y[(n1+1):n]
            res <- NA
            if(any(!is.na(x1)) && any(!is.na(x2))){
                x1 <- x1[!is.na(x1)]
                x2 <- x2[!is.na(x2)]
                xbar1 <- sample(x1,size=1,replace=TRUE)
                xbar2 <- sample(x2,size=1,replace=TRUE)
                if(xbar1==0)
                    res <- -10
                else
                    res <- 1 - xbar2/xbar1
            }
            res
        }
    
    .diff.means <- function(data,indices)
    {
        d <- data[indices,]
        tmp <- apply(d,1,.rmean)
        rd <- as.numeric(tmp)
        dbar <- tapply(rd, d[,1], mean, na.rm=TRUE)
        delta <- diff(dbar)
        as.numeric(c(delta,dbar))
    }
    
    .rmean2 <- function(x)
    {
        mcid <- x[1]
        y <- as.numeric(x[-c(1,2)])
        n <- length(y)
        n1 <- round(n/2)
        x1 <- y[1:n1]
        x2 <- y[(n1+1):n]
        res <- NA
        if(any(!is.na(x1)) && any(!is.na(x2))){
            x1 <- x1[!is.na(x1)]
            x2 <- x2[!is.na(x2)]
            xbar1 <- sample(x1,size=1,replace=TRUE)
            xbar2 <- sample(x2,size=1,replace=TRUE)
            if(xbar1==0)
                tmp <- -10
            else
                tmp <- 1 - xbar2/xbar1
            
            if(tmp > mcid)
                res <- 1
            else
                res <- 0
        }
        res
    }
    
    .diff.means2 <- function(data,indices)
    {
        d <- data[indices,]
        tmp <- apply(d,1,.rmean2)
        rd <- as.numeric(tmp)
        dbar <- tapply(rd, d[,2], mean, na.rm=TRUE)
        delta <- diff(dbar)
        as.numeric(c(delta,dbar))
    }
    
    if(is.null(mcid)){
        out1 <- boot::boot(x, .diff.means, R = round(iter),
                     stype = "i", strata = x[,1])
    }else{
        x2 <- cbind(mcid,x)
        out1 <- boot::boot(x2, .diff.means2, R = round(iter),
                     stype = "i", strata = x2[,2])
    }
    
    out2 <- boot::boot.ci(out1, conf.level=conf.level)
    n <- ncol(x)
    xb1 <- apply(x[,c(2:((1+n)/2))],1,mean,na.rm=TRUE)
    xb2 <- apply(x[,c(((1+n)/2)+1):n],1,mean,na.rm=TRUE)
    d <- 1 - xb2/xb1
    ##Normal <- c(ll=out2$normal[2],ul=out2$normal[3])
    ##Basic <- c(ll=out2$basic[4],ul=out2$basic[5])
    ##Studentized <- c(ll=out2$student[4],ul=out2$student[5])
    ##Percentile <- c(ll=out2$percent[4],ul=out2$percent[5])
    BCa <- c(ll=out2$bca[4],ul=out2$bca[5])
    ##Test <- data.frame(Normal,Basic,Studentized, Percentile,BCa);
    if(is.null(mcid)){
        Mean <- as.numeric(tapply(d,x[,1],mean,na.rm=TRUE))
        Mean <- c(diff(Mean),Mean)
        Bias <- apply(out1$t,2,mean)
        SE <- apply(out1$t,2,sd)
        se0 <- tapply(d,x[,1],sd,na.rm=TRUE)/sqrt(tapply(d,x[,1],length))
        se0 <- c(sqrt(sum(se0^2)),se0)
        Stat <- data.frame(mean.orig=Mean,se.orig=se0,
                           mean.boot=Bias,se.boot=SE) 
        rownames(Stat) <- c("Diff",levels(x[,1]))
        cat("\nRelative changes (1-t1/t0):\n")
        print(round(Stat,3))
        if(nlevels(x[,1])==2){
            pv1 <- t.test(d~x[,1])$p.value
            pv2 <- wilcox.test(d~x[,1])$p.value
            
            Diff <- data.frame(t.test=round(pv1,4),
                               Wilcox.test=round(pv2,4),
                               BCa.ll=round(out2$bca[4],2),
                               BCa.ul=round(out2$bca[5],2))
            rownames(Diff) <- paste(levels(x[,1])[1],"-",levels(x[,1])[2],sep='')
            cat("\nComparisons:\n")
            print(Diff)
        }
    }else{## responder analysis
        d <- 1.0 * (d > mcid)
        Mean <- as.numeric(tapply(d,x[,1],mean,na.rm=TRUE))
        Mean <- c(diff(Mean),Mean)
        Bias <- apply(out1$t,2,mean)
        SE <- apply(out1$t,2,sd)
        Stat <- data.frame(mean.orig=Mean,
                           mean.boot=Bias,se.boot=SE) 
        rownames(Stat) <- c("Diff",levels(x[,1]))
        cat("\nResponder analysis (relative changes 1-t1/t0):\n")
        print(round(Stat,3))
        if(nlevels(x[,1])==2){
            sele <- !is.na(d)
            n <- tapply(d[sele],x[sele,1],length)
            fx <- tapply(d[sele],x[sele,1],sum,na.rm=TRUE)
            ##print(cbind(fx,n))
            pv1 <- prop.test(x=fx,n=n)$p.value
            TeaTasting <-
                matrix(c(fx,n-fx), nrow = 2)
            pv2 <- fisher.test(TeaTasting)$p.value
            
            Diff <- data.frame(t.test=round(pv1,4),
                               Wilcox.test=round(pv2,4),
                               BCa.ll=round(out2$bca[4],2),
                               BCa.ul=round(out2$bca[5],2))
            rownames(Diff) <- paste(levels(x[,1])[1],"-",
                                    levels(x[,1])[2],sep='')
            cat("\nComparisons:\n")
            print(Diff)
        }
    }
                       
    ##cat("\n", out2$bca[1]*100,"% Confidence Intervals:\n",sep="")
    ##print(round(Test,3))
    ##out <- list(Statistics=Stat, CI=Test)
    invisible(NULL)
}

bootPRO <- function(x,type="relative",MCID,iter=999,conf.level=0.95){
    if ((length(conf.level) != 1L) ||
        is.na(conf.level) ||
        (conf.level <= 0) ||
        (conf.level >= 1)) 
        stop("'conf.level' must be a single number between 0 and 1")
    if ((length(iter) != 1L) ||
        is.na(iter) ||
        (iter < 1)) 
        stop("'conf.level' must be a single number between 0 and 1")
    iter <- round(iter)
    type <- match.arg(tolower(type),
                      c("absolute","relative","percent","mean",
                        "average"))
    ## check data format: (1) col1 gives the treatment group; (2) col2
    ## to the end gives the PRO measures from diary. Equal number of
    ## repeated measures are expected for each subject. If not, fill
    ## with NA's.
    x[,1] <- as.factor(x[,1])
    nl2 <- nlevels(x[,1]) 
    if(nl2 != 2)
        stop("there must be two arms (groups)")
    if(missing(MCID)){
        mcid <- NULL
    }else{
        if(!is.numeric(MCID)){
            warning("MCID must be numeric")
            mcid <- NULL
        }else{
            mcid <- MCID[1]
            if(mcid<0){
                warning("MCID must be positive")
                mcid <- NULL
            }
        }
    }
    options(warn=-4)
    if(type=="absolute"||
       type=="mean"||
       type=="average")
        .bootPRO.mean(x,iter=iter,conf.level=conf.level,mcid=mcid)
    else
        .bootPRO.prop(x,iter=iter,conf.level=conf.level,mcid=mcid)

    options(warn=0)

    invisible(NULL)
}
