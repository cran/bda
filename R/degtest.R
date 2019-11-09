## DEG.Test: 2018/07/10

## this test is developed specifically for DEG detection for NGS
## data. We preprocessed the data and feed the algorithm so this
## algorithm works for other types of data.

deg.test <- function(normal.x, normal.y,counts.x, counts.y, gene,
                     weights.x, weights.y){
    ## we assume the counts data are available, and we normalize the
    ## data using certain normalization algorithms. Rank data are
    ## generated based on rank data (or normalized data when counts
    ## data are not available.
    if(missing(counts.x)||missing(counts.y)){
        if(missing(normal.x)||missing(normal.y))
            stop("both normalized and counts data are incomplete")
        ## otherwise if both normalized datasets are available...
        .CheckData(normal.x, normal.y, type=1)
        rank.x <- .Convert2Ranks(normal.x)
        rank.y <- .Convert2Ranks(normal.y)
        counts.x <- NULL
        counts.y <- NULL
        weighs.x <- rep(1, ncol(normal.x))
        weighs.y <- rep(1, ncol(normal.y))
        cat("\nCounts (X):\tnot available")
        cat("\nCounts (Y):\tnot available")
        tmp <- deparse(substitute(normal.x))
        cat("\nNormalized (X):\t", tmp)
        tmp <- deparse(substitute(normal.y))
        cat("\nNormalized (Y):\t", tmp)
    }else{
        tmp <- deparse(substitute(counts.x))
        cat("\nCounts (X):\t", tmp)
        tmp <- deparse(substitute(counts.y))
        cat("\nCounts (Y):\t", tmp)
        .CheckData(counts.x, counts.y, type=1)
        ## normalize counts data using Q3-like normalization
        if(missing(normal.x)){# normalize the data
            normal.x <- .Q3Normalize(counts.x)
            cat("\nNormalized (X):\t upper quartile normalized")
        }else{ # check row names and dim
            .CheckData(normal.x, counts.x, type=2)
            tmp <- deparse(substitute(counts.x))
            cat("\nNormalized (X):\t", tmp)
        }

        ## normalize counts data using Q3-like normalization
        if(missing(normal.y)){# normalize the data
            normal.y <- .Q3Normalize(counts.y)
            cat("\nNormalized (Y):\t upper quartile normalized")
        }else{ # check row names and dim
            .CheckData(normal.y, counts.y, type=2)
            tmp <- deparse(substitute(normal.y))
            cat("\nNormalized (Y):\t", tmp)
        }
        rank.x <- .Convert2Ranks(counts.x)
        rank.y <- .Convert2Ranks(counts.y)

        qx3 <- .upperqtl(counts.x)
        qy3 <- .upperqtl(counts.y)
        qtl3 <- max(qx3, qy3)
        
        if(missing(weights.x)){
            weights.x <- apply(counts.x, 2, quantile, probs=qtl3, names=FALSE)
        }else{
            if(is.logical(weights.x)){
                if(weights.x){
                    weights.x <- apply(counts.x, 2, quantile, probs=qtl3, names=FALSE)
                }else{
                    weights.x <- rep(0, ncol(counts.x))
                }
            }else{
                stopifnot(is.numeric(weights.x))
                ntmp <- length(weights.x)
                stopifnot(ncol(counts.x)==ntmp)
                stopifnot(all(weights.x > 0))
            }
        }
        
        if(missing(weights.y)){
            weights.y <- apply(counts.y, 2, quantile, probs=qtl3, names=FALSE)
        }else{
            if(is.logical(weights.y)){
                if(weights.y){
                    weights.y <- apply(counts.y, 2, quantile, probs=qtl3, names=FALSE)
                }else{
                    weights.y <- rep(0, ncol(counts.y))
                }
            }else{
                stopifnot(is.numeric(weights.y))
                ntmp <- length(weights.y)
                stopifnot(ncol(counts.y)==ntmp)
                stopifnot(all(weights.y > 0))
            }
        }
    }
    
    gnames <- row.names(normal.x)
    
    if(missing(gene)){
        genes <- gnames
        Batch <- TRUE
    }else{
        if(length(gene) > 1){
            genes <- gene
            Batch <- TRUE
        }else{
            Batch <- FALSE
        }
    }
    
    out <- NULL
    if(Batch){
        res <- matrix(0, nrow=length(genes),ncol=12)
        mypv <- NULL;
        method <- NULL
        i <- 0; isele <- NULL
        for(gene in genes){
            i <- i + 1; j <- 1
            irow <- match(gene, gnames, nomatch=0)[1]
            if(irow == 0){
                cat("\nrow=", irow, ", gene=", gene, "not found!")
            }else{
                cat("\nrow=", irow, ", gene=", gene, "processing ...")
                isele <- c(isele, irow)
                tmp <- .geneDEG3(x1=as.numeric(normal.x[irow,]),
                                y1=as.numeric(normal.y[irow,]),
                                x2=as.numeric(counts.x[irow,]),
                                y2=as.numeric(counts.y[irow,]),
                                x3=as.numeric(rank.x[irow,]),
                                y3=as.numeric(rank.y[irow,]),
                                wx=weights.x,
                                wy=weights.y,
                                gene=gene)
                mypv <- c(mypv, tmp$p.value)
                method <- c(method, tmp$method)
                res[i,j] <- tmp$min.x; j <- j + 1;
                res[i,j] <- tmp$max.x; j <- j + 1;
                res[i,j] <- tmp$min.y; j <- j + 1;
                res[i,j] <- tmp$max.y; j <- j + 1;
                res[i,j] <- tmp$mean.x; j <- j + 1;
                res[i,j] <- tmp$sd.x; j <- j + 1;
                res[i,j] <- tmp$mean.y; j <- j + 1;
                res[i,j] <- tmp$sd.y; j <- j + 1;
                res[i,j] <- tmp$mean.rank.x; j <- j + 1;
                res[i,j] <- tmp$sd.rank.x; j <- j + 1;
                res[i,j] <- tmp$mean.rank.y; j <- j + 1;
                res[i,j] <- tmp$sd.rank.y; j <- j + 1;
            }
        }
        cat("\n\nDONE (results in obj$table).\n")
        output <- data.frame(res,
                             p.adj=p.adjust(mypv,method='fdr'),
                             method=method)
        row.names(output) <- gnames[isele]
        nams <- c("Min.X","Max.X",
                  "Min.Y","Max.Y",
                  "Mean.X","SD.X",
                  "Mean.Y","SD.Y",
                  "Mean.Rank.X","SD.Rank.X",
                  "Mean.Rank.Y","SD.Rank.Y",
                  'p.value', 'Method')
        names(output) <- nams
        out$table <- output
        out$gene <- NULL
        out <- structure(out, class='DEG')

    }else{
        irow <- match(gene, gnames, nomatch=0)[1]
        if(irow==0){
            cat("\nrow=", irow, ", gene=", gene, "not found!")
        }else{
            cat("\nrow=", irow, ", gene=", gene, "processing ...")
            out <- .geneDEG3(x1=as.numeric(normal.x[irow,]),
                            y1=as.numeric(normal.y[irow,]),
                            x2=as.numeric(counts.x[irow,]),
                            y2=as.numeric(counts.y[irow,]),
                            x3=as.numeric(rank.x[irow,]),
                            y3=as.numeric(rank.y[irow,]),
                            wx=weights.x,
                            wy=weights.y,
                            gene=gene)
            out2 <- .mlnormtest(gene=gene,
                                x1=as.numeric(normal.x[irow,]),
                                y1=as.numeric(normal.y[irow,]),
                                x2=as.numeric(counts.x[irow,]),
                                y2=as.numeric(counts.y[irow,]),
                                x3=as.numeric(rank.x[irow,]),
                                y3=as.numeric(rank.y[irow,]))
            out$fx <- out2$fx
            out$fy <- out2$fy
            if(out$p.value > out2$p.value){
                out$p.value <- out2$p.value
                out$method <- out2$method
            }
        }
        cat("Done.\n")
    }
    invisible(out)
}

.upperqtl <- function(x){
    tmp <- apply(x==min(x), 2, mean, na.rm=TRUE)
    p0max <- max(tmp)
    if(p0max > 0.99)
        stop("some profile(s) contain more than 99% zeros")
    if(p0max < 0.75){
        pn <- 0.75
    }else{
        pn <- 0.5 * p0max + 0.5
    }
    pn
}

.CheckData <- function(x,y,type=2){
    ## type=1 doesn't check numbers of columns (profiles).  Type=2
    ## check data that were transformed from each other.
    if(any(is.na(x)))
        stop("'x' contains missing values")
    if(any(is.na(y)))
        stop("'y' contains missing values")
    if(type==2){
        gnamesx <- row.names(x)
        gnamesy <- row.names(y)
        if(any(gnamesx != gnamesy))
            stop("gene name not matched")
        if(ncol(x) != ncol(y)){
            stop("profile sizes not matched")
        }
    }else{
        gnamesx <- row.names(x)
        gnamesy <- row.names(y)
        if(any(gnamesx != gnamesy))
            stop("gene name not matched")
    }
    NULL
}

.Convert2Ranks <- function(x){
    ## compute ranks and normalize to [0,1]
    tmp <- apply(x, 2, rank)
    for(i in 1:ncol(x)){
        y <- as.numeric(tmp[,i])
        xrng <- range(y)
        y <- (y-xrng[1])/(diff(xrng))
        x[,i] <- y
    }
    x
}

.deg.test <- function(normal.x, normal.y,counts.x, counts.y, gene){
    ##raw.x=normal.x, raw.y=normal.y,
    ##counts.x=counts.x, counts.y=counts.y,
    ##gene=NULL){
    if(missing(counts.x)||missing(counts.y)){
        out <- .degtest2(normal.x=normal.x,
                         normal.y=normal.y,
                         gene=gene)
    }else{
        out <- .degtest1(counts.x=counts.x,
                         counts.y=counts.y,
                         normal.x=normal.x,
                         normal.y=normal.y,
                         gene=gene)
    }
    
    out
}

.degtest2 <- function(normal.x, normal.y,
                      gene=NULL){
    stopifnot(all(normal.x>=0))
    stopifnot(all(normal.y>=0))
    xnam <- deparse(substitute(normal.x))
    ynam <- deparse(substitute(normal.y))
    nx <- ncol(normal.x)
    ny <- ncol(normal.y)
    ngene <- nrow(normal.x)
    if(nrow(normal.y) != ngene)
        stop("gene names not match between ", xnam, " and ", ynam)
        
    gnames <- row.names(normal.x)
    if(any(row.names(normal.y) != gnames))
        stop("gene names no match between ",xnam, " and ", ynam)

    rank.x <- apply(normal.x, 2, rank)
    rank.y <- apply(normal.y, 2, rank)

    if(is.null(gene)){
        pv1 <- rep(NA, ngene)
        pv2 <- pv3 <- pv1
        max1 <- max2 <- pv1
        zero1 <- zero2 <- pv1
        mean1c <- mean2c <- mean1 <- mean2 <- mean1r <- mean2r <- pv1
        sd1c <- sd2c <- sd1 <- sd2 <- sd1r <- sd2r <- pv1
        ng1 <- ng2 <- pv1
        
        for(x.sele in 1:ngene){
            gene <- gnames[x.sele]
            out <- .degtest2fun(gene,xnam,ynam,
                            x2=as.numeric(normal.x[x.sele,]),
                            y2=as.numeric(normal.y[x.sele,]),
                            x3=as.numeric(rank.x[x.sele,]),
                            y3=as.numeric(rank.y[x.sele,]))
            zero1[x.sele] <- out$stat[1,2]
            zero2[x.sele] <- out$stat[2,2]
            max1[x.sele] <- out$stat[1,3]
            max2[x.sele] <- out$stat[2,3]

            mean1c[x.sele] <- out$test[1,1]
            mean2c[x.sele] <- out$test[1,3]
            sd1c[x.sele] <- out$test[1,2]
            sd2c[x.sele] <- out$test[1,4]
            pv1[x.sele] <- out$test[1,5]

            mean1[x.sele] <- out$test[2,1]
            mean2[x.sele] <- out$test[2,3]
            sd1[x.sele] <- out$test[2,2]
            sd2[x.sele] <- out$test[2,4]
            pv2[x.sele] <- out$test[2,5]

            mean1r[x.sele] <- out$test[3,1]
            mean2r[x.sele] <- out$test[3,3]
            sd1r[x.sele] <- out$test[3,2]
            sd2r[x.sele] <- out$test[3,4]
            pv3[x.sele] <- out$test[3,5]

            ng1[x.sele] <- round(out$fit.x$npar + 1)/3
            ng2[x.sele] <- round(out$fit.y$npar + 1)/3
        }
        
        out <- data.frame(Zeros.x=zero1, Zeros.y=zero2,
                          Max.x=max1, Max.y=max2,
                          MeanCount.x=mean1c,
                          SDCount.x=sd1c,
                          MeanCount.y=mean2c,
                          SDCount.x=sd2c,
                          p.Fisher=pv1,
                          Mean.x=mean1,
                          SD.x=sd1,
                          Mean.y=mean2,
                          SD.y=sd2,
                          p.Perm=pv2,
                          MeanRank.x=mean1r,
                          SDRank.x=sd1r,
                          MeanRank.y=mean2r,
                          SDRank.y=sd2r,
                          p.Perm.Rank=pv3,
                          ng.x=ng1, ng.y=ng2
                          )
        row.names(out) <- gnames
    }else{
        x.sele <- match(gene, gnames)[1]
        if(is.na(x.sele)){
            warning("gene ", gene, "not found")
            out <- NULL
        }else{
            out <- .degtest2fun(gene,xnam,ynam,
                            x2=as.numeric(normal.x[x.sele,]),
                            y2=as.numeric(normal.y[x.sele,]),
                            x3=as.numeric(rank.x[x.sele,]),
                            y3=as.numeric(rank.y[x.sele,]))
        }
    }
    invisible(out)
}

.degtest1fun <- function(gene, xnam1,ynam1,xnam2,ynam2,
                         x1,y1,x2,y2,x3,y3){
    tmp <- table(x1)
    xt1 <- as.numeric(names(tmp))
    xf1 <- as.numeric(tmp)
    nx1 <- max(x1)
    tmp <- table(y1)
    yt1 <- as.numeric(names(tmp))
    yf1 <- as.numeric(tmp)
    ny1 <- max(y1)
    xmax <- max(x1)
    ymax <- max(y1)
    if(xmax > 10){
        xmax <- max(x2)
    }
    if(ymax > 10){
        ymax <- max(y2)
    }
    
    stat <- data.frame(Size=c(sum(xf1), sum(yf1)),
                       Zeros=c(sum(x1==0),sum(y1==0)),
                       Max=c(xmax,ymax)
                       )
    row.names(stat) <- c(xnam2, ynam2)
    test <- NULL
    xytab <- NULL; rnam <- NULL
    if(nx1 < 4 || ny1 < 4){
        for(i in 1:4){
            tmp <- c(sum(x1==i-1), sum(y1==i-1))
            if(any(tmp > 0)){
                xytab <- cbind(xytab, tmp)
                rnam <- c(rnam, as.character(i-1))
            }
        }
        tmp <- c(sum(x1>=4), sum(y1>=4))
        if(any(tmp > 0)){
            xytab <- cbind(xytab, tmp)
            rnam <- c(rnam, "4+")
        }
        xytab <- as.data.frame(xytab)
        names(xytab) <- rnam
        row.names(xytab) <- c(xnam1, ynam1)
    }else if(nx1 <= 10 && ny1 <= 10){
        for(i in 1:4){
            tmp <- c(sum(x1==i-1), sum(y1==i-1))
            if(any(tmp > 0)){
                xytab <- cbind(xytab, tmp)
                rnam <- c(rnam, as.character(i-1))
            }
        }
        tmp <- c(sum(x1>=4), sum(y1>=4))
        if(any(tmp > 0)){
            xytab <- cbind(xytab, tmp)
            rnam <- c(rnam, "4+")
        }
        xytab <- as.data.frame(xytab)
        names(xytab) <- rnam
        row.names(xytab) <- c(xnam1, ynam1)
    }else{
        xytab <- .binxy(x2,y2)
        row.names(xytab) <- c(xnam2, ynam2)
    }
    if(dim(xytab)[2] == 1){
        pv1 <- NA
        pv2 <- NA
        pv3 <- NA
    }else if(dim(xytab)[2] > 3){
        pv1 <- fisher.test(xytab, alternative='two.sided',
                           simulate.p.value=TRUE)$p.value
        pv2 <- perm.test(x2, y2, iter=99999)$p.value
        pv3 <- perm.test(x3, y3, iter=99999)$p.value
    }else{
        pv1 <- fisher.test(xytab, alternative='two.sided')$p.value
        pv2 <- perm.test(x2, y2, iter=99999)$p.value
        pv3 <- perm.test(x3, y3, iter=99999)$p.value
    }

       
    mu.x1 <- mean(x1)
    mu.x2 <- mean(x2)
    mu.x3 <- mean(x3)
    sd.x1 <- sd(x1)
    sd.x2 <- sd(x2)
    sd.x3 <- sd(x3)

    mu.y1 <- mean(y1)
    mu.y2 <- mean(y2)
    mu.y3 <- mean(y3)
    sd.y1 <- sd(y1)
    sd.y2 <- sd(y2)
    sd.y3 <- sd(y3)

    test <- data.frame(
        mu1=round(c(mu.x1,mu.x2,mu.x3),3),
        sd1=round(c(sd.x1,sd.x2,sd.x3),3),
        mu2=round(c(mu.y1,mu.y2,mu.y3),3),
        sd2=round(c(sd.y1,sd.y2,sd.y3),3),
        pv=signif(c(pv1,pv2,pv3),6))
    row.names(test) <- c("Counts","normalized","Rank")
    names(test) <- c(paste(xnam2,".",c("Mean","SD"),sep=''),
                     paste(ynam2,".",c("Mean","SD"),sep=''),
                     "p.value")

    nx1 <- length(xf1)
    if(nx1<3){
        fx <- fit.mlnorm(x1, method='mop')
    }else if(nx1 <= 5){
        fx <- fit.mlnorm(x2, k=1, method='lognormal')
    }else{
        fx <- fit.mlnorm(x2, k=3, method='fnmm')
    }
        
    ny1 <- length(yf1)
    if(ny1<3){
        fy <- fit.mlnorm(y1, method='mop')
    }else if(ny1 <= 5){
        fy <- fit.mlnorm(y2, k=1, method='lognormal')
    }else{
        fy <- fit.mlnorm(y2, k=3, method='fnmm')
    }
    
    list(gene=gene, stat = stat, counts = xytab,
         test = test, fit.x=fx, fit.y=fy
         )
}

.binxy <- function(x,y){
    nx0 <- sum(x==0)
    ny0 <- sum(y==0)
    x <- x[x>0]; y <- y[y>0]; xy <- c(x,y)
    xy <- log(xy); x <- log(x); y <- log(y)
    xhist <- hist(xy, nclass=6, plot=FALSE)
    sele <- xhist$counts == 0
    xbreaks <- xhist$breaks
    if(any(sele)){
        isele <- which(sele) + 1
        xbreaks <- xbreaks[-isele]
    }
    x1 <- hist(x, breaks=xbreaks, plot=FALSE)$counts
    y1 <- hist(y, breaks=xbreaks, plot=FALSE)$counts
    k <- length(xbreaks)
    l1 <- exp(xbreaks[-k])
    l2 <- exp(xbreaks[-1])
    if(nx0 <5 && ny0 <5){
        x1[1] <- x1[1] + nx0
        y1[1] <- y1[1] + ny0
        l1[1] <- 0;
        l2[1] <- 0;
        k <- k - 1
    }else{
        x1 <- c(nx0, x1)
        y1 <- c(ny0, y1)
        l1 <- c(0,l1)
        l2 <- c(0,l2)
    }

    out <- rbind(x1,y1)
    while(k > 2){
        xmin <- min(out)
        if(xmin == 0){
            tmp <- apply(out == xmin,2,sum)
            sele <- tmp > 0
            isele <- which(sele)
            xsum <- apply(out,2,sum)
            if(length(isele) > 1){
                isele2 <- which(xsum == min(xsum[isele]))[1]
            }else{
                isele2 <- isele
            }
            if(isele2==1){
                out[,2] <- out[,2]+out[,1]
                out <- out[,-1]
                l1 <- l1[-2]
                l2 <- l2[-2]
            }else if(isele2==k){
                out[,k-1] <- out[,k-1]+out[,k]
                out <- out[,-k]
                l1 <- l1[-(k-1)]
                l2 <- l2[-(k-1)]
            }else{
                sum1 <- xsum[isele2-1]
                sum2 <- xsum[isele2+1]
                if(sum1 < sum2){
                    out[,isele2-1] <- out[,isele2-1]+out[,isele2]
                    out <- out[,-isele2]
                    l1 <- l1[-isele2]
                    l2 <- l2[-isele2]
                }else{
                    out[,isele2+1] <- out[,isele2+1]+out[,isele2]
                    out <- out[,-isele2]
                    l1 <- l1[-isele2]
                    l2 <- l2[-isele2]
                }
            }
            k <- k - 1
        }else{
            break;
        }
    }
    
    nams <- NULL
    l1 <- round(l1,1)
    l2 <- round(l2,1)
    if(l1[1] == l2[1]){
        nams <- c(nams, l1[1])
    }else{
        tmp <- paste("[",l1[1],",",l2[1],")", sep='')
        nams <- c(nams, tmp)
    }
    for(i in 2:k){
        if(l1[i]==l2[i]){
            nams <- c(nams, l2[i])
        }else{
            tmp <- paste("[",l1[i],"-",l2[i],")",sep='')
            nams <- c(nams, tmp)
        }
    }


    out <- as.data.frame(out)
    names(out) <- nams
    out
}

.degtest1 <- function(counts.x, counts.y,
                     normal.x, normal.y,
                     gene=NULL){
    stopifnot(all(counts.x>=0))
    stopifnot(all(counts.y>=0))
    stopifnot(all(normal.x>=0))
    stopifnot(all(normal.y>=0))
    xnam1 <- deparse(substitute(counts.x))
    ynam1 <- deparse(substitute(counts.y))
    xnam2 <- deparse(substitute(normal.x))
    ynam2 <- deparse(substitute(normal.y))
    nx <- ncol(counts.x)
    ny <- ncol(counts.y)
    if(ncol(normal.x) != nx)
        stop("profile numbers not match between ", xnam1, " and ", xnam2)
    if(ncol(normal.y) != ny)
        stop("profile numbers not match between ", ynam1, " and ", ynam2)
    ngene <- nrow(counts.x)
    if(nrow(counts.y) != ngene)
        stop("gene names not match between ", xnam1, " and ", ynam1)
    if(nrow(normal.x) != ngene)
        stop("gene names not match between ", xnam1, " and ", xnam2)
    if(nrow(normal.y) != ngene)
        stop("gene names not match between ", xnam1, " and ", ynam2)
        
    gnames <- row.names(counts.x)
    if(any(row.names(counts.y) != gnames))
        stop("gene names no match between ", xnam1, " and ", ynam1)
    if(any(row.names(normal.x) != gnames))
        stop("gene names no match between ", xnam1, " and ", xnam2)
    if(any(row.names(normal.y) != gnames))
        stop("gene names no match between ",xnam1, " and ", ynam2)

    rank.x <- apply(normal.x, 2, rank)
    rank.y <- apply(normal.y, 2, rank)

    if(is.null(gene)){
        pv1 <- rep(NA, ngene)
        pv2 <- pv3 <- pv1
        max1 <- max2 <- pv1
        zero1 <- zero2 <- pv1
        mean1c <- mean2c <- mean1 <- mean2 <- mean1r <- mean2r <- pv1
        sd1c <- sd2c <- sd1 <- sd2 <- sd1r <- sd2r <- pv1
        ng1 <- ng2 <- pv1
        
        for(x.sele in 1:ngene){
            gene <- gnames[x.sele]
            cat("\nGene name :=",gene, "Row number =", x.sele,"\n")
            out <- .degtest1fun(gene,xnam1,ynam1,xnam2,ynam2,
                            x1=as.numeric(counts.x[x.sele,]),
                            y1=as.numeric(counts.y[x.sele,]),
                            x2=as.numeric(normal.x[x.sele,]),
                            y2=as.numeric(normal.y[x.sele,]),
                            x3=as.numeric(rank.x[x.sele,]),
                            y3=as.numeric(rank.y[x.sele,]))
            zero1[x.sele] <- out$stat[1,2]
            zero2[x.sele] <- out$stat[2,2]
            max1[x.sele] <- out$stat[1,3]
            max2[x.sele] <- out$stat[2,3]

            mean1c[x.sele] <- out$test[1,1]
            mean2c[x.sele] <- out$test[1,3]
            sd1c[x.sele] <- out$test[1,2]
            sd2c[x.sele] <- out$test[1,4]
            pv1[x.sele] <- out$test[1,5]

            mean1[x.sele] <- out$test[2,1]
            mean2[x.sele] <- out$test[2,3]
            sd1[x.sele] <- out$test[2,2]
            sd2[x.sele] <- out$test[2,4]
            pv2[x.sele] <- out$test[2,5]

            mean1r[x.sele] <- out$test[3,1]
            mean2r[x.sele] <- out$test[3,3]
            sd1r[x.sele] <- out$test[3,2]
            sd2r[x.sele] <- out$test[3,4]
            pv3[x.sele] <- out$test[3,5]

            ng1[x.sele] <- round(out$fit.x$npar + 1)/3
            ng2[x.sele] <- round(out$fit.y$npar + 1)/3
        }
        
        out <- data.frame(Zeros.x=zero1, Zeros.y=zero2,
                          Max.x=max1, Max.y=max2,
                          MeanCount.x=mean1c,
                          SDCount.x=sd1c,
                          MeanCount.y=mean2c,
                          SDCount.x=sd2c,
                          p.Fisher=pv1,
                          Mean.x=mean1,
                          SD.x=sd1,
                          Mean.y=mean2,
                          SD.y=sd2,
                          p.Perm=pv2,
                          MeanRank.x=mean1r,
                          SDRank.x=sd1r,
                          MeanRank.y=mean2r,
                          SDRank.y=sd2r,
                          p.Perm.Rank=pv3,
                          ng.x=ng1, ng.y=ng2
                          )
        row.names(out) <- gnames
    }else{
        x.sele <- match(gene, gnames)[1]
        if(is.na(x.sele)){
            warning("gene ", gene, "not found")
            out <- NULL
        }else{
            out <- .degtest1fun(gene,xnam1,ynam1,xnam2,ynam2,
                            x1=as.numeric(counts.x[x.sele,]),
                            y1=as.numeric(counts.y[x.sele,]),
                            x2=as.numeric(normal.x[x.sele,]),
                            y2=as.numeric(normal.y[x.sele,]),
                            x3=as.numeric(rank.x[x.sele,]),
                            y3=as.numeric(rank.y[x.sele,]))
        }
    }
    invisible(out)
}

.degtest2fun <- function(gene, xnam,ynam,
                         x2,y2,x3,y3){
    
    stat <- data.frame(Size=c(length(x2), length(y2)),
                       Zeros=c(sum(x2==0),sum(y2==0)),
                       Max=c(max(x2),max(y2))
                       )
    row.names(stat) <- c(xnam, ynam)
    px0 <- mean(x2==0)
    py0 <- mean(y2==0)
    nx1 <- sum(x2>0)
    ny1 <- sum(y2>0)
    nx <- length(x2)
    ny <- length(y2)
    
    test <- NULL
    xytab <- NULL;
    if(px0 > 0 || py0 > 0){
        if(nx1 < 10 || ny1 < 10){
            tmp <- c(sum(x2 == 0), sum(y2 == 0))
            tmp1 <- c(sum(x2 > 0), sum(y2 > 0))
            if(any(tmp1 > 0)){
                xytab <- cbind(tmp, tmp1)
                xytab <- as.data.frame(xytab)
                row.names(xytab) <- c(xnam, ynam)
                names(xytab) <- c("Class1", "Class2")
            }else{
                xytab <- matrix(tmp, nrow=2)
                xytab <- as.data.frame(xytab)
                row.names(xytab) <- c(xnam, ynam)
                names(xytab) <- c("Class1")
            }
        }else{
            xy <- c(x2[x2>0],y2[y2>0])
            cut <- median(xy)
            tmp <- c(sum(x2 == 0), sum(y2 == 0))
            tmp1 <- c(sum(x2 <= cut), sum(y2 <= cut))
            if(any(tmp1 > 0))
                tmp <- cbind(tmp, tmp1)
            tmp2 <- c(sum(x2 > cut), sum(y2 > cut))
            if(any(tmp2 > 0))
                tmp <- cbind(tmp, tmp2)
            xytab <- tmp
            xytab <- as.data.frame(xytab)
            row.names(xytab) <- c(xnam, ynam)
            nclass <- ncol(xytab)
            rnam <- paste("Class", 1:nclass, sep='')
            names(xytab) <- rnam
        }
    }else{
        xy <- c(x2,y2)
        cut <- as.numeric(quantile(xy, prob=c(.33,.66)))
        tmp <- NULL
        tmp1 <- c(sum(x2 <= cut[1]), sum(y2 <= cut[1]))
        if(any(tmp1 > 0))
            tmp <- cbind(tmp, tmp1)
        tmp2 <- c(sum(x2 > cut[1] & x2 <= cut[2]),
                  sum(y2 > cut[1] & y2 <= cut[2]))
        if(any(tmp2 > 0))
            tmp <- cbind(tmp, tmp2)
        tmp3 <- c(sum(x2 > cut[2]), sum(y2 > cut[2]))
        if(any(tmp3 > 0))
            tmp <- cbind(tmp, tmp3)
        xytab <- tmp
        xytab <- as.data.frame(xytab)
        row.names(xytab) <- c(xnam, ynam)
        names(xytab) <- c("Class1", "Class2","Class3")
    }
    
    
    if(dim(xytab)[2] == 1){
        pv1 <- NA
        pv2 <- NA
        pv3 <- NA
    }else{
        pv1 <- fisher.test(xytab, alternative='two.sided')$p.value
        pv2 <- perm.test(x2, y2, iter=99999)$p.value
        pv3 <- perm.test(x3, y3, iter=99999)$p.value
    }

       
    mu.x1 <- NA
    mu.x2 <- mean(x2)
    mu.x3 <- mean(x3)
    sd.x1 <- NA
    sd.x2 <- sd(x2)
    sd.x3 <- sd(x3)

    mu.y1 <- NA
    mu.y2 <- mean(y2)
    mu.y3 <- mean(y3)
    sd.y1 <- NA
    sd.y2 <- sd(y2)
    sd.y3 <- sd(y3)

    test <- data.frame(
        mu1=round(c(mu.x1,mu.x2,mu.x3),3),
        sd1=round(c(sd.x1,sd.x2,sd.x3),3),
        mu2=round(c(mu.y1,mu.y2,mu.y3),3),
        sd2=round(c(sd.y1,sd.y2,sd.y3),3),
        pv=signif(c(pv1,pv2,pv3),6))
    row.names(test) <- c("Counts","normalized","Rank")
    names(test) <- c(paste(xnam,".",c("Mean","SD"),sep=''),
                     paste(ynam,".",c("Mean","SD"),sep=''),
                     "p.value")

    if(nx1 < 10){
        fx <- fit.mlnorm(x2, method='mop')
    }else if(nx1 <= 30){
        fx <- fit.mlnorm(x2, k=1, method='lognormal')
    }else{
        fx <- fit.mlnorm(x2, k=3, method='fnmm')
    }
        
    if(ny1 < 10){
        fy <- fit.mlnorm(y2, method='mop')
    }else if(ny1 <= 30){
        fy <- fit.mlnorm(y2, k=1, method='lognormal')
    }else{
        fy <- fit.mlnorm(y2, k=3, method='fnmm')
    }
    
    list(gene=gene, stat = stat, counts = xytab,
         test = test, fit.x=fx, fit.y=fy
         )
}

.stdrank <- function(x){
    tmp <- rank(x)
    rng <- range(tmp)
    (tmp - rng[1])/diff(rng)
}

.q3norm <- function(x,p){
    if(missing(p)) p <- 0.75
    if(p<=0 || p>=1)
        stop("invalid quantile level")
    q3 <- as.numeric(quantile(x,p))
    x/q3
}

.Q3Normalize <- function(x){
    tmp <- apply(x==min(x), 2, mean, na.rm=TRUE)
    p0max <- max(tmp)
    if(p0max > 0.99)
        stop("some profile(s) contain more than 99% zeros")
    pn <- 0.5 * p0max + 0.5
    for(i in 1:ncol(x)){
        x[,i] <- .q3norm(x[,i], pn)
    }
    invisible(x)
}

.bin2Counts <- function(x,y){
    xy <- c(x,y)
    ##if(all(x==0) || all(y==0))
    ##stop("all 'x' or 'y' are zeros")
    if(any(xy == 0)){
        xy1 <- xy[xy > 0]
        nx0 <- sum(x==0)
        ny0 <- sum(y==0)
        x1 <- x[x > 0]
        y1 <- y[y > 0]
    }else{
        x1 <- x
        y1 <- y
        xy1 <- xy
        nx0 <- 0
        ny0 <- 0
    }
    
    xp <- quantile(xy1, prob=c(1:4)/5)
    dxp <- diff(xp)
    sele <- dxp == 0
    if(any(sele)){
        if(all(sele)){
            xbreaks <- c(min(xy1)-1, max(xy1)+1)
        }else{
            xbreaks <- c(min(xy1)-1, xp[-which(sele)], max(xy1)+1)
        }
    }else{
        xbreaks <- c(min(xy1)-1, xp, max(xy1)+1)
    }
    xc <- hist(x1, breaks=xbreaks, plot=FALSE)
    yc <- hist(y1, breaks=xbreaks, plot=FALSE)
    tmp <- rbind(xc$counts,yc$counts)
    xbreaks[1] <- 0
    if(max(nx0, ny0) > 0){
        tmp <- cbind(c(nx0, ny0), tmp)
        xbreaks <- c(0,xbreaks)
    }
    .shrinkTbl(tmp, xbreaks)
}

.glmtest <- function(x,y,wx,wy){
    xy <- c(x,y)
    ## weight the event by sequencing depth for each profile. Make
    ## sure the weight will not be more than [wMax=2].  (We should not give
    ## too much weights due to the sequencing difference?)
    wMax <- 2
    wxy <- as.numeric(c(wx,wy))
    if(any(wxy <= 0)){
        wxy <- rep(0, length(wxy))
    }else if(length(table(wxy)) == 1){
        wxy <- rep(0, length(wxy))
    }else{
        wxy <- wxy / min(wxy)
        if(max(wxy) < wMax) wMax <- max(wxy)
        wxy <- (wxy - min(wxy))/diff(range(wxy)) * wMax
    }

    Y <- xy > 0
    W <- xy
    W[W==0] <- 1
    W[W>15] <- 15
    WTS <- W + wxy - 1
    
    Group <- as.factor(c(rep(0, length(x)), rep(1, length(y))))
    glmout <- glm(Y~Group, weights=round(W), family='binomial')

    pv <- rev(summary(glmout)$coef[2,])[1]
    pv
}

.testCounts <- function(x,y,wx,wy){
    ## works for count data only
    x <- round(x)
    y <- round(y)

    if(all(x==0) && all (y==0)){
        pv <- NA
    }else if(all(x==0)){
        if(max(y)<=1){
            p0 <- mean(y==1)
            pv <- binom.test(x=0,n=length(x),p=p0,
                             alternative='less')$p.value
            ##cat("\nDebugging....\n")
            ##print(table(x))
            ##print(table(y))
            ##print(pv)
        }else{
            ## we do an approximation by adding one event to the all
            ## zero sample
            x[1] <- 1
            pv <- .glmtest(x,y,wx,wy)
        }
    }else if(all(y==0)){
        if(max(x)<=1){
            p0 <- mean(x==1)
            pv <- binom.test(x=0,n=length(y),p=p0,
                             alternative='less')$p.value
        }else{
            ## we do an approximation by adding one event to the all
            ## zero sample
            y[1] <- 1
            pv <- .glmtest(x,y,wx,wy)
        }
    }else{
        pv <- .glmtest(x,y,wx,wy)
    }
    xy <- c(x,y)
    xt <- as.numeric(names(table(xy)))
    tmp <- matrix(0, nrow=2, ncol=length(xt))
    for(i in 1:length(xt)){
        tmp[1,i] <- sum(x==xt[i])
        tmp[2,i] <- sum(y==xt[i])
    }
    tmp <- as.data.frame(tmp)
    row.names(tmp) <- c("X","Y")
    names(tmp) <- paste("count=", xt, sep='')
    
    list(p.value=pv, table=tmp)
}

.shrinkTbl <- function(x,xbrks){
    if(max(xbrks)<10){
        xbrks <- round(xbrks, 4)
    }else{
        xbrks <- round(xbrks, 2)
    }
    stopifnot(is.matrix(x))
    stopifnot(min(dim(x))==2)
    if(nrow(x)==2) x <- t(x)
    nr <- nrow(x)
    while(nr > 2){
        if(any(x < 5)){
            isele <- which(x == min(x))[1]
            if(isele > nr) irow <- isele - nr
            else irow <- isele
            if(irow == 1){
                x[2,] <- x[2,] + x[irow,]
                xbrks <- xbrks[-(irow+1)]
            }else if(irow == nr){
                x[irow-1,] <- x[irow-1,] + x[irow,]
                xbrks <- xbrks[-(irow)]
            }else{
                min1 <- min(x[irow-1,])
                min2 <- min(x[irow+1,])
                if(min1 < min2){
                    x[irow-1,] <- x[irow-1,] + x[irow,]
                    xbrks <- xbrks[-(irow)]
                }else if(min1 > min2){
                    x[irow+1,] <- x[irow+1,] + x[irow,]
                    xbrks <- xbrks[-(irow+1)]
                }else{
                    min1 <- sum(x[irow-1,])
                    min2 <- sum(x[irow+1,])
                    if(min1 < min2){
                        x[irow-1,] <- x[irow-1,] + x[irow,]
                        xbrks <- xbrks[-(irow)]
                    }else{
                        x[irow+1,] <- x[irow+1,] + x[irow,]
                        xbrks <- xbrks[-(irow+1)]
                    }
                }
            }
            x <- x[-irow,]
            nr <- nr - 1
        }else{
            if(nr <= 5){
                break
            }else{
                isele <- which(x == min(x))[1]
                if(isele > nr) irow <- isele - nr
                else irow <- isele
                if(irow == 1){
                    x[2,] <- x[2,] + x[irow,]
                    xbrks <- xbrks[-(irow+1)]
                }else if(irow == nr){
                    x[irow-1,] <- x[irow-1,] + x[irow,]
                    xbrks <- xbrks[-(irow)]
                }else{
                    min1 <- min(x[irow-1,])
                    min2 <- min(x[irow+1,])
                    if(min1 < min2){
                        x[irow-1,] <- x[irow-1,] + x[irow,]
                        xbrks <- xbrks[-(irow)]
                    }else{
                        x[irow+1,] <- x[irow+1,] + x[irow,]
                        xbrks <- xbrks[-(irow+1)]
                    }
                }
                x <- x[-irow,]
                nr <- nr - 1
            }
        }
    }
    out <- as.data.frame(x)
    names(out) <- c("X","Y")
    k <- length(xbrks)
    xbrks[k] <- Inf
    
    if(xbrks[2]==0){
        ll <- xbrks[-c(1,k)];
        ul <- xbrks[-c(1,2)]
        nam <- paste("[",ll,",",ul,")",sep='')
        nam[1] <- paste("(",ll[1],",",ul[1],")",sep='')
        row.names(out) <- c("0",nam)
    }else{
        ll <- xbrks[-k];
        ul <- xbrks[-1]
        nam <- paste("[",ll,",",ul,")",sep='')
        row.names(out) <- nam
    }
    out
}

.testNorm2 <- function(x,y){
    tmp <- .bin2Counts(x,y)
    ##pv <- fisher.test(t(as.matrix(tmp)), simulate.p.value=TRUE)$p.value
    tmp2 <- as.matrix(tmp)
    if(max(dim(tmp2)) < 3){
        pv <- fisher.test(tmp2, simulate.p.value=FALSE)$p.value
    }else if(any(tmp2 < 5)){
        pv <- fisher.test(tmp2, simulate.p.value=TRUE)$p.value
    }else{
        pv <- chisq.test(tmp2)$p.value
    }
    list(p.value=pv, table=tmp)
}


##  We can develop a test to handle the data that the majority are
##  similar but a certain portion of data are larger than the 95%.
##  This is very similar to the testing of X>detection limit. This is
##  very similar to the Wilcoxon test and the detailed measures
##  doesn't matter (which is incorrect). The KS.test won't work as
##  well.  We consider to weight the data

.myGLMtest <- function(x,y, wx, wy){
    ## the data should have obviously different range.  If not,
    ## ignore.
    pv <- 1
    rx <- diff(range(x))
    ry <- diff(range(y))
    if(rx > 0 & ry > 0){
        if(rx/ry > 5 || ry/rx > 5){
            qx99 <- quantile(x, probs=0.95, names=FALSE)
            qy99 <- quantile(y, probs=0.95, names=FALSE)
            q99 <- min(qx99, qy99)
            ##x1 <- sum(x >= q99, na.rm=TRUE)
            ##x2 <- sum(y >= q99, na.rm=TRUE)
            ##n1 <- length(x)
            ##n2 <- length(y)
            ##pv1 <- prop.test(x=c(x1,x2),n=c(n1,n2))$p.value

            wMax <- 2
            wxy <- as.numeric(c(wx,wy))
            if(any(wxy <= 0)){
                wxy <- rep(0, length(wxy))
            }else if(length(table(wxy)) == 1){
                wxy <- rep(0, length(wxy))
            }else{
                wxy <- wxy / min(wxy)
                if(max(wxy) < wMax) wMax <- max(wxy)
                wxy <- (wxy - min(wxy))/diff(range(wxy)) * wMax
            }

            
            xy <- c(x,y)
            Y <- xy >= q99
            group <- factor(c(rep(0, length(x)),rep(1,length(y))))
            W <- ceiling(xy/q99)
            W[W>15] <- 15

            WTS <- round(wxy + W -1)
            glmout <- glm(Y~group, weights=round(W), family='binomial')
            pv <- rev(summary(glmout)$coef[2,])[1]
        }
        ##pv <- c(pv1,pv2)
    }
    pv
}

.testNorm <- function(x,y,wx,wy){
    wxy <- c(wx, wy)
    tmp <- .bin2Counts(x,y)
    tmp2 <- as.matrix(tmp)
    options(warn=-1)
    xy <- c(x,y)
    if(mean(xy==0) > 0.85){
        pv <- .testCounts(x,y,wx,wy)$p.value
    }else{
        pv <- ks.test(x,y)$p.value
    }
    options(warn=0, error=NULL)
    list(p.value=pv, table=tmp)
}

.geneDEG <- function(x1,y1,x2,y2,x3,y3,gene){
    out <- NULL
    out$gene <- gene
    ## Fisher's test based on counts
    if(is.null(x2)||is.null(y2)){
        out$xCstat <- NA
        out$yCstat <- NA
        out$Ctable <- NA
        out$Cpv <- NA
    }else{
        out$xCstat <- c(min(x2),max(x2))
        out$yCstat <- c(min(y2),max(y2))
        tmp <- .testCounts(x2,y2)
        out$Ctable <- tmp$Tab
        out$Cpv <- tmp$p.value
    }
    
    ## Permutation test based on normalized data If the total percent
    ## of zeros is large (say larger than 90%), we test using Fisher's
    ## test.  Otherwise, permutation test.
    out$xstat <- c(mean(x2),sd(x2))
    out$ystat <- c(mean(y2),sd(y2))
    tmp <- .testNorm(x1,y1)
    out$perm.pv <- tmp$p.value
    
    ## Permutation test based on rank based data.
    out$xRstat <- c(mean(x3),sd(x3))
    out$yRstat <- c(mean(y3),sd(y3))
    tmp <- .testNorm(x1,y1)
    out$Rperm.pv <- tmp$p.value

    ## choose the best result/method
    if(max(c(x2, y2)) < 3){
        out$p.value <- out$Cpv
        out$method <- "Fisher's"
    }else if(max(c(x2, y2)) < 20){
        if(min(c(out$xRstat[1],out$yRstat[1])) > 0.90){
            out$p.value <- out$perm.pv
            out$method <- "Permutation"
        }else{
            tmp <- c(out$Cpv, out$perm.pv, out$Rper.pv)
            meth <- c("Fisher's","Permutation","Permutation (rank)")
            sele <- is.na(tmp)
            if(any(sele)){
                if(sum(sele)==3){
                    out$p.value <- NA
                    out$method <- NA
                }else{
                    tmp <- tmp[!sele]
                    meth <- meth[!sele]
                    pv0 <- min(tmp)
                    isele <- which(tmp==pv0)[1]
                    out$p.value <- pv0
                    out$method <- meth[isele]
                }
            }else{
                pv0 <- min(tmp)
                isele <- which(tmp==pv0)[1]
                out$p.value <- pv0
                out$method <- meth[isele]
            }
        }
    }else{
        if(min(c(out$xRstat[1],out$yRstat[1])) > 0.90){
            out$p.value <- out$perm.pv
            out$method <- "Permutation"
        }else{
            tmp <- c(out$perm.pv, out$Rper.pv)
            meth <- c("Permutation","Permutation (rank)")
            sele <- is.na(tmp)
            if(any(sele)){
                if(sum(sele)==2){
                    out$p.value <- NA
                    out$method <- NA
                }else{
                    tmp <- tmp[!sele]
                    meth <- meth[!sele]
                    pv0 <- min(tmp)
                    isele <- which(tmp==pv0)[1]
                    out$p.value <- pv0
                    out$method <- meth[isele]
                }
            }else{
                pv0 <- min(tmp)
                isele <- which(tmp==pv0)[1]
                out$p.value <- pv0
                out$method <- meth[isele]
            }
        }
    }
    

    structure(out, class='DEG')
}

plot.DEG <- function(x,col=c(1,2),
                     xlim=NULL,ylim=NULL,
                     boxplot=FALSE,
                     legend=FALSE,...){
    if(is.null(x$gene)){
        cat("\n view data in obj$table\n")
    }else{
        if(is.null(col)){
            col1 <- 1;
            col2 <- 1
        }else{
            if(length(col)==1){
                col1 <- col;
                col2 <- col
            }else{
                col1 <- col[1];
                col2 <- col[2]
            }
        }
        xy <- c(x$fx$data, x$fy$data)
        if(is.null(xlim)){
            xmin <- min(xy)*.95; xmax <- max(xy)+1.1
            xlim <- c(xmin, xmax)
        }
        if(is.null(ylim)){
            ymax <- max(x$fx$y, x$fy$y) * 1.1
            ylim <- c(0,ymax)
        }
        ## initiate the plot
        plot(0,0,xlim=xlim, ylim=ylim, type='n',
             xlab="gene expression", ylab='density',
             main=x$gene)
        ## draw the first histogram manually
        if(is.null(x$fx$breaks)){
            tmp <- hist(x$fx$data, plot=FALSE,...)
            xbrks <- tmp$breaks
            xcount <- tmp$count
        }else{
            xbrks <- x$fx$breaks
            tmp <- hist(x$fx$data, breaks=xbrks, plot=FALSE,...)
            xcount <- tmp$count
        }
        ## relative frequencies
        xsum <- sum(diff(xbrks) * xcount)
        xcount <- xcount/xsum
        for(i in 1:length(xcount)){
            x0 <- c(xbrks[i], xbrks[i+1]); x0 <- c(x0, rev(x0))
            y0 <- c(0,0,xcount[i],xcount[i])
            polygon(x0,y0,border=col1,...)
        }
        x0 <- seq(xlim[1],xlim[2], length=400)
        f0 <- dmlnorm(x0, x$fx$p, x$fx$meanlog, x$fx$sdlog)
        lines(f0~x0,col=col1,...)

        ## draw another histogram manually
        if(is.null(x$fy$breaks)){
            tmp <- hist(x$fy$data, plot=FALSE,...)
            xbrks <- tmp$breaks
            xcount <- tmp$count
        }else{
            xbrks <- x$fy$breaks
            tmp <- hist(x$fy$data, breaks=xbrks, plot=FALSE,...)
            xcount <- tmp$count
        }
        ## relative frequencies
        xsum <- sum(diff(xbrks) * xcount)
        xcount <- xcount/xsum
        for(i in 1:length(xcount)){
            x0 <- c(xbrks[i], xbrks[i+1]); x0 <- c(x0, rev(x0))
            y0 <- c(0,0,xcount[i],xcount[i])
            polygon(x0,y0,border=col2,...)
        }
        x0 <- seq(xlim[1],xlim[2], length=400)
        f0 <- dmlnorm(x0, x$fy$p, x$fy$meanlog, x$fy$sdlog)
        lines(f0~x0,col=col2,...)
        
        if(legend){
            if(boxplot){
                legend("topleft", legend=c("X","Y"),
                       lty=c(1,1), col=c(col1,col2))
            }else{
                legend("topright", legend=c("X","Y"),
                       lty=c(1,1), col=c(col1,col2))
            }
        }

        if(boxplot){
            op2 <- par(fig = c(0.6,1,.6,1), new = TRUE)
            X <- log(1+x$fx$data)
            Y <- log(1+x$fy$data)
            tmp <- c(col1, col2)
            if(tmp[1]==1) tmp[1] <- 0
            if(tmp[2]==1) tmp[2] <- 0
            boxplot(X, Y, col=tmp,names=c("x","y"))
            par(op2)
        }
    }
}
        
print.DEG <- function(x,...){
    if(is.null(x$gene)){
        cat("\ncomplete results in obj$table:\n")
        pv <- x$table$p.value
        n <- min(nrow(x$table), 20)
        print(x$table[order(pv)[1:n],])
    }else{
        cat("\n\nGene name: ", x$gene,
            "\n  p.value=", x$p.value,
            " [",x$method,"]\n")
        print(x$table)
        
        tmp <- data.frame(
            Counts = paste("[",c(x$min.x, x$min.y),",",
                           c(x$max.x, x$max.y), "]", sep=''),
            Noramlized = paste(round(c(x$mean.x, x$mean.y),3)," (",
                               round(c(x$sd.x, x$sd.y),3), ")", sep=''),
            Rank = paste(round(c(x$mean.rank.x, x$mean.rank.y),3)," (",
                         round(c(x$sd.rank.x, x$sd.rank.y),3), ")", sep='')
        )
        
        row.names(tmp) <- c("X","Y")
        print(tmp)
        
        tmp <- .chkoutlier(x$fx$data)
        if(!is.null(tmp)){
            cat("\nWaning: potential outliers in 'X':", tmp,"\n")
        }
        print(x$fx)

        tmp <- .chkoutlier(x$fy$data)
        if(!is.null(tmp)){
            cat("\nWaning: potential outliers in 'Y':", tmp,"\n")
        }
        print(x$fy)
    }
}

.fitgexpression <- function(x){
    ## if there are enought positive measures, we fit a mixture of 'k'
    ## lognormal components.  If the number of positive measures is
    ## less than 30, fit using MoP (method of percentiles) when the
    ## number of distinct values is two or less; otherwise, fit a
    ## lognormal distribution ('k=1').
    
    np <- sum(x > 0)
    
    if(np > 20){
        fx <- fit.mlnorm(x, k=3, method='fnmm')
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


## this is an old version.  Saved on July 29, 2018.
.deg.test.bak1 <- function(normal.x, normal.y,counts.x, counts.y, gene){
    ## If gene is NULL, perform analysis to all genes.  We can filter
    ## the dataset to perform the test to a small set of genes. If
    ## normalized data are not available, we normalize the counts data
    ## using a 75 percentile (or if the proportion of zeros exceeds
    ## 75%, we use (1+p0)/2.  If normalized data are available, and
    ## counts data are not available, we ingore the analysis using
    ## counts data. If both are available, do a complete test.

    if(missing(counts.x)||missing(counts.y)){
        if(missing(normal.x)||missing(normal.y))
            stop("both normalized and counts data are incomplete")
        ## otherwise if both normalized datasets are available...
        .CheckData(normal.x, normal.y, type=1)
        rank.x <- apply(normal.x, 2, rank)
        rank.y <- apply(normal.y, 2, rank)
        counts.x <- NULL
        counts.y <- NULL
        cat("\nCounts (X):\tnot available")
        cat("\nCounts (Y):\tnot available")
        tmp <- deparse(substitute(normal.x))
        cat("\nNormalized (X):\t", tmp)
        tmp <- deparse(substitute(normal.y))
        cat("\nNormalized (Y):\t", tmp)
    }else{
        tmp <- deparse(substitute(counts.x))
        cat("\nCounts (X):\t", tmp)
        tmp <- deparse(substitute(counts.y))
        cat("\nCounts (Y):\t", tmp)
        .CheckData(counts.x, counts.y, type=1)
        ## normalize counts data using Q3-like normalization
        if(missing(normal.x)){# normalize the data
            normal.x <- .Q3Normalize(counts.x)
            cat("\nNormalized (X):\t upper quartile normalized")
        }else{ # check row names and dim
            .CheckData(normal.x, counts.x, type=2)
            tmp <- deparse(substitute(counts.x))
            cat("\nNormalized (X):\t", tmp)
        }

        ## normalize counts data using Q3-like normalization
        if(missing(normal.y)){# normalize the data
            normal.y <- .Q3Normalize(counts.y)
            cat("\nNormalized (Y):\t upper quartile normalized")
        }else{ # check row names and dim
            .CheckData(normal.y, counts.y, type=2)
            tmp <- deparse(substitute(normal.y))
            cat("\nNormalized (Y):\t", tmp)
        }
        rank.x <- apply(counts.x, 2, .stdrank)
        rank.y <- apply(counts.y, 2, .stdrank)
    }

    gnames <- row.names(normal.x)

    if(missing(gene)){
        genes <- row.names(normal.x)
        Batch <- TRUE
    }else{
        if(length(gene) > 1){
            genes <- gene
            Batch <- TRUE
        }else{
            Batch <- FALSE
        }
    }
    
    out <- NULL
    if(Batch){
        res <- matrix(0, nrow=length(genes),ncol=15)
        mypv <- NULL;
        method <- NULL
        i <- 0; isele <- NULL
        cat("\n...Processing: \n")
        for(gene in genes){
            i <- i + 1; j <- 1
            irow <- match(gene, gnames, nomatch=0)[1]
            cat("\nrow=", irow, ", gene=", gene, "...")
            if(irow > 0){
                isele <- c(isele, irow)
                tmp <- .geneDEG(x1=as.numeric(normal.x[irow,]),
                                y1=as.numeric(normal.y[irow,]),
                                x2=as.numeric(counts.x[irow,]),
                                y2=as.numeric(counts.y[irow,]),
                                x3=as.numeric(rank.x[irow,]),
                                y3=as.numeric(rank.y[irow,]),
                                gene=gene)
                mypv <- c(mypv, tmp$p.value)
                method <- c(method, tmp$method)
                res[i,j] <- tmp$xCstat[1]; j <- j + 1;
                res[i,j] <- tmp$xCstat[2]; j <- j + 1;
                res[i,j] <- tmp$yCstat[1]; j <- j + 1;
                res[i,j] <- tmp$yCstat[2]; j <- j + 1;
                res[i,j] <- tmp$Cpv; j <- j + 1;
                res[i,j] <- tmp$xstat[1]; j <- j + 1;
                res[i,j] <- tmp$xstat[2]; j <- j + 1;
                res[i,j] <- tmp$ystat[1]; j <- j + 1;
                res[i,j] <- tmp$ystat[2]; j <- j + 1;
                res[i,j] <- tmp$perm.pv; j <- j + 1;
                res[i,j] <- tmp$xRstat[1]; j <- j + 1;
                res[i,j] <- tmp$xRstat[2]; j <- j + 1;
                res[i,j] <- tmp$yRstat[1]; j <- j + 1;
                res[i,j] <- tmp$yRstat[2]; j <- j + 1;
                res[i,j] <- tmp$Rperm.pv; j <- j + 1;
            }
        }
        cat("\n\nDONE (results in obj$output).\n")
        output <- data.frame(res,p.adj=p.adjust(mypv),method=method)
        row.names(output) <- gnames[isele]
        nams <- c("Min.Count.x","Max.Count.x",
                  "Min.Count.y","Max.Count.y",
                  "pv1",
                  "Mean.x","SD.x","Mean.y","SD.y",
                  "pv2",
                  "Mean.Rank.x","SD.Rank.x",
                  "Mean.Rank.y","SD.Rank.y",
                  "pv3",
                  'p.value', 'method')
        names(output) <- nams
        out$output <- output
        out$gene <- NULL
        out <- structure(out, class='DEG')

    }else{
        irow <- match(gene, gnames, nomatch=0)[1]
        if(irow==0){
            cat("\nGene being checked: NOT FOUND!\n")
        }else{
            cat("\nGene being checked:\t", gene,"\n")
            out <- .geneDEG(x1=as.numeric(normal.x[irow,]),
                            y1=as.numeric(normal.y[irow,]),
                            x2=as.numeric(counts.x[irow,]),
                            y2=as.numeric(counts.y[irow,]),
                            x3=as.numeric(rank.x[irow,]),
                            y3=as.numeric(rank.y[irow,]),
                            gene=gene)
            
            if(is.null(counts.x) || is.null(counts.y)){
                out$fx <- .fitgexpression(normal.x[irow,])
                out$fy <- .fitgexpression(normal.y[irow,])
            }else{
                if(max(c(counts.x[irow,], counts.y[irow,])) <= 10){
                    out$fx <- .fitgexpression(counts.x[irow,])
                    out$fy <- .fitgexpression(counts.y[irow,])
                }else{
                    out$fx <- .fitgexpression(normal.x[irow,])
                    out$fy <- .fitgexpression(normal.y[irow,])
                }
            }
        }
    }
    invisible(out)
}

.mlnormtest <- function(gene,x1,y1,x2,y2,x3,y3){
    if(is.null(x2)||is.null(y2)){
        fx <- .fitgexpression(x1)
        fy <- .fitgexpression(y1)
        method <- "KS-mlnorm (normalized)"
    }else{
        if(max(c(x2, y2)) <= 10){
            fx <- .fitgexpression(x2)
            fy <- .fitgexpression(y2)
            method <- "KS-mlnorm (counts)"
        }else{
        fx <- .fitgexpression(x1)
            fy <- .fitgexpression(y1)
            method <- "KS-mlnorm (normalized)"
        }
    }
    fx$x.name <- paste("X: ",gene,sep='')
    fy$x.name <- paste("Y: ",gene,sep='')
    ## this test cannot be applied to data that has very few distinct
    ## values (say all zeros).
    if(sum(x1>0)>5 & sum(y1 > 5)){
        x0 <- seq(0, max(c(fx$data,fy$data)), length=401)
        Fx <- pmlnorm(x0, fx$p, fx$meanlog, fx$sdlog)
        Fy <- pmlnorm(x0, fy$p, fy$meanlog, fy$sdlog)
        D <- max(abs(Fy-Fx))
        pv <- .Fortran(.F_pks2, p=as.double(D),
                       as.integer(length(fx$data)),
                       as.integer(length(fy$data)))$p
    }else{
        pv <- 1
    }
    list(fx=fx, fy=fy, p.value=pv, method = method)
}


.geneDEG2 <- function(x1,y1,x2,y2,x3,y3,wx,wy,gene){
    out <- NULL
    out$gene <- gene
    ## Fisher's test based on counts
    if(is.null(x2)||is.null(y2)){
        out$min.x <- NA
        out$max.x <- NA
        out$min.y <- NA
        out$max.y <- NA
        out$fx <- .fitgexpression(x1)
        out$fy <- .fitgexpression(y1)
    }else{
        out$min.x <- min(x2)
        out$max.x <- max(x2)
        out$min.y <- min(y2)
        out$max.y <- max(y2)
        if(max(c(x2, y2)) <= 10){
            out$fx <- .fitgexpression(x2)
            out$fy <- .fitgexpression(y2)
        }else{
            out$fx <- .fitgexpression(x1)
            out$fy <- .fitgexpression(y1)
        }
    }
    out$fx$x.name <- paste("X: ",gene,sep='')
    out$fy$x.name <- paste("Y: ",gene,sep='')
    out$mean.x <- mean(x1)
    out$sd.x <- sd(x1)
    out$mean.y <- mean(y1)
    out$sd.y <- sd(y1)
    
    out$mean.rank.x <- mean(x3)
    out$sd.rank.x <- sd(x3)
    out$mean.rank.y <- mean(y3)
    out$sd.rank.y <- sd(y3)
    ## if rank data are not available, we test using chisq.test or
    ## fisher's exact test
    if(is.null(x2)||is.null(y2)){
        res1 <- .testNorm(x1,y1,wx,wy) 
        res2 <- .testNorm(x3,y3,wx,wy)
        pv1 <- res1$p.value
        pv2 <- res2$p.value

        if(pv1 > pv2){
            out$p.value <- pv1;
            out$table <- res1$table
            out$method <- "normalized"
        }else{
            out$p.value <- pv2;
            out$table <- res2$table
            out$method <- "rank-based"
        }
    }else{
        if(max(c(x2,y2)) <= 10){
            res <- .testCounts(x2,y2,wx,wy);
            out$p.value <- res$p.value
            out$table <- t(res$table)
            out$method <- "counts"
        }else{
            res1 <- .testNorm(x1,y1,wx,wy) 
            res2 <- .testNorm(x3,y3,wx,wy)
            pv1 <- res1$p.value
            pv2 <- res2$p.value

            out$table <- res1$table
            if(pv1 > pv2){
                out$p.value <- pv1;
                out$method <- "normalized"
            }else{
                out$p.value <- pv2;
                out$method <- "rank-based"
            }
        }
    }

    fx <- out$fx
    fy <- out$fy
    ##cvr.x <- fx$SD/fx$Mean
    ##cvr.y <- fy$SD/fy$Mean
    ##rrng <- diff(range(fy$data))/diff(range(fx$data))
    ##cat("\n Ratios = [", cvr.x, ",", cvr.y,
    ##"]\nPerc of zeros = [",
    ##round(mean(fx$data==0)*100,2), "%,",
    ##round(mean(fy$data==0)*100,2),
    ##"%]\nRatio of ranges =",
    ##rrng,"\n")
    ##if(cvr.x > 10 || cvr.y > 10){
    x0 <- seq(0, max(c(fx$data,fy$data)), length=10001)
    Fx <- pmlnorm(x0, fx$p, fx$meanlog, fx$sdlog)
    Fy <- pmlnorm(x0, fy$p, fy$meanlog, fy$sdlog)
    D <- max(abs(Fy-Fx))
    pv <- .Fortran(.F_pks2, p=as.double(D),
                   as.integer(length(fx$data)),
                   as.integer(length(fy$data)))$p
    ##cat("KS-test by lognorm fitting. p-value=", pv,
    ##", D=",D,"\n")
    ##}
    if(pv < out$p.value){
        out$p.value <- pv
        out$method <- "KS-mlnorm"
    }

    pv <- .myGLMtest(x1,y1,wx,wy)
    if(pv < out$p.value){
        out$p.value <- pv
        out$method <- "Logit"
    }


    
    
    structure(out, class='DEG')
}

.chkoutlier <- function(x){
    ol <- NULL
    if(sum(x>0) > 15){
        x <- x + runif(length(x))*.001
        x0 <- x
        x <- sort(x, decreasing = TRUE)
        cv0 <- sd(x)/mean(x)
        cv1 <- sd(x[-1])/mean(x[-1])
        cv2 <- sd(x[-c(1:2)])/mean(x[-c(1:2)])
        cv3 <- sd(x[-c(1:3)])/mean(x[-c(1:3)])
        if(cv0/cv1 > 2){
            isele <- which(x0==x[1])
            ol <- c(ol,isele)
        }
        if(cv1/cv2 > 2.5){
            isele <- which(x0==x[2])
            ol <- c(ol,isele)
        }
        if(cv2/cv3 > 2.5){
            isele <- which(x0==x[3])
            ol <- c(ol,isele)
        }
    }
    ol
}


.geneDEG3 <- function(x1,y1,x2,y2,x3,y3,wx,wy,gene){
    out <- NULL
    out$gene <- gene
    ## Fisher's test based on counts
    if(is.null(x2)||is.null(y2)){
        out$min.x <- NA
        out$max.x <- NA
        out$min.y <- NA
        out$max.y <- NA
    }else{
        out$min.x <- min(x2)
        out$max.x <- max(x2)
        out$min.y <- min(y2)
        out$max.y <- max(y2)
    }
    out$mean.x <- mean(x1)
    out$sd.x <- sd(x1)
    out$mean.y <- mean(y1)
    out$sd.y <- sd(y1)
    
    out$mean.rank.x <- mean(x3)
    out$sd.rank.x <- sd(x3)
    out$mean.rank.y <- mean(y3)
    out$sd.rank.y <- sd(y3)
    ## check counts data only if max(x,y) <= 10.  Otherwise don't.
    ## Also, if counts data are available and max(x,y) <= 10, check
    ## both normalized and rank data and compare results. If max(x,y)
    ## > 10, don't use counts data directly.

    ## test based on counts
    tbl1 <- .tabulateXY(x2,y2)
    out$table.counts <- tbl1
    ##print(tbl1) # <-------  This need to be removed later ########## DELETE

    res <- .testCounts2(x2,y2);
    
    ##cat("\nMethod=", res$method, ", p-value=", res$p.value,"\n") # <--- DELETE--
    out$p.value <- res$p.value
    out$method <- res$method
    out$table <- tbl1
    
    if(max(c(x2,y2)) > 10){
        res1 <- .testNorm(x1,y1,wx,wy) 
        out$table.normalized <- res1$table
        res2 <- .testNorm(x3,y3,wx,wy)
        out$table.rank <- res2$table
        
        pv1 <- res1$p.value
        pv2 <- res2$p.value
        if(pv1 > pv2){
            if(pv1 < out$p.value){
                out$p.value <- pv1;
                out$method <- "normalized"
                out$table <- res1$table
            }
        }else{
            if(pv2 < out$p.value){
                out$p.value <- pv2;
                out$method <- "rank-based"
                out$table <- res2$table
            }
        }
    }
    
    pv <- .myGLMtest2(x1,y1)
    if(pv < out$p.value){
        out$p.value <- pv
        out$method <- "glm.test(binned)"
    }
    structure(out, class='DEG')
}

.tabulateXY <- function(x,y){
    xy <- c(x,y)
    xt <- as.numeric(names(table(xy)))
    k <- length(xt)
    if(k > 10){
        k0 <- 10
        nams <- c(paste("x=",xt[1:9],sep=''),
                  paste("x=", xt[10],"+",sep=''))
    }else{
        k0 <- k
        nams <- paste("x=",xt,sep='')
    }
    X <- rep(0, k0)
    Y <- rep(0,k0)
    for(i in 1:(k0-1)){
        X[i] <- sum(x==xt[i])
        Y[i] <- sum(y==xt[i])
    }
    X[k0] <- sum(x>=xt[k0])
    Y[k0] <- sum(y>=xt[k0])
    tmp <- data.frame(X,Y)
    row.names(tmp) <- nams
    tmp
}

.glmtest2 <- function(x,y){
    xy <- c(x, y)
    Y <- xy > 0
    W <- xy
    W[W<1] <- 1
    W[W>10] <- 10
    ##    wts <- round(log2(W)+1)
    Group <- as.factor(c(rep(0, length(x)), rep(1, length(y))))
    glmout <- glm(Y~Group, weights=W, family='binomial')

    if(!glmout$converged){
        glmout <- glm(Y~Group, weights=round(W), family='binomial')
        if(glmout$converged){
            pv <- rev(summary(glmout)$coef[2,])[1]
        }else{
            pv <- 1
        }
    }else{
        pv <- rev(summary(glmout)$coef[2,])[1]
    }
    pv
}

.testCounts2 <- function(x,y){
    if(is.null(x) || is.null(y)){
        out <- list(p.value=1, method="none")
    }else{ 
        ## works for count data only
        x <- round(x); y <- round(y)

        if(max(x) > 100 && max(y) > 100){
            out <- list(p.value=1, method="none")
        }else{

            if(all(x==0) && all(y==0)){
                out <- list(p.value=1, method="none")
            }else if(max(c(x,y))==1){
                if(all(x==0)){
                    p0 <- mean(y==1)
                    pv <- binom.test(x=0,n=length(x),p=p0,
                                     alternative='less')$p.value
                }else if(all(y==0)){
                    p0 <- mean(x==1)
                    pv <- binom.test(x=0,n=length(y),p=p0,
                                     alternative='less')$p.value
                }else{
                    p0 <- mean(y==1)
                    pv <- binom.test(x=sum(x>0),n=length(x),p=p0,
                                     alternative='less')$p.value
                }
                out <- list(p.value=pv, method="binom.test(counts)")
            }else{
                if(all(x==0)){
                    z <- y
                    tmp <- sum(z)
                    nz <- length(z)
                    if(tmp > nz) tmp <- nz
                    pv <- binom.test(x=0,n=length(x),p=tmp/nz,
                                     alternative='less')$p.value
                    out <- list(p.value=pv, method="binom.test(counts)")
                }else if(all(y==0)){
                    z <- x
                    tmp <- sum(z)
                    nz <- length(z)
                    if(tmp > nz) tmp <- nz
                    pv <- binom.test(x=0,n=length(y),p=tmp/nz,
                                     alternative='less')$p.value
                    out <- list(p.value=pv, method="binom.test(counts)")
                }else if(max(x)==1 && sum(x) <= 3){
                    z <- y
                    tmp <- sum(z)
                    nz <- length(z)
                    if(tmp > nz) tmp <- nz
                    pv <- binom.test(x=sum(x),n=length(x),p=tmp/nz,
                                     alternative='less')$p.value
                    out <- list(p.value=pv, method="binom.test(counts)")
                }else if(max(y)==1 && sum(y) <= 3){
                    z <- x
                    tmp <- sum(z)
                    nz <- length(z)
                    if(tmp > nz) tmp <- nz
                    pv <- binom.test(x=sum(y),n=length(y),p=tmp/nz,
                                     alternative='less')$p.value
                    out <- list(p.value=pv, method="binom.test(counts)")
                }else{
                    x1 <- x; y1 <- y
                    if(all(x1>0)){
                        isele <- which(x1==min(x1))
                        x1[isele[1]] <- 0
                    }
                    if(all(y1>0)){
                        isele <- which(y1==min(y1))
                        y1[isele[1]] <- 0
                    }
                    pv <- .glmtest2(x1,y1)
                    out <- list(p.value=pv, method="glm.test(counts)")
                }
            }
        }
    }
    out
}


.myGLMtest2 <- function(x,y){
    ## the data should have obviously different range.  If not,
    ## ignore.  Also, if the data has small number of distinct values,
    ## ignore!!!
    pv <- 1
    if(sum(x>0) > 10 && sum(y>0) > 10){
        rx <- diff(range(x))
        ry <- diff(range(y))
        if(rx > 0 & ry > 0){
            if(rx/ry > 5 || ry/rx > 5){
                qx99 <- quantile(x, probs=0.95, names=FALSE)
                qy99 <- quantile(y, probs=0.95, names=FALSE)
                q99 <- min(qx99, qy99)
                xy <- c(x,y)
                Y <- xy >= q99
                group <- factor(c(rep(0, length(x)),rep(1,length(y))))
                W <- ceiling(xy/q99)
                W[W>15] <- 15
                glmout <- glm(Y~group, weights=round(W), family='binomial')
                pv <- rev(summary(glmout)$coef[2,])[1]
            }
        }
    }
    pv
}
