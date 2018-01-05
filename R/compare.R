## to develop a function to compare means, proportions when there are
## two or more samples exist.

print.comp <- function(x,...){
    cat("\nComparisons between x[=1] ",x$x.name,
        "and y[=2] ",x$y.name)
    cat("\nNon-parametric method recommended? ",
        ifelse(x$np.method,"yes","no"))
    cat("\nSummary Statistics:\n")
    if(is.matrix(x$stat)){
        cat("Freq Dist:\n\t")
        print(x$stat)
        cat("Perc (cell) Dist:\n\t")
        print(prop.table(x$stat))
        cat("Perc (row) Dist:\n\t")
        print(prop.table(x$stat,1))
        cat("Perc (col) Dist:\n\t")
        print(prop.table(x$stat,2))
    }else{
        print(x$stat,...)
    }
    
    cat("\nTest Results:\n")
    print(x$test,...)
    invisible(x)
}

compare <- function(y, x, var.name,
                    alternative = c("two.sided", "less", "greater"),
                    paired = FALSE, 
                    conf.level = 0.95, ...)
    UseMethod("compare")

compare.matrix <- function(y, x, var.name,
                           alternative = c("two.sided", "less", "greater"),
                           paired = FALSE, 
                           conf.level = 0.95, ...)
{
    warning("The matrix will be analyzed by columns...")
    y <- as.data.frame(y)
    compare.data.frame(y=y, x=y, var.name=var.name,
                       alternative = alternative,
                       paired = paired, conf.level = conf.level, ...)
}

compare.data.frame <- function(y, x, var.name,
                               alternative = c("two.sided", "less", "greater"),
                               paired = FALSE, 
                               conf.level = 0.95, ...)
{
    alpha <- 1-conf.level
    nams <- names(y)
    ## define output variables.  save to two tables for num and cat
    ## vars
    vname1 <- NULL
    stats1 <- NULL
    pv1 <- NULL
    vname2 <- NULL
    stats2 <- NULL
    pv2 <- NULL
    m <- ncol(y)
    if(is.numeric(x)) stop("'x' must be a factor or categorical variable")
    x <- as.character(x)
    lvnames <- levels(as.factor(x))
    k <- length(lvnames)
    if(k<2) stop("Two or more levels are needed for comparisons...")
    for(i in 1:m){
        xi <- y[,i]
        if(is.numeric(xi)){
            vname1 <- c(vname1, nams[i])
            out <- .test.num(xi,x,k,alpha)
            stats1 <- rbind(stats1,out$stats)
            pv1 <- c(pv1,out$pv)
        }else{
            vname2 <- c(vname2, nams[i])
            out <- .test.cat(xi,x)
            vname2 <- c(vname2, out$levels)
            tmp <- rep("",k)
            stats2 <- rbind(stats2,tmp)
            stats2 <- rbind(stats2,out$stats)
            pv2 <- c(pv2,out$pv)
            tmp <- rep(NA, nrow(out$stats))
            pv2 <- c(pv2,tmp)
        }
    }
    if(length(vname1)>0){
        stats1 <- as.data.frame(stats1)
        names(stats1) <- lvnames
        pv1 <- round(pv1,4)
        out1 <- data.frame(Var=vname1,
                           stats1,
                           p.value=pv1)
        row.names(out1) <- NULL
    }else{
        out1 <- NULL
    }
    
    if(length(vname2)>0){
        stats2 <- as.data.frame(stats2)
        names(stats2) <- lvnames
        sele <- is.na(pv2)
        pv2 <- as.character(round(pv2,4))
        pv2[sele] <- ""
        out2 <- data.frame(Var=vname2,
                           stats2,
                           p.value=pv2)
        row.names(out2) <- NULL
    }else{
        out2 <- NULL
    }

    list(Numerical=out1,Categorical=out2)
}

.test.num <- function(y,x,k,alpha){
    mu <- tapply(y, x, mean, na.rm=TRUE)
    s <- tapply(y, x, sd, na.rm=TRUE)
    tmp2 <- paste(round(mu,2), " (",round(s,2),")",sep="")
    stats <- as.character(tmp2)
    pv <- .pv.num(y,x,k,alpha)
    list(stats=stats,pv=pv)
}

.test.cat <- function(y,x){
    y <- as.character(y)
    x <- as.character(x)
    lvnames.y <- levels(as.factor(y))
    lvnames.x <- levels(as.factor(x))
    tbl <- table(y, x)
    pct <- round(prop.table(tbl,2)*100,2)
    out <- NULL
    for(i in 1:nrow(tbl)){
        tmp2 <- paste(tbl[i,], " (",round(pct[i,],2),"%)",sep="")
        out <- rbind(out,tmp2)
    }
    stats <- out
    row.names(stats) <- NULL
    stats <- as.data.frame(stats)
    pv <- .pv.cat(tbl)
    list(levels=lvnames.y,stats=stats,pv=pv)
}

.pv.num <- function(y,x,k,alpha){
    ## x must be categorical.  We won't consider paired cases.  So we
    ## can leave the previous version of "compare(num, num,...)"
    ## function to take care the paired sample cases.
    x <- as.character(x)
    if(k>1){
        ## large samples (CLT applicable?)
        ns <- tapply(y,x, length)
        if(any(ns < 30))
            large.yn <- FALSE
        else
            large.yn <- TRUE
        
        ## normality test
        tmp <- tapply(y,x, shapiro.test)
        normal <- TRUE
        for(j in 1:k){
            if(tmp[[j]]$p.value < alpha) normal <- FALSE
        }
        
        tmp <- bartlett.test(y,x)$p.value
        equ.var <- tmp > alpha
    }
    if(k>2){## anova or kruskal
        if((normal||large.yn) && equ.var){ 
            tmp <- kruskal.test(y~as.factor(x))
            pv0 <- tmp$p.value
        }else{
            tmp <- aov(y~x)
            pv0 <- summary(tmp)[[1]][[5]][1]
        }
    }else if(k==2){
        if(normal||large.yn){#t-test
            if(equ.var){
                tmp <- t.test(y~x, paired=FALSE, var.equal=TRUE)
            }else{
                tmp <- t.test(y~x, paired=FALSE, var.equal=FALSE)
            }
            pv0 <- tmp$p.value
        }else{ ## np method: wilcox
            tmp <- wilcox.test(y~x, paired=FALSE)
            pv0 <- tmp$p.value
        }
    }else{ # not applicable if ONE group ONLY.
        pv0 <- NA
    }
    pv0
}

.pv.cat <- function(x){
    nams <- row.names(x)
    sele <- nams == ""
    if(any(sele)) x <- x[!sele,]
    if(is.null(nrow(x))){
        out <- NA
    }else if(nrow(x)==1){
        out <- NA
    }else{
        if(any(x<30)){
            out <- fisher.test(x)$p.value
        }else{
            out <- chisq.test(x)$p.value
        }
    }
    out
}

compare.default <- function(y, x, var.name,
                            alternative = c("two.sided", "less", "greater"),
                            paired = FALSE, 
                            conf.level = 0.95, ...)
{
    stopifnot(conf.level > 0 & conf.level < 1)
    stopifnot(is.logical(paired))
    if(missing(var.name)){
        x.name <- deparse(substitute(x))
        y.name <- deparse(substitute(y))
    }else{
        if(length(var.name) == 2){
            x.name <- var.name[2]
            y.name <- var.name[1]
        }else{
            x.name <- deparse(substitute(x))
            y.name <- deparse(substitute(y))
        }
    }
    
    npmethod <- FALSE
    alpha <- 1 - conf.level

    if(is.numeric(y)){
        if(is.numeric(x)){
            summ <- .summarize1(y,x)
            if(paired){
                if(length(x) != length(y))
                    stop("'x' and 'y' have different sizes")
                sele <- is.na(x) | is.na(y)
                x <- x[!sele]; y <- y[!sele]
            }else{
                x <- summ$x; y <- summ$y
            }
            
            tmp <- .compare.num(y, x,
                                alternative = alternative,
                                paired = paired, 
                                conf.level = conf.level,
                                ...)
            res <- tmp$res
            npmethod <- tmp$npmethod
        }else{
            x <- as.factor(x)
            summ <- .summarize2(y,x)
            x <- summ$x; y <- summ$y
            tmp <- .compare.fac(y, x,
                                alternative = alternative,
                                paired = paired, 
                                conf.level = conf.level,
                                ...)
            res <- tmp$res
            npmethod <- tmp$npmethod
        }
    }else{
        y <- as.factor(y)
        if(is.numeric(x)){
            summ <- .summarize2(x,y)
            x <- summ$x; y <- summ$y
            tmp <- .compare.fac(y, x,
                                alternative = alternative,
                                paired = paired, 
                                conf.level = conf.level,
                                ...)
            res <- tmp$res
            npmethod <- tmp$npmethod
        }else{
            x <- as.factor(x)
            summ <- .summarize3(y,x)
            tmp <- .compare.cat(summ$out, alternative = alternative)
            res <- tmp$res
            npmethod <- NULL
        }
    }
    
    out <- structure(list(stat=summ$out, test=res,
                          x.name=x.name, y.name=y.name,
                          np.method=npmethod), 
                     class = 'comp')
}

.summarize1 <- function(y,x){
    sele.x <- is.na(x)
    sele.y <- is.na(y)
    na.x <- sum(sele.x)
    na.y <- sum(sele.y)
    if(any(sele.x)){
        x <- x[!sele.x]
    }
    if(any(sele.y)){
        y <- y[!sele.y]
    }
    n.x <- length(x)
    n.y <- length(y)
    mu.x <- mean(x); mu.y <- mean(y)
    s.x <- sd(x); s.y <- sd(y)

    g <- c(rep("1",n.x),rep("2",n.y))
    eq.var <- bartlett.test(c(x,y),g)$p.value

    res2 <- data.frame(n=c(n.x,n.y),
                       n.miss=c(na.x,na.y),
                       Mean=c(mu.x, mu.y), Std.Dev=c(s.x,s.y),
                       Var.equal=c(eq.var,NA))
    list(out=res2,x=x,y=y)
}


.summarize2 <- function(y,x){
    n0 <- tapply(y,x,length)
    if(length(x) != length(y))
        stop("'x' and 'y' have different sizes")
    sele <- is.na(y)
    if(any(sele)){
        y <- y[!sele]
        x <- x[!sele]
    }
    n <- tapply(y,x,length)
    n.miss <- n0 - n
    mu <- tapply(y,x, mean)
    s <- tapply(y,x, sd)
    eq.var <- bartlett.test(y,x)$p.value
    k <- nlevels(x)
    res2 <- data.frame(n=n, n.miss=n.miss,
                       Mean=mu, Std.Dev=s,
                       Var.equal=c(eq.var,rep(NA,k-1)))
    row.names(res2) <- levels(x)
    list(out=res2,x=x,y=y)
}

.summarize3 <- function(y,x){
    M <- table(x,y)
    dimnames(M) <- list(x = levels(x),
                        y = levels(y))
    ## M2 <- M
    ## n <- apply(M,1,sum)
    ## n <- length(x)
    ## M1 <- M/n
    list(out = M, x=x,y=y)
}

.compare.cat <- function(x,
                         alternative = c("two.sided", "less", "greater")
)
{
    npmethod <- FALSE
    pv <- NULL
    method <- c("ChiSq", "Fisher's exact")

    tmp <- chisq.test(x)
    pv <- c(pv,tmp$p.value)
    tmp <- fisher.test(x)
    pv <- c(pv,tmp$p.value)

    if(ncol(x)==2){
        method <- c(method, "Prop.test.2")
        x0 <- x[,2]
        n <- apply(x,1,sum)
        tmp <- prop.test(x0,n,alternative=alternative)
        pv <- c(pv,tmp$p.value)
    }
    res <- data.frame(p.value=pv)
    
    row.names(res) <- method
    list(res=res,npmethod=npmethod)
}

.compare.fac <- function(y, x,
                          alternative = c("two.sided", "less", "greater"),
                          paired = FALSE, 
                          conf.level = 0.95, ...)
{
    paired <- FALSE
    out <- NULL
    if(nlevels(x) == 2){
        sele <- x==levels(x)[1]
        out <- .compare.num(y[sele], y[!sele],
                            alternative = alternative,
                            paired = paired, 
                            conf.level = conf.level, ...)
    }else if(nlevels(x) > 2){
        out <- .compare.aov(y, x,conf.level = conf.level)
    }else{
        stop("'x' has only one level...")
    }
    out
}

.compare.aov <- function(y, x, conf.level=0.95)
{
    pv <- NULL
    ## normality test
    npmethod <- FALSE
    alpha <- 1 - conf.level

    method <- c("ANOVA", "Kruskal-Wallis rank sum")
    tmp <- tapply(y,x, shapiro.test)
    npmethod <- FALSE
    k <- nlevels(x)
    for(i in 1:k){
        if(tmp[[i]]$p.value < alpha) npmethod <- TRUE
    }
        
    tmp <- aov(y~x)
    pv <- c(pv,summary(tmp)[[1]][[5]][1])
    tmp <- kruskal.test(y~x)
    pv <- c(pv,tmp$p.value)
    res <- data.frame(p.value=pv)
    row.names(res) <- method
    list(res=res,npmethod=npmethod)
}

.compare.num <- function(y, x,
                          alternative = c("two.sided", "less", "greater"),
                          paired = FALSE, 
                          conf.level = 0.95, ...)
{
    pv <- NULL
    ll <- NULL
    ul <- NULL
    stats <- NULL
    npmethod <- FALSE
    alpha <- 1 - conf.level

    if(paired){
        tmp <- shapiro.test(y-x)
        if(tmp$p.value < alpha) npmethod <- TRUE
        method <- c("Pearson","Spearman","Kendall","t-test (var.equal)",
                    "t-test (var.unequal)", "Mann-Whitney U")
        tmp <- cor.test(y, x, alternative = alternative,
                        conf.level = conf.level,
                        method="pearson")
        pv <- c(pv,tmp$p.value)
        ll <- c(ll, tmp$conf.int[[1]])
        ul <- c(ul, tmp$conf.int[[2]])
        stats <- c(stats, tmp$est)
        tmp <- cor.test(y, x, alternative = alternative,
                        conf.level = conf.level,
                        method="spearman")
        pv <- c(pv,tmp$p.value)
        ll <- c(ll, NA)
        ul <- c(ul, NA)
        stats <- c(stats, tmp$est)
        tmp <- cor.test(y, x, alternative = alternative,
                        conf.level = conf.level,
                        method="kendall")
        pv <- c(pv,tmp$p.value)
        ll <- c(ll, NA)
        ul <- c(ul, NA)
        stats <- c(stats, tmp$est)
    }else{
        tmp1 <- shapiro.test(y)
        tmp2 <- shapiro.test(x)
        if(tmp1$p.value < alpha || tmp2$p.value < alpha)
            npmethod <- TRUE
        method <- c("t-test (var.equal)","t-test (var.unequal)",
                    "Mann-Whiteney U-test")
    }
    
    tmp <- t.test(y, x, alternative = alternative,
                  conf.level = conf.level,paired=paired,
                  var.equal=TRUE)
    pv <- c(pv,tmp$p.value)
    ll <- c(ll, tmp$conf.int[[1]])
    ul <- c(ul, tmp$conf.int[[2]])
    stats <- c(stats, tmp$stat)
    
    tmp <- t.test(y, x, alternative = alternative,
                  conf.level = conf.level,paired=paired,
                  var.equal=FALSE)
    pv <- c(pv,tmp$p.value)
    ll <- c(ll, tmp$conf.int[[1]])
    ul <- c(ul, tmp$conf.int[[2]])
    stats <- c(stats, tmp$stat)
    
    tmp <- wilcox.test(y, x, alternative = alternative,
                       conf.int=TRUE,
                       conf.level = conf.level,
                       paired=paired)
    pv <- c(pv,tmp$p.value)
    ll <- c(ll, tmp$conf.int[[1]])
    ul <- c(ul, tmp$conf.int[[2]])
    stats <- c(stats, tmp$stat)

    res <- data.frame(stat=stats, p.value=pv, lower.limit=ll,
                      upper.limit=ul)
    row.names(res) <- method
    
    list(res=res,npmethod=npmethod)
}


