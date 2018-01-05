## A graphical method
## migrated from S-Plus code  by Dr. Jiayang Sun.

## 2017/03/31: do eda to qualitative data as well.
eda <- function(x, plot=FALSE)
    UseMethod("eda")

eda.default <- function(x, plot=FALSE)
{
    if(is.numeric(x)){
        out <- .eda.num(x)
    }else if(is.logical(x)){
        out <- .eda.log(x)
    }else{
        x <- as.character(x)
        sele <- x==''
        if(any(sele))
            x[sele] <- "BLANK"
        x <- as.factor(x)
        out <- .eda.char(x)
    }
    out
}

eda.matrix <- function(x, plot=FALSE)
{
    warning("The matrix will be analyzed by columns")
    x <- as.data.frame(x)
    eda.data.frame(x)
}

eda.data.frame <- function(x, plot=FALSE){
    ## we create a table with results in one row for each variable
    ## (numerical).
    m <- ncol(x)
    n <- nrow(x)
    ## separate the data.frame into two parts: "dat1" for numeric
    ## variables and "dat2" for the rest.  2017/06/05: problem if
    ## there is only one column for either data set.  We don't split
    ## the data and do the analysis directly.
    cat("\nCategorical variables:\n") # show categorical data results first

    nams <- names(x)
    res.num <- NULL
    res.cat <- NULL
    sele <- NULL
    for(i in 1:m){
        if(is.numeric(x[,i])){
            sele <- c(sele, TRUE)
            tmp <- .eda.num(x[,i], plot=plot)
            out <- c(tmp$n, tmp$n.na, tmp$mean, tmp$sd, tmp$p.norm)
            res.num <- rbind(res.num, out)
        }else{
           sele <- c(sele, FALSE)
           cat("\nVariable name :=", nams[i],"\n")
           tmp <- as.character(x[,i])
           res <- table(tmp)
           tmp <- round(100*as.numeric(res)/n,2)
           out <- data.frame(Categories=names(res),
                             Counts=as.numeric(res),
                             Percentages=tmp)
           print(out)
        }
    }
    if(sum(sele)>1){
        res.num <- as.matrix(res.num)
        res.num <- as.data.frame(res.num)
        names(res.num) <- c("n", "n.miss","Mean","SD","normality")
        row.names(res.num) <- nams[sele]
        cat("\nNumerical variables:\n")
        print(res.num) 
    }else if(sum(sele)==1){
        cat("\nVariable name :=", nams[sele])
        cat("\n\tn=",res.num[1])
        cat("\n\tn.miss=",res.num[2])
        cat("\n\tMean=",res.num[3])
        cat("\n\tStd.Dev=",res.num[4])
        cat("\n\tNormality p.value=",res.num[5],"\n")
    }
    invisible(res.num)
}
           
.eda.char <- function(x, plot=FALSE)
{
    xtbl <- table(x)
    x.freq <- round(as.numeric(xtbl)/length(x)*100,2)
    nam <- names(xtbl)
    out <- data.frame(Freq=as.numeric(xtbl), Perc=x.freq)
    row.names(out) <- nam
    structure(list(x=out, type='char'),
              class = 'desc')
}

.eda.log <- function(x, plot=FALSE)
{
    n <- length(x)
    x.na <- is.na(x)
    n.na <- sum(x.na)
    if(n.na > 0){
        x <- x[!x.na];
    }
    if(n.na==n){
        out <- NULL
    }else{
        xtbl <- table(x)
        x.freq <- round(as.numeric(xtbl)/length(x)*100,2)
        nam <- names(xtbl)
        out <- data.frame(Freq=as.numeric(xtbl), Perc=x.freq)
        row.names(out) <- nam
    }
    structure(list(x=out, type='char'),
              class = 'desc')
}

.eda.num <- function(x, plot=FALSE)
{
    n <- length(x)
    x.na <- is.na(x)
    n.na <- sum(x.na)
    if(n.na > 0){
        x <- x[!x.na];
    }
    mu <- mean(x)
    s <- sd(x)
    md <- median(x)
    m2 <- mean((x-mu)^2)
    m3 <- mean((x-mu)^3)
    m4 <- mean((x-mu)^4)
    skewness <- m3/m2^1.5
    kurtosis <- m4/m2^2-3
    structure(list(n=n, n.na=n.na,
                   mean=mu,m=md,sd=s,
                   skewness=skewness,
                   kurtosis=kurtosis,
                   p.norm=shapiro.test(x)$p.value,
                   type='num',x=x), 
              class = 'desc')
}

print.desc <- function(x,...){
    if(x$type=='num'){
        cat("\n\tTotal # of obs =",x$n)
        cat("\n\tMissing values =",x$n.na)
        cat("\n\tMean =",x$mean)
        cat("\n\tMedian =",x$m)
        cat("\n\tStd.Dev =",x$sd)
        cat("\n\tSkewness =",x$skewness)
        cat("\n\tKurtosis =",x$kurtosis)
        cat("\n\tNormality (Shapiro-Wilk) =",x$p.norm,"\n")
    }else{
        if(is.null(x$x)){
            cat("\n\tAll missing values...\n")
        }else{
            if(nrow(x$x)<15){
                print(x$x,...)
            }else{
                cat("\n\tNum of levels =",nrow(x$x))
                cat("\n\tToo many levels. Results not shown...\n")
            }
        }
    }
    invisible(x)
}

plot.desc <- function(x,...){
    x0 <- x
    if(x$type=='num'){
        x <- x$x
        par(mfrow = c(2., 3.))
        qqnorm(x)
        qqline(x)
        boxplot(x)
        title("Boxplot")
        hist(x, main = "Histogram")
        iqd <- summary(x)[5.] - summary(x)[2.]
        plot(density(x, width = 2. * iqd),
             main = "Density Plot",
             ylab = "Density", type = "l")
        ts.plot(x)
        title("Time Series Plot")
        acf(x)
    }
    invisible(x0)
}
