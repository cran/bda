## Created on 06/18/2014.

gof.test <- function(x, dist, pars, method='ks')
    UseMethod("gof.test")

gof.test.default <- function(x, dist, pars, method='ks')
{
    if(!is.character(method)) method <- 'ks'
    
    method <- match.arg(tolower(method),
                        c("chisquare","ks","kolmogorov-smirnov"))
    if(missing(dist)||missing(pars)){
        if(method=="chisquare"){
            out <- .gof.chi(x$data, dist=class(x),pars=x$pars)
        }else{
            out <- .gof.ks(x$data, dist=class(x),pars=x$pars)
        }
    }else{
        if(method=='chisquare'){
            out <- .gof.chi(x,dist=dist, pars=pars)
        }else{
            out <- .gof.ks(x,dist=dist, pars=pars)
        }
    }
    
    out
}

gof.test.FMKL <- function(x, dist, pars, method='ks')
{
    xhist <- x$xhist
    gof.test(xhist, dist="fmkl", pars=x$pars, method=method)
}

gof.test.histogram <- function(x, dist, pars, method='ks')
{
    xhist <- binning(x)
    gof.test(xhist, dist=dist, pars=pars, method=method)
}

gof.test.GB <- function(x, dist, pars, method='ks')
{
    xhist <- x$xhist
    gof.test(xhist, dist="gb", pars=x$pars, method=method)
}

gof.test.bdata <- function(x, dist, pars, method='ks')
{
    if(!is.character(method)) method <- 'ks'
    
    method <- match.arg(tolower(method),
                        c("chisquare","ks","kolmogorov-smirnov"))

    if(method=='chisquare'){
        out <- .gof.chi(x,dist=dist, pars=pars)
    }else{
        out <- .gof.ks(x,dist=dist, pars=pars)
    }
    out
}

.gof.chi <- function(x, dist, pars){

    if(missing(dist)) dist <- 'normal'
    dist <- match.arg(tolower(dist),## add RST-GLD later
                      c("normal", "weibull", "fmkl", "gb",
                        "pareto", "lognormal", "gpd",
                        "mixlnorm","mixlnorm", "gld"))
    ##out <- .compact(x)
    if(class(x)=="bdata"||class(x)=="histogram"){
        cts <- x$counts;
        xn <- x$breaks
    }else
        stop("'x' must be binned.")
    
    Fx <- switch(dist,
                 "normal" = pnorm(xn, pars),
                 "weibull" = pweibull(xn, pars),
                 "fmkl" = pgld(xn, pars),
                 "gld" = pgld(xn, pars),
                 "pareto" = pPareto(xn,pars[1],pars[2]),
                 "gpd" = pGPD(xn,pars[1],pars[2],pars[3]),
                 "lnorm" = plnorm(xn,pars[1],pars[2]),
                 "lognormal" = plnorm(xn,pars[1],pars[2]),
                 "mixlnorm" = pmlnorm(xn,pars[1],pars[2],pars[3]),
                 "mixnorm" = pmixnorm(xn,pars[1],pars[2],pars[3]),
                 "gb" = pgbeta(xn, pars)
                 )
    Fx[1] <- 0.0
    Fx[length(Fx)] <- 1.0
    Nhat <- sum(cts) * diff(Fx)
    X2 <- sum((Nhat - cts)^2/Nhat)
    ps <- x$nclass - length(pars) - 1
    parameters <- ifelse(ps>0, ps, 1)
    pval <- pchisq(X2, parameters, lower.tail=FALSE)
    dname <- x$name
    
    names(X2) <- "Chi-square statistic"
    names(parameters) <- "df"
    RVAL <- list(statistic = X2, parameter = parameters,
                 method = "Chi-square Goodness-of-Fit Test",
                 alternative = "two-sided",
                 p.value= pval, observed = cts,
                 expected = Nhat, residuals = sqrt(X2),
                 data.name=dname)
    class(RVAL) <- "htest"
    return(RVAL)
}


.gof.ks <- function(x, dist, pars){

    if(class(x)=="bdata"){
        dname <- x$name
        xn <- x$breaks
        Fn <- c(0,cumsum(x$counts)/sum(x$counts))
        k <- length(xn)
        xn <- xn[-c(1,k)]
        Fn <- Fn[-c(1,k)]
    }else{
        dname <- deparse(substitute(x))
        x.t <- table(x)
        xn <- as.numeric(x.t)
        Fn <- cumsum(x.t)/sum(x.t)
        k <- length(xn)
        xn <- xn[-k]
        Fn <- Fn[-k]
    }

    if(missing(dist))
        stop("'dist' missing")
    if(missing(pars))
        stop("'pars' missing")
    
    dist <- match.arg(tolower(dist),## add RST-GLD later
                      c("normal", "weibull", "fmkl", "gb",
                        "pareto", "lognormal", "gpd",
                        "mixlnorm","mixlnorm", "gld",
                        "pareto"))
    Fx <- switch(dist,
                 "normal" = pnorm(xn, pars),
                 "weibull" = pweibull(xn, pars),
                 "fmkl" = pgld(xn, pars),
                 "gld" = pgld(xn, pars),
                 "pareto" = pPareto(xn,pars[1],pars[2]),
                 "gpd" = pGPD(xn,pars[1],pars[2],pars[3]),
                 "lognormal" = plnorm(xn,pars[1],pars[2]),
                 "mixlnorm" = pmlnorm(xn,pars[1],pars[2],pars[3]),
                 "mixnorm" = pmixnorm(xn,pars[1],pars[2],pars[3]),
                 "gb" = pgbeta(xn, pars),
                 "pareto" = pPareto(xn, xm=pars[1], alpha=pars[2])
                 )

    D <- max(abs(Fx - Fn))
    ##pval <- .Fortran(.F_KSPvalue, pv=as.double(D))$pv
    pval <- .Fortran(.F_KSP2x,
                     pv=as.double(D),
                     as.integer(sum(Fn)))$pv

    ##pval <- 1-.Call(C_pKolmogorov2x, D, sum(Fn))

    names(D) <- "D"
    RVAL <- list(statistic = D,
                 method = "One-sample Kolmogorov-Smirnov Test",
                 alternative = "two-sided",
                 p.value= pval,
                 data.name=dname)
    class(RVAL) <- "htest"
    return(RVAL)
}

