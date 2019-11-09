.aic <- function(llk,npar,size){
    ntotal <- size
    K <- npar
    llk0 <- llk
    AIC = -2.0 * llk0 + 2.0 * K;
    AICc = AIC + 2.0 * K * (K+1.) / (ntotal - K - 1.0);
    BIC = -2.0 * llk0 + log(0.5*ntotal/pi) * K;
    data.frame(AIC=AIC, BIC=BIC, AICc=AICc)
}

.llk <- function(x,dist,pars){
    stopifnot(class(x)=='bdata')
    dist <- match.arg(tolower(dist),
                      c('gld','pareto','gpd'))
    if(dist == 'gld'){
        n <- x$counts
        q <- x$breaks
        tmp <- diff(pgld(q, pars))
        if(any(is.na(tmp))){
            llk <- 999999999999.+runif(1)
        }else if(any(tmp==0)){
            llk <- 999999999999.+runif(1)
        }else{
            llk <- sum(log(tmp)*n)
        }
        
    }else if(dist == 'pareto'){
        n <- x$counts
        q <- x$breaks
        tmp <- diff(pPareto(q, pars[1], pars[2]))
        if(any(tmp==0)){
            llk <- 999999999999.+runif(1)
        }else{
            llk <- sum(log(tmp)*n)
        }
    }else if(dist == 'gpd'){
        n <- x$counts
        q <- x$breaks
        tmp <- diff(pGPD(q, pars[1], pars[2], pars[3]))
        if(any(is.na(tmp))){
            llk <- 999999999999.+runif(1)
        }else if(any(tmp==0)){
            llk <- 999999999999.+runif(1)
        }else{
            llk <- sum(log(tmp)*n)
        }
    }else{
        warning("log-liklihood not computable for this distribution")
        llk <- NA
    }
    llk
}

