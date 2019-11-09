## FSD is fitted using four different methods: Pareto, GPD
## (Generalized Pareto Distribution), Lognormal mixture models, and
## GLD (Generalized Lambda Distribution).  Four main functions (1)
## fit.FSD, (2) print.FSD, (3) lines.FSD, and (4) plot.FSD.

## Firm size data are usually of histogram-type. If not, we bin the
## data first then fit lognormal distributions.
fit.FSD <- function(x, method, k=1, lbound, ubound){
    
    if(class(x) != 'bdata')
        x <- binning(x)
    if(missing(method)){
        out <- fit.Pareto(x, lbound=lbound,method='mle')
        pars <- out$pars
        meth <- "Pareto"
        ModSel <- out$ModSel
        BIC1 <- out$ModSel$BIC
        nams <- meth
        
        out <- fit.GPD(x)
        BIC2 <- out$ModSel$BIC
        ModSel <- rbind(ModSel, out$ModSel)
        if(BIC2 < BIC1 && BIC2 > 0){
            BIC1 <- BIC2
            pars <- out$pars
            meth <- "GPD"
        }
        nams <- c(nams, "GPD")
                  
        out <- fit.GLD(x,lbound=lbound,ubound=ubound)
        BIC2 <- out$ModSel$BIC
        ModSel <- rbind(ModSel, out$ModSel)
        if(BIC2 < BIC1 && BIC2 > 0){
            BIC1 <- BIC2
            pars <- out$pars
            meth <- "GLD"
        }
        nams <- c(nams, "GLD")

        K <- min(5, round(length(x$counts)/2))
        for(i in 1:K){
            out <- fit.lnmix(x,method='lognormal', k=i,
                             lbound=lbound, ubound=ubound)
            BIC2 <- out$ModSel$BIC
            ModSel <- rbind(ModSel, out$ModSel)
            tmp <- paste("Lognormal(k=",i,")",sep='')
            nams <- c(nams, tmp)
            if(BIC2 < BIC1 && BIC2 > 0){
                BIC1 <- BIC2
                pars <- out$pars
                meth <- tmp
            }
        }
        ModSel <- data.frame(ModSel)
        rownames(ModSel) <- nams
    }else{
        method <- match.arg(tolower(method),
                            c("lognormal", "lnmix", "lnmm", "fmkl", "gld",
                              "pareto", "gpd"))
        if(method=="lognormal" || method=="lnmix" || method=="lnmm"){
            out <- fit.lnmix(x,method='lognormal', k=k,
                             lbound=lbound, ubound=ubound)
            pars <- out$pars
            meth <- paste("Lognormal(k=",k,")",sep='')
            ModSel <- out$ModSel
        }else if(method=="fmkl" || method=="gld"){
            out <- fit.GLD(x,lbound=lbound,ubound=ubound)
            pars <- out$pars
            meth <- "GLD"
            ModSel <- out$ModSel
        }else if(method=="pareto"){
            out <- fit.Pareto(x, lbound=lbound,method='mle')
            pars <- out$pars
            meth <- "Pareto"
            ModSel <- out$ModSel
        }else if(method=="gpd"){
            out <- fit.GPD(x)
            pars <- out$pars
            meth <- "GPD"
            ModSel <- out$ModSel
        }
    }
    structure(list(xhist = x,
                   method = meth,
                   pars = pars,
                   ModSel = ModSel
                   ),
              class='FSD')
}


print.FSD <- function(x,...){
    cat("\nFitted distribution:", x$method)
    cat("\nParameter estimates:\n")
    print(round(x$pars,3))
    cat("\nModel selection:\n")
    print(round(x$ModSel,3))
    cat("\n")
}


.plotzpif1 <- function(x,...){
    stopifnot(class(x)=="bdata")
    x0 <- x$breaks
    k <- length(x0)
}



plot.FSD <- function(x,xlim=NULL,ngrid=1000,...){
    if(is.null(xlim)){
        xbrks <- x$xhist$breaks
        k <- length(xbrks)
        if(is.finite(xbrks[1]))
            xmin <- xbrks[1]
        else
            xmin <- xbrks[2]
        if(is.finite(xbrks[k]))
            xmax <- xbrks[k]
        else
            xmax <- xbrks[k-1]
        
        xlim <- c(xmin, xmax)
    }
    
    plot(x$xhist,xlim=xlim,...)
    lines(x,xlim=xlim,ngrid=ngrid,...)
}

lines.FSD <- function(x,xlim=NULL,ngrid=1000,log=FALSE,...){
    ngrid <- round(ngrid)
    stopifnot(ngrid > 1)
    if(is.null(xlim)){
        xbrks <- x$xhist$breaks
        k <- length(xbrks)
        if(is.finite(xbrks[1]))
            xmin <- xbrks[1]
        else
            xmin <- xbrks[2]
        if(is.finite(xbrks[k]))
            xmax <- xbrks[k]
        else
            xmax <- xbrks[k-1]
        xlim <- c(xmin, xmax)
    }

    if(log){
        if(xlim[1] < 0){
            cat("'log-transformation' is not applicable")
            x0 <- seq(xlim[1], xlim[2], length=ngrid)
        }else{
            x0 <- seq(log(xlim[1]), log(xlim[2]), length=ngrid)
        }
    }else{
        x0 <- seq(xlim[1], xlim[2], length=ngrid)
    }

    meth <- tolower(x$method)
    if(meth=='normal'){
        f0 <- dmixnorm(x0,x$pars$p,x$pars$mean,x$pars$sd)
    }else if(meth=='gld'){
        f0 <- dgld(x0, x$pars)
    }else if(meth=='pareto'){
        f0 <- dPareto(x0, x$pars[1], x$pars[2])
    }else if(meth=='gpd'){
        f0 <- dGPD(x0, x$pars[1], x$pars[2], x$pars[3])
    }else{
        f0 <- dmlnorm(x0,x$pars$p,x$pars$mean,x$pars$sd)
    }
    sele <- f0 > 0
    lines(x0[sele],f0[sele],...)
}

