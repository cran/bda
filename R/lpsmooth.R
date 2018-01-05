## To do: for y = r(x) + e when the data is very large, we can group
## the data and then fit lpr based on the grouped data.  Use a
## masurement error model instead.  That says, we fit y.bar = r(x.bar)
## + e, and x = x.bar + e2.  Laplace error can be assumed for e2 to
## speed up the algorithm.



## we use four method to compute the variance of r(x).

## method 1) Larry Wasserman--nearly unbiased.  This method based on
## an lps object;

## method 2) Rice 1984

## method 3) Gasser et al (1986) -- a variation of method 3.

## method 4) For heteroscedastic errors. Need to estimate based on an
## lpr object. Yu and Jones (2004).

lps.variance <- function(y,x, gridsize, bw, method="Rice"){
    method <- match.arg(tolower(method),
                        c("wasserman","rice","gasser",
                          "heteroscedastic"))
    sele <- is.na(y)|is.na(x)
    y <- y[!sele]; x <- x[!sele]
    
    n <- length(y)
    if(method == "wasserman"){
        ## we fit lpr
        if(missing(bw)) stop("'bw' is missing")
        ## we fit an lpr to get temporary results.  sd.y does not matter.
        fhat <- lpsmooth(y=y,x=x,bw=bw,sd.y=1)
        tmp <- .DesignMatrix(x,bw)
        nu1 <- sum(diag(tmp))
        dl2 <- apply(tmp^2,2,sum)
        nu2 <- sum(dl2)
        res <- sum((y-fhat$y)^2)/(n-2*nu1+nu2)
    }else if(method == "rice"){
        ox <- order(x)
        oy <- y[ox];
        ox <- sort(x)
        res <- 0.5* mean((diff(oy))^2)
    }else if(method == "gasser"){
        ox <- order(x)
        oy <- y[ox];
        ox <- sort(x)
        ai <- (ox[3:n] - ox[2:(n-1)])/(ox[3:n] - ox[1:(n-2)])
        bi <- (ox[2:(n-1)] - ox[1:(n-2)])/(ox[3:n] - ox[1:(n-2)])
        ci2 <- 1/(ai^2+bi^2+1)
        di <- ai * oy[1:(n-2)] + bi * oy[3:n] - oy[2:(n-1)]
        res <- mean(ci2*di^2)
    }else{
        if(missing(bw)) stop("'bw' is missing")
        ## we fit an lpr to get temporary results.  sd.y does not matter.
        
        from <- min(x)
        to <- max(x)
        if(missing(gridsize)) gridsize <- 512L
        stopifnot(gridsize > 10)
        gpoints <- seq(from, to, length=gridsize)
        fhat <- .Fortran(.F_lpsmooth,
                         fx = as.double(x), as.integer(n),
                         as.double(x), as.double(y), as.integer(n),
                         bw = as.double(bw), as.integer(0),
                         as.double(c(from, to)), as.integer(0),
                         ellx = double(n), kappa=double(1))
        z <- 2.0 * log(abs(y - fhat$fx))
        qhat <- .Fortran(.F_lpsmooth,
                         fx = as.double(gpoints), as.integer(gridsize),
                         as.double(x), as.double(z), as.integer(n),
                         bw = as.double(bw), as.integer(0),
                         as.double(c(from, to)), as.integer(0),
                         ellx = double(gridsize), kappa=double(1))
        res <- exp(qhat$fx)
    }
    invisible(res)
}

## this function can be used to compute the design (smoothing) matrix
.DesignMatrix <- function(x,bw){
    stopifnot(bw > 0)
    x <- x[!is.na(x)]; # remove missing values
    ndim = length(x)
    
    L = rep(0,ndim*ndim)
    DM = .Fortran(.F_DesignMatrix, as.double(x), as.integer(ndim),
        as.double(bw), dm=as.double(L))$dm
    matrix(DM, ncol = ndim)
}



## 2015/08/3: add two options lscv=T/F and adaptive=T/F.

## If lscv = FALSE, use the given bandwidth to fit lpr directly.

## If lscv = TRUE and adaptive = FALSE, compute lscv bandwidth and fit
## lpr.  Initial bandwidth should be given.

## If lscv = TRUE and adaptive = TURE, compute lscv bandwidth and then
## compute varying smoothing parameter, then fit lpr.  This algorithm
## could be extremeely slow when the sample size is very large.

lpsmooth <-
    function(y,x,bw,sd.y,lscv=FALSE, adaptive=FALSE,
             from,to,gridsize,conf.level=0.95)
    UseMethod("lpsmooth")

lpsmooth.default <-
    function(y,x,bw,sd.y,lscv=FALSE, adaptive=FALSE,
             from,to,gridsize,conf.level=0.95)
        {
            sele <- is.na(y)|is.na(x)
            y <- y[!sele]; x <- x[!sele]

            size.limit <- 1000
            n <- length(y)
            if(missing(from)) from <- min(x)
            if(missing(to)) to <- max(x)
            stopifnot(to > from)
        
            if(missing(bw)){
                if(n < size.limit){
                    lscv <- 1
                    bw <- bw.nrd(x)
                    adaptive = 1  # use adaptive bandwidth selector
                    lscv <- 0
                }else{
                    stop("'bw' is missing.")
                }
            }else{
                stopifnot(bw>0)
            }

            if(n >= size.limit){
                if(lscv || adaptive)
                    stop("sample size is too large. consider binning the data")
                lscv <- 0
                adaptive <- 0
            }else{
                lscv <- ifelse(lscv, 1, 0)
                adaptive <- ifelse(adaptive, 1, 0)  
            }
            
            ##  Compute the variance based on the raw data
            if(missing(sd.y)){
                sd.y <- lps.variance(y=y,x=x,bw=bw)
            }else{
                stopifnot(is.numeric(sd.y))
                stopifnot(!any(sd.y < 0))
            }
        
            if(missing(gridsize)) gridsize <- 512L
            stopifnot(gridsize > 10)


            gpoints <- seq(from, to, length=gridsize)
            n <- length(x)
            stopifnot(length(y) == n)
            if(any(is.na(x)|is.na(y)))
                stop("Missing value(s) in 'x' and/or 'y'")
            if(any(!is.finite(x)|!is.finite(y)))
                stop("Inifite value(s) in 'x' and/or 'y'")
            out <- .Fortran(.F_lpsmooth,
                            fx = as.double(gpoints), as.integer(gridsize),
                            as.double(x), as.double(y), as.integer(n),
                            bw = as.double(bw), as.integer(lscv),
                            as.double(c(from, to)), as.integer(adaptive),
                            ellx = double(gridsize), kappa=double(1))
            
            ## print(out$kappa)
            stopifnot(conf.level<1 & conf.level >0)
            
            cv <-  .Fortran(.F_tubecv,
                            cv=as.double(out$kappa), as.double(conf.level))$cv
            ## print(cv)
            
            y = out$fx
            sele1 = is.na(y) | !is.finite(out$fx)
            if(any(sele1)) y[sele1] = 0.0
            MOE <- cv * out$ellx * sd.y
            ll = y - MOE; ##ll[ll<0] <- 0
            ul = y + MOE;
            
            pars <- list(cv=cv, kappa=out$kappa, bw=out$bw, npar=NULL)
            
            structure(list(y = y, x = gpoints,
                           conf.level = conf.level,
                           type="lpsmooth",
                           pars = pars,
                           ucb=ul, lcb=ll,
                           call = match.call()
                           ),
                      class = 'histosmooth')
        }

lpsmooth.histogram <-
    function(y,x,bw,sd.y,lscv=FALSE, adaptive=FALSE,
             from,to,gridsize,conf.level=0.95)
    {
        f.call <- match.call()
        xhist <- binning(y)
        out <- lpsmooth(y=xhist,from=from,to=to,
                        gridsize=gridsize,
                        conf.level=conf.level)
        out$call <- f.call
        out
   }

lpsmooth.bdata <-
    function(y,x,bw,sd.y,lscv=FALSE, adaptive=FALSE,
             from,to,gridsize,conf.level=0.95)
    {
        f.call <- match.call()
        out <- .histonpr(x=y,from=from,to=to,
                         gridsize=gridsize,
                         conf.level=conf.level)
        out$call <- f.call
        out
   }

## by default, split the histogram into k*nclass bins
.split <- function(x,k)
{
    if(missing(k)) k <- 2
    k <- round(k)
    if(k > 3)
        warning("Number of sub-classes 'f' might be too large")
    if(k <= 1){
        out <- x
    }else{
        a <- min(x$breaks)
        b <- max(x$breaks)
        F <- x$counts
        nbin <- length(F)
        M <- k*nbin + 1
        bw <- diff(x$breaks)/k
        bws <- rep(bw, rep(k, nbin))
        breaks <- a + c(0,cumsum(bws))  # good for unequal bins
        
        out <- .bootkde(x)
        fx <- out$y;
        x0 <- out$x
        Fx <- cumsum(fx); Fx <- Fx/max(Fx)
        Fb <- approx(x0, Fx, breaks)$y
        m <- length(Fb)
        if(is.na(Fb[m])) Fb[m] <- 1
        if(is.na(Fb[1])) Fb[1] <- 0
        mx <- diff(Fb)
        counts <- NULL
        for(i in 1:nbin){
            total <- sum(mx[((i-1)*k+1):(i*k)])
            if(F[i]==0 || total ==0){
                counts <- c(counts, rep(0,k))
            }else{
                tmp2 <- 0; tmp <- 0
                for(j in 1:(k-1)){
                    p <- mx[k*(i-1)+j]/total
                    tmp <- round(p*F[i])
                    counts <- c(counts, tmp)
                    tmp2 <- tmp2 + tmp
                }
                counts <- c(counts, max(0,F[i]-tmp2))
            }
        }
        
        nclass <- length(breaks)-1
        mids <- 0.5 * (breaks[-1]+breaks[-(nclass+1)])
        bw <- diff(breaks)
        
        out <- structure(
            list(breaks = breaks,
                 counts = counts,
                 mids = mids,
                 size = sum(counts),
                 nclass = nclass,
                 equalwidth=x$equalwidth,
                 top.coded=x$top.coded,
                 name = x$name,
                 call = x$call),
            class="bdata")
    }
    out
}

.histonpr <- function(x,from,to,gridsize,conf.level=0.95)
    {
        if(!x$equalwidth)
            stop("Not applicable for unequal width bins")
        if(missing(conf.level)) conf.level <- 0.95
        
        xhist <- x
        mcounts <- mean(x$counts)
        if(mcounts > 50)
            x <- .split(x,k=2)
        else if(mcounts > 100)
            x <- .split(x, k=3)
        
        F <- x$counts;
        X <- x$mids;
        nbin <- length(X);
        n <- sum(F)
        a <- x$breaks[1]; b <- x$breaks[nbin+1];

        ##        if(missing(bw)){
        lscv <- 1
        A <- x$breaks[1:nbin]; B <- x$breaks[-1]
        x0 <- runif(n,rep(A,F),rep(B,F))
        bw <- bw.nrd(x0)
        adaptive = 1  # use adaptive bandwidth selector
        ##        }else{
        ##            lscv <- 0
        ##            adaptive = 0  # don't use adaptive bandwidth selector
        ##        }

        if(missing(from)) from <- a - 3*bw
        if(missing(to)) to <- b + 3*bw
        stopifnot(to > from)

        if(missing(gridsize)) gridsize <- 512L
        stopifnot(gridsize > 10)
        ##  root-transformation
        Y = sqrt(nbin/n*(F+0.25))
        ##print(x$counts)
        ## call npr(...) to estimate the nonparametric regression function
        out <- lpsmooth(y=Y, x=X, bw=bw, sd.y= 0.5*sqrt(nbin/n),
                   from=from, to=to, gridsize=gridsize,
                   conf.level=conf.level)
        y <- out$y; x <- out$x
        ll <- out$lcb; ul <- out$ucb;
        y <- y*y; tot.mass <- sum(y)*(to-from)/gridsize;
        ll[ll<0] <- 0 # density can not be negative
        ll <- ll*ll/tot.mass
        ul <- ul*ul/tot.mass
        
        f0 <- y/tot.mass
        x0 <- x;

        pars <- list(cv=out$cv, kappa=out$kappa, bw=out$bw, npar=NULL)

        structure(list(y = y/tot.mass, x = x,
                       conf.level = out$conf.level,
                       type="lpsmooth",
                       pars = pars,
                       ucb=ul, lcb=ll,
                       call = match.call()
                       ),
                  class = 'histosmooth')

    }


