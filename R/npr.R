npr <-
    function(y,x,sd.x,bw,kernel="decon",
             optimal=FALSE,adaptive=FALSE,
             x0,from,to,gridsize,conf.level=0.95)
{
    if(missing(x0)){
        if(missing(from)) from <- min(x,na.rm=TRUE)
        if(missing(to)) to <- max(x,na.rm=TRUE)
        stopifnot(to > from)
        if(missing(gridsize)) gridsize <- 512L
        stopifnot(gridsize > 10)
        gpoints <- seq(from, to, length=gridsize)
    }else{
        x0 <- x0[!is.na(x0)]
        gpoints <- unique(sort(x0))
        gridsize <- length(gpoints)
    }

    kernel <- match.arg(tolower(kernel),
                        c("normal","gauss","nw","decon","lp",
                          "nadaraya-watson","repeated"))

    if(is.matrix(x)){
        if(is.matrix(y)){
            if(nrow(y)==nrow(x)){
                my <- apply(y,1,mean,na.rm=TRUE)
            }else{
                stop("Different number of records for 'x' and 'y'")
            }
        }else{
            if(length(y)==nrow(x)){
                my <- y
            }else{
                stop("Different number of records for 'x' and 'y'")
            }
        }
        x1 <- NULL
        y1 <- NULL
        w1 <- NULL
        nx <- apply(!is.na(x),1,sum)
        for(i in 1:ncol(x)){
            y1 <- c(y1,my)
            x1 <- c(x1,x[,i])
            w1 <- c(w1, 1/nx)
        }
        sele <- !is.na(x1) & !is.na(y1)
        x1 <- x1[sele]
        y1 <- y1[sele]
        w1 <- w1[sele]
        ##if(any(!is.finite(w1)))
        ##stop("invalid weights. check raw data.")
        out <- .npr2(y=y1, x=x1, w=w1,
                     bw=bw, optimal=optimal,adaptive=adaptive,
                     gpoints=gpoints,gridsize=gridsize,
                     conf.level=conf.level)
    }else{
        if(is.matrix(y)){
            if(nrow(y)==length(x)){
                my <- apply(y,1,mean,na.rm=TRUE)
            }else{
                stop("Different number of records for 'x' and 'y'")
            }
        }else{
            my <- y
        }
        if(kernel=="repeated"){
            sele <- !is.na(x) & !is.na(my)
            x1 <- x[sele]
            y1 <- my[sele]
            w1 <- rep(1,length(x1))
            out <- .npr2(y=y1, x=x1, w=w1,
                         bw=bw, optimal=optimal,adaptive=adaptive,
                         gpoints=gpoints,gridsize=gridsize,
                         conf.level=conf.level)
        }else{
            out <- .npr1(y=my, x=x,sd.x=sd.x,bw=bw,kernel=kernel,
                         optimal=optimal,gpoints=gpoints,
                         gridsize=gridsize,conf.level=conf.level)
        }
    }
    out
}

.ahopt <- function(x){
    k <- 3
    x <- x[!is.na(x)]
    x <- sort(x)
    n <- length(x)
    if(n<10) stop("too few data points")
    l <- c(1:n)-k
    r <- c(1:n)+k
    sele.l <- l < 1
    l[sele.l] <- 1
    r[sele.l] <- 1+2*k

    sele.r <- r > n
    l[sele.r] <- n-2*k
    r[sele.r] <- n
    h <- x[r] - x[l]
}

.npr2 <-
    function(y,x,w,bw,optimal=FALSE,adaptive=FALSE,
             gpoints,gridsize,conf.level=0.95)
{
    nx <- length(x)
    if(missing(bw)){
        hopt <- diff(range(x))/nx*5
    }else{
        if(length(bw) != 1)
            stop("invalid 'bw'.")
        if(bw < 0)
            stop("invalid 'bw'.")
        hopt <- bw
    }
    
    if(adaptive&&nx>10){
        hopt <- bw.nrd0(x) #.ahopt(x)
        out <- .Fortran(.F_awlprNorm,
                        y = as.double(gpoints),
                        as.integer(gridsize),
                        as.double(x),
                        as.double(y),
                        as.double(w),
                        as.integer(nx),
                        h = as.double(hopt))
        hopt <- mean(out$h)
        rx <- out$y
    }else{
        optim <- ifelse(optimal, 1,-1) #<=0 no automatic bandwidth selection
        out <- .Fortran(.F_wlprNorm,
                        y = as.double(gpoints),
                        as.integer(gridsize),
                        as.double(x),
                        as.double(y),
                        as.double(w),
                        as.integer(nx),
                        h = as.double(hopt),
                        as.integer(optim))
        hopt <- out$h
        rx <- out$y
    }
    
    list(y = rx, x = gpoints, bw=hopt)
}

.npr1 <-
    function(y,x,sd.x,bw,kernel="decon",optimal=FALSE,
             gpoints,gridsize,conf.level=0.95)
{
    nx <- length(x)
    if(length(y) != nx)
        stop("'x' and 'y' have different lengths")
    if(missing(sd.x)){
        sd.x <- rep(0, nx)
        homo <- TRUE
    }else{
        if(length(sd.x)>1){
            if(length(sd.x) != nx)
                stop("'x' and 'sd.x' have different lengths")
            if(any(sd.x < 0))
                stop("invalid 'sd.x' value(s)")
            homo <- FALSE
        }else{
            sd.x <- rep(sd.x,nx)
            homo <- TRUE
        }
    }
    
    sele <- is.na(y)|is.na(x)
    y <- y[!sele];y.raw <- y
    x <- x[!sele];x.raw <- x
    sd.x <- sd.x[!sele]
    nx <- length(x)
    loo <- 1
    optim <- ifelse(optimal, 1,-1) #<=0 no automatic bandwidth selection
    if(missing(bw)){
        bw <- .bwdmise(x, sd.x)
    }else if(any(bw <= 0)){
        stop("invalid 'bw'.")
    }
    
    if(kernel=="decon"){# support lp=0 or 1, nothing else yet
        if(homo){
            bw <- mean(bw, na.rm=TRUE) #in case bw has length >1
            out <- .Fortran(.F_DkNpReg,
                            as.double(x),
                            as.double(y),
                            as.double(sd.x),
                            as.integer(nx),
                            h = as.double(bw),
                            y = as.double(gpoints),
                            as.integer(gridsize),
                            as.double(loo),
                            gcv = as.double(1))
            gcv <- out$gcv
            mx <- out$y
            hopt <- bw
        }else{
            if(optimal){
                imax <- 10
                if(length(bw)==1){
                    a <- 0.9*bw
                    b <- 1.2*bw
                }else{
                    a <- min(bw)
                    b <- max(bw)
                }
                h0 <- seq(a, b, length=imax)
                gcv0 <- NULL
                for(h1 in h0){
                    out <- .Fortran(.F_nprHLap,
                                    y = as.double(x),
                                    as.integer(nx),
                                    as.double(x),
                                    as.double(y),
                                    as.double(sd.x),
                                    as.integer(nx),
                                    h = as.double(h1),
                                    gcv = as.double(1))
                    gcv0 <- c(gcv0, out$gcv)
                }
                ##print(gcv0)
                sele <- !is.na(gcv0)
                h0 <- h0[sele]
                gcv0 <- gcv0[sele]
                i <- which(gcv0==min(gcv0))[1]
                bw <- h0[i]
                gcv <- list(x=h0,y=gcv0)
                if(i == 1 || i == length(h0)){
                    warning("optimal bw found on boundary!")
                }else{
                    a <- h0[i-1]
                    b <- h0[i+1]
                    h0 <- seq(a, b, length=10)
                    gcv0 <- NULL
                    for(h1 in h0){
                        out <- .Fortran(.F_nprHLap,
                                        y = as.double(x),
                                        as.integer(nx),
                                        as.double(x),
                                        as.double(y),
                                        as.double(sd.x),
                                        as.integer(nx),
                                        h = as.double(h1),
                                        gcv = as.double(1))
                        gcv0 <- c(gcv0, out$gcv)
                    }
                    ##print(gcv0)
                    sele <- !is.na(gcv0)
                    h0 <- h0[sele]
                    gcv0 <- gcv0[sele]
                    i <- which(gcv0==min(gcv0))[1]
                    bw <- h0[i]
                    gcv <- list(x=h0,y=gcv0)
                }
                out <- .Fortran(.F_nprHLap,
                                y = as.double(gpoints),
                                as.integer(gridsize),
                                as.double(x),
                                as.double(y),
                                as.double(sd.x),
                                as.integer(nx),
                                h = as.double(bw),
                                gcv = as.double(-1))
                hopt <- bw
                mx <- out$y
            }else{
                out <- .Fortran(.F_nprHLap,
                                y = as.double(gpoints),
                                as.integer(gridsize),
                                as.double(x),
                                as.double(y),
                                as.double(sd.x),
                                as.integer(nx),
                                h = as.double(bw),
                                gcv = as.double(-1))
                gcv <- list(x=bw,y=out$gcv)
                hopt <- bw
                mx <- out$y
            }
        }
    }else if(kernel=='lp'){
        ##out <- lpsmooth(y=y,x=x,bw=bw,lscv=TRUE,
        ##from=from,to=to,gridsize=gridsize,
        ##conf.level=conf.level)
        ##mx <- out$y
        ##gpoints <- out$x
        ##hopt <- out$pars
        ##gcv <- NA
        if(homo){
            out <- .Fortran(.F_lprLap,
                            y = as.double(gpoints),
                            as.integer(gridsize),
                            as.double(x),
                            as.double(y),
                            as.double(sd.x[1]),
                            as.integer(nx),
                            h = as.double(bw),
                            gcv = as.double(-1))
            gcv <- list(x=bw,y=out$gcv)
            hopt <- bw
            mx <- out$y
        }else{
            out <- .Fortran(.F_lprHLap,
                            y = as.double(gpoints),
                            as.integer(gridsize),
                            as.double(x),
                            as.double(y),
                            as.double(sd.x),
                            as.integer(nx),
                            h = as.double(bw),
                            gcv = as.double(-1))
            gcv <- list(x=bw,y=out$gcv)
            hopt <- bw
            mx <- out$y
        }
    }else{
        bw <- mean(bw, na.rm=TRUE) #in case bw has length >1
        out <- .Fortran(.F_NWReg,
                        as.double(x),
                        as.double(y),
                        as.integer(nx),
                        bw = as.double(bw),
                        y = as.double(gpoints),
                        as.integer(gridsize),
                        as.double(loo),
                        as.integer(optim),
                        gcv = as.double(loo))
        gcv <- out$gcv
        mx <- out$y
        hopt <- out$bw

    }
    
    list(y = mx, x = gpoints, bw=hopt, gcv=gcv,kernel=kernel)
}

bootsmooth <- function(y,x, iter=100,conf.level=0.95){
  stopifnot(conf.level>0 ||conf.level<1)
  alpha <- (1-conf.level)/2
  y <- as.matrix(y)
  x <- as.matrix(x)
  
  ry <- nrow(y); cy <- ncol(y)
  rx <- nrow(x); cx <- ncol(x)
  stopifnot(rx==ry)
  ## two-stage bootstrapping
  ##require(bda)
  xbar <- apply(x,1,mean,na.rm=TRUE)
  ybar <- apply(y,1,mean,na.rm=TRUE)
  ybar <- 1-ybar/xbar
  tmp <- npr(y=ybar,x=xbar,sd.x=0,kernel='lp')
  x0 <- tmp$x
  y0 <- tmp$y
  ly <- NULL
  for(i in 1:iter){
      selexr <- sample(1:rx,replace=TRUE)
      seleyr <- sample(1:ry,replace=TRUE)
      selexc <- sample(1:cx,size=rx,replace=TRUE)
      seleyc <- sample(1:cy,size=ry,replace=TRUE)
      y1 <- y[cbind(seleyr,seleyc)]
      x1 <- x[cbind(selexr,selexc)]
      y1 <- 1-y1/x1
      tmp <- npr(y=y1,x=x1,sd.x=0,kernel='lp',x0=x0)
      ly <- cbind(ly,tmp$y)
  }
  qy <- apply(ly,1,quantile,prob=c(alpha,.50,1-alpha))
  list(x=x0,y=y0,
       ym=as.numeric(qy[2,]),
       ll=as.numeric(qy[1,]),
       ul=as.numeric(qy[3,]))
}
