##  2013/04/07: -- This file contails some basic functions for
##  survival data analysis.


##  0 or FALSE for censored times, 1 or TRUE for events

Survival <- function(x,censor,method="Nelson",from, to,gridsize){

  stopifnot(is.numeric(x))
  if(any(is.na(x))) stop("Missing value(s) in 'x'")
  n <- length(x)
  xname <- deparse(substitute(x))

  if(missing(censor)) censor <- rep(1, n)
  if(is.logical(censor)) censor <- censor *1
  else if(is.numeric(censor)){
    if(any(censor!=0 & censor != 1))
      stop("Invalid censoring status in 'censor'")
  }else stop("Invalid value(s) in 'censor'")

  if(missing(method)) method <- "Nelson"
  method <- match.arg(tolower(method),
                      c('km', 'pl','product','kaplan',
                        'na','aalen','nelson','wkde',
                        'mae.exp','mae.weibull','weibull'))
  res <- switch(method,
                km = , pl = , product = ,
                kaplan = .fun.S.KM(x,censor),
                na = , aalen = , nelson = .fun.S.NA(x,censor),
                wkde = .fun.S.wkde(x,censor,from=from, to=to,gridsize=gridsize),
                weibull = .fun.S.weibull(x,censor,from=from, to=to,gridsize=gridsize),
                mae.exp = .fun.S.emae(x,censor,from=from, to=to,gridsize=gridsize),
                mae.weibull = .fun.S.wmae(x,censor,from=from, to=to,gridsize=gridsize),
                stop(paste("Method '", method, "'is not supported"))
                )
  
  structure(list(y = res$y, x = res$x,
                 pars = res$pars,
                 var = res$var,
                 method = res$method,
                 n = res$n, nevent = res$nevent,
                 ncensor = res$ncensor,
                 totMass = res$totMass,
                 data.name = xname,
                 type="S(t)"),
            class='surv')
}

print.surv <- function (x, digits = NULL, ...) 
{
  cat("\nData: ", x$data.name, 
      " (", x$n, " obs., ",
      x$nevent, " events, ",
      x$ncensor, " censored)\n", sep = "")
  cat("Type: y = '", x$type, "'\t method = '", x$method,
      "'\t para. = [ ", sep = "")
  if(is.null(x$pars)){
    cat("NA ]\n",sep='')
  }else{
    k <- length(x$pars)
    if(k==1){
      cat(x$pars, " ]\n",sep='')
    }else{
      for(i in 1:(k-1)){
        cat(x$pars[i], ", ",sep='')
      }
      cat(x$pars[k], " ]\n", sep='')
    }
  }
  cat("Total mass = ", x$totMass, "\n", sep='')
  if(length(x$x)<30){
    print(as.data.frame(x[c("x", "y")]),
          digits = digits, ...)
  }else{
    print(summary(as.data.frame(x[c("x", "y")])),
          digits = digits, ...)
  }
  invisible(x)
}

.fun.S.KM <- function(x,w){
  n <- length(x)
  x.censor <- w==0 # 0 censored, 1 events
  x.events <- x[!x.censor]
  n.events <- tapply(w[!x.censor], x.events, length)
  t.events <- as.numeric(names(n.events))
  n.events <- as.numeric(n.events)
  n.atrisk <-   apply(as.matrix(t.events,ncol=1),1,.fun.KM.count,y=x)
  y <- cumprod(1-n.events/n.atrisk)
  vy <- cumsum(n.events/n.atrisk/(n.atrisk-n.events))
  y[n.atrisk <= 0] <- 0.0
  totMass <- 1-((rev(n.atrisk - n.events))[1])/length(x)

  ox <- order(x); wo <- w[ox]
  wsum <- cumsum(rev(wo))
  xpos <- match(0, rev(wsum))
  if(is.na(xpos)) totMass <- 1.0
  else totMass <- (xpos-1)/n
  
  list(x=t.events, y=y, totMass=totMass,
       var = y*y*vy,
       n=n, nevent=sum(w),
       ncensor=n-sum(w), pars=NULL,
       method="Kaplan-Meier")
}

.fun.S.NA <- function(x,w){
  out <- .fun.H.NA(x,w)
  out$y <- exp(-out$y)
  out
}

.fun.KM.count <- function(x,y) sum(y>=x)

.weibull.mle <- function(x,w){
  n <- length(x)
  par1 <- rev(.weibullplot(x)) # bdde.R (reverse order)
  par2 <- .Fortran(.F_WeibullMleNMMIN, as.double(x),
                   as.double(w),as.integer(n),
                   pars = as.double(par1))$pars
  
  khat <- uniroot(.weibullkappa, c(0, 10*par1[1]), x=x, w=w)
  lhat <- (sum(x^khat$root)/sum(w))^(1/khat$root)
  par3 <- c(khat$root,lhat)
  list(WeibullPlot=par1, Nelder=par2, MLE=par3)
}

.weibullkappa <- function(k,x,w)
  1/k+sum(w*log(x))/sum(w) - sum(x^k*log(x))/sum(x^k)

.fun.S.weibull <- function(x,w,from,to,gridsize){
  out <- .fun.S.KM(x=x, w=w)
  pars <- (.weibull.mle(x,w))$MLE
  a <- ifelse(missing(from), 0, from);
  b <- ifelse(missing(to), max(x*w), to);
  if(missing(gridsize)) gridsize <- 512L
  gpoints <- seq(a, b, length = gridsize)
  Sx <- exp(-(gpoints/pars[[2]])^pars[[1]])
  
  list(x=gpoints, y=Sx, totMass=out$totMass,
       n=out$n, nevent=out$nevent,
       ncensor=out$ncensor,
       pars=pars,
       method="Weibull"
       )
}

.fun.S.emae <- function(x, w, from,to,gridsize){
  out <- .fun.S.KM(x=x, w=w)
  if(missing(gridsize)) gridsize <- 512L
  ##  a <- min(x); b <- max(x);
  a <- ifelse(missing(from), 0, from);
  b <- ifelse(missing(to), max(x*w), to);
  gpoints <- seq(a, b, length = gridsize)

  ## sort the data
  ox <- order(x); x <- x[ox]; w <- w[ox]
  n <- length(x)
  y <- .Fortran(.F_expmae,as.double(x),as.double(w),
                as.integer(n), y=as.double(gpoints),
                as.integer(gridsize))$y
  Fx <- cumsum(y)
  Sx <- 1.0 - Fx/max(Fx)* out$totMass
  
  list(x=gpoints, y=Sx, totMass=out$totMass,
       n=out$n, nevent=out$nevent,
       ncensor=out$ncensor,
       pars=NULL,
       method="Weibull(MAE)"
       )
}

.fun.S.wmae <- function(x, w, from,to,gridsize){
  out <- .fun.S.KM(x=x, w=w)
  pars <- (.weibull.mle(x,w))$WeibullPlot
  if(missing(gridsize)) gridsize <- 512L
  ##  a <- min(x); b <- max(x);
  a <- ifelse(missing(from), 0, from);
  b <- ifelse(missing(to), max(x*w), to);
  gpoints <- seq(a, b, length = gridsize)

  ## sort the data
  ox <- order(x); x <- x[ox]; w <- w[ox]
  n <- length(x)
  y <- .Fortran(.F_weibullmae,as.double(x),as.double(w),
                as.integer(n), as.double(pars),
                y=as.double(gpoints),
                as.integer(gridsize))$y
  Fx <- cumsum(y)
  Sx <- 1.0 - Fx/max(Fx)* out$totMass
  
  list(x=gpoints, y=Sx, totMass=out$totMass,
       n=out$n, nevent=out$nevent,
       ncensor=out$ncensor,
       pars=NULL,
       method="Weibull(MAE)"
       )
}

.fun.S.wkde <- function(x,w,from,to,gridsize){
  out <- .fun.S.KM(x=x, w=w)
  x1 <- out$x; y <- out$y
  w1 <- diff(c(0,1-y)); 
  a <- ifelse(missing(from), 0, from);
  b <- ifelse(missing(to), max(x*w), to);
  res <- wkde(x=x1, weights=w1, from=a, to=b, gridsize=gridsize)
  Fx <- cumsum(res$y)
  Sx <- 1.0 - Fx/max(Fx)* out$totMass

  list(x=res$x, y=Sx, totMass=out$totMass,
       n=out$n, nevent=out$nevent,
       ncensor=out$ncensor,
       pars=res$bw,
       method="wkde"
       )
}


Hazard <- function(x,censor,method="Nelson",from, to,gridsize){

  stopifnot(is.numeric(x))
  if(any(is.na(x))) stop("Missing value(s) in 'x'")
  n <- length(x)
  xname <- deparse(substitute(x))

  if(missing(censor)) censor <- rep(1, n)
  if(is.logical(censor)) censor <- censor *1
  else if(is.numeric(censor)){
    if(any(censor!=0 & censor != 1))
      stop("Invalid censoring status in 'censor'")
  }else stop("Invalid value(s) in 'censor'")

  if(missing(method)) method <- "Nelson"
  method <- match.arg(tolower(method),
                      c('km', 'pl','product','kaplan',
                        'na','aalen','nelson','wkde',
                        'mae.exp','mae.weibull','weibull'))
  res <- switch(method,
                km = , pl = , product = ,
                kaplan = .fun.H.KM(x,censor),
                na = , aalen = , nelson = .fun.H.NA(x,censor),
                wkde = .fun.H.wkde(x,censor,from=from, to=to,gridsize=gridsize),
                weibull = .fun.H.weibull(x,censor,from=from, to=to,gridsize=gridsize),
                mae.exp = .fun.H.emae(x,censor,from=from, to=to,gridsize=gridsize),
                mae.weibull = .fun.H.wmae(x,censor,from=from, to=to,gridsize=gridsize),
                stop(paste("Method '", method, "'is not supported"))
                )
  
  structure(list(y = res$y, x = res$x,
                 pars = res$pars,
                 var = res$var,
                 method = res$method,
                 n = res$n, nevent = res$nevent,
                 ncensor = res$ncensor,
                 totMass = res$totMass,
                 data.name = xname,
                 type="H(t)"),
            class='surv')
}


.fun.H.KM <- function(x,w){
  out <- .fun.S.KM(x,w)
  out$y <- -log(out$y)
  out
}

.fun.H.NA <- function(x,w){
  n <- length(x)
  x.censor <- w==0 # 0 censored, 1 events
  x.events <- x[!x.censor]
  n.events <- tapply(w[!x.censor], x.events, length)
  t.events <- as.numeric(names(n.events))
  n.events <- as.numeric(n.events)
  n.atrisk <-   apply(as.matrix(t.events,ncol=1),1,.fun.KM.count,y=x)
  y <- cumsum(n.events/n.atrisk)
  vy <- cumsum(n.events/(n.atrisk*n.atrisk))
  totMass <- 1-((rev(n.atrisk - n.events))[1])/length(x)

  ox <- order(x); wo <- w[ox]
  wsum <- cumsum(rev(wo))
  xpos <- match(0, rev(wsum))
  if(is.na(xpos)) totMass <- 1.0
  else totMass <- (xpos-1)/n

  list(x=t.events,y=y,totMass=totMass,
       var = vy,
       n=n, nevent=sum(w),
       ncensor=n-sum(w),
       pars=NULL,
       method="Nelson-Aalen"
       )
}

.fun.H.weibull <- function(x,w,from,to,gridsize){
  out <- .fun.S.KM(x=x, w=w)
  pars <- (.weibull.mle(x,w))$MLE
  a <- ifelse(missing(from), 0, from);
  b <- ifelse(missing(to), max(x*w), to);
  if(missing(gridsize)) gridsize <- 512L
  gpoints <- seq(a, b, length = gridsize)
  Sx <- (gpoints/pars[[2]])^pars[[1]]
  
  list(x=gpoints, y=Sx, totMass=out$totMass,
       n=out$n, nevent=out$nevent,
       ncensor=out$ncensor,
       pars=pars,
       method="Weibull"
       )
}

.fun.H.emae <- function(x, w, from,to,gridsize){
  out <- .fun.S.emae(x,w,from=from,to=to,gridsize=gridsize)
  out$y <- -log(out$y)
  out
}

.fun.H.wmae <- function(x, w, from,to,gridsize){
  out <- .fun.S.wmae(x,w,from=from,to=to,gridsize=gridsize)
  out$y <- -log(out$y)
  out
}

.fun.H.wkde <- function(x,w,from,to,gridsize){
  out <- .fun.S.wkde(x,w,from=from,to=to,gridsize=gridsize)
  out$y <- -log(out$y)
  out
}


hazard <- function(x,censor,method="Nelson",from, to,gridsize){

  stopifnot(is.numeric(x))
  if(any(is.na(x))) stop("Missing value(s) in 'x'")
  n <- length(x)
  xname <- deparse(substitute(x))

  if(missing(censor)) censor <- rep(1, n)
  if(is.logical(censor)) censor <- censor *1
  else if(is.numeric(censor)){
    if(any(censor!=0 & censor != 1))
      stop("Invalid censoring status in 'censor'")
  }else stop("Invalid value(s) in 'censor'")

  if(missing(method)) method <- "Nelson"
  method <- match.arg(tolower(method),
                      c('km', 'pl','product','kaplan',
                        'na','aalen','nelson','wkde',
                        'lp','localpolynomial','smoothing',
                        'mae.exp','mae.weibull','weibull'))
  res <- switch(method,
                km = , pl = , product = ,
                kaplan = .fun.h.KM(x,censor),
                na = , aalen = ,
                nelson = .fun.h.NA(x,censor,from=from, to=to,
                  gridsize=gridsize,kernel="Epan"),
                smoothing = .fun.h.NA(x,censor,from=from, to=to,
                  gridsize=gridsize,kernel="Biweight"),
                wkde = .fun.h.wkde(x,censor,from=from, to=to,
                  gridsize=gridsize),
                weibull = .fun.h.weibull(x,censor,from=from,
                  to=to,gridsize=gridsize),
                mae.exp = .fun.h.emae(x,censor,from=from,
                  to=to,gridsize=gridsize),
                mae.weibull = .fun.h.wmae(x,censor,from=from,
                  to=to,gridsize=gridsize),
                localpolynomial = ,
                lp = .fun.h.lpr(x,censor,from=from, to=to,
                  gridsize=gridsize),
                stop(paste("Method '", method, "'is not supported"))
                )
  
  structure(list(y = res$y, x = res$x,
                 var = res$var,
                 pars = res$pars,
                 method = res$method,
                 n = res$n, nevent = res$nevent,
                 ncensor = res$ncensor,
                 totMass = res$totMass,
                 data.name = xname,
                 type="h(t)"),
            class='surv')
}

.fun.h.lpr <- function(x, w, from,to,gridsize){
  out <- .fun.H.NA(x,w)
  y <- diff(c(0, out$y)) # delta(H(t))
  y.finite <- is.finite(y)
  res <- npr(y=y[y.finite], x=out$x[y.finite],
             from=from, to=to, gridsize = gridsize)

  list(x=res$x, y=res$y, totMass=out$totMass,
       n=out$n, nevent=out$nevent,
       ncensor=out$ncensor,
       pars=res$bw,
       method="Nelson-Aalen(lp1)"
       )
}

.fun.h.smooth <- function(x, w, from, to, gridsize){
  out <- .fun.H.NA(x,w)
  y <- diff(c(0, out$y)) # delta(H(t))
  y.finite <- is.finite(y)
  res <- npr(y=y[y.finite], x=out$x[y.finite],
             from=from,to=to,gridsize = gridsize)

  list(x=res$x, y=res$y, totMass=out$totMass,
       n=out$n, nevent=out$nevent,
       ncensor=out$ncensor,
       pars=res$bw,
       method="Nelson-Aalen(lpsmmoth)"
       )
}

.fun.h.KM <- function(x,w){
  out <- .fun.S.KM(x,w)
  Sx0 <- out$y
  Sx1 <- c(1, rev((rev(Sx0))[-1]))
  out$y <- 1-Sx0/Sx1
  out
}

.fun.h.NA <- function(x,w,from,to,gridsize,kernel){
  if(missing(from)&&missing(to)&&missing(gridsize)&&
     missing(kernel)){
    out <- .fun.S.NA(x,w)
    Sx0 <- out$y
    Sx1 <- c(1, rev((rev(Sx0))[-1]))
    out$y <- 1-Sx0/Sx1
    res <- out
  }else{
    out <- .fun.H.NA(x,w)
    Hx <- out$y; tx <- out$x; vH <- out$var
    dHx <- diff(c(0,Hx));
    dvH <- diff(c(0, vH));
    out2 <- .smhazard(dHx,tx,dvH,kernel=kernel,
                      from=from,to=to,gridsize=gridsize)
    meth <- paste("Nelson(", out2$type, ")",sep='')
    res <- list(x=out2$x, y=out2$y,
                totMass=out$totMass,
                n=out$n, nevent=out$nevent,
                ncensor=out$ncensor,
                pars = out2$bw,
                method = meth
                )
  }
  res
}

.fun.h.weibull <- function(x,w,from,to,gridsize){
  out <- .fun.S.KM(x=x, w=w)
  pars <- (.weibull.mle(x,w))$MLE
  a <- ifelse(missing(from), 0, from);
  b <- ifelse(missing(to), max(x*w), to);
  if(missing(gridsize)) gridsize <- 512L
  gpoints <- seq(a, b, length = gridsize)
  Sx <- pars[[1]]/pars[[2]]*(gpoints/pars[[2]])^(pars[[1]]-1)
  
  list(x=gpoints, y=Sx, totMass=out$totMass,
       n=out$n, nevent=out$nevent,
       ncensor=out$ncensor,
       pars=pars,
       method="Weibull"
       )
}

.fun.h.emae <- function(x, w, from,to,gridsize){
  out <- .fun.S.KM(x=x, w=w)
  if(missing(gridsize)) gridsize <- 512L
  ##  a <- min(x); b <- max(x);
  a <- ifelse(missing(from), 0, from);
  b <- ifelse(missing(to), max(x*w), to);
  delta <- (b-a)/(gridsize - 1)
  gpoints <- seq(a, b, length = gridsize)

  ## sort the data
  ox <- order(x); x <- x[ox]; w <- w[ox]
  n <- length(x)
  y <- .Fortran(.F_expmae,as.double(x),as.double(w),
                as.integer(n), y=as.double(gpoints),
                as.integer(gridsize))$y
  Fx <- cumsum(y)*delta*out$totMass
  Sx <- 1 - Fx
  
  list(x=gpoints, y=out$totMass * y/Sx, totMass=out$totMass,
       n=out$n, nevent=out$nevent,
       ncensor=out$ncensor,
       pars=NULL,
       method="Weibull(MAE)"
       )
}

.fun.h.wmae <- function(x, w, from,to,gridsize){
  out <- .fun.S.KM(x=x, w=w)
  pars <- (.weibull.mle(x,w))$WeibullPlot
  if(missing(gridsize)) gridsize <- 512L
  ##  a <- min(x); b <- max(x);
  a <- ifelse(missing(from), 0, from);
  b <- ifelse(missing(to), max(x*w), to);
  delta <- (b-a)/(gridsize - 1)
  gpoints <- seq(a, b, length = gridsize)

  ## sort the data
  ox <- order(x); x <- x[ox]; w <- w[ox]
  n <- length(x)
  y <- .Fortran(.F_weibullmae,as.double(x),as.double(w),
                as.integer(n), as.double(pars),
                y=as.double(gpoints),
                as.integer(gridsize))$y
  Fx <- cumsum(y)*delta*out$totMass
  Sx <- 1 - Fx
  
  list(x=gpoints, y=out$totMass * y/Sx, totMass=out$totMass,
       n=out$n, nevent=out$nevent,
       ncensor=out$ncensor,
       pars=NULL,
       method="Weibull(MAE)"
       )
}

.fun.h.wkde <- function(x,w,from,to,gridsize){
  out <- .fun.S.KM(x=x, w=w)
  x1 <- out$x; y <- out$y
  w1 <- diff(c(0,1-y)); 
  res <- wkde(x=x1, weights=w1, from=from,to=to, gridsize=gridsize)
  x0 <- res$x; y0 <- res$y
  delta <- diff(range(x0))/(length(x0)-1)
  Fx <- cumsum(y0)*delta*out$totMass
  Sx <- 1 - Fx

  list(x=x0, y=out$totMass*y0/Sx, totMass=out$totMass,
       n=out$n, nevent=out$nevent,
       ncensor=out$ncensor,
       pars=res$bw,
       method="wkde"
       )
}

.smhazard <- function(y,x,v,bw,kernel='epan',from,to,gridsize){
  if(missing(from)) from <- min(x)
  if(missing(to)) to <- max(x)
  stopifnot(to > from)
  stopifnot(from >= 0)
    
  kernel <- match.arg(tolower(kernel),
                      c('epanechnikov','biweight','uniform'))
  ikernel <- switch(kernel, epanechnikov=0, biweight=1, uniform=2)

  if(missing(bw)){
    lscv <- 1
    bw <- bw.nrd(x)
  }else{
    stopifnot(bw>0)
    lscv <- 0
  }
  if(missing(gridsize)) gridsize <- 512L
  stopifnot(gridsize>10)
  gpoints <- seq(from, to, length=gridsize)
  n <- length(x)
  stopifnot(length(y) == n)
  stopifnot(length(v) == n)
  if(any(is.na(x)|is.na(y)|is.na(v)))
    stop("Missing value(s) in 'x', 'y', or 'v'")
  x.infinite <- !is.finite(x)|!is.finite(y)|!is.finite(v)
  if(any(x.infinite))
    warning("Inifite value(s) in 'x' and/or 'y'")
  if(mean(x.infinite)==1){
    stop("Invalid data")
  }else{
    x <- x[!x.infinite]
    y <- y[!x.infinite]
    v <- v[!x.infinite]
  }
  
  out <- .Fortran(.F_lpshazard,
                  fx = as.double(gpoints),
                  var = as.double(gpoints),
                  as.integer(gridsize),
                  as.double(x), as.double(y),
                  as.double(v), as.integer(n),
                  bw = as.double(bw), as.integer(lscv),
                  as.integer(ikernel))

  structure(list(y = out$fx, x = gpoints,
                 bw=out$bw,
                 ucb=out$var, lcb=out$var,
                 type = kernel,
                 call = match.call()
                 ),
            class = 'smooth')
}

