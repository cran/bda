blinding.est <- function(x, group, guess, type='cpe'){
    if (!(is.numeric(x) && is.null(dim(x))))
        stop("'x' must be a numeric vector, not a matrix or other type.")
    if (any(is.na(x)))
        stop("missing value(s) not allowed in 'x'.")
    if (any(is.na(group)))
        stop("missing value(s) not allowed in 'group'.")
    
    n <- length(x)
    if(length(group) != n)
        stop("'x' and 'group' lengths differ")
    if(length(guess) != n)
        stop("'x' and 'guess' lengths differ")
    
    fg <- as.factor(group)
    if(nlevels(fg) != 2)
        stop("'group' must have two levels")
    
    if (!all(group %in% c(0, 1))) {
        stop("group must be a 0 (ctrl)/1(active) vector (or coercible to 0/1).")
    }

    out <- NULL
    type <- match.arg(
        tolower(type), c("cpe","change-point","changepoint","naive",
                         "classic","simple", "adjusted","unblinding"))
    if(type=="simple"||type=="naive"||type=="classic"){
        lm0 <- lm(x ~ group)
        out$estimate <- summary(lm0)$coef
        ## bootstrapping method not applicable as the estimate does
        ## not depend on 'guess'
    }else if(type=="adjusted"||type=="unblinding"){
        lm0 <- lm(x ~ group + guess)
        out$estimate <- summary(lm0)$coef
        
        mydata = data.frame(y = x, group = group, guess = guess)
        z = matrix(1,nrow=200,ncol=1)
        out0 = apply(z,1,FUN=.subboot,db=mydata)
        res2 = quantile(out0,probs=c(0.0125,0.025,0.5,0.975,0.9875))
        out$boot <- res2

    }else{
        ## assuming "NA" or 2 will be equally likely to be 0/1
        sele <- is.na(guess)
        mydata = data.frame(y=x, group = group, guess = guess)
        if(any(sele)){
            z = matrix(1, nrow=100, ncol=1)
            out2 = apply(z,1,FUN=.subcpe2,db=mydata)
            ##if(nrow(out2) == 3){
            res1 <- .subRubinRule(thetai=out2[1,], ui=out2[2,])
            mu2 <- mean(out2[3,])
            se2 <- sd(out2[3,])
            ME2 <- 1.96 * se2
            res2 <- c(mu2, se2, mu2 - ME2, mu2 + ME2)
            ##}else{
            ##res1 <- .subRubinRule(thetai=out2[1,], ui=out2[2,])
            ##res2 <- .subRubinRule(thetai=out2[3,], ui=out2[4,])
            ##}
            tmp <- rbind(res1, res2)
            out <- as.data.frame(tmp)
            rownames(out) <- c("group","sham.cpe")
            names(out) <- c("Estimate", "Std. Error", "95% LCL (2-Sided)",
                            "95% UCL (2-Sided)")
        }else{
            out <- .subcpe2(0, db=mydata)
        }
    }
    out
}

.subboot <- function(x,db){
    ## db$guess2 <- ave(db$guess, db$group, FUN = function(x) sample(x))
    db$guess2 <- sample(db$guess)
    if(x == 1){ # only this level is used 
        lm0 <- lm(y ~ group + guess2, data=db)
        res <- lm0$coef[2]
    }else{
        lm0 <- lm(y ~ group, data=db)
        res <- lm0$coef[2]
    }
    res
}

.subcpe2 <- function(x,db){
    if(x==1){
        sele <- is.na(db$guess)
        k <- sum(sele)
        guess2 <- db$guess[!sele]
        nl <- nlevels(as.factor(guess2))
        if(nl == 2){
            smp <- c(0:1)
            g2 <- sample(smp, size=k, replace=TRUE)
            db$guess[sele] <- g2
        }else{
            smp <- c(0:2)
            g2 <- sample(smp, size=k, replace=TRUE)
            db$guess[sele] <- g2
            ##out <- .cpe3(db)
            sele <- db$guess == 1
            db$guess[!sele] <- 0
        }
        out2 <- .cpe2(db)
        out <- c(out2$estimate[2,1],out2$estimate[2,2],out2$sham)
    }else{
        out <- .cpe2(db)
    }
    out
}

.cpe2 <- function(db){
    lmout1 <- lm(y ~ group, data = db)
    res1 <- lmout1$coef[[2]]
    lmout2 <- lm(y ~ group + guess, data = db)
    res2 <- lmout2$coef[[2]]
    usham <- max(abs(res1),abs(res2),abs(res1+lmout1$coef[[1]]))

    out <- NULL
    m <- 400
    d <- seq(0,2.5,length=m) * usham
    b1 <- sapply(d, FUN=.cpe2dev, db=db)
    ##out$x <- d
    ##out$y <- b1
    usham2 <- .cpfind(d,b1)[2]
    y2 <- db$y - usham2 * db$guess
    group <- db$group
    lm0 <- lm(y2 ~ group)
    out$estimate <- summary(lm0)$coef
    out$sham.CPE <- usham2
    out
}

.cpe2dev <- function(d,db){
    db$y <- db$y - d * db$guess
    sele <- db$guess == 1
    db$guess[sele] <- 1 
    db$guess[!sele] <- 0 
    options(warn = -1)
    g1 <- glm(guess ~ y, family = binomial(), data=db)
    options(warn = 0)
    p <- fitted(g1)
    ## Clip extreme probabilities to avoid log(0)
    eps <- 1e-12
    p_clipped <- p
    p_clipped[p < eps] <- eps
    p_clipped[p > 1 - eps] <- 1 - eps
    ## Compute stabilized deviance
    dev_stable <- -2 * sum(db$guess * log(p_clipped)
                           + (1 - db$guess) * log(1 - p_clipped))
    dev_stable
}

#cpe1 <- function(y, group, guess) {
#  .Call("cpe1_call", y, group, guess, PACKAGE = "yourpkg")
#}

#cpe2 <- function(y, group, guess, sig) {
#  .Call("cpe2_call", y, group, guess, sig, PACKAGE = "yourpkg")
#}

#.cpe3 <- function(db){
#    guess <- db$guess
#    group <- db$group
#    y <- db$y
#    res1 <- .Fortran(.F_cpe1_call, y=y, group=group, guess=guess)
#    res2 <- .Fortran(.F_cpe2_call, y=y, group=group, guess=guess, sig = res1$sigma)
#    sim1 <- .cpePred(res2$x, res2$y[,1]) # effect size
    ##sim2 <- .cpePred(res2$x, res2$y[,2]) # ctrl
    ##sim3 <- .cpePred(res2$x, res2$y[,3]) # active
#    sim4 <- .cpePred(res2$x, res2$y[,4]) # sham
    ##sim5 <- .cpePred(res2$x, res2$y[,5]) # sigma
#    c(sim1,sim4)
#}

.cpePred <- function(x,y){
    x1 <- x
    x2 <- x1^2
    ##x3 <- x1^3
    ##lm1 <- lm(y~x1)
    lm2 <- summary(lm(y~x1+x2))
    ##lm3 <- lm(y~x1+x2+x3)
    ##tmp <- c(lm1$coef[[1]],lm2$coef[[1]],lm3$coef[[1]])
    c(lm2$coef[1,1], lm2$coef[1,2])
}

.cpfind <- function(x,y){
  n <- length(x)
  xm <- which.max(y)
  k2 = max(30, xm)
  x = x[1:xm]
  y = y[1:xm]

  n <- length(x)
  xm <- which.max(y)
  sele = y <= y[xm]
  if(sum(sele) > 30){
    if(n-xm > 10){
      x = x[1:(xm+10)]
      y = y[1:(xm+10)]
    }else{
      x = x[1:xm]
      y = y[1:xm]
    }
  }else{
    if(n-xm > 20){
      x = x[1:(xm+20)]
      y = y[1:(xm+20)]
    }
  }
  n <- length(x)
  if(length(y) != n) stop("'x' and 'y' length differ")
  if(n<10) stop("too few data points to find changing point")
  ## starting point cannot be < 4
  lm2 <- lm(y~x)
  RSS2 <- sum((lm2$resi)^2)
  k <- max(round(n*0.1),4)
  y2 <- NULL
  x2 <- NULL
  for(i in k:(n-k-1)){
    y0 <- y[1:i]
    x0 <- x[1:i]
    x02 <- x0^2
    lm0 <- lm(y0~x0+x02)
    RSS0 <- sum((lm0$resi)^2)
    y1 <- y[(i+1):n]
    x1 <- x[(1+i):n]
    x12 <- x1^2
    lm1 <- lm(y1~x1+x12)
    RSS1 <- sum((lm1$resi)^2)
    RSS <- RSS0 + RSS1
    x2 <- c(x2,i)
    y2 <- c(y2,RSS)
  }
  xid <- which.min(y2)[1]
  
  c(x2[xid],x[xid])
}
