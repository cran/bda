pro.test <- function(x,y,group,cutoff,class.names,conf.level,x.range,delta){
  options(warn=-2)
  stopifnot(length(x) == length(y))
  stopifnot(length(x) == length(group))
  sele <- is.na(x) | is.na(y)
  if(any(sele)){
    x <- x[!sele]
    y <- y[!sele]
    group <- group[!sele]
  }
  group <- as.factor(as.character(group))

  if(is.numeric(x) && is.numeric(y)){
    out <- .mytestnum(x=x,y=y,group=group,cutoff=cutoff,
                      class.names=class.names,
                      conf.level=conf.level,
                      x.range=x.range,
                      delta=delta)    
  }else{
    x <- as.factor(as.character(x))
    y <- as.factor(as.character(y))
    out <- .mytestchar(x=x,y=y,group=group)
  } 

  options(warn=1)
  out
}

.mytestchar <- function(x,y,group){
  if(nlevels(x)>5 || nlevels(y)>5)
    stop("too many states/levels")
  Outcome <- table(x,y,group)
  out.tbl <- ftable(Outcome)
  out.CMH <- mantelhaen.test(Outcome)
  out <- list(Results=out.tbl, Test=out.CMH)
}

.mytestnum <- function(x,y,group,cutoff,class.names,conf.level,x.range,delta){
  
  if(any(x<0 | y<0))
    stop("negative values found in 'x' or 'y'")

  ## Absolute changes ######################################
  dx <- x - y
  sele <- is.na(dx) | is.na(group)
  dx <- dx[!sele]
  g <- group[!sele]
  ## t-test
  tmp1 <- tapply(dx, g,mean,na.rm=TRUE)
  tmp2 <- t.test(dx~g)$p.value
  tout <- c(tmp1, tmp2,NA)
  ## Mann-Whitney U-test
  tmp1 <- tapply(dx, g,median,na.rm=TRUE)
  tmp2 <- wilcox.test(dx~g)  
  tout <- rbind(tout, c(tmp1,tmp2$p.value,NA))
  ## KS-test
  tmp1 <- .subedf(dx, g)
  tout <- rbind(tout, tmp1)
  ## Relative changes ######################################
  dx <- 1 - y/x
  sele <- is.na(dx) | is.na(group)
  dx <- dx[!sele]
  g <- group[!sele]
  ## t-test
  tmp1 <- tapply(dx, g,mean,na.rm=TRUE)
  tmp2 <- t.test(dx~g)$p.value
  tout <- rbind(tout, c(tmp1, tmp2,NA))
  ## Mann-Whitney U-test
  tmp1 <- tapply(dx, g,median,na.rm=TRUE)
  tmp2 <- wilcox.test(dx~g)  
  tout <- rbind(tout, c(tmp1,tmp2$p.value,NA))
  ## KS-test
  tmp1 <- .subedf(dx, g)
  tout <- rbind(tout, tmp1)
  tout <- as.data.frame(tout)
  names(tout)[3] <- "p-value"
  names(tout)[4] <- "Ref"
  row.names(tout) <- c("t-test (abs)","U-test (abs)","KS-test (abs)",
                       "t-test (rel)", "U-test (rel)","KS-test (rel)")
  #print(tout)
  
  
  out <- .subMyTest(x=x,y=y,group=group,
                    cutoff=cutoff,class.names=class.names)
  
  if(!missing(conf.level)){
    stopifnot(is.numeric(conf.level))  
    stopifnot(length(conf.level) == 1)
    stopifnot(conf.level<1 && conf.level > 0)
    alpha <- min(conf.level,1-conf.level)
    tmp <- .subMyCI(x=x,y=y,group=group,cutoff=cutoff,
                    class.names=class.names,
                    x.range=x.range,delta=delta,
                    alpha=alpha)    
    
      
    out.test <- data.frame(estimate=c(out$Test$p.value,out$RA),
                           t(tmp))
    names(out.test) <- c("estimate", "median", "lcl","ucl")
    rownames(out.test)[1] <- "CMH.test" 
    rownames(out.test)[4] <- "2p.test"
    rownames(out.test)[2] <- paste("resp.rate",rownames(out.test)[2])
    rownames(out.test)[3] <- paste("resp.rate",rownames(out.test)[3])
    out <- list(Multi.State.Analysis=out$Results,conf.level=1-alpha,Test=out.test)
    
  }
  
  tmp <- paste("x <=",cutoff)
  tmp[-1] <- paste(cutoff[-1],"<",tmp[-1])
  tmp <- c(tmp,paste(">", rev(cutoff)[1]))
  
  if(missing(class.names)){
    class.names <- paste("state-",1:length(tmp),sep='')
  }else{
    stopifnot(length(class.names) == length(tmp))    
  }
  
  State.tbl <- data.frame(state=class.names, range=tmp)
  ##print(State.tbl)
  out$states <- State.tbl
  out$Responder.Analysis <- tout
  out
}

.subedf <- function(x, group){
  x0 <- unique(sort(x))
  n <- length(x0)
  g <- as.factor(as.character(group))
  sele <- g==levels(g)[1]
  x1 <- x[sele]; n1 <- length(x1)
  x2 <- x[!sele]; n2 <- length(x2)
  pv <- ks.test(x1,x2)$p.value
  Fx1 <- Fx2 <- rep(0,n)
  for(i in 1:n){
    Fx1[i] <- sum(x1<x0[i])
    Fx2[i] <- sum(x2<x0[i])
  }
  Fx1 <- Fx1/n1
  Fx2 <- Fx2/n2
  dFx <- abs(Fx1-Fx2)
  isele <- which(dFx==max(dFx))
  if(length(isele)>1){
    warning("maximum differences reached at multiple levels")
    #print(x0[isele])
    isele <- isele[1]    
  }
  c(Fx1[isele],Fx2[isele],pv,x0[isele])
}

.subMyCI <- function(x,y,group,cutoff,class.names,x.range,delta,alpha=0.05){
  if(is.character((x))) 
    x0 <- as.numeric(x)
  else if(is.factor(x))
    x0 <- as.numeric(as.character(x))
  else
    x0 <- x
  
  if(is.character((y))) 
    y0 <- as.numeric(y)
  else if(is.factor(y))
    y0 <- as.numeric(as.character(y))
  else y0 <- y
  
  x1 <- c(x0,y0)
  xmin <- min(x1, na.rm=TRUE)
  xmax <- max(x1, na.rm=TRUE)
  if(missing(x.range)){
    dx <- min(diff(c(xmin,sort(cutoff),xmax)))/4
  }else{
    stopifnot(length(x.range) == 2)
    stopifnot(x.range[1] < x.range[2])
    stopifnot(x.range[1] <= xmin)
    stopifnot(x.range[2] >= xmax)
    dx <- min(diff(c(x.range[1],sort(cutoff),x.range[2])))/4
  }
  
  if(missing(delta)){
    d0 <- dx
  }else{
    stopifnot(is.numeric(delta))
    stopifnot(length(delta) == 1)
    if(abs(delta) > dx)
      warning("'delta' might be too large")
    d0 <- abs(delta)
  }
  re <- d0*c(1,.5,0,-.5,-1)
  ## consider limitted cases with length of cutoff no more than 4 
  ## (up to five states)

  l <- length(cutoff)
  if(l==1)
    tmp <- matrix(re, ncol=1)
  if(l==2)
    tmp <- expand.grid(re,re)
  else if(l==3)
    tmp <- expand.grid(re,re,re)
  else if(l==4)
    tmp <- expand.grid(re,re,re,re)
  else
    stop("too many states to compute the CIs")
  
  ctof = sort(cutoff)
  tmp1 <- matrix(ctof, ncol=l, nrow=nrow(tmp),byrow = TRUE)
  cutoffs <- tmp + tmp1
  #print(cutoffs)
  res <- apply(cutoffs,1,.mytest,
               x=x,y=y,group=group,
               class.names=class.names)
  tmp <- apply(res,1,quantile,
               probs=c(0.5,alpha/2,1-alpha/2),
               na.rm=TRUE)
}
  
.mytest <- function(cutoff,x,y,group,class.names){
  out <- .mytestnum(x=x,y=y,group=group,cutoff=cutoff,class.names=class.names)
  c(CMH.test=out$Test$p.value,out$RA)
}


.subMyTest <- function(x,y,group,cutoff,class.names){  
  if(missing(cutoff)){
    cutoff <- median(x)
    if(missing(class.names)){
      class.names <- c("1-low","2-high")    
    }else{
      if(length(class.names) ==2)
        class.names <- paste(1:2,"-", class.names, sep='')
      else
        class.names <- c("1-low","2-high")    
    }
  }else{
    if(any(diff(cutoff)<=0))
      stop("'cutoff' needs to be sorted in increasing order")
    xy <- c(x,y)
    if(min(xy) >= min(cutoff) || max(xy) <= max(cutoff))
      stop("invalid cutoff values")
    k <- length(cutoff)
    if(missing(class.names)){
      class.names <- paste(1:(k+1),"-level", sep='')
    }else{
      if(length(class.names) == k+1){
        class.names <- paste(1:(k+1),"-", class.names, sep='')
      }
      else{
        class.names <- paste(1:(k+1),"-level", sep='')
        warning("incorrect class names")        
      }
    }
  }
  
  .labelclass <- function(x,x0,z){
    k <- length(x0)
    tmp <- rep(z[1],length(x))
    for(i in 1:(k-1)){
      sele <- x>x0[i] & x<=x0[i+1]
      tmp[sele] <- z[i+1]
    }
    sele <- x>x0[k]
    tmp[sele] <- z[k+1]
    tmp
  }
  
  Before <- .labelclass(x,cutoff,class.names)
  After <- .labelclass(y,cutoff,class.names)
  
  ## Show summaries by arm/group
  #require(gmodels)
  #out.tbl2 <- NULL
  #for(i in levels(group)){
  #  sele <- group == i
    ##print(CrossTable(Before[sele], After[sele]))
    #out.tbl2 <- CrossTable(Before[sele], After[sele])
  #}
  
  ## print 3-way table for comparisons
  Outcome <- table(Before,After,group)
  ## print(ftable(Outcome))
  
  out.tbl <- ftable(Outcome)
  ## Cochran-Mantel-Haenszel chi-squared test of the null hypothesis 
  ## that two nominal variables are conditionally independent in each 
  ## stratum, assuming that there is no three-way interaction. 
  ## x is a 3 dimensional contingency table, where the last dimension 
  ## refers to the strata.
  
  out.CMH <- mantelhaen.test(Outcome)
  ##print(out.CMH)

  out <- list(Results=out.tbl, Test=out.CMH)
  
  ## Loglinear Models: https://www.statmethods.net/stats/frequencies.html
  #mydat <- data.frame(group=group,Before=Before,After=After)
  #tmp <- xtabs(~Before+After+group, data=mydat)
  #print(tmp)
  #out.loglin <- loglin(~Before+After+group, tmp) #mutual independend
  #print(out.loglin)
  #group is independent of A*B  
  #out.loglin2 <- loglin(~group+Before+After+Before*After, tmp) 
  #print(out.loglin2)
  
  ## if there are only two arms
  ndim <- dim(Outcome)
  if(ndim[3]==2){
    gnames <- colnames(Outcome[1,,])
    xcount <- rep(0,ndim[3])
    n1 <- sum(Outcome[,,1])
    n2 <- sum(Outcome[,,2])
    for(k in 1:ndim[3]){
      inames <- rownames(Outcome[,,k])
      jnames <- colnames(Outcome[,,k])
      xt <- Outcome[,,k]
      for(i in 1:ndim[1]){
        ilvl <- as.numeric(substr(inames[i],1,1))
        for(j in 1:ndim[2]){
          jlvl <- as.numeric(substr(jnames[j],1,1))
          if(ilvl > jlvl) xcount[k] <- xcount[k] + xt[i,j]
        }
      }
    }
    
    #cat("\nResponder rate:")
    #cat("\n  ", gnames[1], ":", xcount[1], "/", n1,"=",
    #    round(xcount[1]/n1*100,2),"%")
    #cat("\n  ", gnames[2], ":", xcount[2], "/", n2,"=",
    #    round(xcount[2]/n2*100,2),"%")
    ptest <- prop.test(x=xcount,n=c(n1,n2))
    #cat("\n   p-value = ", ptest$p.value)
    #cat("\n\t", ptest$method,"\n")
    
    out.responder <- c(Rate1=xcount[1]/n1,Rate2=xcount[2]/n2,p.value=ptest$p.value)
    names(out.responder) <- c(gnames[1], gnames[2],"p-value")
    out <- list(Results=out.tbl, Test=out.CMH, RA=out.responder)
    
  }
  
  invisible(out)
}
