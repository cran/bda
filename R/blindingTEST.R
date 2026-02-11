

#print.ctest <- function(x,...){
#    cat("\nData summary:\n")
#    print(x$BI.test$data)
#    cat("\nFisher randomization test:\n")
#    print(round(x$BI.test$BI,3))
#    cat("\nLatent-shift test of contamination:\n  p-value = ",
#        x$p.value,"\n\n")
#}

## the followiing function summarize the blinding responses and
## compute bootstrapping CI of the BIs when applicable: (1) with
## missing responses, or (2) with IDK responses.

## make 'iter' internal as number of imputation datasets. Even if the
## number of missing values is not large, say n<3 or p<5%, we can keep
## iter=100 for now. This can be updated later to permutation test (to
## increase speed and accuracy).

## combine the results using Rubin's method.


blinding.BI <- function(group, guess){
    out <- NULL
    iter <- 100 #this might be too large, can be reduced to 50
    
    if (!all(group %in% c(0, 1))) {
        stop("group must be a 0 (ctrl)/1(active) vector (or coercible to 0/1).")
    }
    
    n <- length(group)
    if(length(guess) != n)
        stop("'group' and 'guess' must have the same length.")

    sele <- is.na(guess)
    k <- sum(sele)
    if(k > 0){
        xtbl <- .BIEst3tbl(group=group,guess=guess)
        out$tbl <- xtbl
        mydata = data.frame(group = group, guess = guess)
        z = matrix(1,nrow=iter,ncol=1)
        out2 = apply(z,1,FUN=.subptbBI3,db=mydata)

        res1 <- .subRubinRule(thetai=out2[1,], ui=out2[2,])
        ## JBL range [0,+1]
        res1[3] <- max(res1[3],0)
        res1[4] <- min(res1[4],1)
        
        res2c <- .subRubinRule(thetai=out2[3,], ui=out2[4,])
        res2t <- .subRubinRule(thetai=out2[5,], ui=out2[6,])
        ## BBL range [-1,+1]
        res2c[3] <- max(res2c[3],-1)
        res2c[4] <- min(res2c[4],1)
        res2t[3] <- max(res2t[3],-1)
        res2t[4] <- min(res2t[4],1)

        tmp <- rbind(NULL,res1)
        res0 <- as.data.frame(tmp)
        rownames(res0) <- "Overall"
        names(res0) <- c("Estimate", "Std. Error", "95% LCL (2-Sided)",
                         "95% UCL (2-Sided)")
        JBI <- res0
        
        tmp <- rbind(res2t, res2c)
        res0 <- as.data.frame(tmp)
        rownames(res0) <- c("Treatment", "Control")
        names(res0) <- c("Estimate", "Std. Error", "95% LCL (2-Sided)",
                         "95% UCL (2-Sided)")
        BBI <- res0
        tmp <- list(JamesBI = JBI, BangBI = BBI)
        out$BI <- tmp
    }else{
        out <- .BIEst3(group=group,guess=guess)
    }
    out
}

.BIEst3tbl <- function(group,guess){
    sele <- is.na(guess)
    xtbl0 <- table(sele, group)
    group2 <- group[!sele];#group2
    guess2 <- guess[!sele];#guess2
    if (all(guess2 %in% c(0, 1))) {
        arm <- group2
        arm[arm==0] = "B:Placebo"
        arm[arm==1] = "A:Treatment"
        blind = guess2
        blind[blind!=1] <- "2:Placebo"
        blind[blind==1] <- "1:Treatment"
        xtbl1 <- table(blind, arm)
    }else{
        arm = group2; #arm
        arm[arm==0] = "B:Placebo"
        arm[arm==1] = "A:Treatment"
        blind = rep("3:DNK",length(guess2))
        blind[guess2==0] <- "2:Placebo"
        blind[guess2==1] <- "1:Treatment"
        xtbl1 = table(blind, arm)
    }
    xtbl <- rbind(xtbl1, "NA"=xtbl0[2,])
}

.BIEst3 <- function(group,guess){
    if (all(guess %in% c(0, 1))) {
        arm <- group
        arm[arm==0] = "B:Placebo"
        arm[arm==1] = "A:Treatment"
        blind = guess
        blind[blind!=1] <- "2:Placebo"
        blind[blind==1] <- "1:Treatment"
        xtbl <- table(blind, arm)
        xtbl2 <- rbind(xtbl, "3:DNK"=c(0,0))
        res <- BI::BI(xtbl2)
    }else{
        arm <- group
        arm[arm==0] = "B:Placebo"
        arm[arm==1] = "A:Treatment"
        blind = rep("3:DNK",length(guess))
        blind[guess==0] <- "2:Placebo"
        blind[guess==1] <- "1:Treatment"
        xtbl <- table(blind, arm)
        res <- BI::BI(xtbl)
    }
    list(tbl=xtbl,BI=res)
}

.subptbBI3 <- function(x,db){
  if(x!=1) #redundant for now
      stop("imputing method not defined")

  sele <- is.na(db$guess)
  k <- sum(sele)
  gnna <- db$guess[!sele]
  nl <- nlevels(as.factor(gnna))
  guess2 <- sample(c(0:(nl-1)), size=k,replace=TRUE)
  db$guess[sele] <- guess2
  BIres = .BIEst3(group=db$group,guess=db$guess)
  ##  EXPORT SELECTED BI:BI RESULTS to use Rubin's method

  res <- c(BIres$BI$JamesBI[1,1],BIres$BI$JamesBI[1,2],
           BIres$BI$BangBI[2,1],BIres$BI$BangBI[2,2],
           BIres$BI$BangBI[1,1],BIres$BI$BangBI[1,2])
  return(res)
}


.subRubinRule <- function(thetai,ui){
    m <- length(thetai)
    stopifnot(length(ui) == m)
    ui <- ui^2
    theta <- mean(thetai)
    u <- mean(ui)
    B <- var(thetai)
    T <- u + (1+1/m)*B
    SE <- sqrt(T)
    nu <- (m-1)*(1+u/(1+1/m)/B)^2
    ME <- abs(qt(0.975,df=nu))
    ll <- theta - ME * SE
    ul <- theta + ME * SE
    c(theta,SE,ll,ul)
}
