% This is a template R document file
% Distributed under GPL 3 or later

\name{fit.FSD}
\alias{fit.FSD}
\alias{fit.GLD}
\alias{fit.Copula}
\alias{print.FSD}
\alias{plot.FSD}
\alias{lines.FSD}

\title{Fitting firm size-age distributions}
\description{
 To fit firm size and/or age distributions based on binned data.
}
\usage{
fit.FSD(x,breaks,dist)
}
\arguments{
  \item{x}{'x' can be a vector or a matrix, or a 'histogram'.}
  \item{breaks}{a matrix with two columns if 'x' is a matrix. Otherwise, it is a vector. Can be missing if 'x' is a vector and 'x' will be grouped uisng the default parameters with \code{hist}. }
  \item{dist}{distribution type for 'x' and/or 'y'. Options include \code{Weibull}, \code{gpd} -- generalized Pareto distribution, \code{pd} or
 \code{Pareto}, \code{EWD}-- exponentiated Weibull distribtion.}
}

\value{
 \item{x.fit, y.fit}{Fitted marginal distribution for row- and column-data when 'x' is a matrix.}
 \item{Psi}{Plackett estimate of the Psi. ONLY for 'fit.FSD2'}

 If each of \code{x.fit} and \code{y.fit}, the following values are available:

 \item{xhist}{histogram.}
 \item{dist}{distribution type. Options include \code{Weibull},
   \code{gpd} -- generalized Pareto distribution, \code{pd} or
   \code{Pareto}, \code{EWD}-- exponentiated Weibull distribtion.}
 \item{size}{Total number of observations (sample size)}
 \item{pars}{Estimates of parameters.}
 \item{y,y2,x}{the pdf (\code{y}) and cdf (\code{y2}) values evaluated on a grid \code{x}}
}
\examples{
x <- rweibull(1000,2,1)
(out <- fit.FSD(x))

data(FSD)
b <- c(0,1:5,6,11,16,21,26,38,Inf);
x <- as.numeric(FirmAge[38,]);
## not run
##(out <- fit.FSD(x))#treated as raw data
xh <- binning(counts=x, breaks=b)
(out <- fit.FSD(xh))
#(out <- fit.FSD(xh,dist="ewd"))
(out <- fit.FSD(xh,dist="pd"))
(out <- fit.FSD(xh,dist="gpd"))


x <- as.numeric(FirmSize[nrow(FirmSize),])
brks.size <- c(0,4.5,9.5,19.5,49.5,99.5,249.5,499.5,
               999.5,2499.5,4999.5, 9999.5,Inf)
xh <- binning(counts=x,breaks=brks.size)
Fn <- cumsum(x)/sum(x);Fn
k <- length(Fn)
i <- c((k-5):(k-1))
out1 <- fit.GLD(xh, qtl=brks.size[i+1],
               qtl.levels=Fn[i],lbound=0)

i <- c(2,3,4,10,11)
out2 <- fit.GLD(xh, qtl=brks.size[i+1],
               qtl.levels=Fn[i],lbound=0)

i <- c(1,2,8,9,10)
out3 <- fit.GLD(xh, qtl=brks.size[i+1],
               qtl.levels=Fn[i],lbound=0)

plot(xh,xlim=c(0,120))
lines(out1, col=2)
lines(out2, col=3)
lines(out3, col=4)

ZipfPlot(xh,plot=TRUE)
lines(log(1-out1$y2)~log(out1$x), col=2)
lines(log(1-out2$y2)~log(out2$x), col=3)
lines(log(1-out3$y2)~log(out3$x), col=4)

## sample codes for the figures and tables in the PLoS ONE manuscript

## Table 1 *****************************************************
#xtable(Firm2)

## Figure 1 ****************************************************

#rm(list=ls())
#require(bda)
#data(FSD)

##postscript(file='fig1.eps',paper='letter')
#par(mfrow=c(2,2))

#tmp <- ImportFSD(FirmAge, year=2014,type="age");
#xh <- tmp$age
#plot(xh, xlab="X", main="(a) Histogram of 2014 firm age data")
#fit0 <- fit.FSD(xh, dist='exp')
#lines(fit0$x.fit$ly~fit0$x.fit$lx, lty=2)

#ZipfPlot(fit0, lty=2,
#         xlab="log(b)",ylab="log(r)",
#         main="(b) Zipf plot of firm age (2014)")

#tmp <- ImportFSD(FirmAge, year=1988,type="age");
#xh2 <- tmp$age
#fit1 <- fit.FSD(xh2, dist='exp')
#ZipfPlot(fit1, lty=2,
#         xlab="log(b)",ylab="log(r)",
#         main="(c) Zipf plot of firm age (1988)")

#ZipfPlot(fit0, lty=2, type='l',
#         xlab="log(b)",ylab="log(r)",
#         main="(d) Zipf plots of firm age (1979-2014)")
#for(year in 1979:2014){
#    tmp <- ImportFSD(FirmAge, year=year,type="age")
#    tmp2 <- ZipfPlot(tmp$age, plot=FALSE)
#    lines(tmp2)
#}

#dev.off()

## Figure 2 ****************************************************

#rm(list=ls())
#require(bda)
#data(FSD)

#postscript(file='fig2.eps',paper='letter')
#par(mfrow=c(2,2))
## plot (a)
#tmp <- ImportFSD(FirmSize, year=2014,type="size");
#yh <- tmp$size
#plot(yh, xlab="Y", main="(a) Histogram of firm size (2014)")

## plot (b)
#lbrks <- log(yh$breaks); lbrks[1] <- log(1)
#cnts <- yh$freq
#xhist2 <- binning(counts=cnts,breaks=lbrks)
#plot(xhist2,xlab="log(Y)",
#     main="(b) Histogram of firm size (2014, log-scale)")

## plot (c)
#ZipfPlot(yh, plot=TRUE,
#         main="(c) Zipf plot firm size (2014)",
#         xlab="log(r)",ylab="log(b)")

## plot (d)

#ZipfPlot(yh, plot=TRUE,type='l',
#         main="(d) Zipf plots of firm size (1977-2014)",
#         xlab="log(r)",ylab="log(b)")

#res <- NULL
#for(year in 1977:2014){
#    tmp <- ImportFSD(FirmSize,year=year,type="size");
#    yh0 <- tmp$size
#    zipf1 <- ZipfPlot(yh0,plot.new=FALSE)
#    lines(zipf1)
#    res <- c(res, zipf1$slope)
#}

#dev.off()

#mean(res); sd(res)
#quantile(res, prob=c(0.025,0.097))


## Figure 3a & 3b ****************************************************
#rm(list=ls())
#require(bda)
#data(FSD)

#postscript(file='fig3a.eps',paper='letter')

tmp <- ImportFSD(FirmAge, year=2014,type="age");
xh <- tmp$age
(fit1 <- fit.FSD(xh, dist='exp'));
(fit2 <- fit.FSD(xh, dist='weibull'));
#(fit3 <- fit.FSD(xh, dist='ewd')); #slow
#(fit0 <- fit.FSD(xh, dist='gpd')); # test, not used
(fit4 <- fit.FSD(xh, dist='gld'));

#plot(xh, xlab="X", main="(a) Fitted Distribtions")
#lines(fit1, lty=1, col=1,lwd=3)
#lines(fit2, lty=2, col=1,lwd=3)
#lines(fit3, lty=3, col=1,lwd=3)
#lines(fit4, lty=4, col=1,lwd=3)
#lines(fit0, lty=4, col=2,lwd=3)

legend("topright",cex=2,
       legend=c("EXP","Weibull","EWD","GLD"),
       lty=c(1:4),lwd=rep(3,4),col=rep(1,4))
dev.off()

## Table 3
#r2 <- c(fit1$x.fit$Dn.Zipf, fit2$x.fit$Dn.Zipf,
#        fit3$x.fit$Dn.Zipf, fit4$x.fit$Dn.Zipf)
#dn <- c(fit1$x.fit$Dn, fit2$x.fit$Dn,
#        fit3$x.fit$Dn, fit4$x.fit$Dn)
#aic <- c(fit1$x.fit$AIC, fit2$x.fit$AIC,
#         fit3$x.fit$AIC, fit4$x.fit$AIC)
#bic <- c(fit1$x.fit$BIC, fit2$x.fit$BIC,
#         fit3$x.fit$BIC, fit4$x.fit$BIC)
#aicc <- c(fit1$x.fit$AICc, fit2$x.fit$AICc,
#          fit3$x.fit$AICc, fit4$x.fit$AICc)
#tbl3 <- data.frame(
#    Dn.Zipf = r2, Dn=dn, AIC=aic, BIC=bic,AICc=aicc)
#rownames(tbl3) <- c("EXP","WD","EWD","GLD")
##save(tbl3, file='tbl3.Rdata')
#tbl3
#require(xtable)
#xtable(tbl3,digits=c(0,4,4,0,0,0))


#postscript(file='fig3b.eps',paper='letter')
#par(mfrow=c(1,1))
#ZipfPlot(fit1, plot.new=TRUE,col=1,lwd=3,lty=1,
#         xlab="log(x)",ylab="log(S(x))",
#         main="(b) Zipf Plots of Fitted Distribtions")
#ZipfPlot(fit2, plot.new=FALSE,col=1,lty=2,lwd=3)
#ZipfPlot(fit3, plot.new=FALSE,col=1,lty=3,lwd=3)
#ZipfPlot(fit4, plot.new=FALSE,col=1,lty=4,lwd=3)

#legend("bottomleft",cex=2,
#       legend=c("EXP","Weibull","EWD","GLD"),
#       lty=c(1:4),lwd=rep(3,4),col=rep(1,4))
#dev.off()

## TABLE 4 ############################################
## More Results **************************************

#mysummary <- function(x,dist){
#    .winner <- function(x,DIST=dist){
#        isele <- which(x==min(x))
#        if(length(isele)>1)
#            warning("multiple winners, only the first is used")
#        DIST[isele[1]]
#    }
#    .mysum <- function(x,y) sum(y==x)
#    mu1 <- tapply(x$Dn.Zipf,x$Dist, mean, na.rm=TRUE)
#    mu2 <- tapply(x$Dn,x$Dist, mean, na.rm=TRUE)
#    sd1 <- tapply(x$Dn.Zipf,x$Dist, sd, na.rm=TRUE)
#    sd2 <- tapply(x$Dn,x$Dist, sd, na.rm=TRUE)
#    win1 <- tapply(x$Dn.Zipf,x$Year, .winner,DIST=dist)
#    win2 <- tapply(x$Dn,x$Year, .winner,DIST=dist)
#    win3 <- tapply(x$BIC,x$Year, .winner,DIST=dist)
#    dist0 <- levels(as.factor(x$Dist))
    
#    out1 <- sapply(dist0, .mysum,y=win1)
#    out2 <- sapply(dist0, .mysum,y=win2)
#    out3 <- sapply(dist0, .mysum,y=win3)
#    #sele1 <- match(names(out1), dist0)
#    #sele2 <- match(names(out2), dist0)
#    #sele3 <- match(names(out3), dist0)
    
#    out <- data.frame(Mean.Dn.Zipf=mu1, SD.Dn.Zipf=sd1,
#                      Mean.Dn=mu2, SD.Dn=sd2,
#                      Win.Zipf=out1,
#                      Win.Dn=out2,
#                      Win.BIC=out3)
#    out
#}

#require(bda) #version 14.3.11+
#data(FSD)

#Dn.Zipf <- NULL
#Dn <- NULL
#BIC <- NULL
#Year <- NULL
#DIST <- NULL; dist0 <- c("EXP","WD","EWD","GLD")
#for(year in 1983:2014){
#    tmp <- ImportFSD(FirmAge, year=year,type="age");
#    xh <- tmp$age
#    fit1 <- fit.FSD(xh, dist='exp');fit1
#    DIST <- c(DIST,"EXP")
#    Year <- c(Year, year)
#    Dn.Zipf <- c(Dn.Zipf, fit1$x.fit$Dn.Zipf)
#    Dn <- c(Dn, fit1$x.fit$Dn)
#    BIC <- c(BIC, fit1$x.fit$BIC)
#    fit1 <- fit.FSD(xh, dist='weibull');fit1
#    DIST <- c(DIST,"WD")
#    Year <- c(Year, year)
#    Dn.Zipf <- c(Dn.Zipf, fit1$x.fit$Dn.Zipf)
#    Dn <- c(Dn, fit1$x.fit$Dn)
#    BIC <- c(BIC, fit1$x.fit$BIC)
#    fit1 <- fit.FSD(xh, dist='ewd');fit1
#    DIST <- c(DIST,"EWD")
#    Year <- c(Year, year)
#    Dn.Zipf <- c(Dn.Zipf, fit1$x.fit$Dn.Zipf)
#    Dn <- c(Dn, fit1$x.fit$Dn)
#    BIC <- c(BIC, fit1$x.fit$BIC)
#    fit1 <- fit.FSD(xh, dist='gld');fit1
#    DIST <- c(DIST,"GLD")
#    Year <- c(Year, year)
#    Dn.Zipf <- c(Dn.Zipf, fit1$x.fit$Dn.Zipf)
#    Dn <- c(Dn, fit1$x.fit$Dn)
#    BIC <- c(BIC, fit1$x.fit$BIC)
#}

#RES <- data.frame(Year=Year,Dist=DIST,Dn.Zipf=Dn.Zipf,Dn=Dn,BIC=BIC)
#save(RES, file='tbl4.Rdata')
#load(file='tbl4.Rdata')

#DIST <- c("EXP","WD","EWD","GLD")
#sele <- RES$Year>=1983 & RES$Year<=1987;sum(sele)
#(out <- mysummary(RES[sele,],dist=DIST))
#xtable(out[c(2,4,1,3),],digits=c(0,4,4,4,4,0,0,0))

#sele <- RES$Year>=1988 & RES$Year<=1992;sum(sele)
#(out <- mysummary(RES[sele,],dist=DIST))
#xtable(out[c(2,4,1,3),],digits=c(0,4,4,4,4,0,0,0))

#sele <- RES$Year>=1993 & RES$Year<=1997;sum(sele)
#(out <- mysummary(RES[sele,],dist=DIST))
#xtable(out[c(2,4,1,3),],digits=c(0,4,4,4,4,0,0,0))

#sele <- RES$Year>=1998 & RES$Year<=2002;sum(sele)
#(out <- mysummary(RES[sele,],dist=DIST))
#xtable(out[c(2,4,1,3),],digits=c(0,4,4,4,4,0,0,0))

#sele <- RES$Year>=2003 & RES$Year<=2014;sum(sele)
#(out <- mysummary(RES[sele,],dist=DIST))
#xtable(out[c(2,4,1,3),],digits=c(0,4,4,4,4,0,0,0))

#(out <- mysummary(RES,dist=DIST))
#xtable(out[c(2,4,1,3),],digits=c(0,4,4,4,4,0,0,0))

## Figure 4a & 4b ******************************************
#rm(list=ls())
#require(bda)
#data(FSD)


tmp <- ImportFSD(FirmSize, year=2014,type="size");
yh <- tmp$size

(fit1 <- fit.FSD(yh, dist='lognormal'));
(fit2 <- fit.FSD(yh, dist='pareto'));
(fit3 <- fit.FSD(yh, dist='gpd'));
(fit4 <- fit.FSD(yh, dist='gld'));

#postscript(file='fig4a.eps',paper='letter')
#par(mfrow=c(1,1))
#plot(yh, xlab="Y", main="(a) Fitted Distributions",
#     xlim=c(0,50))
#lines(fit1, lty=1, col=1,lwd=3)
#lines(fit2, lty=2, col=1,lwd=3)
#lines(fit3, lty=3, col=1,lwd=3)
#lines(fit4, lty=4, col=1,lwd=3)


#legend("topright",cex=2,
#       legend=c("LN","PD","GPD","GLD"),
#       lty=c(1:4),lwd=rep(3,4),col=rep(1,4))
#dev.off()

#postscript(file='fig4b.eps',paper='letter')
#par(mfrow=c(1,1))
#ZipfPlot(fit1, plot.new=TRUE,col=1,lwd=3,lty=1,
#         xlab="log(y)",ylab="log(S(y))",
#         main="(b) Zipf Plots of Fitted Distribtions")
#ZipfPlot(fit2, plot.new=FALSE,col=1,lty=2,lwd=3)
#ZipfPlot(fit3, plot.new=FALSE,col=1,lty=3,lwd=3)
#ZipfPlot(fit4, plot.new=FALSE,col=1,lty=4,lwd=3)

#legend("bottomleft",cex=2,
#       legend=c("LN","PD","GPD","GLD"),
#       lty=c(1:4),lwd=rep(3,4),col=rep(1,4))
#dev.off()

## TABLE 5 ############################################
## More information about firm size


#Dn.Zipf <- NULL
#Dn <- NULL
#BIC <- NULL
#Year <- NULL
#DIST <- NULL;

#dist0 <- c("LN","PD","GPD","GLD")
#for(year in 1977:2014){
#    DIST <- c(DIST,dist0)
#    Year <- c(Year, rep(year,length(dist0)))
#    tmp <- ImportFSD(FirmSize, year=year,type="size");
#    xh <- tmp$size
#    fit1 <- fit.FSD(xh, dist='lognormal');
#    Dn.Zipf <- c(Dn.Zipf, fit1$x.fit$Dn.Zipf)
#    Dn <- c(Dn, fit1$x.fit$Dn)
#    BIC <- c(BIC, fit1$x.fit$BIC)
#    fit1 <- fit.FSD(xh, dist='pd');
#    Dn.Zipf <- c(Dn.Zipf, fit1$x.fit$Dn.Zipf)
#    Dn <- c(Dn, fit1$x.fit$Dn)
#    BIC <- c(BIC, fit1$x.fit$BIC)
#    fit1 <- fit.FSD(xh, dist='gpd');
#    Dn.Zipf <- c(Dn.Zipf, fit1$x.fit$Dn.Zipf)
#    Dn <- c(Dn, fit1$x.fit$Dn)
#    BIC <- c(BIC, fit1$x.fit$BIC)
#    fit1 <- fit.FSD(xh, dist='gld');
#    Dn.Zipf <- c(Dn.Zipf, fit1$x.fit$Dn.Zipf)
#    Dn <- c(Dn, fit1$x.fit$Dn)
#    BIC <- c(BIC, fit1$x.fit$BIC)
#}

#RES <- data.frame(Year=Year,Dist=DIST,Dn.Zipf=Dn.Zipf,Dn=Dn,BIC=BIC)
#save(RES, file='tbl5.Rdata')
#load(file='tbl5.Rdata')

#dist0 <- c("LN","PD","GPD","GLD")
#(out <- mysummary(RES,dist=dist0))
#require(xtable)
#xtable(out[c(3,2,1),],digits=c(0,4,4,4,4,0,0,0))


## Figure 5  ****************************************************

#xy <- tmp <- ImportFSD(Firm2);
#out1 <- fit.FSD(xy$xy, breaks=xy$breaks,dist=c("EWD","GPD"));out1
#out2 <- fit.FSD(xy$xy, breaks=xy$breaks,dist=c("EWD","GLD"));out2
#out3 <- fit.FSD(xy$xy, breaks=xy$breaks,dist=c("GLD","GPD"));out3
#out4 <- fit.FSD(xy$xy, breaks=xy$breaks,dist=c("GLD","GLD"));out4

#postscript(file="fig5.eps",paper='letter')

#par(mfrow=c(2,2))
#out <- out1
#res2 <- plot(out,grid.size=40,nlevels=30,
#             ylim=c(0,15),xlim=c(0,30),
#             xlab="Firm Age",ylab="Firm Size",
#             main="(a) Contour Plot -- (EWD+GPD)")

#out <- out2
#res2 <- plot(out,grid.size=40,nlevels=30,
#             ylim=c(0,15),xlim=c(0,30),
#             xlab="Firm Age",ylab="Firm Size",
#             main="(b) Contour Plot -- (EWD+GLD)")
#out <- out3
#res2 <- plot(out,grid.size=40,nlevels=30,
#             ylim=c(0,15),xlim=c(0,30),
#             xlab="Firm Age",ylab="Firm Size",
#             main="(c) Contour Plot -- (GLD+GPD)")
#
#out <- out4
#res2 <- plot(out,grid.size=40,nlevels=30,
#             ylim=c(0,15),xlim=c(0,30),
#             xlab="Firm Age",ylab="Firm Size",
#             main="(d) Contour Plot -- (GLD+GLD)")
##dev.off()


## Figure 6  ****************************************************
#xy2 <- tmp <- ImportFSD(Firm2);

## this is an example showing how to use partial data to fit FSD.
## get marginal frequency distribution for firm age:
#(X <- apply(xy2$xy,1,sum));
#brks.age <- c(0,1,2,3,4,5,6,11,16,21,26,38,Inf)
## get marginal frequency distribution for firm size:
#(Y <- apply(xy2$xy,2,sum));
#(Y <- Y[-1])
#brks.size <- c(5,10,20,50,100,250,500,1000,2500,5000,10000,Inf)
#mxy2 <- xy2$xy[,-1]

#(fitx1 <- fit.FSD(X, breaks=brks.age, dist="ewd"))
#(fitx2 <- fit.FSD(X, breaks=brks.age, dist="gld"))
#(fity1 <- fit.FSD(Y, breaks=brks.size, dist="gld"))
#(fity2 <- fit.FSD(Y, breaks=brks.size, dist="pd"))
#(fity3 <- fit.FSD(Y, breaks=brks.size, dist="gpd"))

#(out11 <- fit.Copula(fitx1, fity1, mxy2))
#(out12 <- fit.Copula(fitx1, fity2, mxy2))
#(out21 <- fit.Copula(fitx2, fity1, mxy2))
#(out22 <- fit.Copula(fitx2, fity2, mxy2))
#(out13 <- fit.Copula(fitx1, fity3, mxy2))
#(out23 <- fit.Copula(fitx2, fity3, mxy2))

##postscript(file="fig6.eps",paper='letter')

#par(mfrow=c(1,1))
#plot(out11,grid.size=40,nlevels=20,lty=2,
#     ylim=c(4.8,16),xlim=c(0,18),
#     xlab="Firm Age",ylab="Firm Size",
#     main="(d) Contour Plot -- (GLD+GLD)")
#plot(out12,grid.size=40,nlevels=30, col=2, plot.new=FALSE)
## or use the command below
## plot(out2,grid.size=50,nlevels=50, col=4, add=TRUE)
##dev.off()

}

\keyword{distribution}

