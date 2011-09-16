

bootkde <- function(x, method='z.score',scale=1, rounding = 'nearest',
                    from, to, alpha=0.05, gridsize=512L,na.rm=TRUE)
{
  method <- match.arg(tolower(method),c("z.score","quantile","cdf"))
  name <- deparse(substitute(x))
  sele = is.na(x)
  if(any(sele)){
    if(na.rm) x = x[!sele]
    else stop("'x' contains missing values.")
  }
  n = length(x)
  out = binning(x, scale=scale, method=rounding)

  X = out$x; F=out$counts; B=out$width/2; A=-B; N=length(X)

  if(missing(alpha))alpha=0.05
  if(alpha<=0|alpha>=1)stop("Invalid confidence level.")
  alpha = min(alpha,1-alpha)

  ucb=NULL; lcb=NULL;
  ## Set default bandwidth
  h0 = bw.nrd(rep(X,F)+runif(sum(F),-B,B))
  h = .Fortran(.F_hbmise,as.double(X), as.double(F),as.double(2*B),
    as.integer(N), hopt=as.double(h0))$hopt

  if(missing(from)) from = min(X) - 3*h
  if(missing(to)) to = max(X) + 3*h
  a <- from;  b <- to
  ## set grid to evaluate the pdf/cdf
  M <- 2^(ceiling(log2(gridsize)))
  gpoints <- seq(a, b, length = M)

  
  y = .Fortran(.F_ofcpdf,
    as.double(X), as.double(F),as.double(-B),as.double(B),
    as.integer(N), y=as.double(gpoints), as.integer(M),
    para=as.double(h))$y
  y=cumsum(y)
  
  tmp = as.matrix(rep(sum(F),1000), ncol=1);
  ## simulate the rounding errors by a pilot estimate of F(x) based on a BME
  out = apply(tmp,1, .bootemp,x=gpoints,y=X,f=F,lb=A,ub=B,
    Fx=y,from=from,to=to)

  if(method=="z.score"){
    sigs = apply(out,1,sd);
    z0 = qnorm(1-alpha/2.);
    ym = apply(out,1,median);
    lcb = ym-sigs*z0; lcb[lcb<0]=0;
    ucb = ym+sigs*z0; ucb[ucb>1]=1
    y=out[,1];
    type="Density"
  }else if(method=="quantile"){
    ym = apply(out,1,median);
    lcb = as.numeric(apply(out,1,quantile, alpha/2)); lcb[lcb<0]=0;
    ucb = as.numeric(apply(out,1,quantile, 1-alpha/2)); ucb[ucb>1]=1
    y=out[,1];
    type="Density"
  }else{
    out = apply(out,2,cumsum)*(b-a)/M
    ym = apply(out,1,median);
    lcb = as.numeric(apply(out,1,quantile, alpha/2)); lcb[lcb<0]=0;
    ucb = as.numeric(apply(out,1,quantile, 1-alpha/2)); ucb[ucb>1]=1
    y=out[,1];
    type="Probability"
  }
  
  return(structure(list(y=y,x=gpoints,ym=ym,n = sum(F),
                        ucb=ucb,lcb=lcb,alpha=alpha,
                        call = match.call(), data.name = name
                        ), class = "bcb"))
}

.bootemp <- function(n,x,y,f,lb,ub,Fx,from,to){
  tmp = cbind(y,f,lb,ub);M=length(x);u=runif(n);N=length(f);
  smpl = .Fortran(.F_remp,
    as.integer(N),as.double(y),as.double(f),as.double(lb), as.double(ub),
    as.integer(M), as.double(Fx),as.double(x),smpl=as.double(u))$smpl
  density(smpl,from=from,to=to)$y
}

plot.bcb <-
  function (x, main = NULL, xlab = "x", ylab = "Probability",
            lwd=1, col=1, lty=1, zero.line = TRUE,
            bgcol='gray',scb=FALSE, ...) 
{
  if (is.null(main)) main <- deparse(x$data.name)
  plot.default(x$x,x$ucb, main = main,
               xlab = xlab, ylab = ylab, type='n',...)
  if(scb){
    cord.x = c(x$x,rev(x$x))
    cord.y = c(x$lcb, rev(x$ucb))
    polygon(cord.x,cord.y,border=col,col=bgcol)
  }
  lines(x$x,x$y, lwd=lwd, col=col,lty=lty,...)
  
  if (zero.line) 
    abline(h = 0, lwd = 0.1, col = "gray")
  invisible(NULL)
}

lines.bcb <-
  function (x, lwd=1, col=1, lty=1, bgcol='gray', scb=FALSE, ...) 
{
  if(scb){
    cord.x = c(x$x,rev(x$x))
    cord.y = c(x$lcb, rev(x$ucb))
    polygon(cord.x,cord.y,border=col,col=bgcol)
  }
  lines(x$x,x$y, lwd=lwd, col=col,lty=lty,...)
  
  invisible(NULL)
}

print.bcb <- function (x, digits = NULL, ...) 
{
  cat("\nData: ", x$data.name, 
      " (", x$n, " obs.);\n", sep = "")
  print(summary(as.data.frame(x[c("x", "y","lcb","ucb")])),
        digits = digits, 
        ...)
  invisible(x)
}
