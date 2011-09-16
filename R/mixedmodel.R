.bfmm <- function(x,m=2,mu,type='gaussian',range.x, ...)
{  
  type = match.arg(tolower(type),
    c('gaussian', 'normal','beta','gamma','weibull'))
  idist = switch(type,
    gaussian = 0,
    normal = 0,
    beta = 1,
    gamma = 2,
    weibull = 3
    )
  
  n = length(x$x)
  if (missing(range.x))
    range.x <- c(x$x[1] - 0.5 * x$widths[1], x$x[n] + 0.5 * x$widths[n])
  xmin <- range.x[1L]
  xmax <- range.x[2L]
  xspan = xmax - xmin;

  x0 = (x$x - xmin)/xspan;
  x0wd = x$widths/xspan

  llk = 0 # compute llk.
  if(missing(mu)){
    m = round(m)
    if(m<1) stop("Invalid number of components.")
    if(n < m)
      warning("Number of mixing components is larger than the distinct bin centers!")
    m = min(m, n - 2)
    type = ifelse(m>1,'normal',type)
    idist = ifelse(m>1,0,idist)
    par = rep(0,3*m)  ## works only for two parameter distribution families.
    Iter = min(100, 2*choose(n-2,m))
    x0s = x0[-c(1,n)]
    par[3*(1:m)-1] = sample(x0s, m)
    out = .Fortran(.F_fitmm,
      as.double(x0), as.double(x$counts), as.double(x0wd),
      as.integer(n), as.integer(idist), 
      as.integer(m), par = as.double(par), llk=as.double(llk)
      )
    for(i in 1:Iter){
      par[3*(1:m)-1] = sample(x0s, m)
      out2 = .Fortran(.F_fitmm,
        as.double(x0), as.double(x$counts), as.double(x0wd),
        as.integer(n), as.integer(idist), 
        as.integer(m), par = as.double(par), llk=as.double(llk)
        )
      if(out2$llk[1] > out$llk[1]) out = out2
    }
  }else{
    m = length(mu);
    type = ifelse(m>1,'normal',type)
    idist = ifelse(m>1,0,idist)
    par = rep(0,3*m)  ## works only for two parameter distribution families.
    mu = (mu-xmin)/xspan
    par[3*(1:m)-1] = mu

    out =  .Fortran(.F_fitmm,
      as.double(x0), as.double(x$counts), as.double(x0wd),
      as.integer(n), as.integer(idist), 
      as.integer(m), par = as.double(par), llk=as.double(llk))
  }

  pars =  matrix(out$par, nrow=3, ncol=m)
  sele = pars[1,] > 0.00001
  pars = pars[,sele]
  pars = as.data.frame(pars)
  row.names(pars) <- c("Proportion","Para.1","Para.2")
  m = sum(sele)
  
  llk0 = out$llk
  K = ifelse(m == 1, 2., 3. * m - 1.);
  N = sum(x$counts);
  llktransform = 2.0 * sum(x$counts) * log(xspan)
  AIC = -2.0 * llk0 + 2.0 * K + llktransform;
  AICc = AIC + 2.0 * K * (K+1.) / (N - K - 1.0);
  BIC = -2.0 * llk0 + log(N) * K + llktransform;
  structure(list(para = pars, type=type, m=m,
                 range.x=c(xmin,xmax), data = x,
                 llk = c("AIC"=AIC, "BIC"=BIC,"AICc"=AICc)),
            class='mm')
}

bfmm <- function(x,m=2,mu,type='gaussian',method='nelder',range.x, ...)
UseMethod("bfmm")

bfmm.default <- function(x,m=2,mu,type='gaussian',method='nelder',range.x, ...)
{  
  stopifnot(class(x) == 'bdata')
  method = match.arg(tolower(method), c('nelder', 'newton'))
  switch(method,
         nelder = .bfmm(x,m=m,mu=mu,type=type,range.x=range.x,...),
         newton = .binem(x,m=m,mu=mu,...)
         )
}

density.mm <- function(x,x0,gridsize=500,...)
  {
    if(is.null(x$range.x)){
      xmin = min(x$data$x); xmax = max(x$data$x);
      xspan = xmax - xmin;
      if(missing(x0)){
        x.ext = 1.5 * x$data$s; 
        x0 = seq(xmin-x.ext,xmax+x.ext, length=gridsize)
      }else{
        gridsize = length(x0)
      }
      xpts = x0
      fx = rep(0,gridsize)
      for(i in 1:x$m){
        fx0 = dnorm(x0,x$para[2,i],x$para[3,i])
        fx = fx + x$para[1,i] * fx0;
      }
    }else{
      xmin = x$range.x[1]; xmax = x$range.x[2];
      xspan = xmax - xmin;
      if(missing(x0)){
        s = x$data$s; z = x$data$x
        x.ext = 1.5*s/(max(z)-min(z)+1) 
        x0 = switch(x$type,
          gaussian = seq(0-x.ext,1+x.ext, length=gridsize),
          normal = seq(0,1+x.ext, length=gridsize),
          beta = seq(0,1, length=gridsize),
          gamma = seq(0,1+x.ext, length=gridsize),
          weibull = seq(0,1+x.ext, length=gridsize)
          )
        xpts = x0*xspan + xmin
      }else{
        gridsize = length(x0)
        xpts = x0 #save raw x0
        x0 = (x0 - xmin)/xspan
      }
      fx = rep(0,gridsize)
      for(i in 1:x$m){
        fx0 = switch(x$type,
          gaussian = dnorm(x0,x$para[2,i],x$para[3,i]),
          normal = dnorm(x0,x$para[2,i],x$para[3,i]),
          beta = dbeta(x0,x$para[2,i],x$para[3,i]),
          gamma = dgamma(x0,shape=x$para[2,i],scale=x$para[3,i]),
          weibull = dweibull(x0,shape=x$para[2,i],scale=x$para[3,i])
          )
        fx = fx + x$para[1,i] * fx0;
      }
      fx = fx/xspan
    }
    list(x = xpts, y = fx)
  }

bfmm.data.frame <- function(x,m=2,mu,type='gaussian',method='nelder',range.x,...){
  name <- deparse(substitute(x))
  if(is.null(x$x)||is.null(x$widths)||is.null(x$counts))
    stop("Information missing in data 'x'. Try 'help(rounding)'.")
  
  out = structure(list(x = x$x, widths = x$widths, counts = x$counts,
    scale = 1.0, data.name = name), class='bdata')
  bfmm(out,m=m,mu=mu,type=type,method=method,...)
}

bfmm.numeric <- function(x,m=2,mu,type='gaussian',method='nelder',
                         range.x,scale=1, rounding = 'nearest',...){
  out = rounding(x, scale=scale, method=rounding)
  bfmm(out,m=m,mu=mu,type=type,method=method,...)  
}

print.mm <- function(x,...)
  {
    if(!is.null(x$range.x))
      cat("\nData were rescaled from [",
          x$range.x[1], ",", x$range.x[2], "] to [0,1].",sep = "")
    
    cat("\nFitted a mixed model of ( ", x$m, " ) ", x$type,
        " component(s).\n", sep = "")
    print(x$llk)
    cat("\n\tParamters:\n")
    print(x$para)
    cat("\n\nUse density(x) to compute the fitted density function.\n\n")
    invisible(x)

  }

