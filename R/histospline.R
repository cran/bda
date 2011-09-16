##  Part of R package BDA
##  Copyright (C) 2009-2010 Bin Wang
##
##  Unlimited use and distribution (see LICENCE).

histospline <-
  function(x,f,gridsize=512L,na.rm=TRUE,just="center",binned=FALSE, scale=1.0){
    name <- deparse(substitute(x))
    if(scale <= 0) stop("Wrong 'scale' level.")
    x = x/scale
    ##  bin data if x is not binned.
    if(binned){
      if(missing(f))stop("Frequencies missing!")
      ijust = match.arg(tolower(just),
        c("center","left","right"))
      a = switch(ijust,center=0,left=0.5,right=-0.5)
      x1 = data.frame(x=x+a,f=f,b=0.5)
      X = x1$x; F = f;
    }else{
      if(any(is.na(x))) x = x[!is.na(x)]
      x0 = .discretize(x,na.rm=na.rm,just=just)
      X=x0$x; F=x0$f;
    }
    rf = F/sum(F)
    out = spline(X,rf,n=gridsize)
    f0 = out$y
    f0[f0<0]=0
    gpoints = out$x
    return(structure(list(y=f0/scale,x=gpoints*scale,bw=NULL,scale=scale,
                          call = match.call(), data.name = name),
                     class='bde'))
}
