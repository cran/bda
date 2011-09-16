perm.test <- function(x,y,fun,alternative = "two.sided", 
                      trials = 1000,...)
  UseMethod("perm")

perm <- function(x,y,fun,alternative = "two.sided", 
                 trials = 1000, ...)
  UseMethod("perm")

##  if fun is missing, define it as comparing distributional
##  difference; otherwise, user-defined.
perm.default <- function(x,y,fun,alternative = "two.sided", 
                      trials = 1000,...)
{
  xnam <- deparse(substitute(x));
  ynam <- deparse(substitute(y));
  nam = paste(xnam, '(',length(x), ') vs ', ynam, '(',length(y),")")
  if(trials<500) warning("Trial number might be too small.")
  ##  check data (x,y)
  if (!is.numeric(x) && !is.complex(x) && !is.logical(x)) {
    warning("'x' is not numeric or logical: returning NA")
    return(NA_real_)
  }
  x <- x[!is.na(x)]
  if (!is.numeric(y) && !is.complex(y) && !is.logical(y)) {
    warning("'y' is not numeric or logical: returning NA")
    return(NA_real_)
  }
  y <- y[!is.na(y)]

  alternative = match.arg(tolower(alternative),
      c("one.sided","two.sided"))

    
  ## use of replicate() with parameters:
  if(missing(fun)) fun = .edfperm

  D = fun(x,y)
  rfun <- function(x,y){
    nx = length(x); ny = length(y); n = nx + ny;
    n0 = sample(n)
    xy = c(x,y)
    fun(xy[n0[1:nx]], xy[n0[(nx+1):n]])
  }
  bar <- function(n, x,y) replicate(n, rfun(x=x,y=y))
  z = bar(trials, x=x,y=y)

  ##  z = replicate(trials, fun(x=x,y=y))
  ##  z = apply(as.matrix(rep(1,trials),ncol=1),1,fun,x=x,y=y)
  pv = switch(alternative,
    two.sided = mean(abs(D) <= abs(z)),
    one.sided = mean(D <= z),
    stop("Wrong test type!")
    )

  RVAL <- list(statistic = c(D = D), p.value = pv,
               method = "Permutation test to compare two samples/populations.", 
               data.name = nam)
  class(RVAL) <- "htest"
  return(RVAL)

}

.edfperm <- function(x,y){
  ## no randomization within this function.  It could be defined by
  ## users. The radomization should be a step in the main program. 
  rngx = range(c(x,y))
  x0 = seq(rngx[1],rngx[2],length=256);
  Fx1 = edf(x,x0);
  Fx2 = edf(y,x0);
  max(abs(Fx1$y-Fx2$y))
}

