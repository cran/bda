dsimpson <- function(x){
  j=0:4
  out = 0.5*dnorm(x)
  for(j in 0:4){
    out = out +0.1*dnorm(x,j/2-1,0.1)
  }
  list(y=out,x=x)
}

rsimpson <- function(n){
  x1 = rnorm(n)
  x2 = rnorm(n,0,0.1)
  g = runif(n)
  sele1 = g<0.5
  sele2 = g>=0.5 & g<0.6
  sele3 = g>=0.6 & g<0.7
  sele4 = g>=0.7 & g<0.8
  sele5 = g>=0.8 & g<0.9
  sele6 = g>=0.9
  if(sum(sele1)>0) out = x1[sele1]
  if(sum(sele2)>0) out = c(out, x2[sele2]-1)
  if(sum(sele3)>0) out = c(out, x2[sele3]-.5)
  if(sum(sele4)>0) out = c(out, x2[sele4])
  if(sum(sele5)>0) out = c(out, x2[sele5]+.5)
  if(sum(sele6)>0) out = c(out, x2[sele6]+1)
  out
}
