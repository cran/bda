\name{blinding.BI}
\alias{blinding.BI}

\title{Computing Blinding Index}
\usage{
  blinding.BI(group, guess)
}
\arguments{
\item{guess}{guess response from blinding survey}
\item{group}{group assignments. 0='control',1='treatment' or 'active',
2='I don't know', or 'DNK'. Missing values are allowed.}
}
\value{
\itemize{
\item blinding index (with bootstrapping CI in the presence of missing responses.
}
}
\description{
To be updated.
}
\examples{
\donttest{
u1      = 5.5 # trt
u2      = 2.0 # ctrl
theta   = 3.2 # sham
sigma2  = 2.5   # v(rij)
ntreat  = 500      
nsham   = 500

beta0 = 1.0
beta1 = 2.0
beta2 = 1.0 # no contamination

Tind  = c(rep(1, ntreat), rep(0,nsham))  #treatment group indicator
u1v   = rep(u1,ntreat)
u2v   = rep(u2,nsham)
uv    = c(u1v,u2v)
tauv  = uv - rep(u2, ntreat+nsham)
r = rnorm(ntreat + nsham, mean = 0, sd = sqrt(sigma2))
q = 1/(1 + exp(-(beta0 + beta1*Tind + beta2*(tauv+r))))
bernGen = function(qq){rbinom(1,1,qq)}
I = sapply(q,bernGen)
x = uv + theta*I + r   # fixed sham effect
## I have concerns about the error term(s). x.sham~N(theta,sigma.sham)?
sigma.sham = 1.5
r2 = rnorm(ntreat + nsham, mean = 0, sd = sqrt(sigma.sham))
x = (uv + r) + theta*I #+ r2   # fixed sham effect

out1 <- blinding.BI(group=Tind,guess=I);
out1

}
}
