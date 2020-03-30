
bootkde <- function(x,freq,a,b,from, to, gridsize=512L,
                    method='boot')
{
    method <- match.arg(tolower(method),
                        c("maximum","mle",
                          "berkson","convolution",
                          "boot","bootstrap"))
    ## assume raw data have rounded in the same way with additive
    ## uniform errors following Uniform(a,b). Data can also be given
    ## as (x,F). If 'F' is missing, all with frequency 1.
    X <- unique(sort(x))
    N <- length(X)
    n <- length(x)
    if(missing(freq)){
        xt <- table(x)
        F <- as.numeric(xt)
    }else{
        F <- round(freq)
        if(any(F <= 0))
            stop("invalid 'freq' found")
        if(length(F) != n)
            stop("'x' and 'freq' have different length")
        if(N != n)
            stop("'x' values are not unique")
    }
    
    if(missing(a)){
        if(missing(b)){
            stop("'a' and/or 'b' needs to be specified")
        }else{
            stopifnot(b>0)
            a <- -b
        }
    }else{
        if(missing(b)){
            stopifnot(a<0)
            b <- -a
        }else{
            stopifnot(b>a)
            d <- 0.5 * (a+b) 
            if(d != 0){
                a <- -abs(d)
                b <- -a
                X <- X + d
            }
        }
    }
    B <- rep(b,N)
    A <- rep(a,N)

    a <- ifelse(missing(from), X[1]+A[1], from)
    b <- ifelse(missing(to), X[N] + B[N], to)
    gpoints <- seq(a, b, length = gridsize)

    if(method=="berkson"||method=="convolution"){
        ## Set default bandwidth
        h0 = bw.nrd(rep(X,F)+runif(sum(F),A,B))
        h = .Fortran(.F_hbmise,as.double(X), as.double(F),
            as.double(2*B),as.integer(N),
            hopt=as.double(h0))$hopt
        
        y <- .Fortran(.F_ofcpdf,
                      as.double(X),
                      as.double(F),
                      as.double(A),
                      as.double(B),
                      as.integer(N),
                      y=as.double(gpoints),
                      as.integer(gridsize),
                      para=as.double(h))$y
        pout <- h
    }else if(method=="mle"||method=="maximum"){
        xbar <- sum(X*F)/sum(F)
        s2 <- sum((X-xbar)^2*F)/(sum(F)-1)
        theta <- c(xbar,sqrt(s2))
        Theta <- .Fortran(.F_mlensimp,
                          as.double(X),
                          as.double(F),
                          as.double(A),
                          as.double(B),
                          as.integer(N),
                          para=as.double(theta))$para
        y <- dnorm(gpoints, Theta[1], Theta[2])
        pout <- Theta
    }else{
        ## Set default bandwidth
        h0 = bw.nrd(rep(X,F)+runif(sum(F),A,B))
        h = .Fortran(.F_hbmise,as.double(X), as.double(F),
            as.double(2*B),as.integer(N),
            hopt=as.double(h0))$hopt
        
        y <- .Fortran(.F_ofcpdf,
                      as.double(X),
                      as.double(F),
                      as.double(A),
                      as.double(B),
                      as.integer(N),
                      y=as.double(gpoints),
                      as.integer(gridsize),
                      para=as.double(h))$y
        pout <- h
        Fx <- cumsum(y)
        ##iter <- 100
        ##tmp <- as.matrix(rep(sum(F),iter), ncol=1);
        ## simulate the rounding errors by a pilot estimate of F(x)
        ## based on a BME
        ##out <- apply(tmp, 1, .bootemp,
        ##             x=gpoints,y=X,f=F,lb=A,ub=B,
        ##             Fx=Fx,from=from,to=to)
        ##y <- apply(out,1,median);
        y <- .bootemp(sum(F),gpoints,X,F,A,B,Fx,from,to)
    }
    
    
    list(y=y,x=gpoints,pars=pout)
}

.bootemp <- function(n,x,y,f,lb,ub,Fx,from,to){
    M=length(x);u=runif(n);N=length(f);
    
    smpl = .Fortran(.F_remp, as.integer(N), as.double(y), as.double(f),
        as.double(lb), as.double(ub), as.integer(M), as.double(Fx),
        as.double(x), smpl=as.double(u), as.integer(n))$smpl
    
    density(smpl,n=M,from=from,to=to)$y
}

