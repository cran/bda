### modified based on function density(...)
## 2012/07/10 by Bin Wang

tkde <-
function (x, func.W, bw = "nrd0", adjust = 1, kernel = c("gaussian", 
    "epanechnikov", "rectangular", "triangular", "biweight", 
    "cosine", "optcosine"), window = kernel, 
     width, n = 512, from, to, cut = 3, na.rm = FALSE, 
    ...) 
{
    if (length(list(...))) 
        warning("non-matched further arguments are disregarded")

    stopifnot(!missing(func.W))
    
    if (!is.numeric(x)) 
        stop("argument 'x' must be numeric")
    name <- deparse(substitute(x))
    x <- as.vector(x)
    n.user <- n
    n <- max(n, 512)
    if (n > 512) 
      n <- 2^ceiling(log2(n))

    xt = func.W(x)

    bw0 = bw.nrd(x)
    
    if(missing(from)) from = min(x) - cut * bw0
    if(missing(to)) to = max(x) + cut * bw0
    from.t = func.W(from)
    to.t = func.W(to)

    out = density(xt, bw=bw, adjust=adjust, kernel=kernel, weights=NULL,
      window=window, width=width, give.Rkern = FALSE, n=n,
      from = from.t, to = to.t, cut = cut, na.rm = na.rm)

    x <- seq.int(from, to, length.out = n.user)
    x0 = func.W(x)
    y =  approx(out$x, out$y, x0)$y
    kappa = (to - from)/n.user * sum(y)
    list(x = x, y = y/kappa, bw = out$bw,
         data.name = name)
  }
