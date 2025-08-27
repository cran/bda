
blinding <- function(x, y)
{
    n <- length(x)
    if(length(y) != n) stop("'x' and 'y' differ in length")

    if(any(is.na(x) | is.na(y))) stop("missing values not allowed")
    ox <- order(x)
    x <- sort(x)
    y <- y[ox]

    k <- round(0.1*n)

    out <- .Fortran(.F_chgpt,
                    as.double(x),
                    as.double(y),
                    as.integer(n),
                    as.integer(k),
                    ans = as.double(1))
}

