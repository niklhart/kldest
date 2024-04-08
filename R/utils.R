# utility functions

#' Matrix trace operator
#'
#' @param M A square matrix
#' @return The matrix trace (a scalar)
tr <- function(M) sum(diag(M))


#' Constant plus diagonal matrix
#'
#' Specify a matrix with constant values on the diagonal and on the off-diagonals.
#' Such matrices can be used to vary the degree of dependency in covariate matrices,
#' for example when evaluating accuracy of KL-divergence estimation algorithms.
#'
#' @param dim Dimension
#' @param diag Value at the diagonal
#' @param offDiag Value at off-diagonals
#' @return A \code{dim}-by-\code{dim} matrix
#' @examples
#' constDiagMatrix(dim = 3, diag = 1, offDiag = 0.9)
#'
#' @export
constDiagMatrix <- function(dim = 1, diag = 1, offDiag = 0) {
    offDiag * matrix(1, nrow = dim, ncol = dim) + (diag - offDiag) * diag(dim)
}


#' Trapezoidal integration in 1 or 2 dimensions
#'
#' @param h A length `d` numeric vector of grid widths.
#' @param fx A `d`-dimensional array (or a vector, if `d=1`).
#' @returns The trapezoidal approximation of the integral.
#' @examples
#' # 1D example
#' trapz(h = 1, fx = 1:10)
#' # 2D example
#' trapz(h = c(1,1), fx = matrix(1:10, nrow = 2))
#' @export
trapz <- function(h, fx) {

    d    <- length(h)
    fx   <- as.array(fx)
    dims <- dim(fx)
    stopifnot(length(dims) == d)

    prod(h) * switch(d,
           "1" = {
               sum(fx) - 0.5*(fx[1] + fx[dims])
           },
           "2" = {
               sum(fx) - 0.5*(sum(fx[c(1,dims[1]),]) + sum(fx[,c(1,dims[2])])) +
                   0.25 * sum(fx[c(1,dims[1]),c(1,dims[2])])
           },
           stop("Case d>2 not implemented yet.")
    )

}


#' Combinations of input arguments
#'
#' @param ... Any number of atomic vectors.
#' @return A data frame with columns named as the inputs, containing all input
#'      combinations.
#' @examples
#' combinations(a = 1:2, b = letters[1:3], c = LETTERS[1:2])
#' @export
combinations <- function(...) {
    args  <- list(...)
    largs <- vapply(args, length, 1)
    ncomb <- prod(largs)
    ntarg <- ncomb / largs
    neach  <- cumprod(largs) / largs
    ntimes <- ntarg / neach

    comb <- mapply(function(x, nt, ne) rep(x, times = nt, each = ne),
                   args, ntimes, neach, SIMPLIFY = FALSE)
    do.call(data.frame,comb)
}



#' Probability density function of multivariate Gaussian distribution
#'
#' @param x A vector of length `d` at which Gaussian density is evaluated.
#' @param mu A vector of length `d`, mean of Gaussian distribution.
#' @param Sigma A `d`-by-`d` matrix, covariance matrix of Gaussian distribution.
#' @returns The probability density of \eqn{N(\mu,\Sigma)} evaluated at `x`.
#' @examples
#' # 1D example
#' mvdnorm(x = 2, mu = 1, Sigma = 2)
#' dnorm(x = 2, mean = 1, sd = sqrt(2))
#' # Independent 2D example
#' mvdnorm(x = c(2,2), mu = c(1,1), Sigma = diag(1:2))
#' prod(dnorm(x = c(2,2), mean = c(1,1), sd = sqrt(1:2)))
#' # Correlated 2D example
#' mvdnorm(x = c(2,2), mu = c(1,1), Sigma = matrix(c(2,1,1,2),nrow=2))
#' @export
mvdnorm <- function(x, mu, Sigma) {

    d <- length(mu)
    iSigma <- solve(Sigma)

    (2*pi)^(-d/2)*sqrt(det(iSigma))*exp(-0.5*as.vector(t(x-mu) %*% iSigma %*% (x-mu)))
}



#' Detect if a one- or two-sample problem is specified
#'
#' @param Y A vector, matrix, data frame or `NULL`
#' @param q A function or `NULL`.
#' @return `TRUE` for a two-sample problem (i.e., `Y` non-null and `q = NULL`)
#'     and `FALSE` for a one-sample problem (i.e., `Y = NULL` and `q` non-null).
is_two_sample <- function(Y, q) {

    two.sample <- !is.null(Y)
    one.sample <- !is.null(q)

    if (!xor(two.sample,one.sample)) stop("Either input `Y` or `q` must be specified.")
    if (two.sample && is.function(Y)) stop("Input `Y` is a function -- argument `q` may have been provided by position rather than by name.")

    two.sample
}

