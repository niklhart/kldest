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


#' Trapezoidal integration in `d` dimensions
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
               sum(fx) - 0.5*(sum(fx[c(1,dims[1]),]) + sum(fx[,c(1,dims[2])])) + 0.25 * sum(fx[c(1,dims[1]),c(1,dims[2])])
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


