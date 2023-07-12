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
