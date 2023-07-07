# utility functions

#' Matrix trace operator
#'
#' @param M A square matrix
#' @return The matrix trace (a scalar)
tr <- function(M) sum(diag(M))
