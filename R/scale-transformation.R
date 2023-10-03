# scale transformation function

#' Transform samples to uniform scale
#'
#' Since Kullback-Leibler divergence is scale-invariant, its sample-based
#' approximations can be computed on a conveniently chosen scale. This helper
#' functions transforms each variable in a way that all marginal distributions
#' of the joint dataset \eqn{(X,Y)} are uniform. In this way, the scales of
#' different variables are rendered comparable, with the idea of a better
#' performance of neighbour-based methods in this situation.
#'
#' @inherit kld_est_nn params
#' @returns A list with fields `X` and `Y`, containing the transformed samples.
#' @examples
#' # 2D example
#' n <- 10L
#' X <- cbind(rnorm(n, mean = 0, sd = 3),
#'            rnorm(n, mean = 1, sd = 2))
#' Y <- cbind(rnorm(n, mean = 1, sd = 2),
#'            rnorm(n, mean = 0, sd = 2))
#' to_uniform_scale(X, Y)
#' @export
to_uniform_scale <- function(X, Y) {

    XY <- rbind(X,Y)
    F_XY <- apply(X = XY, MARGIN = 2, FUN = ecdf)

    d <- ncol(X)

    FX <- matrix(NA, nrow = nrow(X), ncol = d)
    FY <- matrix(NA, nrow = nrow(Y), ncol = d)

    for (i in 1:d) {

        FX[ ,i] <- F_XY[[i]](X[,i])
        FY[ ,i] <- F_XY[[i]](Y[,i])

    }

    list(X = FX, Y = FY)
}
