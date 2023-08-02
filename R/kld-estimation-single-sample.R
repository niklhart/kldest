# KL-D estimation based on a single sample from p and a density function q.

#' Single-sample version of 1-NN divergence estimator from Wang et al. (2009).
#'
#' This function experiments with partial nearest-neighbour estimation, i.e.
#' only for p, but not for q, which may be available analytically (the estimated
#' approximate model). For now, it really doesn't work well.
#'
#' @param X An `n`-by-`d` matrix, representing `n` samples from the true
#'    distribution \eqn{P} in `d` dimensions. Vector input is treated as a
#'    column matrix.
#' @param q The density of the approximate model \eqn{Q}.
#' @examples
#' X <- rnorm(100)
#' Y <- rnorm(100, mean = 1, sd = 2)
#' q <- function(x) dnorm(x, mean = 1, sd = 2)
#' kld_gaussian(0,1,1,2^2)
#' kld_est_1nn(X, Y)
#' kld_est_1nn_1sample(X, q)
#' @returns A scalar, the estimated Kullback-Leibler divergence \eqn{D_{KL}(P||Q)}.
#' @export
kld_est_1nn_1sample <- function(X, q) {

    # get important dimensions
    X <- as.matrix(X)
    d <- ncol(X) # number of dimensions, must be the same in X and Y
    n <- nrow(X) # number of samples in X

    # get distances to nearest neighbors from kdTree using nn2 from the RANN package
    r <- RANN::nn2(X,X, k=2, eps=.01)[[2]][,2] # get 2 closest neighbors, then take the second (the closest is the point itself) to get n x 1 matrix

    log_phat_X <- -log(n-1) + lgamma(0.5*d+1) - 0.5*d*log(pi) - d*log(r)
    log_q_X    <- log(apply(X, MARGIN = 1, FUN = q))

    return (mean(log_phat_X-log_q_X))
}
