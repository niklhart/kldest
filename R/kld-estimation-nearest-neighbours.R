# Different variants of nearest neighbour Kullback-Leibler divergence estimation

#' k-nearest neighbour KL divergence estimator
#'
#' This function estimates Kullback-Leibler divergence \eqn{D_{KL}(P||Q)} between
#' two continuous distributions \eqn{P} and \eqn{Q} using nearest-neighbour (NN)
#' density estimation in a Monte Carlo approximation of \eqn{D_{KL}(P||Q)}.
#'
#' Input for estimation is a sample `X` from \eqn{P} and either the density
#' function `q` of \eqn{Q} (one-sample problem) or a sample `Y` of \eqn{Q}
#' (two-sample problem). In the two-sample problem, it is the estimator in Eq.(5)
#' of Wang et al. (2009). In the one-sample problem, the asymptotic bias (the
#' expectation of a Gamma distribution) is substracted, see Pérez-Cruz (2008),
#' Eq.(18).
#'
#' References:
#'
#' Wang, Kulkarni and Verdú, "Divergence Estimation for Multidimensional
#' Densities Via k-Nearest-Neighbor Distances", IEEE Transactions on Information
#' Theory, Vol. 55, No. 5 (2009).
#'
#' Pérez-Cruz, "Kullback-Leibler Divergence Estimation of Continuous
#' Distributions", IEEE International Symposium on Information Theory (2008).
#'
#' @param X,Y  `n`-by-`d` and `m`-by-`d` matrices, representing `n` samples from
#'    the true distribution \eqn{P} and `m` samples from the approximate distribution
#'    \eqn{Q}, both in `d` dimensions. Vector input is treated as a column matrix.
#'    `Y` can be left blank if `q` is specified (see below).
#' @param q The density function of the approximate distribution \eqn{Q}. Either
#'    `Y` or `q` must be specified.
#' @param k The number of nearest neighbours to consider for NN density estimation.
#'    Larger values for `k` generally increase bias, but decrease variance of the
#'    estimator. Defaults to `k = 1`.
#' @param eps Error bound in the nearest neighbour search. A value of `eps = 0`
#'    implies exact nearest neighbour search, otherwise approximate nearest
#'    neighbours are sought. Defaults to `eps = 0.01`.
#' @param log.q If `TRUE`, function `q` is the log-density rather than the density
#'    of the approximate distribution \eqn{Q} (default: `log.q = FALSE`).
#' @returns A scalar, the estimated Kullback-Leibler divergence \eqn{D_{KL}(P||Q)}.
#' @example examples/nn-estimators.R
#' @export
kld_est_nn <- function(X, Y = NULL, q = NULL, k = 1L, eps = 0.01, log.q = FALSE) {

    # get important dimensions
    X <- as.matrix(X)
    d <- ncol(X) # number of dimensions
    n <- nrow(X) # number of samples in X

    # get distances to nearest neighbors from a kd-tree using nn2 from the RANN package
    # (the closest is the point itself, hence we consider one more neighbour)
    nnXX <- RANN::nn2(X, X, k = k+1, eps = eps)$nn.dists[ ,k+1]

    # check validity of input: one- or two-sample problem?
    if (!xor(is.null(Y),is.null(q))) stop("Either input Y or q must be provided.")

    if (!is.null(Y)) {
        # two-sample problem
        Y <- as.matrix(Y) # number of dimensions must be the same in X and Y
        m <- nrow(Y)      # number of samples in Y

        nnYX <- RANN::nn2(Y, X, k = k, eps = eps)$nn.dists[ ,k]

        return (log(m/(n-1)) - d*mean(log(nnXX/nnYX)))

    } else {
        # one-sample problem
        log_phat_X <- -log(n-1) + lgamma(0.5*d+1) - 0.5*d*log(pi) - d*log(nnXX) + digamma(k)
        log_q_X    <- if (log.q) {
            apply(X, MARGIN = 1, FUN = q)
        } else {
            log(apply(X, MARGIN = 1, FUN = q))
        }

        return (mean(log_phat_X-log_q_X))
    }
}



#' Generalized k-nearest neighbour KL divergence estimator.
#'
#' This function implements the generalized k-nearest neighbour estimator in
#' Wang et al. (2009), Eq.(17). In this estimator, the number of nearest
#' neighbours to consider may differ between samples and sample points.
#' Currently, this estimator is only implemented for the two-sample problem.
#'
#' Reference:
#' Wang, Kulkarni and Verdú, "Divergence Estimation for Multidimensional
#' Densities Via k-Nearest-Neighbor Distances", IEEE Transactions on Information
#' Theory, Vol. 55, No. 5 (2009).
#'
#' @inherit kld_est_nn params return
#' @param l,k Scalars or numeric vectors of length `n` and `m`, representing the
#'   number of nearest neighbours to use for nearest neighbour density estimation
#'   of P and Q, respectively, for each of the data points \code{X[i,]}, with
#'   \code{i} ranging from \code{1} to \code{n}, and \code{Y[j,]}, with
#'   \code{j} ranging from \code{1} to \code{m}. The default is `l = k = 1`.
#'   In the special case that `l = k` and `k` is scalar, the estimator coincides
#'   with `kld_est_nn(X, Y, k).
#'
kld_est_gnn <- function(X, Y = NULL, q = NULL, l = k, k = 1, eps = 0.01, log.q = FALSE) {

    # get important dimensions
    X <- as.matrix(X)
    d <- ncol(X) # number of dimensions, must be the same in X and Y
    n <- nrow(X) # number of samples in X

    # check validity of input l
    if (length(l) == 1) l <- rep(l, n)
    stopifnot(length(l) == n, max(l) < n)

    # nearest neighbours within X (delete the NN in X to X, which is X itself)
    nnXX <- RANN::nn2(X, X, k = max(l)+1, eps = eps)$nn.dists[ ,-1, drop = FALSE]

    # distance to the l-th nearest neighbour
    rho_l <- vapply(1:n, function(i) nnXX[i,l[i]], 1)

    # check validity of input: one- or two-sample problem?
    if (!xor(is.null(Y),is.null(q))) stop("Either input Y or q must be provided.")

    if (!is.null(Y)) {
        # two-sample problem
        Y <- as.matrix(Y)
        m <- nrow(Y) # number of samples in Y

        # check validity of input k
        if (length(k) == 1) k <- rep(k, n)
        stopifnot(ncol(Y) == d, length(k) == n, max(k) <= m)

        # nearest neighbours to X within Y
        nnYX <- RANN::nn2(Y, X, k = max(k), eps = eps)$nn.dists

        # distance to the k-th nearest neighbour
        nu_k <- vapply(1:n, function(i) nnYX[i,k[i]], 1)

        # equation (17) from Wang et al.
        return(d*mean(log(nu_k/rho_l)) + mean(digamma(l) - digamma(k)) + log(m/(n-1)))

    } else {
        # one-sample problem
        log_phat_X <- -log(n-1) + lgamma(0.5*d+1) - 0.5*d*log(pi) - d*log(rho_l) + digamma(l)
        log_q_X    <- if (log.q) {
            apply(X, MARGIN = 1, FUN = q)
        } else {
            log(apply(X, MARGIN = 1, FUN = q))
        }
        return(mean(log_phat_X-log_q_X))
    }

}


#' Bias-reduced generalized k-nearest-neighbour KL divergence estimation
#'
#' This is the bias-reduced generalized k-NN based KL divergence estimator from
#' Wang et al. (2009) specified in Eq.(29).
#'
#' `TODO`: more details on the algorithm!
#'
#' Reference:
#' Wang, Kulkarni and Verdú, "Divergence Estimation for Multidimensional
#' Densities Via k-Nearest-Neighbor Distances", IEEE Transactions on Information
#' Theory, Vol. 55, No. 5 (2009). DOI: https://doi.org/10.1109/TIT.2009.2016060
#'
#' @inherit kld_est_nn params return examples
#' @param max.k Maximum numbers of nearest neighbours to compute (default: `50`).
#'    A larger `max.k` may yield a more accurate KL-D estimate (see `warn.max.k`),
#'    but will always increase the computational cost.
#' @param warn.max.k If `TRUE` (the default), warn if `max.k` is such that more
#'    than `max.k` neighbours are within the neighbourhood for some data
#'    points. In this case, only the first `max.k` neighbours will be counted.
#'    As a consequence, `max.k` may required to be increased.
#' @export
kld_est_brnn <- function(X, Y, max.k = 50, warn.max.k = TRUE, eps = 0.01) {

    # get important dimensions
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    d <- ncol(X) # number of dimensions, must be the same in X and Y
    n <- nrow(X) # number of samples in X
    m <- nrow(Y) # number of samples in Y

    # get distances to nearest neighbours from a kd-tree using nn2 from the RANN package
    nnXX <- RANN::nn2(X, X, k = min(max.k+1,n), eps = eps)$nn.dists[,-1, drop = FALSE]    # delete the NN in X to X, which is X itself
    nnYX <- RANN::nn2(Y, X, k = min(max.k,m), eps = eps)$nn.dists

    # maximum of nearest neighbours in X and Y
    delta <- pmax(nnXX[,1], nnYX[,1])

    # how many data points are closer than delta?
    l <- rowSums(nnXX <= delta)
    k <- rowSums(nnYX <= delta)

    # check that max.k was chosen large enough and produce a warning otherwise
    if (warn.max.k && max.k %in% c(l,k)) {
        warning(paste('The "max.k"-th neighbour is inside the "delta" environment',
                      'for some datapoint, and further neighbours are ignored.',
                      'Increase input "max.k" to investigate the impact of this.'))
    }

    # distance to the k-th / l-th nearest neighbours
    rho_l <- vapply(1:n, function(i) nnXX[i,l[i]], 1)
    nu_k <- vapply(1:n, function(i) nnYX[i,k[i]], 1)

    # equation (17) from Wang et al.
    d/n*sum(log(nu_k/rho_l)) + mean(digamma(l) - digamma(k)) + log(m/(n-1))
}

