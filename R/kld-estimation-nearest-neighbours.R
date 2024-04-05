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
#'    (the default) implies an exact nearest neighbour search, for `eps > 0`
#'    approximate nearest neighbours are sought, which may be somewhat faster for
#'    high-dimensional problems.
#' @param log.q If `TRUE`, function `q` is the log-density rather than the density
#'    of the approximate distribution \eqn{Q} (default: `log.q = FALSE`).
#' @returns A scalar, the estimated Kullback-Leibler divergence \eqn{\hat D_{KL}(P||Q)}.
#' @example examples/nn-estimators.R
#' @export
kld_est_nn <- function(X, Y = NULL, q = NULL, k = 1L, eps = 0, log.q = FALSE) {

    # get important dimensions
    X <- as.matrix(X)
    d <- ncol(X) # number of dimensions
    n <- nrow(X) # number of samples in X

    # early return if there aren't enough data points in X
    if (k >= n) return(NA_real_)

    # get distances to nearest neighbors from a kd-tree using nn2 from the RANN package
    # (the closest is the point itself, hence we consider one more neighbour)
    nnXX <- RANN::nn2(X, X, k = k+1, eps = eps)$nn.dists[ ,k+1]

    # check validity of input: one- or two-sample problem?
    two.sample <- is_two_sample(Y, q)

    if (two.sample) {
        # two-sample problem
        Y <- as.matrix(Y) # number of dimensions must be the same in X and Y
        m <- nrow(Y)      # number of samples in Y

        # early return if there aren't enough data points in Y
        if (k > m) return(NA_real_)

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




#' Bias-reduced generalized k-nearest-neighbour KL divergence estimation
#'
#' This is the bias-reduced generalized k-NN based KL divergence estimator from
#' Wang et al. (2009) specified in Eq.(29).
#'
#' Finite sample bias reduction is achieved by an adaptive choice of the number
#' of nearest neighbours. Fixing the number of nearest neighbours upfront, as
#' done in [kld_est_nn()], may result in very different distances
#' \eqn{\rho^l_i,\nu^k_i} of a datapoint \eqn{x_i} to its \eqn{l}-th nearest
#' neighbours in \eqn{X} and \eqn{k}-th nearest neighbours in \eqn{Y},
#' respectively, which may lead to unequal biases in NN density estimation,
#' especially in a high-dimensional setting.
#' To overcome this issue, the number of neighbours \eqn{l,k} are here chosen
#' in a way to render \eqn{\rho^l_i,\nu^k_i} comparable, by taking the largest
#' possible number of neighbours \eqn{l_i,k_i} smaller than
#' \eqn{\delta_i:=\max(\rho^1_i,\nu^1_i)}.
#'
#' Since the bias reduction explicitly uses both samples `X` and `Y`, one-sample
#' estimation is not possible using this method.
#'
#' Reference:
#' Wang, Kulkarni and Verdú, "Divergence Estimation for Multidimensional
#' Densities Via k-Nearest-Neighbor Distances", IEEE Transactions on Information
#' Theory, Vol. 55, No. 5 (2009). DOI: https://doi.org/10.1109/TIT.2009.2016060
#'
#' @inherit kld_est_nn params return examples
#' @param max.k Maximum numbers of nearest neighbours to compute (default: `100`).
#'    A larger `max.k` may yield a more accurate KL-D estimate (see `warn.max.k`),
#'    but will always increase the computational cost.
#' @param warn.max.k If `TRUE` (the default), warns if `max.k` is such that more
#'    than `max.k` neighbours are within the neighbourhood \eqn{\delta} for some
#'    data point(s). In this case, only the first `max.k` neighbours are counted.
#'    As a consequence, `max.k` may required to be increased.
#' @export
kld_est_brnn <- function(X, Y, max.k = 100, warn.max.k = TRUE, eps = 0) {

    # get important dimensions
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    d <- ncol(X) # number of dimensions, must be the same in X and Y
    n <- nrow(X) # number of samples in X
    m <- nrow(Y) # number of samples in Y


    # early return if there aren't enough data points in X or Y
    if (n <= 1 || m == 0) return(NA_real_)

    # get distances to nearest neighbours from a kd-tree using nn2 from the RANN package
    nnXX <- RANN::nn2(X, X, k = max(min(max.k+1,n),2), eps = eps)$nn.dists[,-1, drop = FALSE]    # delete the NN in X to X, which is X itself
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

    # distance to the k-th / l-th nearest neighbours (similar order of magnitude)
    rho_l <- vapply(1:n, function(i) nnXX[i,l[i]], 1)
    nu_k <- vapply(1:n, function(i) nnYX[i,k[i]], 1)

    # equation (17) from Wang et al.
    d/n*sum(log(nu_k/rho_l)) + mean(digamma(l) - digamma(k)) + log(m/(n-1))
}

