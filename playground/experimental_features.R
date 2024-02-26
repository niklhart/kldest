# experimental features not included into the current version.


#' Nearest-neighbour density ratio based KL divergence estimator.
#'
#' Currently, this estimator is not part of the kldest package. Its performance
#' for small samples is poor, since a single point in `Y` for which no point in
#' `X` is closer than the `k`-th nearest neighbour in `Y` suffices for the KL
#' divergence estimate to be `-Inf`.
#'
#' @inherit kld_est_kde return
#' @inheritParams kld_est_nn
#' @param k Number of nearest neighbours to consider for NN density estimation.
#'    Defaults to `k = neff^(1/(d+1))`, where `neff = min(n,m)`. The choice of
#'    `k` is quite important for accuracy and the default may not be optimal.
#' @examples
#' # 1D example
#' set.seed(0)
#' X <- rnorm(100)
#' Y <- rnorm(100, mean = 1)
#' kld_est_nndr(X,Y)
kld_est_nndr <- function(X, Y, k = NULL, eps = 0) {

    # get important dimensions
    X <- as.matrix(X)
    Y <- as.matrix(Y)  # number of dimensions must be the same in X and Y
    d <- ncol(X)       # number of dimensions
    n <- nrow(X)       # number of samples in X
    m <- nrow(Y)       # number of samples in Y

    # default number of nearest neighbours (Noshad et al., 2017)
    if (is.null(k)) k <- ceiling(min(n,m)^(1/(d+1)))

    # combined dataset
    Z <- rbind(X,Y)

    idx.ZY <- RANN::nn2(Z, Y, k = k, eps = eps)$nn.idx
    Mi <- rowSums(idx.ZY > n)
    Ni <- k - Mi

    # return
    -1/n * sum(log(Ni / (Mi + 1)))

}




#' Zhang/Grabchak KL divergence estimator for samples from discrete distributions
#'
#' CAVE: not implemented yet!
#'
#' The estimator is that from Eq. (1.3) from Zhang and Grabchak (2014).
#'
#' Reference:
#' Zhang and Grabchak, "Nonparametric Estimation of Kullback-Leibler Divergence",
#' Neural Computation 26, 2570–2593 (2014).
#'
#' @inherit kld_est_kde return
#' @param X,Y Two samples from discrete distributions, specified as vectors,
#'    matrices or data frames.
#' @examples
#' # 1D example
#' X <- c(rep('M',5),rep('F',5))
#' Y <- c(rep('M',6),rep('F',4))
#' kld_est_discrete(X, Y)
#' @export
kld_est_zg2014 <- function(X, Y) {

    # get important dimensions
    X <- as.data.frame(X)
    Y <- as.data.frame(Y)
    d <- ncol(X) # number of dimensions, must be the same in X and Y
    n <- nrow(X) # number of samples in X
    m <- nrow(Y) # number of samples in Y

    # check consistency of dimensions
    stopifnot(ncol(Y) == d)

    # ensure consistency of factor levels
    uXY <- if (d == 1) {
        list(unique(c(X[[1]],Y[[1]])))
    } else {
        mapply(function(x,y) unique(c(x,y)),
               lapply(X, unique),
               lapply(Y, unique))
    }

    # convert all columns to factors to correctly account for missing levels
    X[] <- lapply(1:d, function(k) factor(X[[k]], levels = uXY[[k]]))
    Y[] <- lapply(1:d, function(k) factor(Y[[k]], levels = uXY[[k]]))

    # frequency counts via tables
    tX <- table(X)
    tY <- table(Y)

    pX <- stop("Implement this!")

    # estimate KL divergencce via relative frequencies
    kl_div_discrete(tX/n, tY/m)
}


#' Generalized k-nearest neighbour KL divergence estimator.
#'
#' This function implements the generalized k-nearest neighbour estimator in
#' Wang et al. (2009), Eq.(17). In this estimator, the number of nearest
#' neighbours to consider may differ between samples and sample points.
#'
#' The method is fully functional, but since it is not really useful in practice,
#' I haven't included it in the package.
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
#'   with `kld_est_nn(X, ..., k = k).
#'
kld_est_gnn <- function(X, Y = NULL, q = NULL, l = k, k = 1, eps = 0, log.q = FALSE) {

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
    two.sample <- is_two_sample(Y,q)

    if (two.sample) {
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
