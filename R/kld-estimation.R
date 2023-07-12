# Estimation algorithms for Kullback-Leibler divergence

#' Density-based estimation of Kullback-Leibler divergence
#'
#' @description
#' This estimation method approximates the densities of the unknown distributions
#' \eqn{P} and \eqn{Q} by a kernel density estimate, using a sample size- and
#' dimension-dependent bandwidth parameter and a Gaussian kernel.
#'
#' @param X,Y  `n`-by-`d` and `m`-by-`d` matrices, representing `n` samples from
#'    the true distribution \eqn{P} and `m` samples from the approximate distribution
#'    \eqn{Q}, both in `d` dimensions. Vector input is treated as a column matrix.
#' @param hX,hY Positive scalars or length `d` vectors, representing bandwidth
#'    parameters (possibly different in each component) for the density estimates
#'    of P and Q, respectively. The default `NULL` uses Scott's rule,
#'    \deqn{h_i = \frac{\sigma_i}{n^{1/(d+4)}}}
#' @returns A scalar, the estimated Kullback-Leibler divergence \eqn{D_{KL}(P||Q)}.
#' @examples
#' # KL-D between two samples from 1D Gaussians:
#' X <- rnorm(100)
#' Y <- rnorm(100, mean = 1, sd = 2)
#' kl_div_gaussian(mu1 = 0, sigma1 = 1, mu2 = 1, sigma2 = 2^2)
#' kldest_density(X,Y)
#' # KL-D between two samples from 2D Gaussians:
#' X1 <- rnorm(100)
#' X2 <- rnorm(100)
#' Y1 <- rnorm(100)
#' Y2 <- Y1 + rnorm(100)
#' X <- cbind(X1,X2)
#' Y <- cbind(Y1,Y2)
#' kl_div_gaussian(mu1 = rep(0,2), sigma1 = diag(2), mu2 = rep(0,2),
#'                 sigma2 = matrix(c(1,1,1,2),nrow=2))
#' kldest_density(X,Y)
#' @export
kldest_density <- function(X, Y, hX = NULL, hY = NULL, rule = c("Silverman","Scott")) {

    # input processing
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    n <- nrow(X)
    d <- ncol(X)
    m <- nrow(Y)
    stopifnot(ncol(Y) == d)

    # Use heuristics rule for bandwidth parameters, if unspecified
    rule <- match.arg(rule)
    if (is.null(hX)) {
        sd_X <- apply(X, MARGIN = 2, FUN = sd)
        hX <- switch(rule,
                     Silverman = sd_X*(4/((d+2)*n))^(1/(d+4)),
                     Scott     =  sd_X/n^(1/(d+4))
        )
    }
    if (is.null(hY)) {
        sd_Y <- apply(Y, MARGIN = 2, FUN = sd)
        hY <- switch(rule,
                     Silverman = sd_Y*(4/((d+2)*n))^(1/(d+4)),
                     Scott     =  sd_Y/n^(1/(d+4))
        )
    }

    p_X <- numeric(n)
    q_X <- numeric(n)

    for (i in 1:n) {
        kXi <- vapply(seq_len(d),
                     FUN = function(k) dnorm(X[i,k], mean = X[ ,k], sd = hX[k]),
                     FUN.VALUE = numeric(n))
        kYi <- vapply(seq_len(d),
                     FUN = function(k) dnorm(X[i,k], mean = Y[ ,k], sd = hY[k]),
                     FUN.VALUE = numeric(m))

        p_X[i] <- mean(apply(kXi, MARGIN = 1, FUN = prod))
        q_X[i] <- mean(apply(kYi, MARGIN = 1, FUN = prod))
    }

    # KL-divergence using estimated densities
    mean(log(p_X/q_X))

}


#' Vectorized version of kldest_density, currently under development
#' @export
kldest_density_vectorized <- function(X, Y, hX = NULL, hY = NULL, rule = c("Silverman","Scott")) {

    # input processing
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    n <- nrow(X)
    d <- ncol(X)
    m <- nrow(Y)
    stopifnot(ncol(Y) == d)

    # Use heuristics rule for bandwidth parameters, if unspecified
    rule <- match.arg(rule)
    if (is.null(hX)) {
        sd_X <- apply(X, MARGIN = 2, FUN = sd)
        hX <- switch(rule,
                     Silverman = sd_X*(4/((d+2)*n))^(1/(d+4)),
                     Scott     =  sd_X/n^(1/(d+4))
        )
    }
    if (is.null(hY)) {
        sd_Y <- apply(Y, MARGIN = 2, FUN = sd)
        hY <- switch(rule,
                     Silverman = sd_Y*(4/((d+2)*n))^(1/(d+4)),
                     Scott     =  sd_Y/n^(1/(d+4))
        )
    }

    kXX <- vapply(1:d,
                  FUN = function(k) dnorm(outer(X[ ,k], X[ ,k], FUN = "-"),
                                          sd = hX[k]),
                  FUN.VALUE = matrix(0,nrow=n,ncol=n))
    kXY <- vapply(1:d,
                  FUN = function(k) dnorm(outer(X[ ,k], Y[ ,k], FUN = "-"),
                                          sd = hY[k]),
                  FUN.VALUE = matrix(0,nrow=n,ncol=m))

    p_X <- rowMeans(apply(kXX, MARGIN = 1:2, FUN = prod))
    q_X <- rowMeans(apply(kXY, MARGIN = 1:2, FUN = prod))

    # KL-divergence using estimated densities
    mean(log(p_X/q_X))

}



#' 1-D density-based estimation of Kullback-Leibler divergence
#'
#' @description
#' This estimation method approximates the densities of the unknown distributions
#' \eqn{P} and \eqn{Q} by a kernel density estimate using function 'density'.
#'
#' @param X,Y  Numeric vectors or single-column matrices, representing samples
#'    from the true distribution \eqn{P} and the approximate distribution
#'    \eqn{Q}, respectively.
#' @returns A scalar, the estimated Kullback-Leibler divergence \eqn{D_{KL}(P||Q)}.
#' @examples
#' # KL-D between two samples from 1D Gaussians:
#' X <- rnorm(100)
#' Y <- rnorm(100, mean = 1, sd = 2)
#' kl_div_gaussian(mu1 = 0, sigma1 = 1, mu2 = 1, sigma2 = 2^2)
#' kldest_density1(X,Y)
#' @export
kldest_density1 <- function(X, Y) {

    # input processing
    X <- as.vector(X)
    Y <- as.vector(Y)

    # using base R's density function
    pX <- approxfun(density(X))
    pY <- approxfun(density(Y))

    # KL-divergence using estimated densities
    mean(log(pX(X)/pY(X)))
}

#' Kernel density estimation in 1-6 dimensions using package 'ks'.
#'
#' NOTE: seems to be extremely slow...
kldest_density6 <- function(X, Y) {

    # input processing
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    n <- nrow(X)
    d <- ncol(X)
    m <- nrow(Y)
    stopifnot(ncol(Y) == d)

    # KL-divergence using estimated densities
    mean(log(pX$estimate/predict(pY,x=X)))
}

#' Universal 1-nearest neighbour divergence estimator from Wang et al. (2009).
#'
#' Direct implementation of the universal 1-NN divergence estimator from eq. (5)
#' in Wang et al. (2009). The estimator is identical to the one by Perez-Cruz
#' (2008).
#'
#' Code source: https://gist.github.com/lars-von-buchholtz/636f542ce8d93d5a14ae52a6c538ced5
#' Adapted for R from python code at https://gist.github.com/atabakd/ed0f7581f8510c8587bc2f41a094b518
#'
#' @inherit kldest_density params return
#' @export
#' @example examples/all-estimators.R
kl_universal_1nn <- function(X,Y) {

    # get important dimensions
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    d <- ncol(X) # number of dimensions, must be the same in X and Y
    n <- nrow(X) # number of samples in X
    m <- nrow(Y) # number of samples in Y

    # get distances to nearest neighbors from kdTree using nn2 from the RANN package
    r <- RANN::nn2(X,X, k=2, eps=.01)[[2]][,2] # get 2 closest neighbors, then take the second (the closest is the point itself) to get n x 1 matrix
    s <- RANN::nn2(Y,X, k=1, eps=.01)[[2]] # also n x 1 matrix

    # There is a mistake in the paper. In Eq. 14, the right side misses a negative sign
    # on the first term of the right hand side.
    return (- sum(log(r/s)) * d / n + log(m / (n - 1.)))
}

# Generalized nearest-neighbour (Wang)
#' Generalized KL divergence estimation method from Wang et al. (2009)
#'
#' This is the generalized k-NN based KL divergence estimator, eq. (29).
#' `TODO`: more details on the algorithm!
#'
#' Reference: https://doi.org/10.1109/TIT.2009.2016060
#'
#' @inherit kldest_density params return
#' @param max.k Maximum numbers of nearest neighbours to compute (default: `50`)
#' @param warn.max.k If `TRUE` (the default), warn if `max.k` is such that more
#'    than `max.k` neighbours are within the `eps` neighbourhood for some data
#'    points. In this case, only the first `max.k` neighbours will be counted.
#'    As a consequence, `max.k` may required to be increased.
#' @export
kl_generalized_knn_eps <- function(X, Y, max.k = 50, warn.max.k = TRUE) {

    # get important dimensions
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    d <- ncol(X) # number of dimensions, must be the same in X and Y
    n <- nrow(X) # number of samples in X
    m <- nrow(Y) # number of samples in Y

    # check consistency of dimensions
    stopifnot(ncol(Y) == d, max.k < min(n,m))

    # get distances to nearest neighbors from kdTree using nn2 from the RANN package
    nnXX <- RANN::nn2(X,X, k=max.k+1, eps=.01)$nn.dist[,-1, drop = FALSE]    # delete the NN in X to X, which is X itself
    nnYX <- RANN::nn2(Y,X, k=max.k, eps=.01)$nn.dist

    # maximum of nearest neighbours in X and Y
    eps <- pmax(nnXX[,1], nnYX[,1])
    #    print(eps)

    # how many data points are closer than eps?
    l <- rowSums(nnXX <= eps)
    k <- rowSums(nnYX <= eps)

    # check that max.k was chosen large enough and produce a warning otherwise
    if (warn.max.k) {
        tooSmallMaxK <- max.k %in% c(l,k)
        if (tooSmallMaxK) warning('Parameter "eps" too large, increase input "max.k".')
    }

    # distance to the k-th / l-th nearest neighbours
    rho_l <- vapply(1:n, function(i) nnXX[i,l[i]], 1)
    nu_k <- vapply(1:n, function(i) nnYX[i,k[i]], 1)

    # equation (17) from Wang et al.
    d/n*sum(log(nu_k/rho_l)) + mean(digamma(l) - digamma(k)) + log(m/(n-1))
}


# Monte Carlo (using density q)?


