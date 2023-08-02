# Kullback-Leibler divergence estimators between continuous distributions

#' Kernel density-based estimation of Kullback-Leibler divergence
#'
#' Disclaimer: this function doesn't use binning / FFT and hence, it is very slow!
#'
#' @description
#' This estimation method approximates the densities of the unknown distributions
#' \eqn{P} and \eqn{Q} by a kernel density estimate, using a sample size- and
#' dimension-dependent bandwidth parameter and a Gaussian kernel. It works for
#' any number of dimensions but is very slow.
#'
#' @param X,Y  `n`-by-`d` and `m`-by-`d` matrices, representing `n` samples from
#'    the true distribution \eqn{P} and `m` samples from the approximate distribution
#'    \eqn{Q}, both in `d` dimensions. Vector input is treated as a column matrix.
#' @param hX,hY Positive scalars or length `d` vectors, representing bandwidth
#'    parameters (possibly different in each component) for the density estimates
#'    of P and Q, respectively. If unspecified, a heurestic specified via the `rule`
#'    argument is used.
#' @param rule A heuristic for computing arguments `hX` and/or `hY`. The default
#'    `"silverman"` is Silverman's rule
#'    \deqn{h_i = \sigma_i\left(\frac{4}{(2+d)n}\right)^{1/(d+4)}.}
#'    As an alternative, Scott's rule `"scott"` can be used,
#'    \deqn{h_i = \frac{\sigma_i}{n^{1/(d+4)}}.}
#' @returns A scalar, the estimated Kullback-Leibler divergence \eqn{D_{KL}(P||Q)}.
#' @example examples/continuous-estimators.R
#' @export
kld_est_kde <- function(X, Y, hX = NULL, hY = NULL, rule = c("Silverman","Scott")) {

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


#' 1-D kernel density-based estimation of Kullback-Leibler divergence
#'
#' @description
#' This estimation method approximates the densities of the unknown distributions
#' \eqn{P} and \eqn{Q} by a kernel density estimate using function 'density' from
#' package 'stats'.
#'
#' @inherit kld_est_kde return
#' @param X,Y  Numeric vectors or single-column matrices, representing samples
#'    from the true distribution \eqn{P} and the approximate distribution
#'    \eqn{Q}, respectively.
#' @param MC A boolean: use a Monte Carlo approximation instead of numerical
#'    integration via the trapezoidal rule (default: `FALSE`)?
#' @param ... Further parameters to passed on to function `stats::density` (e.g.,
#'    argument `kernel`)
#' @examples
#' # KL-D between two samples from 1D Gaussians:
#' X <- rnorm(100)
#' Y <- rnorm(100, mean = 1, sd = 2)
#' kl_div_gaussian(mu1 = 0, sigma1 = 1, mu2 = 1, sigma2 = 2^2)
#' kldest_density1(X,Y)
#' kldest_density1(X,Y, MC = TRUE)
#' @export
kld_est_kde1 <- function(X, Y, MC = FALSE, ...) {

    # input processing
    X <- as.vector(X)
    Y <- as.vector(Y)

    # using base R's density function
    dX <- density(X, ...)
    dY <- density(Y, from = head(dX$x,1), to = tail(dX$x,1), ...)

    if (MC) {
        pX <- approxfun(dX)
        pY <- approxfun(dY)

        # MC average
        mean(log(pX(X)/pY(X)))

    } else {
        # trapezoidal integration
        trapz(h  = dX$x[2] - dX$x[1],
              fx = dX$y*log(dX$y/dY$y))
    }
}

#' Density-ratio-based estimation of Kullback-Leibler divergence
#'
#' @description
#' This estimation method approximates the ratio of densities of the unknown
#' distributions \eqn{P} and \eqn{Q}  using package 'densratio'.
#'
#' @inherit kld_est_kde params return
#' @param ... Additional options to be passed to function `densratio`.
#' @examples
#' # KL-D between two samples from 1D Gaussians:
#' X <- rnorm(100)
#' Y <- rnorm(100, mean = 1, sd = 2)
#' kl_div_gaussian(mu1 = 0, sigma1 = 1, mu2 = 1, sigma2 = 2^2)
#' kldest_densratio(X,Y)
#' kldest_densratio(X,Y, method = "KLIEP")
#' kldest_densratio(X,Y, method = "RuLSIF")
#' @export
kld_est_densratio <- function(X, Y, ...) {

    # input processing
    X <- as.matrix(X)
    Y <- as.matrix(Y)

    # dimensionality check
    stopifnot(ncol(X) == ncol(Y))

    # using package 'densratio'
    r <- densratio::densratio(x1 = X, x2 = Y, verbose = FALSE,...)

    rX <- r$compute_density_ratio(X)

    # MC average
    mean(log(rX))

}

#' 2-D density-based estimation of Kullback-Leibler divergence
#'
#' @description
#' This estimation method approximates the densities of the unknown distributions
#' \eqn{P} and \eqn{Q} by a kernel density estimate using function 'bkde' from
#' package 'KernSmooth'.
#'
#' @param X,Y  Two-column matrices, representing samples from the bivariate true
#'    distribution \eqn{P} and approximate distribution \eqn{Q}, respectively.
#' @param MC A boolean: use a Monte Carlo approximation instead of numerical
#'    integration via the trapezoidal rule (default: `FALSE`)? Currently, this
#'    option is not implemented, i.e. a value of `TRUE` results in an error.
#' @param hX,hY Bandwidths for the kernel density estimates of \eqn{P} and \eqn{Q},
#'    respectively. The default `NULL` means they are determined by argument `rule`.
#' @param rule A heuristic to derive parameters `hX` and `hY`, default is
#'    `"Silverman", which means that
#'    `\deqn{h_i = \sigma_i\left(\frac{4}{(2+d)n}\right)^{1/(d+4)}.}
#' @param eps A nonnegative scalar; if `eps > 0`, \eqn{Q} is estimated as a mixture
#'    between the kernel density estimate and a uniform distribution on the computational
#'    grid. The weight of the uniform component is `eps` times the maximum density
#'    estimate of \eqn{Q}. Defaults to `eps = 1e-5`.
#' @returns A scalar, the estimated Kullback-Leibler divergence \eqn{D_{KL}(P||Q)}.
#' @examples
#' # KL-D between two samples from 2-D Gaussians:
#' X1 <- rnorm(1000)
#' X2 <- rnorm(1000)
#' Y1 <- rnorm(1000)
#' Y2 <- Y1 + rnorm(1000)
#' X <- cbind(X1,X2)
#' Y <- cbind(Y1,Y2)
#' kl_div_gaussian(mu1 = rep(0,2), sigma1 = diag(2), mu2 = rep(0,2),
#'                 sigma2 = matrix(c(1,1,1,2),nrow=2))
#' kldest_density2(X,Y)
#' # kldest_density2(X,Y, MC = TRUE)
#' @export
kld_est_kde2 <- function(X, Y, MC = FALSE, hX = NULL, hY = NULL,
                        rule = c("Silverman","Scott"), eps = 1e-5) {

    # input processing
    X <- as.matrix(X)
    Y <- as.matrix(Y)

    # check that we are indeed in the 2D case
    d <- 2L
    stopifnot(ncol(X) == d, ncol(Y) == d)

    n <- nrow(X)
    m <- nrow(Y)

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
                     Silverman = sd_Y*(4/((d+2)*m))^(1/(d+4)),
                     Scott     =  sd_Y/m^(1/(d+4))
        )
    }

    # using 2D kernel density estimation implemented in package 'KernSmooth'
    range.x <- list(range(X[,1]),range(X[,2]))
    dX <- KernSmooth::bkde2D(x = X, bandwidth = hX, range.x = range.x)
    dY <- KernSmooth::bkde2D(x = Y, bandwidth = hY, range.x = range.x, truncate = FALSE)

    if (MC) {
        stop("Not implemented yet...")
#        pX <- approxfun(dX)
#        pY <- approxfun(dY)
#        mean(log(pX(X)/pY(X)))

    } else {
        # numerical integration
        h <- c(dX$x1[2]-dX$x1[1], dX$x2[2] - dX$x2[1])
        if (eps > 0) {
            dY$fhat[dY$fhat == 0] <- eps * max(dY$fhat)
            normY <- trapz(h = h, fx = dY$fhat)
            dY$fhat <- dY$fhat/normY
        }
        fX <- dX$fhat * log(dX$fhat/dY$fhat)
        fX[dX$fhat == 0] <- 0
        trapz(h = h, fx = fX)
    }

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
#' @inherit kld_est_kde params return examples
#' @export
kld_est_1nn <- function(X,Y) {

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


#' Generalized k-nearest neighbour KL divergence estimator by Wang et al. (2009)
#' Equation (17)
#'
#' Reference:
#' Wang, Kulkarni and Verdú, "Divergence Estimation for Multidimensional
#' Densities Via k-Nearest-Neighbor Distances", IEEE Transactions on Information
#' Theory, Vol. 55, No. 5 (2009).
#'
#' @inherit kld_est_kde params return examples
#' @param l,k Scalars or numeric vectors of length `n` and `m`, representing the
#'   number of nearest neighbours to use for nearest neighbour density estimation
#'   of P and Q, respectively, for each of the data points \code{X[i,]}, with
#'   \code{i} ranging from \code{1} to \code{n}, and \code{Y[j,]}, with
#'   \code{j} ranging from \code{1} to \code{m}. The default is `l = k = 1`.
#' @returns A scalar, the estimated Kullback-Leibler divergence D(P||Q).
kld_est_knn <- function(X, Y, l = k, k = 1) {

    # get important dimensions
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    d <- ncol(X) # number of dimensions, must be the same in X and Y
    n <- nrow(X) # number of samples in X
    m <- nrow(Y) # number of samples in Y

    # check consistency of dimensions
    if (length(l) == 1) l <- rep(l, n)
    if (length(k) == 1) k <- rep(k, m)
    stopifnot(ncol(Y) == d, length(l) == n, length(k) == n, max(l) < m, max(k) <= n)

    # get distances to nearest neighbors from kdTree using nn2 from the RANN package
    nnXX <- RANN::nn2(X, X, k = max(l)+1, eps = .01)$nn.dists[ ,-1, drop = FALSE]    # delete the NN in X to X, which is X itself
    nnYX <- RANN::nn2(Y, X, k = max(k), eps = .01)$nn.dists

    # distance to the k-th / l-th nearest neighbours
    rho_l <- vapply(1:n, function(i) nnXX[i,l[i]], 1)
    nu_k <- vapply(1:n, function(i) nnYX[i,k[i]], 1)

    # equation (17) from Wang et al.
    d/n*sum(log(nu_k/rho_l)) + mean(digamma(l) - digamma(k)) + log(m/(n-1))
}


#' Generalized nearest-neighbour KL divergence estimation from Wang et al. (2009)
#'
#' This is the generalized k-NN based KL divergence estimator, eq. (29).
#' `TODO`: more details on the algorithm!
#'
#' Reference:
#' Wang, Kulkarni and Verdú, "Divergence Estimation for Multidimensional
#' Densities Via k-Nearest-Neighbor Distances", IEEE Transactions on Information
#' Theory, Vol. 55, No. 5 (2009). DOI: https://doi.org/10.1109/TIT.2009.2016060
#'
#' @inherit kld_est_kde params return examples
#' @param max.k Maximum numbers of nearest neighbours to compute (default: `50`).
#'    A larger `max.k` may yield a more accurate KL-D estimate (see `warn.max.k`),
#'    but will always increase the computational cost.
#' @param warn.max.k If `TRUE` (the default), warn if `max.k` is such that more
#'    than `max.k` neighbours are within the `eps` neighbourhood for some data
#'    points. In this case, only the first `max.k` neighbours will be counted.
#'    As a consequence, `max.k` may required to be increased.
#' @export
kld_est_gknn <- function(X, Y, max.k = 50, warn.max.k = TRUE) {

    # get important dimensions
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    d <- ncol(X) # number of dimensions, must be the same in X and Y
    n <- nrow(X) # number of samples in X
    m <- nrow(Y) # number of samples in Y

    # check consistency of dimensions
    stopifnot(ncol(Y) == d)

    # get distances to nearest neighbors from kdTree using nn2 from the RANN package
    nnXX <- RANN::nn2(X, X, k = min(max.k+1,n), eps = .01)$nn.dists[,-1, drop = FALSE]    # delete the NN in X to X, which is X itself
    nnYX <- RANN::nn2(Y, X, k = min(max.k,m), eps = .01)$nn.dists

    # maximum of nearest neighbours in X and Y
    eps <- pmax(nnXX[,1], nnYX[,1])

    # how many data points are closer than eps?
    l <- rowSums(nnXX <= eps)
    k <- rowSums(nnYX <= eps)

    # check that max.k was chosen large enough and produce a warning otherwise
    if (warn.max.k && max.k %in% c(l,k)) {
        warning(paste('The "max.k"-th neighbour is inside the "eps" environment',
                      'for some datapoint, and further neighbours are ignored.',
                      'Increase input "max.k" to investigate the impact of this.'))
    }

    # distance to the k-th / l-th nearest neighbours
    rho_l <- vapply(1:n, function(i) nnXX[i,l[i]], 1)
    nu_k <- vapply(1:n, function(i) nnYX[i,k[i]], 1)

    # equation (17) from Wang et al.
    d/n*sum(log(nu_k/rho_l)) + mean(digamma(l) - digamma(k)) + log(m/(n-1))
}





