# Estimation algorithms for Kullback-Leibler divergence

# Density-based algorithm

# + Scott's rule for bandwidth selection: h = sqrt(V(p_j))/n^(1/(d+4)),
#   with variance V of marginal p_j, n number of data points and d the dimension.
# + choice of kernel function less important... --> Gaussian
# + estimate f and g by f~(x) = (1/(n*h^d)) * sum(K((x-xi)/h))
# + if available, use density q

#' Density-based estimation of Kullback-Leibler divergence
#'
#' @description
#' This estimation method approximates the unknown densities p and q by a kernel
#' density estimate, using a sample size- and dimension-dependent bandwidth parameter
#' and a Gaussian kernel.
#'
#' @param X A numeric vector or matrix, representing the sample from the true
#'    distribution.
#' @param Y A numeric vector or matrix, representing the sample from the true
#'    distribution. If \code{X} and \code{Y} matrices, they must have the same
#'    number of columns \code{d}.
#' @param h A positive scalar or length \code{d} vector, representing the bandwidth
#'    parameter (possibly different in each component). The default \code{NULL} uses
#'    Scott's rule (EMBED LATEX?)
#' @param K A density function representing the kernel to be used. Defaults to a Gaussian kernel.
#' @returns A scalar, the estimated Kullback-Leibler divergence
#' @examples
#' # KL-D between two samples from 1D Gaussians:
#' X <- rnorm(100)
#' Y <- rnorm(100, mean = 1, sd = 2)
#' kl_divergence_gaussian(mu1 = 0, sigma1 = 1, mu2 = 1, sigma2 = 2^2)
#' kldest_density(X,Y)
#' @export
kldest_density <- function(X, Y, h = NULL, K = dnorm) {

    X <- as.matrix(X)
    Y <- as.matrix(Y)
    n <- nrow(X)
    d <- ncol(X)
    m <- nrow(Y)
    stopifnot(d == ncol(Y))

    # Scott's rule for bandwidth parameter
    sd_X <- apply(X, MARGIN = 2, FUN = sd)
    sd_Y <- apply(Y, MARGIN = 2, FUN = sd)
    hX <- sd_X/n^(1/(d+4))     # Scott's rule
    hY <- sd_Y/m^(1/(d+4))     # Scott's rule

    # Gaussian kernel
    kX <- function(x) t(dnorm(t(x-X), sd = hX))
    kY <- function(x) t(dnorm(t(x-Y), sd = hY))

    # density estimate
    pX <- Vectorize(function(x) colMeans(kX(x)))
    pY <- Vectorize(function(y) colMeans(kY(y)))

    # KL-divergence using estimated densities
    mean(log(pX(X)/pY(X)))
}

# Nearest-neighbour algorithm (plain)


# Generalized nearest-neighbour (Wang)


# Monte Carlo (using density q)?


